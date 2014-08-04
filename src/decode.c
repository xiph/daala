/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "decint.h"
#include "generic_code.h"
#include "laplace_code.h"
#include "filter.h"
#include "dct.h"
#include "intra.h"
#include "partition.h"
#include "pvq.h"
#include "pvq_code.h"
#include "block_size.h"
#include "block_size_dec.h"
#include "tf.h"
#include "state.h"

static int od_dec_init(od_dec_ctx *dec, const daala_info *info,
 const daala_setup_info *setup) {
  int ret;
  (void)setup;
  ret = od_state_init(&dec->state, info);
  if (ret < 0) return ret;
  dec->packet_state = OD_PACKET_DATA;
  return 0;
}

static void od_dec_clear(od_dec_ctx *dec) {
  od_state_clear(&dec->state);
}

daala_dec_ctx *daala_decode_alloc(const daala_info *info,
 const daala_setup_info *setup) {
  od_dec_ctx *dec;
  if (info == NULL) return NULL;
  dec = (od_dec_ctx *)_ogg_malloc(sizeof(*dec));
  if (od_dec_init(dec, info, setup) < 0) {
    _ogg_free(dec);
    return NULL;
  }
  return dec;
}

void daala_decode_free(daala_dec_ctx *dec) {
  if (dec != NULL) {
    od_dec_clear(dec);
    _ogg_free(dec);
  }
}

int daala_decode_ctl(daala_dec_ctx *dec, int req, void *buf, size_t buf_sz) {
  (void)dec;
  (void)buf;
  (void)buf_sz;
  switch (req) {
    default: return OD_EIMPL;
  }
}

static void od_decode_mv(daala_dec_ctx *dec, od_mv_grid_pt *mvg, int vx,
 int vy, int level, int mv_res, int width, int height) {
  int ex;
  int ey;
  generic_encoder *model;
  int *mv_ex;
  int *mv_ey;
  int pred[2];
  int ox;
  int oy;
  int id;
  mv_ex = dec->adapt.mv_ex;
  mv_ey = dec->adapt.mv_ey;
  model = &dec->adapt.mv_model;
  ex = mv_ex[level] >> mv_res;
  ey = mv_ex[level] >> mv_res;
  id = od_decode_cdf_adapt(&dec->ec, dec->adapt.mv_small_cdf, 16,
   dec->adapt.mv_small_increment);
  oy = id >> 2;
  ox = id & 0x3;
  if (ox == 3) ox += generic_decode(&dec->ec, model, width << (3 - mv_res), &ex, 2);
  if (oy == 3) oy += generic_decode(&dec->ec, model, height << (3 - mv_res), &ey, 2);
  mv_ex[level] -= (mv_ex[level] - (ox << mv_res << 16)) >> 6;
  mv_ey[level] -= (mv_ey[level] - (oy << mv_res << 16)) >> 6;
  if (ox && od_ec_dec_bits(&dec->ec, 1)) ox = -ox;
  if (oy && od_ec_dec_bits(&dec->ec, 1)) oy = -oy;
  od_state_get_predictor(&dec->state, pred, vx, vy, level, mv_res);
  mvg->mv[0] = (pred[0] + ox) << mv_res;
  mvg->mv[1] = (pred[1] + oy) << mv_res;
}

struct od_mb_dec_ctx {
  signed char *modes[OD_NPLANES_MAX];
  od_coeff *c;
  od_coeff **d;
  /* holds a TF'd copy of the transform coefficients in 4x4 blocks */
  od_coeff *tf[OD_NPLANES_MAX];
  od_coeff *md;
  od_coeff *mc;
  od_coeff *l;
  int run_pvq[OD_NPLANES_MAX];
  int is_keyframe;
  int nk;
  int k_total;
  int sum_ex_total_q8;
  int ncount;
  int count_total_q8;
  int count_ex_total_q8;
};
typedef struct od_mb_dec_ctx od_mb_dec_ctx;

static void od_decode_compute_pred(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, od_coeff *pred,
  int ln, int pli, int bx, int by, int has_ur) {
  int n;
  int n2;
  int xdec;
  int ydec;
  int w;
  int frame_width;
  signed char *modes;
  od_coeff *d;
  od_coeff *tf;
  od_coeff *md;
  od_coeff *l;
  int x;
  int y;
  int zzi;
  OD_ASSERT(ln >= 0 && ln <= 2);
  n = 1 << (ln + 2);
  n2 = n*n;
  xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
  frame_width = dec->state.frame_width;
  w = frame_width >> xdec;
  modes = ctx->modes[OD_DISABLE_CFL ? pli : 0];
  d = ctx->d[pli];
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  tf = ctx->tf[pli];
  md = ctx->md;
  l = ctx->l;
  if (ctx->is_keyframe) {
    if (bx > 0 && by > 0) {
      if (pli == 0 || OD_DISABLE_CFL) {
        ogg_uint16_t mode_cdf[OD_INTRA_NMODES];
        int m_l;
        int m_ul;
        int m_u;
        int mode;
        od_coeff *coeffs[4];
        int strides[4];
        /*Calculate the intra-prediction.*/
        coeffs[0] = tf + ((by - (1 << ln)) << 2)*w + ((bx - (1 << ln)) << 2);
        coeffs[1] = tf + ((by - (1 << ln)) << 2)*w + ((bx - (0 << ln)) << 2);
        coeffs[2] = tf + ((by - (1 << ln)) << 2)*w + ((bx + (1 << ln)) << 2);
        coeffs[3] = tf + ((by - (0 << ln)) << 2)*w + ((bx - (1 << ln)) << 2);
        if (!has_ur) {
          coeffs[2] = coeffs[1];
        }
        strides[0] = w;
        strides[1] = w;
        strides[2] = w;
        strides[3] = w;
        m_l = modes[by*(w >> 2) + bx - 1];
        m_ul = modes[(by - 1)*(w >> 2) + bx - 1];
        m_u = modes[(by - 1)*(w >> 2) + bx];
        od_intra_pred_cdf(mode_cdf, dec->adapt.mode_probs[pli],
         OD_INTRA_NMODES, m_l, m_ul, m_u);
#if OD_DISABLE_INTRA
        mode = 0;
#else
        mode = od_ec_decode_cdf_unscaled(&dec->ec, mode_cdf, OD_INTRA_NMODES);
#endif
        (*OD_INTRA_GET[ln])(pred, coeffs, strides, mode);
#if OD_DISABLE_INTRA
        OD_CLEAR(pred+1, n2-1);
#endif
        for (y = 0; y < (1 << ln); y++) {
          for (x = 0; x < (1 << ln); x++) {
            modes[(by + y)*(w >> 2) + bx + x] = mode;
          }
        }
        od_intra_pred_update(dec->adapt.mode_probs[pli], OD_INTRA_NMODES, mode,
         m_l, m_ul, m_u);
      }
      else {
        int mode;
        mode = modes[(by << ydec)*(frame_width >> 2) + (bx << xdec)];
        od_chroma_pred(pred, d, l, w, bx, by, ln, xdec, ydec,
         dec->state.bsize, dec->state.bstride,
         OD_INTRA_CHROMA_WEIGHTS_Q8[mode]);
      }
    }
    else {
      int nsize;
      for (zzi = 0; zzi < n2; zzi++) pred[zzi] = 0;
      nsize = ln;
      /*444/420 only right now.*/
      OD_ASSERT(xdec == ydec);
      if (bx > 0) {
        int noff;
        nsize = OD_BLOCK_SIZE4x4(dec->state.bsize, dec->state.bstride,
         (bx - 1) << xdec, by << ydec);
        nsize = OD_MAXI(nsize - xdec, 0);
        noff = 1 << nsize;
        /*Because of the quad-tree structure we can always find our neighbors
           starting offset by rounding to a multiple of his size.*/
        OD_ASSERT(!(bx & (noff - 1)));
        pred[0] = d[((by & ~(noff - 1)) << 2)*w + ((bx - noff) << 2)];
      }
      else if (by > 0) {
        int noff;
        nsize = OD_BLOCK_SIZE4x4(dec->state.bsize, dec->state.bstride,
         bx << xdec, (by - 1) << ydec);
        nsize = OD_MAXI(nsize - xdec, 0);
        noff = 1 << nsize;
        OD_ASSERT(!(by & (noff - 1)));
        pred[0] = d[((by - noff) << 2)*w + ((bx & ~(noff - 1)) << 2)];
      }
      /*Rescale DC for correct transform size.*/
      if (nsize > ln) pred[0] >>= (nsize - ln);
      else if (nsize < ln) pred[0] <<= (ln - nsize);
      if (pli == 0) {
        for (y = 0; y < (1 << ln); y++) {
          for (x = 0; x < (1 << ln); x++) {
            modes[(by + y)*(w >> 2) + bx + x] = 0;
          }
        }
      }
    }
  }
  else {
    int ci;
    ci = 0;
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        pred[ci++] = md[(y + (by << 2))*w + (x + (bx << 2))];
      }
    }
  }
}

static void od_single_band_scalar_decode(daala_dec_ctx *dec, int ln,
 od_coeff *pred, const od_coeff *predt, int q, int pli) {
  int *adapt;
  int vk;
  int zzi;
  int n2;
  ogg_int32_t adapt_curr[OD_NSB_ADAPT_CTXS];
  adapt = dec->adapt.pvq_adapt;
  n2 = 1 << (2*ln + 4);
  vk = generic_decode(&dec->ec, &dec->adapt.model_g[pli], -1,
   &dec->adapt.ex_g[pli][ln], 0);
  laplace_decode_vector(&dec->ec, pred + 1, n2 - 1, vk, adapt_curr, adapt);
  for (zzi = 1; zzi < n2; zzi++) {
    pred[zzi] = pred[zzi]*q + predt[zzi];
  }
  if (adapt_curr[OD_ADAPT_K_Q8] > 0) {
    adapt[OD_ADAPT_K_Q8] += 256*adapt_curr[OD_ADAPT_K_Q8] -
     adapt[OD_ADAPT_K_Q8] >> OD_SCALAR_ADAPT_SPEED;
    adapt[OD_ADAPT_SUM_EX_Q8] += adapt_curr[OD_ADAPT_SUM_EX_Q8] -
     adapt[OD_ADAPT_SUM_EX_Q8] >> OD_SCALAR_ADAPT_SPEED;
  }
  if (adapt_curr[OD_ADAPT_COUNT_Q8] > 0) {
    adapt[OD_ADAPT_COUNT_Q8] += adapt_curr[OD_ADAPT_COUNT_Q8]-
     adapt[OD_ADAPT_COUNT_Q8] >> OD_SCALAR_ADAPT_SPEED;
    adapt[OD_ADAPT_COUNT_EX_Q8] += adapt_curr[OD_ADAPT_COUNT_EX_Q8]-
     adapt[OD_ADAPT_COUNT_EX_Q8] >> OD_SCALAR_ADAPT_SPEED;
  }
}

void od_single_band_decode(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int ln,
 int pli, int bx, int by, int has_ur) {
  int n;
  int xdec;
  int w;
  int frame_width;
  od_coeff *c;
  od_coeff *d;
  od_coeff *tf;
  od_coeff *md;
  od_coeff *mc;
  od_coeff pred[16*16];
  od_coeff predt[16*16];
  int run_pvq;
  int quant;
  int dc_quant;
#ifndef USE_BAND_PARTITIONS
  unsigned char const *zig;
#endif
  OD_ASSERT(ln >= 0 && ln <= 2);
  n = 1 << (ln + 2);
  run_pvq = ctx->run_pvq[pli];
  bx <<= ln;
  by <<= ln;
#ifndef USE_BAND_PARTITIONS
  zig = OD_DCT_ZIGS[ln];
#endif
  xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = dec->state.frame_width;
  w = frame_width >> xdec;
  c = ctx->c;
  d = ctx->d[pli];
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  tf = ctx->tf[pli];
  md = ctx->md;
  mc = ctx->mc;
  /*Apply forward transform to MC predictor.*/
  if (!ctx->is_keyframe) {
    (*OD_FDCT_2D[ln])(md + (by << 2)*w + (bx << 2), w,
     mc + (by << 2)*w + (bx << 2), w);
  }
  od_decode_compute_pred(dec, ctx, pred, ln, pli, bx, by, has_ur);
#ifdef USE_BAND_PARTITIONS
  od_raster_to_coding_order(predt,  n, &pred[0], n, !run_pvq);
#else
  /*Zig-zag*/
  for (y = 0; y < n; y++) {
    for (x = 0; x < n; x++) {
      predt[zig[y*n + x]] = pred[y*n + x];
    }
  }
#endif
  quant = OD_MAXI(1, dec->quantizer[pli]);
  if (run_pvq)
    dc_quant = OD_MAXI(1, quant*OD_PVQ_QM_Q4[pli][ln][0] >> 4);
  else
    dc_quant = (pli==0 || dec->quantizer[pli]==0) ? quant : (quant + 1) >> 1;
  if (run_pvq) {
    pvq_decode(dec, predt, pred, quant, pli, ln,
     OD_PVQ_QM_Q4[pli][ln], OD_PVQ_BETA[pli][ln],
     OD_PVQ_INTER_BAND_MASKING[ln], ctx->is_keyframe);
  }
  else {
    od_single_band_scalar_decode(dec, ln, pred, predt, quant, pli);
  }
  if (OD_DISABLE_HAAR_DC || !ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && run_pvq;
    if (!has_dc_skip || pred[0]) {
      pred[0] = has_dc_skip + generic_decode(&dec->ec,
       &dec->adapt.model_dc[pli], -1, &dec->adapt.ex_dc[pli][ln][0], 2);
      if (pred[0]) pred[0] *= od_ec_dec_bits(&dec->ec, 1) ? -1 : 1;
    }
    pred[0] = pred[0]*dc_quant + predt[0];
  }
  else {
    pred[0] = d[((by << 2))*w + ((bx << 2))];
  }
#ifdef USE_BAND_PARTITIONS
  od_coding_order_to_raster(&d[((by << 2))*w + (bx << 2)], w, pred, n, !run_pvq);
#else
  /*De-zigzag*/
  for (y = 0; y < n; y++) {
    for (x = 0; x < n; x++) {
      d[((by << 2) + y)*w + (bx << 2) + x] = pred[zig[y*n + x]];
    }
  }
#endif
  /*Update the TF'd luma plane with CfL, or all the planes without CfL.*/
  if (ctx->is_keyframe && (pli == 0 || OD_DISABLE_CFL)) {
    od_convert_block_down(tf + (by << 2)*w + (bx << 2), w,
     d + (by << 2)*w + (bx << 2), w, ln, 0, 0);
  }
  /*Apply the inverse transform.*/
  (*OD_IDCT_2D[ln])(c + (by << 2)*w + (bx << 2), w,
   d + (by << 2)*w + (bx << 2), w);
}

static void od_32x32_decode(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int ln,
 int pli, int bx, int by, int has_ur) {
  bx <<= 1;
  by <<= 1;
  od_single_band_decode(dec, ctx, ln - 1, pli, bx + 0, by + 0, 1);
  od_single_band_decode(dec, ctx, ln - 1, pli, bx + 1, by + 0, has_ur);
  od_single_band_decode(dec, ctx, ln - 1, pli, bx + 0, by + 1, 1);
  od_single_band_decode(dec, ctx, ln - 1, pli, bx + 1, by + 1, 0);
}

typedef void (*od_dec_func)(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int ln,
 int pli, int bx, int by, int has_ur);

const od_dec_func OD_DECODE_BLOCK[OD_NBSIZES + 2] = {
  od_single_band_decode,
  od_single_band_decode,
  od_single_band_decode,
  od_32x32_decode
};

#if !OD_DISABLE_HAAR_DC
static void od_decode_haar_dc(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int pli,
 int bx, int by, int l, int xdec, int ydec, od_coeff hgrad, od_coeff vgrad,
 int has_ur) {
  int od;
  int d;
  int w;
  int i;
  int dc_quant;
  od_coeff *c;
  c = ctx->d[pli];
  w = dec->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  od = OD_BLOCK_SIZE4x4(dec->state.bsize,
   dec->state.bstride, bx << l, by << l);
  d = OD_MAXI(od, xdec);
  OD_ASSERT(d <= l);
  if (dec->quantizer[pli] == 0) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, dec->quantizer[pli]*OD_DC_RES[pli] >> 4);
  }
  if (l == 3) {
    int nhsb;
    int quant;
    int l2;
    od_coeff sb_dc_pred;
    od_coeff sb_dc_curr;
    od_coeff *sb_dc_mem;
    nhsb = dec->state.nhsb;
    sb_dc_mem = dec->state.sb_dc_mem[pli];
    l2 = l - xdec + 2;
    if (by > 0 && bx > 0) {
      /* These coeffs were LS-optimized on subset 1. */
      if (has_ur) {
        sb_dc_pred = (22*sb_dc_mem[by*nhsb + bx - 1]
         - 9*sb_dc_mem[(by - 1)*nhsb + bx - 1]
         + 15*sb_dc_mem[(by - 1)*nhsb + bx]
         + 4*sb_dc_mem[(by - 1)*nhsb + bx + 1] + 16) >> 5;
      }
      else {
        sb_dc_pred = (23*sb_dc_mem[by*nhsb + bx - 1]
         - 10*sb_dc_mem[(by - 1)*nhsb + bx - 1]
         + 19*sb_dc_mem[(by - 1)*nhsb + bx] + 16) >> 5;
      }
    }
    else if (by > 0) sb_dc_pred = sb_dc_mem[(by - 1)*nhsb + bx];
    else if (bx > 0) sb_dc_pred = sb_dc_mem[by*nhsb + bx - 1];
    else sb_dc_pred = 0;
    quant = generic_decode(&dec->ec, &dec->adapt.model_dc[pli], -1,
     &dec->adapt.ex_sb_dc[pli], 2);
    if (quant) {
      if (od_ec_dec_bits(&dec->ec, 1)) quant = -quant;
    }
    sb_dc_curr = quant*dc_quant + sb_dc_pred;
    c[(by << l2)*w + (bx << l2)] = sb_dc_curr;
    sb_dc_mem[by*nhsb + bx] = sb_dc_curr;
    if (by > 0) vgrad = sb_dc_mem[(by - 1)*nhsb + bx] - sb_dc_curr;
    if (bx > 0) hgrad = sb_dc_mem[by*nhsb + bx - 1] - sb_dc_curr;
  }
  if (l > d) {
    od_coeff x[4];
    int l2;
    l--;
    bx <<= 1;
    by <<= 1;
    l2 = l - xdec + 2;
    x[0] = c[(by << l2)*w + (bx << l2)];
    x[1] = c[(by << l2)*w + ((bx + 1) << l2)];
    x[2] = c[((by + 1) << l2)*w + (bx << l2)];
    x[3] = c[((by + 1) << l2)*w + ((bx + 1) << l2)];
    for (i = 1; i < 4; i++) {
      int quant;
      quant = generic_decode(&dec->ec, &dec->adapt.model_dc[pli], -1,
       &dec->adapt.ex_dc[pli][l][i-1], 2);
      if (quant) {
        if (od_ec_dec_bits(&dec->ec, 1)) quant = -quant;
      }
      x[i] = quant*dc_quant;
    }
    /* Gives best results for subset1, more conservative than the
       theoretical /4 of a pure gradient. */
    x[1] += hgrad/5;
    x[2] += vgrad/5;
    hgrad = x[1];
    vgrad = x[2];
    OD_HAAR_KERNEL(x[0], x[1], x[2], x[3]);
    c[(by << l2)*w + (bx << l2)] = x[0];
    c[(by << l2)*w + ((bx + 1) << l2)] = x[1];
    c[((by + 1) << l2)*w + (bx << l2)] = x[2];
    c[((by + 1) << l2)*w + ((bx + 1) << l2)] = x[3];
    od_decode_haar_dc(dec, ctx, pli, bx + 0, by + 0, l, xdec, ydec, hgrad,
     vgrad, 0);
    od_decode_haar_dc(dec, ctx, pli, bx + 1, by + 0, l, xdec, ydec, hgrad,
     vgrad, 0);
    od_decode_haar_dc(dec, ctx, pli, bx + 0, by + 1, l, xdec, ydec, hgrad,
     vgrad, 0);
    od_decode_haar_dc(dec, ctx, pli, bx + 1, by + 1, l, xdec, ydec, hgrad,
     vgrad, 0);
  }
}
#endif

static void od_decode_block(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int pli,
 int bx, int by, int l, int xdec, int ydec, int has_ur) {
  int od;
  int d;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  od = OD_BLOCK_SIZE4x4(dec->state.bsize,
   dec->state.bstride, bx << l, by << l);
  d = OD_MAXI(od, xdec);
  OD_ASSERT(d <= l);
  if (d == l) {
    d -= xdec;
    /*Construct the luma predictors for chroma planes.*/
    if (ctx->l != NULL) {
      int w;
      int frame_width;
      OD_ASSERT(pli > 0);
      frame_width = dec->state.frame_width;
      w = frame_width >> xdec;
      od_resample_luma_coeffs(ctx->l + (by << (2 + d))*w + (bx << (2 + d)), w,
       ctx->d[0] + (by << (2 + l))*frame_width + (bx << (2 + l)),
       frame_width, xdec, ydec, d, od);
    }
    (*OD_DECODE_BLOCK[d])(dec, ctx, d, pli, bx, by, has_ur);
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_decode_block(dec, ctx, pli, bx + 0, by + 0, l, xdec, ydec, 1);
    od_decode_block(dec, ctx, pli, bx + 1, by + 0, l, xdec, ydec, has_ur);
    od_decode_block(dec, ctx, pli, bx + 0, by + 1, l, xdec, ydec, 1);
    od_decode_block(dec, ctx, pli, bx + 1, by + 1, l, xdec, ydec, 0);
  }
}

static void od_dec_mv_unpack(daala_dec_ctx *dec) {
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  od_img *img;
  int width;
  int height;
  int mv_res;
  od_mv_grid_pt *mvp;
  od_mv_grid_pt **grid;
  OD_ASSERT(dec->state.ref_imgi[OD_FRAME_PREV] >= 0);
  od_state_mvs_clear(&dec->state);
  nhmvbs = (dec->state.nhmbs + 1) << 2;
  nvmvbs = (dec->state.nvmbs + 1) << 2;
  img = dec->state.io_imgs + OD_FRAME_REC;
  mv_res = dec->state.mv_res = od_ec_dec_uint(&dec->ec, 3);
  width = (img->width + 32) << (3 - mv_res);
  height = (img->height + 32) << (3 - mv_res);
  grid = dec->state.mv_grid;
  /*Level 0.*/
  for (vy = 0; vy <= nvmvbs; vy += 4) {
    for (vx = 0; vx <= nhmvbs; vx += 4) {
      mvp = &grid[vy][vx];
      mvp->valid = 1;
      od_decode_mv(dec, mvp, vx, vy, 0, mv_res, width, height);
    }
  }
  /*Level 1.*/
  for (vy = 2; vy <= nvmvbs; vy += 4) {
    for (vx = 2; vx <= nhmvbs; vx += 4) {
      int p_invalid;
      p_invalid = od_mv_level1_prob(grid, vx, vy);
      mvp = &grid[vy][vx];
      mvp->valid = od_ec_decode_bool_q15(&dec->ec, p_invalid);
      if (mvp->valid) {
        od_decode_mv(dec, mvp, vx, vy, 1, mv_res, width, height);
      }
    }
  }
  /*Level 2.*/
  for (vy = 0; vy <= nvmvbs; vy += 2) {
    for (vx = 2*((vy & 3) == 0); vx <= nhmvbs; vx += 4) {
      mvp = &grid[vy][vx];
      if ((vy-2 < 0 || grid[vy-2][vx].valid)
       && (vx-2 < 0 || grid[vy][vx-2].valid)
       && (vy+2 > nvmvbs || grid[vy+2][vx].valid)
       && (vx+2 > nhmvbs || grid[vy][vx+2].valid)) {
        mvp->valid = od_ec_decode_bool_q15(&dec->ec, 13684);
        if (mvp->valid) {
          od_decode_mv(dec, mvp, vx, vy, 2, mv_res, width, height);
        }
      }
    }
  }
  /*Level 3.*/
  for (vy = 1; vy <= nvmvbs; vy += 2) {
    for (vx = 1; vx <= nhmvbs; vx += 2) {
      mvp = &grid[vy][vx];
      if (grid[vy-1][vx-1].valid && grid[vy-1][vx+1].valid
       && grid[vy+1][vx+1].valid && grid[vy+1][vx-1].valid) {
        mvp->valid = od_ec_decode_bool_q15(&dec->ec, 16384);
        if (mvp->valid) {
          od_decode_mv(dec, mvp, vx, vy, 3, mv_res, width, height);
        }
      }
    }
  }
  /*Level 4.*/
  for (vy = 2; vy <= nvmvbs - 2; vy += 1) {
    for (vx = 3 - (vy & 1); vx <= nhmvbs - 2; vx += 2) {
      mvp = &grid[vy][vx];
      if (grid[vy-1][vx].valid && grid[vy][vx-1].valid
       && grid[vy+1][vx].valid && grid[vy][vx+1].valid) {
        mvp->valid = od_ec_decode_bool_q15(&dec->ec, 16384);
        if (mvp->valid) {
          od_decode_mv(dec, mvp, vx, vy, 4, mv_res, width, height);
        }
      }
    }
  }
}

/*Decode the blocks sizes into their corresponding form in dec->state.bsize.
  See the comment for the `bsize` member of `od_state` for more information
   about the data layout and meaning.*/
static void od_decode_block_sizes(od_dec_ctx *dec) {
  int i;
  int j;
  int nhsb;
  int nvsb;
  od_state_init_border(&dec->state);
  nhsb = dec->state.nhsb;
  nvsb = dec->state.nvsb;
  if (OD_LIMIT_LOG_BSIZE_MIN != OD_LIMIT_LOG_BSIZE_MAX) {
    for (i = 0; i < nvsb; i++) {
      for (j = 0; j < nhsb; j++) {
        od_block_size_decode(&dec->ec, &dec->adapt,
         &dec->state.bsize[4*dec->state.bstride*i + 4*j], dec->state.bstride);
      }
    }
  } else {
    for (i = 0; i < nvsb*4; i++) {
      for (j = 0; j < nhsb*4; j++) {
        dec->state.bsize[dec->state.bstride*i + j] =
         OD_LIMIT_LOG_BSIZE_MIN - OD_LOG_BSIZE0;
      }
    }
  }
}

int daala_decode_packet_in(daala_dec_ctx *dec, od_img *img,
 const ogg_packet *op) {
  int nplanes;
  int pli;
  int frame_width;
  int frame_height;
  int nvsb;
  int nhsb;
  int refi;
  od_mb_dec_ctx mbctx;
  if (dec == NULL || img == NULL || op == NULL) return OD_EFAULT;
  if (dec->packet_state != OD_PACKET_DATA) return OD_EINVAL;
  if (op->e_o_s) dec->packet_state = OD_PACKET_DONE;
  od_ec_dec_init(&dec->ec, op->packet, op->bytes);
  /*Read the packet type bit.*/
  if (od_ec_decode_bool_q15(&dec->ec, 16384)) return OD_EBADPACKET;
  mbctx.is_keyframe = od_ec_decode_bool_q15(&dec->ec, 16384);
  /*Update the buffer state.*/
  if (dec->state.ref_imgi[OD_FRAME_SELF] >= 0) {
    dec->state.ref_imgi[OD_FRAME_PREV] =
     dec->state.ref_imgi[OD_FRAME_SELF];
    /*TODO: Update golden frame.*/
    if (dec->state.ref_imgi[OD_FRAME_GOLD] < 0) {
      dec->state.ref_imgi[OD_FRAME_GOLD] =
       dec->state.ref_imgi[OD_FRAME_SELF];
      /*TODO: Mark keyframe timebase.*/
    }
  }
  /*Select a free buffer to use for this reference frame.*/
  for (refi = 0; refi == dec->state.ref_imgi[OD_FRAME_GOLD]
   || refi == dec->state.ref_imgi[OD_FRAME_PREV]
   || refi == dec->state.ref_imgi[OD_FRAME_NEXT]; refi++);
  dec->state.ref_imgi[OD_FRAME_SELF] = refi;
  nhsb = dec->state.nhsb;
  nvsb = dec->state.nvsb;
  od_adapt_ctx_reset(&dec->adapt, mbctx.is_keyframe);
  od_decode_block_sizes(dec);
  if (!mbctx.is_keyframe) {
    od_dec_mv_unpack(dec);
    od_state_mc_predict(&dec->state, OD_FRAME_PREV);
  }
  frame_width = dec->state.frame_width;
  frame_height = dec->state.frame_height;
  {
    od_coeff *ctmp[OD_NPLANES_MAX];
    od_coeff *dtmp[OD_NPLANES_MAX];
    od_coeff *mctmp[OD_NPLANES_MAX];
    od_coeff *mdtmp[OD_NPLANES_MAX];
    od_coeff *ltmp[OD_NPLANES_MAX];
    od_coeff *lbuf[OD_NPLANES_MAX];
    int xdec;
    int ydec;
    int sby;
    int sbx;
    int h;
    int w;
    int y;
    int x;
    /*Initialize the data needed for each plane.*/
    nplanes = dec->state.info.nplanes;
    for (pli = 0; pli < (OD_DISABLE_CFL ? nplanes : 1); pli++) {
      xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      mbctx.modes[pli] = _ogg_calloc((frame_width >> (2 + xdec))
       *(frame_height >> (2 + ydec)), sizeof(*mbctx.modes[pli]));
    }
    if (mbctx.is_keyframe) {
      for (pli = 0; pli < (OD_DISABLE_CFL ? nplanes : 1); pli++) {
        xdec = dec->state.io_imgs[OD_FRAME_REC].planes[pli].xdec;
        ydec = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ydec;
        w = frame_width >> xdec;
        h = frame_height >> ydec;
        mbctx.tf[pli] = _ogg_calloc(w*h, sizeof(*mbctx.tf[pli]));
      }
    }
    /*Apply the prefilter to the motion-compensated reference.*/
    if (!mbctx.is_keyframe) {
      for (pli = 0; pli < nplanes; pli++) {
        xdec = dec->state.io_imgs[OD_FRAME_REC].planes[pli].xdec;
        ydec = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ydec;
        w = frame_width >> xdec;
        h = frame_height >> ydec;
        mctmp[pli] = _ogg_calloc(w*h, sizeof(*mctmp[pli]));
        mdtmp[pli] = _ogg_calloc(w*h, sizeof(*mdtmp[pli]));
        /*Collect the image data needed for this plane.*/
        {
          unsigned char *mdata;
          int ystride;
          int coeff_shift;
          coeff_shift = dec->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
          mdata = dec->state.io_imgs[OD_FRAME_REC].planes[pli].data;
          ystride = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ystride;
          for (y = 0; y < h; y++) {
            for (x = 0; x < w; x++) {
              mctmp[pli][y*w + x] = (mdata[ystride*y + x] - 128)
               << coeff_shift;
            }
          }
        }
        /*Apply the prefilter across the entire image.*/
        for (sby = 0; sby < nvsb; sby++) {
          for (sbx = 0; sbx < nhsb; sbx++) {
            od_apply_prefilter(mctmp[pli], w, sbx, sby, 3, dec->state.bsize,
             dec->state.bstride, xdec, ydec, (sbx > 0 ? OD_LEFT_EDGE : 0) |
             (sby < nvsb - 1 ? OD_BOTTOM_EDGE : 0));
          }
        }
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /* TODO: We shouldn't be encoding the full, linear quantizer range. */
      dec->quantizer[pli] = od_ec_dec_uint(&dec->ec, 512 << OD_COEFF_SHIFT);
      mbctx.run_pvq[pli] = dec->quantizer[pli] &&
       od_ec_decode_bool_q15(&dec->ec, 16384);
      ctmp[pli] = _ogg_calloc(w*h, sizeof(*ctmp[pli]));
      dtmp[pli] = _ogg_calloc(w*h, sizeof(*dtmp[pli]));
      /*We predict chroma planes from the luma plane.
        Since chroma can be subsampled, we cache subsampled versions of the
         luma plane in the frequency domain.
        We can share buffers with the same subsampling.*/
      if (pli > 0) {
        int plj;
        if (xdec || ydec) {
          for (plj = 1; plj < pli; plj++) {
            if (xdec == dec->state.io_imgs[OD_FRAME_INPUT].planes[plj].xdec
             && ydec == dec->state.io_imgs[OD_FRAME_INPUT].planes[plj].ydec) {
              ltmp[pli] = NULL;
              lbuf[pli] = ltmp[plj];
            }
          }
          if (plj >= pli) {
            lbuf[pli] = ltmp[pli] = _ogg_calloc(w*h, sizeof(*ltmp[pli]));
          }
        }
        else {
          ltmp[pli] = NULL;
          lbuf[pli] = ctmp[pli];
        }
      }
      else lbuf[pli] = ltmp[pli] = NULL;
    }
    for (sby = 0; sby < nvsb; sby++) {
      for (sbx = 0; sbx < nhsb; sbx++) {
        for (pli = 0; pli < nplanes; pli++) {
          mbctx.c = ctmp[pli];
          mbctx.d = dtmp;
          mbctx.mc = mctmp[pli];
          mbctx.md = mdtmp[pli];
          mbctx.l = lbuf[pli];
          xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
          ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
          mbctx.nk = mbctx.k_total = mbctx.sum_ex_total_q8 = 0;
          mbctx.ncount = mbctx.count_total_q8 = mbctx.count_ex_total_q8 = 0;
          if (!OD_DISABLE_HAAR_DC && mbctx.is_keyframe) {
            od_decode_haar_dc(dec, &mbctx, pli, sbx, sby, 3, xdec, ydec, 0, 0,
             sby > 0 && sbx < nhsb - 1);
          }
          od_decode_block(dec, &mbctx, pli, sbx, sby, 3, xdec, ydec,
           sby > 0 && sbx < nhsb - 1);
        }
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /*Apply the postfilter across the entire image.*/
      for (sby = 0; sby < nvsb; sby++) {
        for (sbx = 0; sbx < nhsb; sbx++) {
          od_apply_postfilter(ctmp[pli], w, sbx, sby, 3, dec->state.bsize,
           dec->state.bstride, xdec, ydec, (sby > 0 ? OD_TOP_EDGE : 0) |
           (sbx < nhsb - 1 ? OD_RIGHT_EDGE : 0));
        }
      }
      {
        unsigned char *data;
        int ystride;
        data = dec->state.io_imgs[OD_FRAME_REC].planes[pli].data;
        ystride = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ystride;
        for (y = 0; y < h; y++) {
          for (x = 0; x < w; x++) {
            int coeff_shift;
            coeff_shift = dec->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
            data[ystride*y + x] = OD_CLAMP255(((ctmp[pli][y*w + x]
             + (1 << coeff_shift >> 1)) >> coeff_shift) + 128);
          }
        }
      }
    }
    for (pli = nplanes; pli-- > 0;) {
      _ogg_free(ltmp[pli]);
      _ogg_free(dtmp[pli]);
      _ogg_free(ctmp[pli]);
      if (!mbctx.is_keyframe) {
        _ogg_free(mdtmp[pli]);
        _ogg_free(mctmp[pli]);
      }
    }
    for (pli = 0; pli < (OD_DISABLE_CFL ? nplanes : 1); pli++) {
      _ogg_free(mbctx.modes[pli]);
    }
    if (mbctx.is_keyframe) {
      for (pli = 0; pli < (OD_DISABLE_CFL ? nplanes : 1); pli++) {
        _ogg_free(mbctx.tf[pli]);
      }
    }
  }
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  /*Dump YUV*/
  od_state_dump_yuv(&dec->state, dec->state.io_imgs + OD_FRAME_REC, "out");
#endif
  od_state_upsample8(&dec->state,
   dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF],
   dec->state.io_imgs + OD_FRAME_REC);
  /*Return decoded frame.*/
  *img = dec->state.io_imgs[OD_FRAME_REC];
  img->width = dec->state.info.pic_width;
  img->height = dec->state.info.pic_height;
  dec->state.cur_time++;
  return 0;
}
