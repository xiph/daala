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
  generic_encoder *model;
  int pred[2];
  int ox;
  int oy;
  int id;
  int equal_mvs;
  equal_mvs = od_state_get_predictor(&dec->state, pred, vx, vy, level, mv_res);
  model = &dec->state.adapt.mv_model;
  id = od_decode_cdf_adapt(&dec->ec, dec->state.adapt.mv_small_cdf[equal_mvs],
   16, dec->state.adapt.mv_small_increment);
  oy = id >> 2;
  ox = id & 0x3;
  if (ox == 3) {
    ox += generic_decode(&dec->ec, model, width << (3 - mv_res),
     &dec->state.adapt.mv_ex[level], 6);
  }
  if (oy == 3) {
    oy += generic_decode(&dec->ec, model, height << (3 - mv_res),
     &dec->state.adapt.mv_ey[level], 6);
  }
  if (ox && od_ec_dec_bits(&dec->ec, 1)) ox = -ox;
  if (oy && od_ec_dec_bits(&dec->ec, 1)) oy = -oy;
  mvg->mv[0] = (pred[0] + ox) << mv_res;
  mvg->mv[1] = (pred[1] + oy) << mv_res;
}

struct od_mb_dec_ctx {
  od_coeff *c;
  od_coeff **d;
  od_coeff *md;
  od_coeff *mc;
  od_coeff *l;
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
  int ln, int pli, int bx, int by) {
  int n;
  int n2;
  int xdec;
  int w;
  int frame_width;
  od_coeff *md;
  od_coeff *l;
  int x;
  int y;
  OD_ASSERT(ln >= 0 && ln <= 3);
  n = 1 << (ln + 2);
  n2 = n*n;
  xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = dec->state.frame_width;
  w = frame_width >> xdec;
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  md = ctx->md;
  l = ctx->l;
  if (ctx->is_keyframe) {
    if (pli == 0 || OD_DISABLE_CFL) {
      OD_CLEAR(pred, n2);
    }
    else {
      int i;
      int j;
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          pred[i*n + j] = l[((by << 2) + i)*w + (bx << 2) + j];
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

static void od_block_lossless_decode(daala_dec_ctx *dec, int ln,
 od_coeff *pred, const od_coeff *predt, int pli) {
  int *adapt;
  int vk;
  int zzi;
  int n2;
  ogg_int32_t adapt_curr[OD_NSB_ADAPT_CTXS];
  adapt = dec->state.adapt.pvq_adapt;
  n2 = 1 << (2*ln + 4);
  vk = generic_decode(&dec->ec, &dec->state.adapt.model_g[pli], -1,
   &dec->state.adapt.ex_g[pli][ln], 0);
  laplace_decode_vector(&dec->ec, pred + 1, n2 - 1, vk, adapt_curr, adapt);
  for (zzi = 1; zzi < n2; zzi++) {
    pred[zzi] = pred[zzi] + predt[zzi];
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

static void od_block_decode(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int ln,
 int pli, int bx, int by) {
  int n;
  int xdec;
  int w;
  int frame_width;
  od_coeff *c;
  od_coeff *d;
  od_coeff *md;
  od_coeff *mc;
  od_coeff pred[OD_BSIZE_MAX*OD_BSIZE_MAX];
  od_coeff predt[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int lossless;
  int quant;
  int dc_quant;
  OD_ASSERT(ln >= 0 && ln <= 3);
  n = 1 << (ln + 2);
  lossless = (dec->quantizer[pli] == 0);
  bx <<= ln;
  by <<= ln;
  xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = dec->state.frame_width;
  w = frame_width >> xdec;
  c = ctx->c;
  d = ctx->d[pli];
  md = ctx->md;
  mc = ctx->mc;
  /*Apply forward transform to MC predictor.*/
  if (!ctx->is_keyframe) {
    (*dec->state.opt_vtbl.fdct_2d[ln])(md + (by << 2)*w + (bx << 2), w,
     mc + (by << 2)*w + (bx << 2), w);
  }
  od_decode_compute_pred(dec, ctx, pred, ln, pli, bx, by);
  if (ctx->is_keyframe && pli == 0) {
    od_hv_intra_pred(pred, d, w, bx, by, dec->state.bsize,
     dec->state.bstride, ln);
  }
  od_raster_to_coding_order(predt,  n, &pred[0], n, lossless);
  quant = OD_MAXI(1, dec->quantizer[pli]);
  if (lossless) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, quant*
     dec->state.pvq_qm_q4[pli][od_qm_get_index(ln, 0)] >> 4);
  }
  if (lossless) {
    od_block_lossless_decode(dec, ln, pred, predt, pli);
  }
  else {
    pvq_decode(dec, predt, pred, quant, pli, ln,
     OD_PVQ_BETA[pli][ln], OD_ROBUST_STREAM, ctx->is_keyframe);
  }
  if (OD_DISABLE_HAAR_DC || !ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && !lossless;
    if (!has_dc_skip || pred[0]) {
      pred[0] = has_dc_skip + generic_decode(&dec->ec,
       &dec->state.adapt.model_dc[pli], -1, &dec->state.adapt.ex_dc[pli][ln][0], 2);
      if (pred[0]) pred[0] *= od_ec_dec_bits(&dec->ec, 1) ? -1 : 1;
    }
    pred[0] = pred[0]*dc_quant + predt[0];
  }
  else {
    pred[0] = d[((by << 2))*w + ((bx << 2))];
  }
  od_coding_order_to_raster(&d[((by << 2))*w + (bx << 2)], w, pred, n,
   lossless);
  /*Apply the inverse transform.*/
  (*dec->state.opt_vtbl.idct_2d[ln])(c + (by << 2)*w + (bx << 2), w,
   d + (by << 2)*w + (bx << 2), w);
}

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
    if (dec->quantizer[pli] != 0 && d - xdec == 3) {
      dc_quant = OD_MAXI(1, dec->quantizer[pli]*12*OD_DC_RES[pli] >> 8);
    }
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
    quant = generic_decode(&dec->ec, &dec->state.adapt.model_dc[pli], -1,
     &dec->state.adapt.ex_sb_dc[pli], 2);
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
      quant = generic_decode(&dec->ec, &dec->state.adapt.model_dc[pli], -1,
       &dec->state.adapt.ex_dc[pli][l][i-1], 2);
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

static void od_decode_recursive(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int pli,
 int bx, int by, int l, int xdec, int ydec) {
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
    od_block_decode(dec, ctx, d, pli, bx, by);
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_decode_recursive(dec, ctx, pli, bx + 0, by + 0, l, xdec, ydec);
    od_decode_recursive(dec, ctx, pli, bx + 1, by + 0, l, xdec, ydec);
    od_decode_recursive(dec, ctx, pli, bx + 0, by + 1, l, xdec, ydec);
    od_decode_recursive(dec, ctx, pli, bx + 1, by + 1, l, xdec, ydec);
  }
}

static int od_dec_mv_in_frame(int vx, int vy, int nhmvbs, int nvmvbs) {
  return vx >= 2 && vy >= 2 && vx <= nhmvbs - 2 && vy <= nvmvbs - 2;
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
  mv_res = od_ec_dec_uint(&dec->ec, 3);
  od_state_set_mv_res(&dec->state, mv_res);
  width = (img->width + 32) << (3 - mv_res);
  height = (img->height + 32) << (3 - mv_res);
  grid = dec->state.mv_grid;
  /*Motion vectors outside the frame are always zero.*/
  /*Level 0.*/
  /*We don't modify the loop indices as in the encoder because we need to
    set all level 0 MVs valid.*/
  for (vy = 0; vy <= nvmvbs; vy += 4) {
    for (vx = 0; vx <= nhmvbs; vx += 4) {
      mvp = &grid[vy][vx];
      mvp->valid = 1;
      if (od_dec_mv_in_frame(vx, vy, nhmvbs, nvmvbs)) {
        od_decode_mv(dec, mvp, vx, vy, 0, mv_res, width, height);
      }
    }
  }
  /*Level 1.*/
  for (vy = 2; vy <= nvmvbs; vy += 4) {
    for (vx = 2; vx <= nhmvbs; vx += 4) {
      int p_invalid;
      p_invalid = od_mv_level1_probz(grid, vx, vy);
      mvp = &grid[vy][vx];
      if (p_invalid >= 16384) {
        mvp->valid = od_ec_decode_bool_q15(&dec->ec, p_invalid);
      }
      else {
        mvp->valid = !od_ec_decode_bool_q15(&dec->ec, 32768 - p_invalid);
      }
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
        int p_invalid;
        p_invalid = od_mv_level2_probz(grid, vx, vy);
        if (p_invalid >= 16834) {
          mvp->valid = od_ec_decode_bool_q15(&dec->ec, p_invalid);
        }
        else {
          mvp->valid = !od_ec_decode_bool_q15(&dec->ec, 32768 - p_invalid);
        }
        if (mvp->valid && od_dec_mv_in_frame(vx, vy, nhmvbs, nvmvbs)) {
          od_decode_mv(dec, mvp, vx, vy, 2, mv_res, width, height);
        }
      }
    }
  }
  /*Level 3.*/
  /*Level 3 motion vector flags are complicated on the edges. See the comments
    in encode.c for why this code is complicated.*/
  for (vy = 1; vy <= nvmvbs; vy += 2) {
    for (vx = 1; vx <= nhmvbs; vx += 2) {
      mvp = &grid[vy][vx];
      if (vy < 2 || vy > nvmvbs - 2) {
        if ((vx == 3 && grid[vy == 1 ? vy - 1 : vy + 1][vx - 1].valid)
         || (vx == nhmvbs - 3
         && grid[vy == 1 ? vy - 1 : vy + 1][vx + 1].valid)) {
          mvp->valid = 1;
          /*MV is outside frame and will be zero.*/
        }
        else if (vx > 3 && vx < nhmvbs - 3) {
          if (!(vx & 2) && grid[vy == 1 ? vy - 1 : vy + 1][vx + 1].valid) {
            /*0 = both valid, 1 = only this one, 2 = other one valid*/
            int s;
            s = od_ec_decode_cdf_q15(&dec->ec, OD_UNIFORM_CDF_Q15(3), 3);
            mvp->valid = s != 2;
            grid[vy][vx + 2].valid = s != 1;
            /*MV is outside frame and will be zero.*/
          }
        }
      }
      else if (vx < 2 || vx > nhmvbs - 2) {
        if ((vy == 3 && grid[vy - 1][vx == 1 ? vx - 1 : vx + 1].valid)
         || (vy == nvmvbs - 3
         && grid[vy + 1][vx == 1 ? vx - 1 : vx + 1].valid)) {
          mvp->valid = 1;
          /*MV is outside frame and will be zero.*/
        }
        else if (!(vy & 2) && grid[vy + 1][vx == 1 ? vx - 1 : vx + 1].valid) {
          int s;
          s = od_ec_decode_cdf_q15(&dec->ec, OD_UNIFORM_CDF_Q15(3), 3);
          mvp->valid = s != 2;
          grid[vy + 2][vx].valid = s != 1;
          /*MVs are valid but will be zero.*/
        }
      }
      else if (grid[vy - 1][vx - 1].valid && grid[vy - 1][vx + 1].valid
       && grid[vy + 1][vx + 1].valid && grid[vy + 1][vx - 1].valid) {
        int p_invalid;
        p_invalid = od_mv_level3_probz(grid, vx, vy);
        if (p_invalid >= 16384) {
          mvp->valid = od_ec_decode_bool_q15(&dec->ec, p_invalid);
        }
        else {
          mvp->valid = !od_ec_decode_bool_q15(&dec->ec, 32768 - p_invalid);
        }
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
        int p_invalid;
        p_invalid = od_mv_level4_probz(grid, vx, vy);
        if (p_invalid >= 16384) {
          mvp->valid = od_ec_decode_bool_q15(&dec->ec, p_invalid);
        }
        else {
          mvp->valid = !od_ec_decode_bool_q15(&dec->ec, 32768 - p_invalid);
        }
        if (mvp->valid && od_dec_mv_in_frame(vx, vy, nhmvbs, nvmvbs)) {
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
  if (OD_LIMIT_BSIZE_MIN != OD_LIMIT_BSIZE_MAX) {
    for (i = 0; i < nvsb; i++) {
      for (j = 0; j < nhsb; j++) {
        od_block_size_decode(&dec->ec, &dec->state.adapt,
         &dec->state.bsize[4*dec->state.bstride*i + 4*j], dec->state.bstride);
      }
    }
  }
  else {
    for (i = 0; i < nvsb*4; i++) {
      for (j = 0; j < nhsb*4; j++) {
        dec->state.bsize[dec->state.bstride*i + j] = OD_LIMIT_BSIZE_MIN;
      }
    }
  }
}

static void od_decode_residual(od_dec_ctx *dec, od_mb_dec_ctx *mbctx) {
  int nplanes;
  int pli;
  int xdec;
  int ydec;
  int sby;
  int sbx;
  int h;
  int w;
  int y;
  int x;
  int frame_width;
  int frame_height;
  int nvsb;
  int nhsb;
  od_state *state;
  state = &dec->state;
  /*Initialize the data needed for each plane.*/
  nplanes = state->info.nplanes;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  frame_width = state->frame_width;
  frame_height = state->frame_height;
  /*Apply the prefilter to the motion-compensated reference.*/
  if (!mbctx->is_keyframe) {
    for (pli = 0; pli < nplanes; pli++) {
      xdec = state->io_imgs[OD_FRAME_REC].planes[pli].xdec;
      ydec = state->io_imgs[OD_FRAME_REC].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /*Collect the image data needed for this plane.*/
      {
        unsigned char *mdata;
        int ystride;
        int coeff_shift;
        coeff_shift = dec->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
        mdata = state->io_imgs[OD_FRAME_REC].planes[pli].data;
        ystride = state->io_imgs[OD_FRAME_REC].planes[pli].ystride;
        for (y = 0; y < h; y++) {
          for (x = 0; x < w; x++) {
            state->mctmp[pli][y*w + x] = (mdata[ystride*y + x] - 128)
             << coeff_shift;
          }
        }
      }
      /*Apply the prefilter across the entire image.*/
      for (sby = 0; sby < nvsb; sby++) {
        for (sbx = 0; sbx < nhsb; sbx++) {
          od_apply_prefilter(state->mctmp[pli], w, sbx, sby, 3,
           state->bsize, state->bstride, xdec, ydec,
           (sbx > 0 ? OD_LEFT_EDGE : 0) |
           (sby < nvsb - 1 ? OD_BOTTOM_EDGE : 0));
        }
      }
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    /* TODO: We shouldn't be encoding the full, linear quantizer range. */
    dec->quantizer[pli] = od_ec_dec_uint(&dec->ec, 512 << OD_COEFF_SHIFT);
  }
  for (sby = 0; sby < nvsb; sby++) {
    for (sbx = 0; sbx < nhsb; sbx++) {
      for (pli = 0; pli < nplanes; pli++) {
        mbctx->c = state->ctmp[pli];
        mbctx->d = state->dtmp;
        mbctx->mc = state->mctmp[pli];
        mbctx->md = state->mdtmp[pli];
        mbctx->l = state->lbuf[pli];
        xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
        ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
        mbctx->nk = mbctx->k_total = mbctx->sum_ex_total_q8 = 0;
        mbctx->ncount = mbctx->count_total_q8 = mbctx->count_ex_total_q8 = 0;
        if (!OD_DISABLE_HAAR_DC && mbctx->is_keyframe) {
          od_decode_haar_dc(dec, mbctx, pli, sbx, sby, 3, xdec, ydec, 0, 0,
           sby > 0 && sbx < nhsb - 1);
        }
        od_decode_recursive(dec, mbctx, pli, sbx, sby, 3, xdec, ydec);
      }
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    /*Apply the postfilter across the entire image.*/
    for (sby = 0; sby < nvsb; sby++) {
      for (sbx = 0; sbx < nhsb; sbx++) {
        od_apply_postfilter(state->ctmp[pli], w, sbx, sby, 3, state->bsize,
         state->bstride, xdec, ydec, (sby > 0 ? OD_TOP_EDGE : 0) |
         (sbx < nhsb - 1 ? OD_RIGHT_EDGE : 0));
      }
    }
    {
      unsigned char *data;
      od_coeff *ctmp;
      int ystride;
      int coeff_shift;
      coeff_shift = dec->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
      data = state->io_imgs[OD_FRAME_REC].planes[pli].data;
      ctmp = state->ctmp[pli];
      ystride = state->io_imgs[OD_FRAME_REC].planes[pli].ystride;
      for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
          data[ystride*y + x] = OD_CLAMP255(((ctmp[y*w + x]
           + (1 << coeff_shift >> 1)) >> coeff_shift) + 128);
        }
      }
    }
  }
}

int daala_decode_packet_in(daala_dec_ctx *dec, od_img *img,
 const ogg_packet *op) {
  int refi;
  od_mb_dec_ctx mbctx;
  if (dec == NULL || img == NULL || op == NULL) return OD_EFAULT;
  if (dec->packet_state != OD_PACKET_DATA) return OD_EINVAL;
  if (op->e_o_s) dec->packet_state = OD_PACKET_DONE;
  od_ec_dec_init(&dec->ec, op->packet, op->bytes);
  /*Read the packet type bit.*/
  if (od_ec_decode_bool_q15(&dec->ec, 16384)) return OD_EBADPACKET;
  mbctx.is_keyframe = od_ec_decode_bool_q15(&dec->ec, 16384);
  if (mbctx.is_keyframe) {
    int nplanes;
    int pli;
    nplanes = dec->state.info.nplanes;
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      for (i = 0; i < OD_QM_SIZE; i++) {
        dec->state.pvq_qm_q4[pli][i] = od_ec_dec_bits(&dec->ec, 8);
      }
    }
  }
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
  od_adapt_ctx_reset(&dec->state.adapt, mbctx.is_keyframe);
  od_decode_block_sizes(dec);
  if (!mbctx.is_keyframe) {
    od_dec_mv_unpack(dec);
    od_state_mc_predict(&dec->state, OD_FRAME_PREV);
  }
  od_decode_residual(dec, &mbctx);
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
