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
#include "tf.h"
#include "state.h"
#include "quantizer.h"

static int od_dec_init(od_dec_ctx *dec, const daala_info *info,
 const daala_setup_info *setup) {
  int ret;
  (void)setup;
  ret = od_state_init(&dec->state, info);
  if (ret < 0) return ret;
  dec->packet_state = OD_PACKET_DATA;
  dec->user_bsize = NULL;
  dec->user_flags = NULL;
  dec->user_mv_grid = NULL;
  dec->user_mc_img = NULL;
  return 0;
}

static void od_dec_clear(od_dec_ctx *dec) {
  od_state_clear(&dec->state);
}

daala_dec_ctx *daala_decode_alloc(const daala_info *info,
 const daala_setup_info *setup) {
  od_dec_ctx *dec;
  if (info == NULL) return NULL;
  dec = (od_dec_ctx *)malloc(sizeof(*dec));
  if (od_dec_init(dec, info, setup) < 0) {
    free(dec);
    return NULL;
  }
  return dec;
}

void daala_decode_free(daala_dec_ctx *dec) {
  if (dec != NULL) {
    od_dec_clear(dec);
    free(dec);
  }
}

int daala_decode_ctl(daala_dec_ctx *dec, int req, void *buf, size_t buf_sz) {
  (void)dec;
  (void)buf;
  (void)buf_sz;
  switch (req) {
    case OD_DECCTL_SET_BSIZE_BUFFER : {
      if(dec == NULL || buf == NULL) return OD_EFAULT;
      /*Check that buf is large enough to hold the block sizes for a frame.*/
      if (buf_sz != sizeof(unsigned char)*dec->state.nvsb*dec->state.nhsb*16) {
        return OD_EINVAL;
      }
      dec->user_bsize = (unsigned char *)buf;
      dec->user_bstride = dec->state.nhsb*4;
      return 0;
    }
    case OD_DECCTL_SET_FLAGS_BUFFER : {
      if (dec == NULL || buf == NULL) return OD_EFAULT;
      /*Check that buf is large enough to hold the band flags for a frame.*/
      if (buf_sz != sizeof(unsigned int)*dec->state.nvsb*8*dec->state.nhsb*8) {
        return OD_EINVAL;
      }
      dec->user_flags = (unsigned int *)buf;
      dec->user_fstride = dec->state.nhsb*8;
      return 0;
    }
    case OD_DECCTL_SET_MV_BUFFER : {
      if (dec== NULL || buf == NULL) return OD_EFAULT;
      if (buf_sz !=
       sizeof(od_mv_grid_pt)*(dec->state.nhmvbs + 1)*(dec->state.nvmvbs + 1)) {
        return OD_EINVAL;
      }
      dec->user_mv_grid = buf;
      return 0;
    }
    case OD_DECCTL_SET_MC_IMG : {
      if (dec == NULL || buf == NULL) return OD_EFAULT;
      if (buf_sz != sizeof(od_img)) return OD_EINVAL;
      dec->user_mc_img = buf;
      return 0;
    }
    default: return OD_EIMPL;
  }
}

static void od_dec_blank_img(od_img *img) {
  int pli;
  int frame_buf_width;
  int frame_buf_height;
  int plane_buf_width;
  int plane_buf_height;
  frame_buf_width = img->width + (OD_UMV_PADDING << 1);
  frame_buf_height = img->height + (OD_UMV_PADDING << 1);
  for (pli = 0; pli < img->nplanes; pli++) {
    plane_buf_width = (frame_buf_width << 1) >> img->planes[pli].xdec;
    plane_buf_height = (frame_buf_height << 1) >> img->planes[pli].ydec;
    memset(img->planes[pli].data, 128, plane_buf_width*plane_buf_height);
  }
}

/*We're decoding an INTER frame, but have no initialized reference
   buffers (i.e., decoding did not start on a key frame).
  We initialize them to a solid gray here.*/
static void od_dec_init_dummy_frame(daala_dec_ctx *dec) {
  dec->state.ref_imgi[OD_FRAME_GOLD] =
   dec->state.ref_imgi[OD_FRAME_PREV] =
   dec->state.ref_imgi[OD_FRAME_SELF] = 0;
  od_dec_blank_img(dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF]);
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
  int use_activity_masking;
  int qm;
  int use_haar_wavelet;
};
typedef struct od_mb_dec_ctx od_mb_dec_ctx;

static void od_decode_compute_pred(daala_dec_ctx *dec, od_mb_dec_ctx *ctx,
 od_coeff *pred, int bs, int pli, int bx, int by) {
  int n;
  int xdec;
  int w;
  int bo;
  int y;
  int x;
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  n = 1 << bs + OD_LOG_BSIZE0;
  xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  w = dec->state.frame_width >> xdec;
  bo = (by << OD_LOG_BSIZE0)*w + (bx << OD_LOG_BSIZE0);
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  if (ctx->is_keyframe) {
    if (pli == 0 || OD_DISABLE_CFL || ctx->use_haar_wavelet) {
      OD_CLEAR(pred, n*n);
    }
    else {
      od_coeff *l;
      l = ctx->l;
      for (y = 0; y < n; y++) {
        for (x = 0; x < n; x++) {
          pred[n*y + x] = l[bo + y*w + x];
        }
      }
    }
  }
  else {
    od_coeff *md;
    md = ctx->md;
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        pred[n*y + x] = md[bo + y*w + x];
      }
    }
  }
}

static int od_ec_dec_unary(od_ec_dec *ec) {
  int ret;
  ret = 0;
  while (od_ec_dec_bits(ec, 1) == 0) ret++;
  return ret;
}

static int od_decode_coeff_split(daala_dec_ctx *dec, int sum, int ctx) {
  int shift;
  int a;
  a = 0;
  if (sum == 0) return 0;
  shift = OD_MAXI(0, OD_ILOG(sum) - 4);
  if (shift) {
    a = od_ec_dec_bits(&dec->ec, shift);
  }
  a += od_decode_cdf_adapt(&dec->ec, dec->state.adapt.haar_coeff_cdf[15*ctx
   + (sum >> shift) - 1], (sum >> shift) + 1,
   dec->state.adapt.haar_coeff_increment) << shift;
  if (a > sum) {
    a = sum;
    dec->ec.error = 1;
  }
  return a;
}

static int od_decode_tree_split(daala_dec_ctx *dec, int sum, int ctx) {
  int shift;
  int a;
  a = 0;
  if (sum == 0) return 0;
  shift = OD_MAXI(0, OD_ILOG(sum) - 4);
  if (shift) {
    a = od_ec_dec_bits(&dec->ec, shift);
  }
  a += od_decode_cdf_adapt(&dec->ec, dec->state.adapt.haar_split_cdf[15*(2*ctx
   + OD_MINI(shift, 1)) + (sum >> shift) - 1], (sum >> shift) + 1,
   dec->state.adapt.haar_split_increment) << shift;
  if (a > sum) {
    a = sum;
    dec->ec.error = 1;
  }
  return a;
}

static void od_decode_sum_tree(daala_dec_ctx *dec, od_coeff *c, int ln,
 od_coeff tree_sum, int x, int y, int dir, int pli) {
  int n;
  int coeff_mag;
  od_coeff children_sum;
  od_coeff children[2][2];
  n = 1 << ln;
  if (tree_sum == 0) return;
  coeff_mag = od_decode_coeff_split(dec, tree_sum, dir
   + 3*(OD_ILOG(OD_MAXI(x,y)) - 1));
  c[y*n + x] = coeff_mag;
  children_sum = tree_sum - coeff_mag;
  /* Decode sum of each four children relative to tree. */
  if (children_sum) {
    int sum1;
    if (dir == 0) {
      sum1 = od_decode_tree_split(dec, children_sum, 0);
      children[0][0] = od_decode_tree_split(dec, sum1, 2);
      children[0][1] = sum1 - children[0][0];
      children[1][0] = od_decode_tree_split(dec, children_sum - sum1, 2);
      children[1][1] = children_sum - sum1 - children[1][0];
    }
    else {
      sum1 = od_decode_tree_split(dec, children_sum, 1);
      children[0][0] = od_decode_tree_split(dec, sum1, 2);
      children[1][0] = sum1 - children[0][0];
      children[0][1] = od_decode_tree_split(dec, children_sum - sum1, 2);
      children[1][1] = children_sum - sum1 - children[0][1];
    }
  }
  else {
    children[0][0] = children[0][1] = children[1][0] = children[1][1] = 0;
  }
  if (4*x < n && 4*y < n) {
    /* Recursive calls. */
    od_decode_sum_tree(dec, c, ln, children[0][0], 2*x, 2*y, dir, pli);
    od_decode_sum_tree(dec, c, ln, children[0][1], 2*x + 1, 2*y, dir, pli);
    od_decode_sum_tree(dec, c, ln, children[1][0], 2*x, 2*y + 1, dir, pli);
    od_decode_sum_tree(dec, c, ln, children[1][1], 2*x + 1, 2*y + 1, dir, pli);
  }
  else {
    c[2*y*n + 2*x] = children[0][0];
    c[2*y*n + 2*x + 1] = children[0][1];
    c[(2*y + 1)*n + 2*x] = children[1][0];
    c[(2*y + 1)*n + 2*x + 1] = children[1][1];
  }
}

static void od_wavelet_unquantize(daala_dec_ctx *dec, int ln, od_coeff *pred,
 const od_coeff *predt, int quant, int pli) {
  int n;
  int i;
  int dir;
  od_coeff tree_sum[2][2];
  n = 1 << ln;
  for (i = 0; i < n; i++) {
    int j;
    for (j = 0; j < n; j++) if (i + j) pred[i*n + j] = 0;
  }
  {
    int bits;
    bits = od_decode_cdf_adapt(&dec->ec, dec->state.adapt.haar_bits_cdf[pli],
     16, dec->state.adapt.haar_bits_increment);
    if (bits == 15) bits += od_ec_dec_unary(&dec->ec);
    /* Theoretical maximum sum is around 2^7 * 2^OD_COEFF_SHIFT * 32x32,
       so 2^21, but let's play safe. */
    if (bits > 24) {
      /* This can't happen, must be a bit-stream desync/corruption. */
      dec->ec.error = 1;
      return;
    }
    else if (bits > 1) {
      tree_sum[0][0] = (1 << (bits - 1)) | od_ec_dec_bits(&dec->ec, bits - 1);
    }
    else tree_sum[0][0] = bits;
    /* Handle diagonal first to make H/V symmetric. */
    tree_sum[1][1] = od_decode_tree_split(dec, tree_sum[0][0], 3);
    tree_sum[0][1] = od_decode_tree_split(dec, tree_sum[0][0] - tree_sum[1][1],
     4);
    tree_sum[1][0] = tree_sum[0][0] - tree_sum[1][1] - tree_sum[0][1];
  }
  od_decode_sum_tree(dec, pred, ln, tree_sum[0][1], 1, 0, 0, pli);
  od_decode_sum_tree(dec, pred, ln, tree_sum[1][0], 0, 1, 1, pli);
  od_decode_sum_tree(dec, pred, ln, tree_sum[1][1], 1, 1, 2, pli);
  for (i = 0; i < n; i++) {
    int j;
    for (j = (i == 0); j < n; j++) {
      int sign;
      od_coeff in;
      in = pred[i*n + j];
      if (in) {
        sign = od_ec_dec_bits(&dec->ec, 1);
        if (sign) in = -in;
      }
      pred[i*n + j] = in;
    }
  }
  for (dir = 0; dir < 3; dir++) {
    int level;
    for (level = 0; level < ln; level++) {
      int bo;
      int q;
      bo = (((dir + 1) >> 1) << level)*n + (((dir + 1) & 1) << level);
      if (quant == 0) q = 1;
      else q = quant*OD_HAAR_QM[dir == 2][level] >> 4;
      for (i = 0; i < 1 << level; i++) {
        int j;
        for (j = 0; j < 1 << level; j++)
          pred[bo + i*n + j] = q*pred[bo + i*n + j] + predt[bo + i*n + j];
      }
    }
  }
}

static void od_block_decode(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int bs,
 int pli, int bx, int by, int skip) {
  int n;
  int xdec;
  int w;
  int bo;
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
  int use_masking;
  const int *qm;
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  n = 1 << (bs + 2);
  lossless = (dec->quantizer[pli] == 0);
  use_masking = ctx->use_activity_masking;
  qm = ctx->qm == OD_HVS_QM ? OD_QM8_Q4_QM_HVS : OD_QM8_Q4_QM_FLAT;
  bx <<= bs;
  by <<= bs;
  xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = dec->state.frame_width;
  w = frame_width >> xdec;
  bo = (by << 2)*w + (bx << 2);
  c = ctx->c;
  d = ctx->d[pli];
  md = ctx->md;
  mc = ctx->mc;
  /*Apply forward transform to MC predictor.*/
  if (!ctx->is_keyframe) {
    if (ctx->use_haar_wavelet) {
      od_haar(md + bo, w, mc + bo, w, bs + 2);
    }
    else {
      (*dec->state.opt_vtbl.fdct_2d[bs])(md + bo, w, mc + bo, w);
      od_apply_qm(md + bo, w, md + bo, w, bs, xdec, 0, qm);
    }
  }
  od_decode_compute_pred(dec, ctx, pred, bs, pli, bx, by);
  if (ctx->is_keyframe && pli == 0 && !ctx->use_haar_wavelet) {
    od_hv_intra_pred(pred, d, w, bx, by, dec->state.bsize,
     dec->state.bstride, bs);
  }
  if (ctx->use_haar_wavelet) {
    int i;
    int j;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        predt[i*n + j] = pred[i*n + j];
      }
    }
  }
  else {
    od_raster_to_coding_order(predt,  n, &pred[0], n);
  }
  quant = OD_MAXI(1, dec->quantizer[pli]);
  if (lossless) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, quant*
     dec->state.pvq_qm_q4[pli][od_qm_get_index(bs, 0)] >> 4);
  }
  if (ctx->use_haar_wavelet) {
    od_wavelet_unquantize(dec, bs + 2, pred, predt, dec->quantizer[pli], pli);
  }
  else {
    unsigned int flags;
    od_pvq_decode(dec, predt, pred, quant, pli, bs,
     OD_PVQ_BETA[use_masking][pli][bs], OD_ROBUST_STREAM, ctx->is_keyframe,
     &flags, skip);
    if (pli == 0 && dec->user_flags != NULL) {
      dec->user_flags[by*dec->user_fstride + bx] = flags;
    }
  }
  if (!ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && !ctx->use_haar_wavelet;
    if (!has_dc_skip || pred[0]) {
      pred[0] = has_dc_skip + generic_decode(&dec->ec,
       &dec->state.adapt.model_dc[pli], -1, &dec->state.adapt.ex_dc[pli][bs][0], 2);
      if (pred[0]) pred[0] *= od_ec_dec_bits(&dec->ec, 1) ? -1 : 1;
    }
    pred[0] = pred[0]*dc_quant + predt[0];
  }
  else {
    pred[0] = d[bo];
  }
  if (ctx->use_haar_wavelet) {
    int i;
    int j;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        d[bo + i*w + j] = pred[i*n + j];
      }
    }
  }
  else {
    od_coding_order_to_raster(&d[bo], w, pred, n);
  }
  if (ctx->use_haar_wavelet) {
    od_haar_inv(c + bo, w, d + bo, w, bs + 2);
  }
  else {
    od_apply_qm(d + bo, w, d + bo, w, bs, xdec, 1, qm);
    /*Apply the inverse transform.*/
    (*dec->state.opt_vtbl.idct_2d[bs])(c + bo, w, d + bo, w);
  }
}

#if !OD_DISABLE_HAAR_DC
static void od_decode_haar_dc_sb(daala_dec_ctx *dec, od_mb_dec_ctx *ctx,
 int pli, int bx, int by, int xdec, int ydec, int has_ur, od_coeff *ohgrad,
 od_coeff *ovgrad) {
  int w;
  int dc_quant;
  od_coeff *c;
  int nhsb;
  int quant;
  int ln;
  od_coeff sb_dc_pred;
  od_coeff sb_dc_curr;
  od_coeff *sb_dc_mem;
  c = ctx->d[pli];
  w = dec->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  if (dec->quantizer[pli] == 0) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, dec->quantizer[pli]*OD_DC_RES[pli] >> 4);
  }
  nhsb = dec->state.nhsb;
  sb_dc_mem = dec->state.sb_dc_mem[pli];
  ln = OD_LOG_BSIZE_MAX - xdec;
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
  c[(by << ln)*w + (bx << ln)] = sb_dc_curr;
  sb_dc_mem[by*nhsb + bx] = sb_dc_curr;
  if (by > 0) *ovgrad = sb_dc_mem[(by - 1)*nhsb + bx] - sb_dc_curr;
  if (bx > 0) *ohgrad = sb_dc_mem[by*nhsb + bx - 1] - sb_dc_curr;
}
#endif

static void od_decode_haar_dc_level(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int pli,
  int bx, int by, int bsi, int xdec, od_coeff *hgrad, od_coeff *vgrad) {
  int i;
  od_coeff x[4];
  int ln;
  int ac_quant[2];
  int dc_quant;
  int w;
  w = dec->state.frame_width >> xdec;
  if (dec->quantizer[pli] == 0) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, dec->quantizer[pli]*OD_DC_RES[pli] >> 4);
  }
  if (dec->quantizer[pli] == 0) ac_quant[0] = ac_quant[1] = 1;
  else {
    ac_quant[0] = dc_quant*OD_DC_QM[xdec][bsi - xdec][0] >> 4;
    ac_quant[1] = dc_quant*OD_DC_QM[xdec][bsi - xdec][1] >> 4;
  }
  ln = bsi - xdec + 2;
  x[0] = ctx->d[pli][(by << ln)*w + (bx << ln)];
  x[1] = ctx->d[pli][(by << ln)*w + ((bx + 1) << ln)];
  x[2] = ctx->d[pli][((by + 1) << ln)*w + (bx << ln)];
  x[3] = ctx->d[pli][((by + 1) << ln)*w + ((bx + 1) << ln)];
  for (i = 1; i < 4; i++) {
    int quant;
    quant = generic_decode(&dec->ec, &dec->state.adapt.model_dc[pli], -1,
     &dec->state.adapt.ex_dc[pli][bsi][i-1], 2);
    if (quant) {
      if (od_ec_dec_bits(&dec->ec, 1)) quant = -quant;
    }
    x[i] = quant*ac_quant[i == 3];
  }
  /* Gives best results for subset1, more conservative than the
     theoretical /4 of a pure gradient. */
  x[1] += *hgrad/5;
  x[2] += *vgrad/5;
  *hgrad = x[1];
  *vgrad = x[2];
  OD_HAAR_KERNEL(x[0], x[1], x[2], x[3]);
  ctx->d[pli][(by << ln)*w + (bx << ln)] = x[0];
  ctx->d[pli][(by << ln)*w + ((bx + 1) << ln)] = x[1];
  ctx->d[pli][((by + 1) << ln)*w + (bx << ln)] = x[2];
  ctx->d[pli][((by + 1) << ln)*w + ((bx + 1) << ln)] = x[3];
}

static void od_decode_recursive(daala_dec_ctx *dec, od_mb_dec_ctx *ctx, int pli,
 int bx, int by, int bsi, int xdec, int ydec, od_coeff hgrad, od_coeff vgrad) {
  int obs;
  int bs;
  int w;
  int skip;
  int frame_width;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  obs = OD_BLOCK_SIZE4x4(dec->state.bsize,
   dec->state.bstride, bx << bsi, by << bsi);
  frame_width = dec->state.frame_width;
  w = frame_width >> xdec;
  skip = 0;
  /* Read the luma skip symbol. A value of 4 means "split the block", while < 4
     means that we code the block. In the latter case, we need to forward
     the skip value to the PVQ decoder. */
  if (ctx->use_haar_wavelet) obs = bsi;
  else if (pli == 0) {
    skip = od_decode_cdf_adapt(&dec->ec,
     dec->state.adapt.skip_cdf[2*bsi + (pli != 0)], 5,
     dec->state.adapt.skip_increment);
    if (skip < 4) obs = bsi;
    else obs = -1;
  }
  bs = OD_MAXI(obs, xdec);
  OD_ASSERT(bs <= bsi);
  if (bs == bsi) {
    bs -= xdec;
    if (pli == 0) {
      int i;
      int j;
      int n4;
      n4 = 1 << bsi;
      /* Save the block size decision so that chroma can reuse it. */
      for (i = 0; i < n4; i++) {
        for (j = 0; j < n4; j++) {
          OD_BLOCK_SIZE4x4(dec->state.bsize, dec->state.bstride, (bx << bsi) + i,
           (by << bsi) + j) = bsi;
        }
      }
    }
    /*Construct the luma predictors for chroma planes.*/
    if (ctx->l != NULL) {
      OD_ASSERT(pli > 0);
      od_resample_luma_coeffs(ctx->l + (by << (2 + bs))*w + (bx << (2 + bs)), w,
       ctx->d[0] + (by << (2 + bsi))*frame_width + (bx << (2 + bsi)),
       frame_width, xdec, ydec, bs, obs);
    }
    if (pli > 0 && !ctx->use_haar_wavelet) {
      /* Decode the skip for chroma. */
      skip = od_decode_cdf_adapt(&dec->ec,
       dec->state.adapt.skip_cdf[2*bsi + (pli != 0)], 5,
       dec->state.adapt.skip_increment);
    }
    od_block_decode(dec, ctx, bs, pli, bx, by, skip);
  }
  else {
    int f;
    int bo;
    bs = bsi - xdec;
    f = OD_FILT_SIZE(bs - 1, xdec);
    bo = (by << (OD_LOG_BSIZE0 + bs))*w + (bx << (OD_LOG_BSIZE0 + bs));
    if (!ctx->is_keyframe) {
      od_prefilter_split(ctx->mc + bo, w, bs, f);
    }
    bsi--;
    bx <<= 1;
    by <<= 1;
    if (ctx->is_keyframe) {
      od_decode_haar_dc_level(dec, ctx, pli, bx, by, bsi, xdec, &hgrad, &vgrad);
    }
    od_decode_recursive(dec, ctx, pli, bx + 0, by + 0, bsi, xdec, ydec, hgrad,
     vgrad);
    od_decode_recursive(dec, ctx, pli, bx + 1, by + 0, bsi, xdec, ydec, hgrad,
     vgrad);
    od_decode_recursive(dec, ctx, pli, bx + 0, by + 1, bsi, xdec, ydec, hgrad,
     vgrad);
    od_decode_recursive(dec, ctx, pli, bx + 1, by + 1, bsi, xdec, ydec, hgrad,
     vgrad);
    bs = bsi - xdec;
    bo = (by << (OD_LOG_BSIZE0 + bs))*w + (bx << (OD_LOG_BSIZE0 + bs));
    od_postfilter_split(ctx->c + bo, w, bs + 1, f);
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
  int log_mvb_sz;
  int level;
  od_mv_grid_pt *mvp;
  od_mv_grid_pt **grid;
  uint16_t *cdf;
  OD_ASSERT(dec->state.ref_imgi[OD_FRAME_PREV] >= 0);
  od_state_mvs_clear(&dec->state);
  nhmvbs = dec->state.nhmvbs;
  nvmvbs = dec->state.nvmvbs;
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
  for (vy = 0; vy <= nvmvbs; vy += OD_MVB_DELTA0) {
    for (vx = 0; vx <= nhmvbs; vx += OD_MVB_DELTA0) {
      mvp = grid[vy] + vx;
      mvp->valid = 1;
      od_decode_mv(dec, mvp, vx, vy, 0, mv_res, width, height);
    }
  }
  for (log_mvb_sz = OD_LOG_MVB_DELTA0, level = 1; log_mvb_sz-- > 0; level++) {
    int mvb_sz;
    mvb_sz = 1 << log_mvb_sz;
    /*Odd levels.*/
    for (vy = mvb_sz; vy <= nvmvbs; vy += 2*mvb_sz) {
      for (vx = mvb_sz; vx <= nhmvbs; vx += 2*mvb_sz) {
        if (grid[vy - mvb_sz][vx - mvb_sz].valid
         && grid[vy - mvb_sz][vx + mvb_sz].valid
         && grid[vy + mvb_sz][vx + mvb_sz].valid
         && grid[vy + mvb_sz][vx - mvb_sz].valid) {
          cdf = od_mv_split_flag_cdf(&dec->state, vx, vy, level);
          mvp = grid[vy] + vx;
          mvp->valid = od_decode_cdf_adapt(&dec->ec,
           cdf, 2, dec->state.adapt.split_flag_increment);
          if (mvp->valid) {
            od_decode_mv(dec, mvp, vx, vy, level, mv_res, width, height);
          }
        }
      }
    }
    level++;
    /*Even Levels.*/
    for (vy = 0; vy <= nvmvbs; vy += mvb_sz) {
      for (vx = mvb_sz*!(vy & mvb_sz); vx <= nhmvbs; vx += 2*mvb_sz) {
        if ((vy - mvb_sz < 0 || grid[vy - mvb_sz][vx].valid)
         && (vx - mvb_sz < 0 || grid[vy][vx - mvb_sz].valid)
         && (vy + mvb_sz > nvmvbs || grid[vy + mvb_sz][vx].valid)
         && (vx + mvb_sz > nhmvbs || grid[vy][vx + mvb_sz].valid)) {
          cdf = od_mv_split_flag_cdf(&dec->state, vx, vy, level);
          mvp = grid[vy] + vx;
          mvp->valid = od_decode_cdf_adapt(&dec->ec,
           cdf, 2, dec->state.adapt.split_flag_increment);
          if (mvp->valid) {
            od_decode_mv(dec, mvp, vx, vy, level, mv_res, width, height);
          }
        }
      }
    }
  }
  if (dec->user_mv_grid != NULL) {
    for (vy = 0; vy <= nvmvbs; vy++) {
      for (vx = 0; vx <= nhmvbs; vx++) {
        memcpy(&dec->user_mv_grid[vy*(nhmvbs + 1) + vx], &grid[vy][vx],
         sizeof(od_mv_grid_pt));
      }
    }
  }
}

static void od_decode_coefficients(od_dec_ctx *dec, od_mb_dec_ctx *mbctx) {
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
      if (!mbctx->is_keyframe && !mbctx->use_haar_wavelet) {
        od_apply_prefilter_frame_sbs(state->mctmp[pli], w, nhsb, nvsb, xdec,
         ydec);
      }
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    dec->quantizer[pli] =
     od_codedquantizer_to_quantizer(od_ec_dec_uint(&dec->ec,
     OD_N_CODED_QUANTIZERS));
  }
  for (sby = 0; sby < nvsb; sby++) {
    for (sbx = 0; sbx < nhsb; sbx++) {
      for (pli = 0; pli < nplanes; pli++) {
        od_coeff hgrad;
        od_coeff vgrad;
        hgrad = vgrad = 0;
        mbctx->c = state->ctmp[pli];
        mbctx->d = state->dtmp;
        mbctx->mc = state->mctmp[pli];
        mbctx->md = state->mdtmp[pli];
        mbctx->l = state->lbuf[pli];
        xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
        ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
        if (mbctx->is_keyframe) {
          od_decode_haar_dc_sb(dec, mbctx, pli, sbx, sby, xdec, ydec,
           sby > 0 && sbx < nhsb - 1, &hgrad, &vgrad);
        }
        od_decode_recursive(dec, mbctx, pli, sbx, sby, OD_NBSIZES - 1, xdec,
         ydec, hgrad, vgrad);
      }
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    if (!mbctx->use_haar_wavelet) {
      od_apply_postfilter_frame_sbs(state->ctmp[pli], w, nhsb, nvsb, xdec,
       ydec);
    }
    for (sby = 0; sby < nvsb; sby++) {
      for (sbx = 0; sbx < nhsb; sbx++) {
        if (mbctx->is_keyframe &&
         OD_BLOCK_SIZE4x4(dec->state.bsize, dec->state.bstride,
         sbx << (OD_NBSIZES - 1), sby << (OD_NBSIZES - 1)) == OD_NBSIZES - 1) {
          int ln;
          ln = OD_LOG_BSIZE_MAX - xdec;
          od_bilinear_smooth(&state->ctmp[pli][(sby << ln)*w + (sbx << ln)],
           ln, w, dec->quantizer[pli], pli);
        }
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
  od_img *ref_img;
  if (dec == NULL || img == NULL || op == NULL) return OD_EFAULT;
  if (dec->packet_state != OD_PACKET_DATA) return OD_EINVAL;
  if (op->e_o_s) dec->packet_state = OD_PACKET_DONE;
  od_ec_dec_init(&dec->ec, op->packet, op->bytes);
  /*Read the packet type bit.*/
  if (od_ec_decode_bool_q15(&dec->ec, 16384)) return OD_EBADPACKET;
  mbctx.is_keyframe = od_ec_decode_bool_q15(&dec->ec, 16384);
  mbctx.use_activity_masking = od_ec_decode_bool_q15(&dec->ec, 16384);
  mbctx.qm = od_ec_decode_bool_q15(&dec->ec, 16384);
  mbctx.use_haar_wavelet = od_ec_decode_bool_q15(&dec->ec, 16384);
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
  if (!mbctx.is_keyframe) {
    /*If there have been no reference frames, and we need one,
       initialize one.*/
    if (dec->state.ref_imgi[OD_FRAME_GOLD] < 0 ||
     dec->state.ref_imgi[OD_FRAME_PREV] < 0 ) {
      od_dec_init_dummy_frame(dec);
    }
  }
  /*Select a free buffer to use for this reference frame.*/
  for (refi = 0; refi == dec->state.ref_imgi[OD_FRAME_GOLD]
   || refi == dec->state.ref_imgi[OD_FRAME_PREV]
   || refi == dec->state.ref_imgi[OD_FRAME_NEXT]; refi++);
  dec->state.ref_imgi[OD_FRAME_SELF] = refi;
  od_adapt_ctx_reset(&dec->state.adapt, mbctx.is_keyframe);
  if (!mbctx.is_keyframe) {
    od_dec_mv_unpack(dec);
    od_state_mc_predict(&dec->state, OD_FRAME_PREV);
    if (dec->user_mc_img != NULL) {
      od_img_copy(dec->user_mc_img, &dec->state.io_imgs[OD_FRAME_REC]);
    }
  }
  od_decode_coefficients(dec, &mbctx);
  if (dec->user_bsize != NULL) {
    int j;
    int nhsb;
    int nvsb;
    nhsb = dec->state.nhsb;
    nvsb = dec->state.nvsb;
    for (j = 0; j < nvsb*4; j++) {
      memcpy(&dec->user_bsize[dec->user_bstride*j],
       &dec->state.bsize[dec->state.bstride*j], nhsb*4);
    }
  }
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  /*Dump YUV*/
  od_state_dump_yuv(&dec->state, dec->state.io_imgs + OD_FRAME_REC, "out");
#endif
  ref_img = dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF];
  OD_ASSERT(ref_img);
  od_img_copy(ref_img, dec->state.io_imgs + OD_FRAME_REC);
  od_img_edge_ext(ref_img);
  /*Return decoded frame.*/
  *img = dec->state.io_imgs[OD_FRAME_REC];
  img->width = dec->state.info.pic_width;
  img->height = dec->state.info.pic_height;
  dec->state.cur_time++;
  return 0;
}
