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
#include "filter.h"
#include "dct.h"
#include "intra.h"
#include "partition.h"
#include "pvq_decoder.h"
#include "block_size.h"
#include "tf.h"
#include "state.h"
#include "quantizer.h"
#include "accounting.h"

static int od_dec_init(od_dec_ctx *dec, const daala_info *info,
 const daala_setup_info *setup) {
  od_img *img;
  od_img_plane *iplane;
  size_t data_sz;
  unsigned char *output_img_data;
  int frame_buf_width;
  int frame_buf_height;
  int plane_buf_width;
  int plane_buf_height;
  int output_bytes;
  int output_bits;
  int ret;
  int pli;
  (void)setup;
  ret = od_state_init(&dec->state, info);
  if (ret < 0) return ret;
  dec->packet_state = OD_PACKET_DATA;
  dec->user_bsize = NULL;
  dec->user_flags = NULL;
  dec->user_mv_grid = NULL;
  dec->user_mc_img = NULL;
  dec->user_dering = NULL;
  data_sz = 0;
  output_bits = 8 + (info->bitdepth_mode - OD_BITDEPTH_MODE_8)*2;
  output_bytes = output_bits > 8 ? 2 : 1;
  /*TODO: Check for overflow before allocating.*/
  frame_buf_width = dec->state.frame_width + (OD_BUFFER_PADDING << 1);
  frame_buf_height = dec->state.frame_height + (OD_BUFFER_PADDING << 1);
  for (pli = 0; pli < info->nplanes; pli++) {
    plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
    plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
    data_sz += plane_buf_width*plane_buf_height*output_bytes;
  }
  dec->output_img_data = output_img_data =
    (unsigned char *)od_aligned_malloc(data_sz, 32);
  if (OD_UNLIKELY(!dec->output_img_data)) {
    return OD_EFAULT;
  }
  img = &dec->output_img;
  img->nplanes = info->nplanes;
  img->width = dec->state.info.pic_width;
  img->height = dec->state.info.pic_height;
  for (pli = 0; pli < img->nplanes; pli++) {
    plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
    plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
    iplane = img->planes + pli;
    iplane->xdec = info->plane_info[pli].xdec;
    iplane->ydec = info->plane_info[pli].ydec;
    iplane->bitdepth = output_bits;
    /*At this moment, our output is always planar.*/
    iplane->xstride = output_bytes;
    iplane->ystride = plane_buf_width*iplane->xstride;
    iplane->data = output_img_data
      + iplane->xstride*(OD_BUFFER_PADDING >> info->plane_info[pli].xdec)
      + iplane->ystride*(OD_BUFFER_PADDING >> info->plane_info[pli].ydec);
    output_img_data += plane_buf_height*iplane->ystride;
  }
#if OD_ACCOUNTING
  od_accounting_init(&dec->acct);
  dec->acct_enabled = 0;
#endif
  return 0;
}

static void od_dec_clear(od_dec_ctx *dec) {
#if OD_ACCOUNTING
  od_accounting_clear(&dec->acct);
#endif
  od_aligned_free(dec->output_img_data);
  od_state_clear(&dec->state);
}

daala_dec_ctx *daala_decode_create(const daala_info *info,
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
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      /*Check that buf is large enough to hold the block sizes for a frame.*/
      OD_RETURN_CHECK(
       buf_sz == sizeof(unsigned char)*dec->state.nvsb*dec->state.nhsb*16,
       OD_EINVAL);
      dec->user_bsize = (unsigned char *)buf;
      dec->user_bstride = dec->state.nhsb*4;
      return OD_SUCCESS;
    }
    case OD_DECCTL_SET_FLAGS_BUFFER : {
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      /*Check that buf is large enough to hold the band flags for a frame.*/
      OD_RETURN_CHECK(
       buf_sz == sizeof(unsigned int)*dec->state.nvsb*8*dec->state.nhsb*8,
       OD_EINVAL);
      dec->user_flags = (unsigned int *)buf;
      dec->user_fstride = dec->state.nhsb*8;
      return OD_SUCCESS;
    }
    case OD_DECCTL_SET_MV_BUFFER : {
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz ==
       sizeof(od_mv_grid_pt)*(dec->state.nhmvbs + 1)*(dec->state.nvmvbs + 1),
       OD_EINVAL);
      dec->user_mv_grid = buf;
      return OD_SUCCESS;
    }
    case OD_DECCTL_SET_MC_IMG : {
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK((buf_sz == sizeof(od_img)), OD_EINVAL);
      dec->user_mc_img = buf;
      return OD_SUCCESS;
    }
#if OD_ACCOUNTING
    case OD_DECCTL_SET_ACCOUNTING_ENABLED: {
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(int), OD_EINVAL);
      dec->acct_enabled = *(int*)buf != 0;
      return OD_SUCCESS;
    }
    case OD_DECCTL_GET_ACCOUNTING : {
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(dec->acct_enabled, OD_EINVAL);
      OD_RETURN_CHECK(buf_sz == sizeof(od_accounting *), OD_EINVAL);
      *(od_accounting **)buf = &dec->acct.acct;
      return OD_SUCCESS;
    }
#endif
    case OD_DECCTL_SET_DERING_BUFFER : {
      int nhdr;
      int nvdr;
      OD_RETURN_CHECK(dec, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      nhdr = dec->state.frame_width >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
      nvdr = dec->state.frame_height >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
      OD_RETURN_CHECK(
       buf_sz == sizeof(unsigned char)*nvdr*nhdr, OD_EINVAL);
      dec->user_dering = (unsigned char *)buf;
      return OD_SUCCESS;
    }
    default: return OD_EIMPL;
  }
}

static void od_dec_blank_img(od_img *img) {
  int pli;
  int frame_buf_height;
  frame_buf_height = img->height + (OD_BUFFER_PADDING << 1);
  for (pli = 0; pli < img->nplanes; pli++) {
    int plane_buf_size;
    int plane_buf_offset;
    plane_buf_size =
     (frame_buf_height >> img->planes[pli].ydec)*img->planes[pli].ystride;
    plane_buf_offset =
     (OD_BUFFER_PADDING >> img->planes[pli].ydec)*img->planes[pli].ystride
     + (OD_BUFFER_PADDING >> img->planes[pli].xdec)*img->planes[pli].xstride;
    memset(img->planes[pli].data - plane_buf_offset, 128, plane_buf_size);
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

static void od_decode_mv(daala_dec_ctx *dec, int num_refs, od_mv_grid_pt *mvg,
 int vx, int vy, int level, int mv_res, int width, int height) {
  generic_encoder *model;
  int pred[2];
  int ox;
  int oy;
  int id;
  int equal_mvs;
  int ref_pred;
  if (num_refs > 1) {
    ref_pred = od_mc_get_ref_predictor(&dec->state, vx, vy, level);
    OD_ASSERT(ref_pred >= 0);
    OD_ASSERT(ref_pred < num_refs);
    mvg->ref = od_decode_cdf_adapt(&dec->ec,
     dec->state.adapt.mv_ref_cdf[ref_pred], num_refs, 256,
     "mv:ref");
  } else {
    mvg->ref = OD_FRAME_PREV;
  }
  equal_mvs = od_state_get_predictor(&dec->state, pred, vx, vy, level,
   mv_res, mvg->ref);
  model = &dec->state.adapt.mv_model;
  id = od_decode_cdf_adapt(&dec->ec, dec->state.adapt.mv_small_cdf[equal_mvs],
   16, dec->state.adapt.mv_small_increment, "mv:low");
  oy = id >> 2;
  ox = id & 0x3;
  if (ox == 3) {
    ox += generic_decode(&dec->ec, model, width << (3 - mv_res),
     &dec->state.adapt.mv_ex[level], 6, "mv:high:x");
  }
  if (oy == 3) {
    oy += generic_decode(&dec->ec, model, height << (3 - mv_res),
     &dec->state.adapt.mv_ey[level], 6, "mv:high:y");
  }
  if (ox && od_ec_dec_bits(&dec->ec, 1, "mv:sign:x")) ox = -ox;
  if (oy && od_ec_dec_bits(&dec->ec, 1, "mv:sign:y")) oy = -oy;
  mvg->mv[0] = (pred[0] + ox)*(1 << mv_res);
  mvg->mv[1] = (pred[1] + oy)*(1 << mv_res);
}

/*Block-level decoder context information.
  Global decoder context information is in od_dec_ctx.*/
struct od_mb_dec_ctx {
  od_coeff *c;
  od_coeff **d;
  od_coeff *md;
  od_coeff *mc;
  od_coeff *l;
  int is_keyframe;
  int num_refs;
  int use_activity_masking;
  int qm;
  int use_haar_wavelet;
  int is_golden_frame;
};
typedef struct od_mb_dec_ctx od_mb_dec_ctx;

static void od_decode_compute_pred(daala_dec_ctx *dec, od_mb_dec_ctx *ctx,
 od_coeff *pred, const od_coeff *d, int bs, int pli, int bx, int by) {
  int n;
  int xdec;
  int w;
  int bo;
  int y;
  int x;
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  n = 1 << (bs + OD_LOG_BSIZE0);
  xdec = dec->output_img.planes[pli].xdec;
  w = dec->state.frame_width >> xdec;
  bo = (by << OD_LOG_BSIZE0)*w + (bx << OD_LOG_BSIZE0);
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  if (ctx->is_keyframe) {
    if (pli == 0 || OD_DISABLE_CFL || ctx->use_haar_wavelet) {
      OD_CLEAR(pred, n*n);
      if (pli == 0 && !ctx->use_haar_wavelet) {
        od_hv_intra_pred(pred, d, w, bx, by, dec->state.bsize,
         dec->state.bstride, bs);
      }
    }
    else {
      od_coeff *l;
      l = ctx->l;
      for (y = 0; y < n; y++) {
        for (x = 0; x < n; x++) {
          pred[n*y + x] = l[n*y + x];
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

#if OD_ACCOUNTING
# define od_ec_dec_unary(ec, str) od_ec_dec_unary_(ec, str)
# define od_decode_coeff_split(dec, sum, ctx, str) od_decode_coeff_split_(dec, sum, ctx, str)
# define od_decode_tree_split(dec, sum, ctx, str) od_decode_tree_split_(dec, sum, ctx, str)
#else
# define od_ec_dec_unary(ec, str) od_ec_dec_unary_(ec)
# define od_decode_coeff_split(dec, sum, ctx, str) od_decode_coeff_split_(dec, sum, ctx)
# define od_decode_tree_split(dec, sum, ctx, str) od_decode_tree_split_(dec, sum, ctx)
#endif

static int od_ec_dec_unary_(od_ec_dec *ec OD_ACC_STR) {
  int ret;
  ret = 0;
  while (od_ec_dec_bits(ec, 1, acc_str) == 0) ret++;
  return ret;
}

static int od_decode_coeff_split_(daala_dec_ctx *dec, int sum, int ctx OD_ACC_STR) {
  int shift;
  int a;
  a = 0;
  if (sum == 0) return 0;
  shift = OD_MAXI(0, OD_ILOG(sum) - 4);
  if (shift) {
    a = od_ec_dec_bits(&dec->ec, shift, acc_str);
  }
  a += od_decode_cdf_adapt(&dec->ec, dec->state.adapt.haar_coeff_cdf[15*ctx
   + (sum >> shift) - 1], (sum >> shift) + 1,
   dec->state.adapt.haar_coeff_increment, acc_str) << shift;
  if (a > sum) {
    a = sum;
    dec->ec.error = 1;
  }
  return a;
}

static int od_decode_tree_split_(daala_dec_ctx *dec, int sum, int ctx OD_ACC_STR) {
  int shift;
  int a;
  a = 0;
  if (sum == 0) return 0;
  shift = OD_MAXI(0, OD_ILOG(sum) - 4);
  if (shift) {
    a = od_ec_dec_bits(&dec->ec, shift, acc_str);
  }
  a += od_decode_cdf_adapt(&dec->ec, dec->state.adapt.haar_split_cdf[15*(2*ctx
   + OD_MINI(shift, 1)) + (sum >> shift) - 1], (sum >> shift) + 1,
   dec->state.adapt.haar_split_increment, acc_str) << shift;
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
   + 3*(OD_ILOG(OD_MAXI(x,y)) - 1), "haar:coeffsplit");
  c[y*n + x] = coeff_mag;
  children_sum = tree_sum - coeff_mag;
  /* Decode sum of each four children relative to tree. */
  if (children_sum) {
    int sum1;
    if (dir == 0) {
      sum1 = od_decode_tree_split(dec, children_sum, 0, "haar:split");
      children[0][0] = od_decode_tree_split(dec, sum1, 2, "haar:split");
      children[0][1] = sum1 - children[0][0];
      children[1][0] = od_decode_tree_split(dec, children_sum - sum1, 2, "haar:split");
      children[1][1] = children_sum - sum1 - children[1][0];
    }
    else {
      sum1 = od_decode_tree_split(dec, children_sum, 1, "haar:split");
      children[0][0] = od_decode_tree_split(dec, sum1, 2, "haar:split");
      children[1][0] = sum1 - children[0][0];
      children[0][1] = od_decode_tree_split(dec, children_sum - sum1, 2, "haar:split");
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
     16, dec->state.adapt.haar_bits_increment, "haar:top");
    if (bits == 15) bits += od_ec_dec_unary(&dec->ec, "haar:top");
    /* Theoretical maximum sum is around 2^7 * 2^OD_COEFF_SHIFT * 32x32,
       so 2^21, but let's play safe. */
    if (bits > 24) {
      /* This can't happen, must be a bit-stream desync/corruption. */
      dec->ec.error = 1;
      return;
    }
    else if (bits > 1) {
      tree_sum[0][0] = (1 << (bits - 1)) | od_ec_dec_bits(&dec->ec, bits - 1,
       "haar:top");
    }
    else tree_sum[0][0] = bits;
    /* Handle diagonal first to make H/V symmetric. */
    tree_sum[1][1] = od_decode_tree_split(dec, tree_sum[0][0], 3, "haar:top");
    tree_sum[0][1] = od_decode_tree_split(dec, tree_sum[0][0] - tree_sum[1][1],
     4, "haar:top");
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
        sign = od_ec_dec_bits(&dec->ec, 1, "haar:sign");
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
  int use_activity_masking;
  const int *qm;
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  n = 1 << (bs + 2);
  lossless = OD_LOSSLESS(dec, pli);
  use_activity_masking = ctx->use_activity_masking;
  qm = ctx->qm == OD_HVS_QM ? OD_QM8_Q4_HVS : OD_QM8_Q4_FLAT;
  bx <<= bs;
  by <<= bs;
  xdec = dec->output_img.planes[pli].xdec;
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
  od_decode_compute_pred(dec, ctx, pred, d, bs, pli, bx, by);
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
    /*Safely initialize d since some coeffs are skipped by PVQ.*/
    od_init_skipped_coeffs(d, pred, ctx->is_keyframe, bo, n, w);
    od_raster_to_coding_order(predt,  n, &pred[0], n);
  }
  quant = OD_MAXI(1, dec->state.quantizer[pli]);
  if (lossless) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, quant*
     dec->state.pvq_qm_q4[pli][od_qm_get_index(bs, 0)] >> 4);
  }
  if (ctx->use_haar_wavelet) {
    od_wavelet_unquantize(dec, bs + 2, pred, predt,
     dec->state.quantizer[pli], pli);
  }
  else {
    unsigned int flags;
    od_pvq_decode(dec, predt, pred, quant, pli, bs,
     OD_PVQ_BETA[use_activity_masking][pli][bs], OD_ROBUST_STREAM,
     ctx->is_keyframe, &flags, skip);
    if (pli == 0 && dec->user_flags != NULL) {
      dec->user_flags[by*dec->user_fstride + bx] = flags;
    }
  }
  if (!ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && !ctx->use_haar_wavelet;
    if (!has_dc_skip || pred[0]) {
      pred[0] = has_dc_skip + generic_decode(&dec->ec,
       &dec->state.adapt.model_dc[pli], -1,
       &dec->state.adapt.ex_dc[pli][bs][0], 2, "dc:mag");
      if (pred[0]) pred[0] *= od_ec_dec_bits(&dec->ec, 1, "dc:sign") ? -1 : 1;
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
  od_coeff *d;
  int nhsb;
  int quant;
  int ln;
  od_coeff sb_dc_pred;
  od_coeff sb_dc_curr;
  od_coeff *sb_dc_mem;
  (void)ydec;
  d = ctx->d[pli];
  w = dec->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  if (OD_LOSSLESS(dec, pli)) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, dec->state.quantizer[pli]*
     dec->state.pvq_qm_q4[pli][od_qm_get_index(OD_NBSIZES - 1, 0)] >> 4);
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
   &dec->state.adapt.ex_sb_dc[pli], 2, "haardc:mag:top");
  if (quant) {
    if (od_ec_dec_bits(&dec->ec, 1, "haardc:sign:top")) quant = -quant;
  }
  sb_dc_curr = quant*dc_quant + sb_dc_pred;
  d[(by << ln)*w + (bx << ln)] = sb_dc_curr;
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
  if (OD_LOSSLESS(dec, pli)) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, dec->state.quantizer[pli]*
     dec->state.pvq_qm_q4[pli][od_qm_get_index(OD_NBSIZES - 1, 0)] >> 4);
  }
  if (OD_LOSSLESS(dec, pli)) ac_quant[0] = ac_quant[1] = 1;
  else {
    ac_quant[0] = (dc_quant*OD_DC_QM[bsi - xdec][0] + 8) >> 4;
    ac_quant[1] = (dc_quant*OD_DC_QM[bsi - xdec][1] + 8) >> 4;
  }
  ln = bsi - xdec + 2;
  x[0] = ctx->d[pli][(by << ln)*w + (bx << ln)];
  for (i = 1; i < 4; i++) {
    int quant;
    quant = generic_decode(&dec->ec, &dec->state.adapt.model_dc[pli], -1,
     &dec->state.adapt.ex_dc[pli][bsi][i-1], 2, "haardc:mag:level");
    if (quant) {
      if (od_ec_dec_bits(&dec->ec, 1, "haardc:sign:level")) quant = -quant;
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

#if OD_SIGNAL_Q_SCALING
static void od_decode_quantizer_scaling(daala_dec_ctx *dec, int bx, int by,
 int skip) {
  int sbx;
  int sby;
  int q_scaling;
  OD_ASSERT(skip == !!skip);
  sbx = bx;
  sby = by;
  if (!skip) {
    int above;
    int left;
    above = sby > 0 ? dec->state.sb_q_scaling[(sby - 1)*dec->state.nhsb + sbx]
     : 0;
    left = sbx > 0 ? dec->state.sb_q_scaling[sby*dec->state.nhsb + (sbx - 1)]
     : 0;
    q_scaling = od_decode_cdf_adapt(&dec->ec,
     dec->state.adapt.q_cdf[above + left*4], 4,
     dec->state.adapt.q_increment, "quant");
  }
  else {
    q_scaling = 0;
  }
  dec->state.sb_q_scaling[sby*dec->state.nhsb + sbx] = q_scaling;
}
#endif

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
  OD_ACCOUNTING_SET_LOCATION(dec, pli, bsi, bx << bsi, by << bsi);
  /* Read the luma skip symbol. A value of 4 means "split the block", while < 4
     means that we code the block. In the latter case, we need to forward
     the skip value to the PVQ decoder. */
  if (ctx->use_haar_wavelet) obs = bsi;
  else if (pli == 0) {
    skip = od_decode_cdf_adapt(&dec->ec,
     dec->state.adapt.skip_cdf[2*bsi + (pli != 0)], 4 + (bsi > 0),
     dec->state.adapt.skip_increment, "skip");
#if OD_SIGNAL_Q_SCALING
    if (bsi == OD_NBSIZES - 1) {
      od_decode_quantizer_scaling(dec, bx, by, skip == 2);
    }
#endif
    if (skip < 4) obs = bsi;
    else obs = -1;
  }
  bs = OD_MAXI(obs, xdec);
  OD_ASSERT(bs <= bsi);
  if (bs == bsi) {
    int i;
    int j;
    bs -= xdec;
    if (pli == 0) {
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
      od_resample_luma_coeffs(ctx->l, 1 << (bs + OD_LOG_BSIZE0),
       ctx->d[0] + (by << (2 + bsi))*frame_width + (bx << (2 + bsi)),
       frame_width, xdec, ydec, bs, obs);
    }
    if (pli > 0 && !ctx->use_haar_wavelet) {
      /* Decode the skip for chroma. */
      skip = od_decode_cdf_adapt(&dec->ec,
       dec->state.adapt.skip_cdf[2*bsi + (pli != 0)], 4,
       dec->state.adapt.skip_increment, "skip");
    }
    od_block_decode(dec, ctx, bs, pli, bx, by, skip);
    for (i = 0; i < 1 << bs; i++) {
      for (j = 0; j < 1 << bs; j++) {
        dec->state.bskip[pli][((by << bs) + i)*dec->state.skip_stride
         + (bx << bs) + j] = (skip == 2) && !ctx->is_keyframe;
      }
    }

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
    if (ctx->is_keyframe) {
      od_decode_haar_dc_level(dec, ctx, pli, 2*bx, 2*by, bsi - 1, xdec, &hgrad,
       &vgrad);
    }
    od_decode_recursive(dec, ctx, pli, 2*bx + 0, 2*by + 0, bsi - 1, xdec, ydec,
     hgrad, vgrad);
    od_decode_recursive(dec, ctx, pli, 2*bx + 1, 2*by + 0, bsi - 1, xdec, ydec,
     hgrad, vgrad);
    od_decode_recursive(dec, ctx, pli, 2*bx + 0, 2*by + 1, bsi - 1, xdec, ydec,
     hgrad, vgrad);
    od_decode_recursive(dec, ctx, pli, 2*bx + 1, 2*by + 1, bsi - 1, xdec, ydec,
     hgrad, vgrad);
    bs = bsi - xdec;
    bo = (by << (OD_LOG_BSIZE0 + bs))*w + (bx << (OD_LOG_BSIZE0 + bs));
    od_postfilter_split(ctx->c + bo, w, bs, f, dec->state.coded_quantizer[pli],
     &dec->state.bskip[pli][(by << bs)*dec->state.skip_stride + (bx << bs)],
     dec->state.skip_stride);
  }
}

static void od_dec_mv_unpack(daala_dec_ctx *dec, int num_refs) {
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
  img = dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF];
  mv_res = od_ec_dec_uint(&dec->ec, 3, "mv:res");
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
      OD_ACCOUNTING_SET_LOCATION(dec, OD_ACCT_MV, 0, vx, vy);
      mvp = grid[vy] + vx;
      mvp->valid = 1;
      od_decode_mv(dec, num_refs, mvp, vx, vy, 0, mv_res, width, height);
    }
  }
  for (log_mvb_sz = OD_LOG_MVB_DELTA0, level = 1; log_mvb_sz-- > 0; level++) {
    int mvb_sz;
    mvb_sz = 1 << log_mvb_sz;
    /*Odd levels.*/
    for (vy = mvb_sz; vy <= nvmvbs; vy += 2*mvb_sz) {
      for (vx = mvb_sz; vx <= nhmvbs; vx += 2*mvb_sz) {
        OD_ACCOUNTING_SET_LOCATION(dec, OD_ACCT_MV, level, vx, vy);
        if (grid[vy - mvb_sz][vx - mvb_sz].valid
         && grid[vy - mvb_sz][vx + mvb_sz].valid
         && grid[vy + mvb_sz][vx + mvb_sz].valid
         && grid[vy + mvb_sz][vx - mvb_sz].valid) {
          cdf = od_mv_split_flag_cdf(&dec->state, vx, vy, level);
          mvp = grid[vy] + vx;
          mvp->valid = od_decode_cdf_adapt(&dec->ec,
           cdf, 2, dec->state.adapt.split_flag_increment, "mv:valid");
          if (mvp->valid) {
            od_decode_mv(dec, num_refs, mvp, vx, vy, level, mv_res,
             width, height);
          }
        }
      }
    }
    level++;
    /*Even Levels.*/
    for (vy = 0; vy <= nvmvbs; vy += mvb_sz) {
      for (vx = mvb_sz*!(vy & mvb_sz); vx <= nhmvbs; vx += 2*mvb_sz) {
        OD_ACCOUNTING_SET_LOCATION(dec, OD_ACCT_MV, level, vx, vy);
        if ((vy - mvb_sz < 0 || grid[vy - mvb_sz][vx].valid)
         && (vx - mvb_sz < 0 || grid[vy][vx - mvb_sz].valid)
         && (vy + mvb_sz > nvmvbs || grid[vy + mvb_sz][vx].valid)
         && (vx + mvb_sz > nhmvbs || grid[vy][vx + mvb_sz].valid)) {
          cdf = od_mv_split_flag_cdf(&dec->state, vx, vy, level);
          mvp = grid[vy] + vx;
          mvp->valid = od_decode_cdf_adapt(&dec->ec,
           cdf, 2, dec->state.adapt.split_flag_increment, "mv:valid");
          if (mvp->valid) {
            od_decode_mv(dec, num_refs, mvp, vx, vy, level, mv_res,
             width, height);
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
  int w;
  int y;
  int x;
  int frame_width;
  int nvsb;
  int nhsb;
  int nhdr;
  int nvdr;
  od_state *state;
  od_img *rec;
  state = &dec->state;
  /*Initialize the data needed for each plane.*/
  nplanes = state->info.nplanes;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  frame_width = state->frame_width;
  rec = state->ref_imgs + state->ref_imgi[OD_FRAME_SELF];
  OD_ACCOUNTING_SET_LOCATION(dec, OD_ACCT_FRAME, 0, 0, 0);
  /* Map our quantizers; we potentially need them to know what reference
     resolution we're working at. */
  for (pli = 0; pli < nplanes; pli++) {
    dec->state.coded_quantizer[pli] = od_ec_dec_uint(&dec->ec,
     OD_N_CODED_QUANTIZERS, "quantizer");
    dec->state.quantizer[pli] =
     od_codedquantizer_to_quantizer(dec->state.coded_quantizer[pli]);
  }
  /*Apply the prefilter to the motion-compensated reference.*/
  if (!mbctx->is_keyframe) {
    for (pli = 0; pli < nplanes; pli++) {
      xdec = rec->planes[pli].xdec;
      ydec = rec->planes[pli].ydec;
      w = frame_width >> xdec;
      /*Collect the image data needed for this plane.*/
      od_ref_plane_to_coeff(state,
       state->mctmp[pli], OD_LOSSLESS(dec, pli), rec, pli);
      if (!mbctx->use_haar_wavelet) {
        od_apply_prefilter_frame_sbs(state->mctmp[pli], w, nhsb, nvsb, xdec,
         ydec);
      }
    }
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
        xdec = dec->output_img.planes[pli].xdec;
        ydec = dec->output_img.planes[pli].ydec;
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
    xdec = dec->output_img.planes[pli].xdec;
    ydec = dec->output_img.planes[pli].ydec;
    w = frame_width >> xdec;
    if (!mbctx->use_haar_wavelet) {
      od_apply_postfilter_frame_sbs(state->ctmp[pli], w, nhsb, nvsb, xdec,
       ydec, dec->state.coded_quantizer[pli], &dec->state.bskip[pli][0],
       dec->state.skip_stride);
    }
  }
  nhdr = state->frame_width >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
  nvdr = state->frame_height >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
  if (dec->state.quantizer[0] > 0) {
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      int size;
      xdec = dec->output_img.planes[pli].xdec;
      ydec = dec->output_img.planes[pli].ydec;
      size = nvsb*nhsb*OD_BSIZE_MAX*OD_BSIZE_MAX >> xdec >> ydec;
      for (i = 0; i < size; i++) {
        state->etmp[pli][i] = state->ctmp[pli][i];
      }
    }
    for (sby = 0; sby < nvdr; sby++) {
      for (sbx = 0; sbx < nhdr; sbx++) {
        int filtered;
        int c;
        int up;
        int left;
        int i;
        int j;
        unsigned char *bskip;
        state->dering_flags[sby*nhdr + sbx] = 0;
        bskip = dec->state.bskip[0] +
         (sby << OD_LOG_DERING_GRID)*dec->state.skip_stride +
         (sbx << OD_LOG_DERING_GRID);
        for (j = 0; j < 1 << OD_LOG_DERING_GRID; j++) {
          for (i = 0; i < 1 << OD_LOG_DERING_GRID; i++) {
            if (!bskip[j*dec->state.skip_stride + i]) {
              state->dering_flags[sby*nhdr + sbx] = 1;
            }
          }
        }
        if (!state->dering_flags[sby*nhdr + sbx]) {
          continue;
        }
        up = 0;
        if (sby > 0) {
          up = state->dering_flags[(sby - 1)*nhdr + sbx];
        }
        left = 0;
        if (sbx > 0) {
          left = state->dering_flags[sby*nhdr + (sbx - 1)];
        }
        c = (up << 1) + left;
        filtered = od_decode_cdf_adapt(&dec->ec, state->adapt.clpf_cdf[c], 2,
         state->adapt.clpf_increment, "clp");
        state->dering_flags[sby*nhdr + sbx] = filtered;
        if (filtered) {
          for (pli = 0; pli < nplanes; pli++) {
            int16_t buf[OD_BSIZE_MAX*OD_BSIZE_MAX];
            od_coeff *output;
            int ln;
            int n;
            int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
            xdec = dec->output_img.planes[pli].xdec;
            ydec = dec->output_img.planes[pli].ydec;
            w = frame_width >> xdec;
            ln = OD_LOG_DERING_GRID + OD_LOG_BSIZE0 - xdec;
            n = 1 << ln;
            OD_ASSERT(xdec == ydec);
            /*buf is used for output so that we don't use filtered pixels in
              the input to the filter, but because we look past block edges,
              we do this anyway on the edge pixels. Unfortunately, this limits
              potential parallelism.*/
            od_dering(state, buf, n,
             &state->etmp[pli][(sby << ln)*w +
             (sbx << ln)], w, ln, sbx, sby, nhdr, nvdr,
             dec->state.quantizer[pli], xdec, dir, pli, &dec->state.bskip[pli]
             [(sby << (OD_LOG_DERING_GRID - ydec))*dec->state.skip_stride
             + (sbx << (OD_LOG_DERING_GRID - xdec))], dec->state.skip_stride);
            output = &state->ctmp[pli][(sby << ln)*w + (sbx << ln)];
            for (y = 0; y < n; y++) {
              for (x = 0; x < n; x++) {
                output[y*w + x] = buf[y*n+ x];
              }
            }
          }
        }
      }
    }
    if (dec->user_dering != NULL) {
      for (sby = 0; sby < nvdr; sby++) {
        for (sbx = 0; sbx < nhdr; sbx++) {
          dec->user_dering[sby*nhdr + sbx] =
           state->dering_flags[sby*nhdr + sbx];
        }
      }
    }
  }
  else {
    if (dec->user_dering != NULL) {
      OD_CLEAR(dec->user_dering, nhdr*nvdr);
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    xdec = dec->output_img.planes[pli].xdec;
    ydec = dec->output_img.planes[pli].ydec;
    w = frame_width >> xdec;
    if (dec->state.quantizer[0] > 0) {
      for (sby = 0; sby < nvsb; sby++) {
        for (sbx = 0; sbx < nhsb; sbx++) {
          if (mbctx->is_keyframe) {
            od_smooth_recursive(state->ctmp[pli], dec->state.bsize,
             dec->state.bstride, sbx, sby, OD_NBSIZES - 1, w, xdec, ydec,
             OD_BLOCK_32X32, dec->state.quantizer[pli], pli);
          }
        }
      }
    }
    /*Move/scale/shift reconstructed data values from transform
      storage back into the SELF reference frame.*/
    od_coeff_to_ref_plane(state, rec, pli,
     state->ctmp[pli], OD_LOSSLESS(dec, pli));
  }
}

int daala_decode_packet_in(daala_dec_ctx *dec, od_img *img,
 const daala_packet *op) {
  int refi;
  od_mb_dec_ctx mbctx;
  od_img *ref_img;
  if (dec == NULL || img == NULL || op == NULL) return OD_EFAULT;
  if (dec->packet_state != OD_PACKET_DATA) return OD_EINVAL;
  if (op->e_o_s) dec->packet_state = OD_PACKET_DONE;
  od_ec_dec_init(&dec->ec, op->packet, op->bytes);
#if OD_ACCOUNTING
  if (dec->acct_enabled) {
    od_accounting_reset(&dec->acct);
    dec->ec.acct = &dec->acct;
  }
  else dec->ec.acct = NULL;
#endif
  OD_ACCOUNTING_SET_LOCATION(dec, OD_ACCT_FRAME, 0, 0, 0);
  /*Read the packet type bit.*/
  if (od_ec_decode_bool_q15(&dec->ec, 16384, "flags")) return OD_EBADPACKET;
  mbctx.is_keyframe = od_ec_decode_bool_q15(&dec->ec, 16384, "flags");
  if (!mbctx.is_keyframe) {
    mbctx.num_refs = od_ec_dec_uint(&dec->ec, OD_MAX_CODED_REFS, "flags") + 1;
  } else {
    mbctx.num_refs = 0;
  }
  mbctx.use_activity_masking = od_ec_decode_bool_q15(&dec->ec, 16384, "flags");
  mbctx.qm = od_ec_decode_bool_q15(&dec->ec, 16384, "flags");
  mbctx.use_haar_wavelet = od_ec_decode_bool_q15(&dec->ec, 16384, "flags");
  mbctx.is_golden_frame = od_ec_decode_bool_q15(&dec->ec, 16384, "flags");
  if (mbctx.is_keyframe) {
    int nplanes;
    int pli;
    nplanes = dec->state.info.nplanes;
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      for (i = 0; i < OD_QM_SIZE; i++) {
        dec->state.pvq_qm_q4[pli][i] = od_ec_dec_bits(&dec->ec, 8, "qm");
      }
    }
  }
  /*Update the buffer state.*/
  if (dec->state.ref_imgi[OD_FRAME_SELF] >= 0) {
    dec->state.ref_imgi[OD_FRAME_PREV] =
     dec->state.ref_imgi[OD_FRAME_SELF];
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
    od_dec_mv_unpack(dec, mbctx.num_refs);
    od_state_mc_predict(&dec->state,
     dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF]);
    if (dec->user_mc_img != NULL) {
      od_img_copy(dec->user_mc_img,
       dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF]);
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
  ref_img = dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF];
  OD_ASSERT(ref_img);
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  /*Dump YUV*/
  od_state_dump_yuv(&dec->state, ref_img, "out");
#endif
  od_img_edge_ext(ref_img);
  /*Return decoded frame.*/
  od_img_copy(&dec->output_img, ref_img);
  *img = dec->output_img;
  dec->state.cur_time++;
  if (mbctx.is_golden_frame) {
    dec->state.ref_imgi[OD_FRAME_GOLD] =
     dec->state.ref_imgi[OD_FRAME_SELF];
  }
  return 0;
}
