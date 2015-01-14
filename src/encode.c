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
#include "encint.h"
#if defined(OD_ENCODER_CHECK)
# include "decint.h"
#endif
#include "generic_code.h"
#include "filter.h"
#include "dct.h"
#include "intra.h"
#include "logging.h"
#include "partition.h"
#include "pvq.h"
#include "pvq_code.h"
#include "block_size.h"
#include "logging.h"
#include "tf.h"
#include "accounting.h"
#include "state.h"
#include "mcenc.h"
#if defined(OD_X86ASM)
# include "x86/x86int.h"
#endif

#if defined(OD_ACCOUNTING)
# define OD_ENC_ACCT_UPDATE(enc, cat, value) \
 OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec), cat, value)
#else
# define OD_ENC_ACCT_UPDATE(enc, cat, value)
#endif

static int od_quantizer_from_quality(int quality) {
  return quality == 0 ? 0 :
   (quality << OD_COEFF_SHIFT >> OD_QUALITY_SHIFT) +
   (1 << OD_COEFF_SHIFT >> 1);
}

void od_enc_opt_vtbl_init_c(od_enc_ctx *enc) {
  enc->opt_vtbl.mc_compute_sad_4x4_xstride_1 =
   od_mc_compute_sad_4x4_xstride_1_c;
  enc->opt_vtbl.mc_compute_sad_8x8_xstride_1 =
   od_mc_compute_sad_8x8_xstride_1_c;
  enc->opt_vtbl.mc_compute_sad_16x16_xstride_1 =
   od_mc_compute_sad_16x16_xstride_1_c;
}

static void od_enc_opt_vtbl_init(od_enc_ctx *enc) {
#if defined(OD_X86ASM)
  od_enc_opt_vtbl_init_x86(enc);
#else
  od_enc_opt_vtbl_init_c(enc);
#endif
}

static int od_enc_init(od_enc_ctx *enc, const daala_info *info) {
  int i;
  int ret;
  ret = od_state_init(&enc->state, info);
  if (ret < 0) return ret;
  od_enc_opt_vtbl_init(enc);
  oggbyte_writeinit(&enc->obb);
  od_ec_enc_init(&enc->ec, 65025);
  enc->packet_state = OD_PACKET_INFO_HDR;
  for (i = 0; i < OD_NPLANES_MAX; i++){
    enc->quality[i] = 10;
  }
  enc->mvest = od_mv_est_alloc(enc);
  if (OD_UNLIKELY(!enc->mvest)) {
    return OD_EFAULT;
  }
  enc->params.mv_level_min = 0;
  enc->params.mv_level_max = 4;
#if defined(OD_ACCOUNTING)
  od_acct_init(&enc->acct);
#endif
  enc->bs = (od_block_size_comp *)_ogg_malloc(sizeof(*enc->bs));
#if defined(OD_ENCODER_CHECK)
  enc->dec = daala_decode_alloc(info, NULL);
#endif
  return 0;
}

static void od_enc_clear(od_enc_ctx *enc) {
  od_mv_est_free(enc->mvest);
  od_ec_enc_clear(&enc->ec);
  oggbyte_writeclear(&enc->obb);
  od_state_clear(&enc->state);
#if defined(OD_ACCOUNTING)
  od_acct_clear(&enc->acct);
#endif
}

daala_enc_ctx *daala_encode_create(const daala_info *info) {
  od_enc_ctx *enc;
  if (info == NULL) return NULL;
  enc = (od_enc_ctx *)_ogg_malloc(sizeof(*enc));
  if (od_enc_init(enc, info) < 0) {
    _ogg_free(enc);
    return NULL;
  }
  return enc;
}

void daala_encode_free(daala_enc_ctx *enc) {
  if (enc != NULL) {
#if defined(OD_ENCODER_CHECK)
    if (enc->dec != NULL) {
      daala_decode_free(enc->dec);
    }
#endif
    _ogg_free(enc->bs);
    od_enc_clear(enc);
    _ogg_free(enc);
  }
}

int daala_encode_ctl(daala_enc_ctx *enc, int req, void *buf, size_t buf_sz) {
  (void)buf;
  (void)buf_sz;
  switch (req) {
    case OD_SET_QUANT:
    {
      int i;
      OD_ASSERT(enc);
      OD_ASSERT(buf);
      OD_ASSERT(buf_sz == sizeof(*enc->quality));
      for (i = 0; i < OD_NPLANES_MAX; i++){
        int tmp = *(int *)buf;
        enc->quality[i] = tmp > 0 ? (tmp << OD_QUALITY_SHIFT) - 8 : 0;
      }
      return OD_SUCCESS;
    }
    case OD_SET_MC_USE_CHROMA:
    {
      int mc_use_chroma;
      OD_ASSERT(enc);
      OD_ASSERT(buf);
      OD_ASSERT(buf_sz == sizeof(mc_use_chroma));
      mc_use_chroma = *(int *)buf;
      if (mc_use_chroma) {
        enc->mvest->flags |= OD_MC_USE_CHROMA;
      }
      else {
        enc->mvest->flags &= ~OD_MC_USE_CHROMA;
      }
      return OD_SUCCESS;
    }
    case OD_SET_MV_RES_MIN:
    {
      int mv_res_min;
      OD_ASSERT(enc);
      OD_ASSERT(buf);
      OD_ASSERT(buf_sz == sizeof(mv_res_min));
      mv_res_min = *(int *)buf;
      if (mv_res_min < 0 || mv_res_min > 2) {
        return OD_EINVAL;
      }
      enc->mvest->mv_res_min = mv_res_min;
      return OD_SUCCESS;
    }
    case OD_SET_MV_LEVEL_MIN:
    {
      int mv_level_min;
      OD_ASSERT(enc);
      OD_ASSERT(buf);
      OD_ASSERT(buf_sz == sizeof(mv_level_min));
      mv_level_min = *(int *)buf;
      if (mv_level_min < 0 || mv_level_min > 4) {
        return OD_EINVAL;
      }
      enc->params.mv_level_min = mv_level_min;
      return OD_SUCCESS;
    }
    case OD_SET_MV_LEVEL_MAX:
    {
      int mv_level_max;
      OD_ASSERT(enc);
      OD_ASSERT(buf);
      OD_ASSERT(buf_sz == sizeof(mv_level_max));
      mv_level_max = *(int *)buf;
      if (mv_level_max < 0 || mv_level_max > 4) {
        return OD_EINVAL;
      }
      enc->params.mv_level_max = mv_level_max;
      return OD_SUCCESS;
    }
    default: return OD_EIMPL;
  }
}

void od_encode_checkpoint(const daala_enc_ctx *enc, od_rollback_buffer *rbuf) {
  od_ec_enc_checkpoint(&rbuf->ec, &enc->ec);
  OD_COPY(&rbuf->adapt, &enc->state.adapt, 1);
}

void od_encode_rollback(daala_enc_ctx *enc, const od_rollback_buffer *rbuf) {
  od_ec_enc_rollback(&enc->ec, &rbuf->ec);
  OD_COPY(&enc->state.adapt, &rbuf->adapt, 1);
}

static void od_img_plane_copy_pad8(od_img_plane *dst_p,
 int plane_width, int plane_height, od_img_plane *src_p,
 int pic_width, int pic_height) {
  unsigned char *dst_data;
  ptrdiff_t dstride;
  int y;
  dstride = dst_p->ystride;
  /*If we have _no_ data, just encode a dull green.*/
  if (pic_width == 0 || pic_height == 0) {
    dst_data = dst_p->data;
    for (y = 0; y < plane_height; y++) {
      OD_CLEAR(dst_data, plane_width);
      dst_data += dstride;
    }
  }
  /*Otherwise, copy what we do have, and add our own padding.*/
  else {
    unsigned char *src_data;
    unsigned char *dst;
    ptrdiff_t sxstride;
    ptrdiff_t systride;
    int x;
    /*Step 1: Copy the data we do have.*/
    sxstride = src_p->xstride;
    systride = src_p->ystride;
    dst_data = dst_p->data;
    src_data = src_p->data;
    dst = dst_data;
    for (y = 0; y < pic_height; y++) {
      if (sxstride == 1) OD_COPY(dst, src_data, pic_width);
      else for (x = 0; x < pic_width; x++) dst[x] = *(src_data + sxstride*x);
      dst += dstride;
      src_data += systride;
    }
    /*Step 2: Perform a low-pass extension into the padding region.*/
    /*Right side.*/
    for (x = pic_width; x < plane_width; x++) {
      dst = dst_data + x - 1;
      for (y = 0; y < pic_height; y++) {
        dst[1] = (2*dst[0] + (dst - (dstride & -(y > 0)))[0]
         + (dst + (dstride & -(y + 1 < pic_height)))[0] + 2) >> 2;
        dst += dstride;
      }
    }
    /*Bottom.*/
    dst = dst_data + dstride*pic_height;
    for (y = pic_height; y < plane_height; y++) {
      for (x = 0; x < plane_width; x++) {
        dst[x] = (2*(dst - dstride)[x] + (dst - dstride)[x - (x > 0)]
         + (dst - dstride)[x + (x + 1 < plane_width)] + 2) >> 2;
      }
      dst += dstride;
    }
  }
}

/*Extend the edge into the padding.*/
static void od_img_plane_edge_ext8(od_img_plane *dst_p,
 int plane_width, int plane_height, int horz_padding, int vert_padding) {
  ptrdiff_t dstride;
  unsigned char *dst_data;
  unsigned char *dst;
  int x;
  int y;
  dstride = dst_p->ystride;
  dst_data = dst_p->data;
  /*Left side.*/
  for (y = 0; y < plane_height; y++) {
    dst = dst_data + dstride * y;
    for (x = 1; x <= horz_padding; x++) {
      (dst-x)[0] = dst[0];
    }
  }
  /*Right side.*/
  for (y = 0; y < plane_height; y++) {
    dst = dst_data + plane_width - 1 + dstride * y;
    for (x = 1; x <= horz_padding; x++) {
      dst[x] = dst[0];
    }
  }
  /*Top.*/
  dst = dst_data - horz_padding;
  for (y = 0; y < vert_padding; y++) {
    for (x = 0; x < plane_width + 2 * horz_padding; x++) {
      (dst - dstride)[x] = dst[x];
    }
    dst -= dstride;
  }
  /*Bottom.*/
  dst = dst_data - horz_padding + plane_height * dstride;
  for (y = 0; y < vert_padding; y++) {
    for (x = 0; x < plane_width + 2 * horz_padding; x++) {
      dst[x] = (dst - dstride)[x];
    }
    dst += dstride;
  }
}

struct od_mb_enc_ctx {
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
typedef struct od_mb_enc_ctx od_mb_enc_ctx;

#if !defined(OD_DUMP_COEFFS)
static void od_encode_compute_pred(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, od_coeff *pred,
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
  xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = enc->state.frame_width;
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

static void od_single_band_lossless_encode(daala_enc_ctx *enc, int ln,
 od_coeff *scalar_out, const od_coeff *cblock, const od_coeff *predt,
 int pli) {
  int *adapt;
  int vk;
  int zzi;
  int n2;
  ogg_int32_t adapt_curr[OD_NSB_ADAPT_CTXS];
  adapt = enc->state.adapt.pvq_adapt;
  vk = 0;
  n2 = 1 << (2*ln + 4);
  for (zzi = 1; zzi < n2; zzi++) {
    scalar_out[zzi] = cblock[zzi] - predt[zzi];
    vk += abs(scalar_out[zzi]);
  }
  generic_encode(&enc->ec, &enc->state.adapt.model_g[pli], vk, -1,
   &enc->state.adapt.ex_g[pli][ln], 0);
  laplace_encode_vector(&enc->ec, scalar_out + 1, n2 - 1, vk, adapt_curr,
   adapt);
  for (zzi = 1; zzi < n2; zzi++) {
    scalar_out[zzi] = scalar_out[zzi] + predt[zzi];
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

static void od_block_encode(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int ln,
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
  od_coeff cblock[OD_BSIZE_MAX*OD_BSIZE_MAX];
  od_coeff scalar_out[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int quant;
  int dc_quant;
  int lossless;
#if defined(OD_OUTPUT_PRED)
  od_coeff preds[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int zzi;
#endif
  OD_ASSERT(ln >= 0 && ln <= 3);
  n = 1 << (ln + 2);
  bx <<= ln;
  by <<= ln;
  xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = enc->state.frame_width;
  w = frame_width >> xdec;
  c = ctx->c;
  d = ctx->d[pli];
  md = ctx->md;
  mc = ctx->mc;
  /* Apply forward transform. */
  if (!ctx->is_keyframe) {
    (*enc->state.opt_vtbl.fdct_2d[ln])(md + (by << 2)*w + (bx << 2), w,
     mc + (by << 2)*w + (bx << 2), w);
  }
  od_encode_compute_pred(enc, ctx, pred, ln, pli, bx, by);
  if (ctx->is_keyframe && pli == 0) {
    od_hv_intra_pred(pred, d, w, bx, by, enc->state.bsize,
     enc->state.bstride, ln);
  }
  lossless = (enc->quantizer[pli] == 0);
#if defined(OD_OUTPUT_PRED)
  for (zzi = 0; zzi < (n*n); zzi++) preds[zzi] = pred[zzi];
#endif
  /* Change ordering for encoding. */
  od_raster_to_coding_order(cblock,  n, &d[((by << 2))*w + (bx << 2)], w,
   lossless);
  od_raster_to_coding_order(predt,  n, &pred[0], n, lossless);
  /* Lossless encoding uses an actual quantizer of 1, but is signalled
     with a 'quantizer' of 0. */
  quant = OD_MAXI(1, enc->quantizer[pli]);
  if (lossless) dc_quant = quant;
  else {
    dc_quant = OD_MAXI(1, quant*
     enc->state.pvq_qm_q4[pli][od_qm_get_index(ln, 0)] >> 4);
  }
  /* This quantization may be overridden in the PVQ code for full RDO. */
  if (OD_DISABLE_HAAR_DC || !ctx->is_keyframe) {
    if (abs(cblock[0] - predt[0]) < dc_quant * 141 / 256) { /* 0.55 */
      scalar_out[0] = 0;
    }
    else {
      scalar_out[0] = OD_DIV_R0(cblock[0] - predt[0], dc_quant);
    }
  }
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_AC_COEFFS);
  if (lossless) {
    od_single_band_lossless_encode(enc, ln, scalar_out, cblock, predt, pli);
  }
  else {
    pvq_encode(enc, predt, cblock, scalar_out, quant, pli, ln,
     OD_PVQ_BETA[pli][ln], OD_ROBUST_STREAM, ctx->is_keyframe);
  }
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  if (OD_DISABLE_HAAR_DC || !ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && !lossless;
    OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_DC_COEFF);
    if (!has_dc_skip || scalar_out[0]) {
      generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli],
       abs(scalar_out[0]) - has_dc_skip, -1, &enc->state.adapt.ex_dc[pli][ln][0], 2);
    }
    if (scalar_out[0]) od_ec_enc_bits(&enc->ec, scalar_out[0] < 0, 1);
    OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
    scalar_out[0] = scalar_out[0]*dc_quant;
    scalar_out[0] += predt[0];
  }
  else {
    scalar_out[0] = cblock[0];
  }
  od_coding_order_to_raster(&d[((by << 2))*w + (bx << 2)], w, scalar_out, n,
   lossless);
  /*Apply the inverse transform.*/
#if !defined(OD_OUTPUT_PRED)
  (*enc->state.opt_vtbl.idct_2d[ln])(c + (by << 2)*w + (bx << 2), w,
   d + (by << 2)*w + (bx << 2), w);
#else
# if 0
  /*Output the resampled luma plane.*/
  if (pli != 0) {
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        preds[y*n + x] = l[((by << 2) + y)*w + (bx << 2) + x] >> xdec;
      }
    }
  }
# endif
  (*enc->state.opt_vtbl.idct_2d[ln])(c + (by << 2)*w + (bx << 2), w, preds, n);
#endif
}
#endif

static void od_compute_dcts(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int pli,
  int bx, int by, int l, int xdec, int ydec) {
  int od;
  int d;
  int w;
  int bo;
  od_coeff *c;
  c = ctx->d[pli];
  w = enc->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  od = OD_BLOCK_SIZE4x4(enc->state.bsize,
   enc->state.bstride, bx << l, by << l);
  d = OD_MAXI(od, xdec);
  OD_ASSERT(d <= l);
  if (d == l) {
    d -= xdec;
    bo = (by << (OD_LOG_BSIZE0 + d))*w + (bx << (OD_LOG_BSIZE0 + d));
    (*enc->state.opt_vtbl.fdct_2d[d])(c + bo, w, ctx->c + bo, w);
#if defined(OD_DUMP_COEFFS)
    {
      int i;
      int j;
      int n;
      n = 1 << (OD_LOG_BSIZE0 + d);
      printf("%d ", n);
      for (j = 0; j < n; j++) for (i = 0; i < n; i++) {
        printf("%d ", c[bo + j*w + i]);
      }
      printf("\n");
    }
#endif
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_compute_dcts(enc, ctx, pli, bx + 0, by + 0, l, xdec, ydec);
    od_compute_dcts(enc, ctx, pli, bx + 1, by + 0, l, xdec, ydec);
    od_compute_dcts(enc, ctx, pli, bx + 0, by + 1, l, xdec, ydec);
    od_compute_dcts(enc, ctx, pli, bx + 1, by + 1, l, xdec, ydec);
    if (!OD_DISABLE_HAAR_DC && ctx->is_keyframe) {
      od_coeff x[4];
      int l2;
      l2 = l - xdec + 2;
      x[0] = c[(by << l2)*w + (bx << l2)];
      x[1] = c[(by << l2)*w + ((bx + 1) << l2)];
      x[2] = c[((by + 1) << l2)*w + (bx << l2)];
      x[3] = c[((by + 1) << l2)*w + ((bx + 1) << l2)];
      OD_HAAR_KERNEL(x[0], x[2], x[1], x[3]);
      c[(by << l2)*w + (bx << l2)] = x[0];
      c[(by << l2)*w + ((bx + 1) << l2)] = x[1];
      c[((by + 1) << l2)*w + (bx << l2)] = x[2];
      c[((by + 1) << l2)*w + ((bx + 1) << l2)] = x[3];
    }
  }
}

#if !OD_DISABLE_HAAR_DC
static void od_quantize_haar_dc(daala_enc_ctx *enc, od_mb_enc_ctx *ctx,
 int pli, int bx, int by, int l, int xdec, int ydec, od_coeff hgrad,
 od_coeff vgrad, int has_ur) {
  int od;
  int d;
  int w;
  int i;
  int dc_quant;
  od_coeff *c;
  c = ctx->d[pli];
  w = enc->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  od = OD_BLOCK_SIZE4x4(enc->state.bsize,
   enc->state.bstride, bx << l, by << l);
  d = OD_MAXI(od, xdec);
  OD_ASSERT(d <= l);
  if (enc->quantizer[pli] == 0) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, enc->quantizer[pli]*OD_DC_RES[pli] >> 4);
  }
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_DC_COEFF);
  if (l == 3) {
    int nhsb;
    int quant;
    int dc0;
    int l2;
    od_coeff sb_dc_pred;
    od_coeff sb_dc_curr;
    od_coeff *sb_dc_mem;
    /* Giving a small resolution boost to 32x32 because its reduced overlap
       means a larger synthesis magnitude. */
    if (enc->quantizer[pli] != 0 && d - xdec == 3) {
      dc_quant = OD_MAXI(1, enc->quantizer[pli]*12*OD_DC_RES[pli] >> 8);
    }
    nhsb = enc->state.nhsb;
    sb_dc_mem = enc->state.sb_dc_mem[pli];
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
    dc0 = c[(by << l2)*w + (bx << l2)] - sb_dc_pred;
    quant = OD_DIV_R0(dc0, dc_quant);
    generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli], abs(quant), -1,
     &enc->state.adapt.ex_sb_dc[pli], 2);
    if (quant) od_ec_enc_bits(&enc->ec, quant < 0, 1);
    sb_dc_curr = quant*dc_quant + sb_dc_pred;
    c[(by << l2)*w + (bx << l2)] = sb_dc_curr;
    sb_dc_mem[by*nhsb + bx] = sb_dc_curr;
    if (by > 0) vgrad = sb_dc_mem[(by - 1)*nhsb + bx] - sb_dc_curr;
    if (bx > 0) hgrad = sb_dc_mem[by*nhsb + bx - 1]- sb_dc_curr;
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
    x[1] -= hgrad/5;
    x[2] -= vgrad/5;
    for (i = 1; i < 4; i++) {
      int quant;
      quant = OD_DIV_R0(x[i], dc_quant);
      generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli], abs(quant), -1,
       &enc->state.adapt.ex_dc[pli][l][i-1], 2);
      if (quant) od_ec_enc_bits(&enc->ec, quant < 0, 1);
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
    od_quantize_haar_dc(enc, ctx, pli, bx + 0, by + 0, l, xdec, ydec, hgrad,
     vgrad, 0);
    od_quantize_haar_dc(enc, ctx, pli, bx + 1, by + 0, l, xdec, ydec, hgrad,
     vgrad, 0);
    od_quantize_haar_dc(enc, ctx, pli, bx + 0, by + 1, l, xdec, ydec, hgrad,
     vgrad, 0);
    od_quantize_haar_dc(enc, ctx, pli, bx + 1, by + 1, l, xdec, ydec, hgrad,
     vgrad, 0);
  }
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
}
#endif

#if !defined(OD_DUMP_COEFFS)
static void od_encode_recursive(daala_enc_ctx *enc, od_mb_enc_ctx *ctx,
 int pli, int bx, int by, int l, int xdec, int ydec) {
  int od;
  int d;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  od = OD_BLOCK_SIZE4x4(enc->state.bsize,
   enc->state.bstride, bx << l, by << l);
  d = OD_MAXI(od, xdec);
  OD_ASSERT(d <= l);
  if (d == l) {
    d -= xdec;
    /*Construct the luma predictors for chroma planes.*/
    if (ctx->l != NULL) {
      int w;
      int frame_width;
      OD_ASSERT(pli > 0);
      frame_width = enc->state.frame_width;
      w = frame_width >> xdec;
      od_resample_luma_coeffs(ctx->l + (by << (2 + d))*w + (bx << (2 + d)), w,
       ctx->d[0] + (by << (2 + l))*frame_width + (bx << (2 + l)),
       frame_width, xdec, ydec, d, od);
    }
    od_block_encode(enc, ctx, d, pli, bx, by);
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_encode_recursive(enc, ctx, pli, bx + 0, by + 0, l, xdec, ydec);
    od_encode_recursive(enc, ctx, pli, bx + 1, by + 0, l, xdec, ydec);
    od_encode_recursive(enc, ctx, pli, bx + 0, by + 1, l, xdec, ydec);
    od_encode_recursive(enc, ctx, pli, bx + 1, by + 1, l, xdec, ydec);
  }
}
#endif

static void od_encode_mv(daala_enc_ctx *enc, od_mv_grid_pt *mvg, int vx,
 int vy, int level, int mv_res, int width, int height) {
  generic_encoder *model;
  int pred[2];
  int ox;
  int oy;
  int id;
  int equal_mvs;
  equal_mvs = od_state_get_predictor(&enc->state, pred, vx, vy, level, mv_res);
  ox = (mvg->mv[0] >> mv_res) - pred[0];
  oy = (mvg->mv[1] >> mv_res) - pred[1];
  /*Interleave positive and negative values.*/
  model = &enc->state.adapt.mv_model;
  id = OD_MINI(abs(oy), 3)*4 + OD_MINI(abs(ox), 3);
  od_encode_cdf_adapt(&enc->ec, id, enc->state.adapt.mv_small_cdf[equal_mvs],
   16, enc->state.adapt.mv_small_increment);
  if (abs(ox) >= 3) {
    generic_encode(&enc->ec, model, abs(ox) - 3, width << (3 - mv_res),
     &enc->state.adapt.mv_ex[level], 6);
  }
  if (abs(oy) >= 3) {
    generic_encode(&enc->ec, model, abs(oy) - 3, height << (3 - mv_res),
     &enc->state.adapt.mv_ey[level], 6);
  }
  if (abs(ox)) od_ec_enc_bits(&enc->ec, ox < 0, 1);
  if (abs(oy)) od_ec_enc_bits(&enc->ec, oy < 0, 1);
}

static void od_img_copy_pad(od_state *state, od_img *img) {
  int pli;
  int nplanes;
  nplanes = img->nplanes;
  /* Copy and pad the image. */
  for (pli = 0; pli < nplanes; pli++) {
    od_img_plane plane;
    int plane_width;
    int plane_height;
    int xdec;
    int ydec;
    *&plane = *(img->planes + pli);
    xdec = plane.xdec;
    ydec = plane.ydec;
    plane_width = ((state->info.pic_width + (1 << xdec) - 1) >> xdec);
    plane_height = ((state->info.pic_height + (1 << ydec) - 1) >> ydec);
    od_img_plane_copy_pad8(&state->io_imgs[OD_FRAME_INPUT].planes[pli],
     state->frame_width >> xdec, state->frame_height >> ydec,
     &plane, plane_width, plane_height);
    od_img_plane_edge_ext8(&state->io_imgs[OD_FRAME_INPUT].planes[pli],
     state->frame_width >> xdec, state->frame_height >> ydec,
     OD_UMV_PADDING >> xdec, OD_UMV_PADDING >> ydec);
  }
}

#if defined(OD_DUMP_IMAGES)
static void od_img_dump_padded(od_state *state) {
  daala_info *info;
  od_img img;
  int nplanes;
  int pli;
  info = &state->info;
  nplanes = info->nplanes;
  /*Modify the image offsets to include the padding.*/
  *&img = *(state->io_imgs+OD_FRAME_INPUT);
  for (pli = 0; pli < nplanes; pli++) {
    img.planes[pli].data -= (OD_UMV_PADDING>>info->plane_info[pli].xdec)
        +img.planes[pli].ystride*(OD_UMV_PADDING>>info->plane_info[pli].ydec);
  }
  img.width += OD_UMV_PADDING<<1;
  img.height += OD_UMV_PADDING<<1;
  od_state_dump_img(state, &img, "pad");
}
#endif

static void od_predict_frame(daala_enc_ctx *enc) {
  int nplanes;
  int pli;
  int frame_width;
  int frame_height;
  nplanes = enc->state.info.nplanes;
  frame_width = enc->state.frame_width;
  frame_height = enc->state.frame_height;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  enc->state.ani_iter = 0;
#endif
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "Predicting frame %i:",
   (int)daala_granule_basetime(enc, enc->state.cur_time)));
  /*2851196 ~= sqrt(ln(2)/6) in Q23.
   The lower bound of 56 is there because we do not yet consider PVQ noref
    flags during the motion search, so we waste far too many bits trying to
    predict unpredictable areas when lamba is too small.
   Hopefully when we fix that, we can remove the limit.*/
  od_mv_est(enc->mvest, OD_FRAME_PREV,
   OD_MAXI((2851196 + (((1 << OD_COEFF_SHIFT) - 1) >> 1) >> OD_COEFF_SHIFT)*
   enc->quantizer[0] >> (23 - OD_LAMBDA_SCALE), 56));
  od_state_mc_predict(&enc->state, OD_FRAME_PREV);
  /*Do edge extension here because the block-size analysis needs to read
    outside the frame, but otherwise isn't read from.*/
  for (pli = 0; pli < nplanes; pli++) {
    od_img_plane plane;
    *&plane = *(enc->state.io_imgs[OD_FRAME_REC].planes + pli);
    od_img_plane_edge_ext8(&plane, frame_width >> plane.xdec,
     frame_height >> plane.ydec, OD_UMV_PADDING >> plane.xdec,
     OD_UMV_PADDING >> plane.ydec);
  }
#if defined(OD_DUMP_IMAGES)
  /*Dump reconstructed frame.*/
  /*od_state_dump_img(&enc->state,enc->state.io_imgs + OD_FRAME_REC,"rec");*/
  od_state_fill_vis(&enc->state);
  od_state_dump_img(&enc->state, &enc->state.vis_img, "vis");
#endif
}

#if defined(OD_DUMP_IMAGES)
static void od_bsize_dump_img(od_state *state, int nvsb, int nhsb) {
  od_img copy;
  int size;
  unsigned char *img;
  int i;
  int j;
  int stride;
  copy = state->io_imgs[OD_FRAME_INPUT];
  stride = copy.planes[0].ystride;
  size = stride*copy.height;
  /* Pad by one pixel to make code simpler. */
  img = copy.planes[0].data = _ogg_malloc(size + stride + 1);
  OD_COPY(img, state->io_imgs[OD_FRAME_INPUT].planes[0].data, size);
  /* Draw black borders around the blocks. */
  for(i = 0; i < 4*nvsb; i++) {
    for(j = 0; j < 4*nhsb; j++) {
      if ((i & 3) == 0 && (j & 3) == 0) {
        int k;
        for(k = 0; k < 32; k++) img[i*stride*8 + j*8 + k] = 0;
        for(k = 0; k < 32; k++) img[(8*i + k)*stride + j*8] = 0;
      }
      if ((i & 1) == 0 && (j & 1) == 0 && state->bsize[i*state->bstride + j]
       == OD_BLOCK_16X16) {
        int k;
        for(k = 0; k < 16; k++) img[i*stride*8 + j*8 + k] = 0;
        for(k = 0; k < 16; k++) img[(8*i + k)*stride + j*8] = 0;
      }
      if (state->bsize[i*state->bstride + j] <= OD_BLOCK_8X8) {
        int k;
        for(k = 0; k < 8; k++) img[i*stride*8 + j*8 + k] = 0;
        for(k = 0; k < 8; k++) img[(8*i + k)*stride + j*8] = 0;
        if (state->bsize[i*state->bstride + j] == OD_BLOCK_4X4) {
          img[(8*i + 4)*stride + j*8 + 3] = 0;
          img[(8*i + 4)*stride + j*8 + 4] = 0;
          img[(8*i + 4)*stride + j*8 + 5] = 0;
          img[(8*i + 3)*stride + j*8 + 4] = 0;
          img[(8*i + 5)*stride + j*8 + 4] = 0;
        }
      }
    }
  }
  od_state_dump_img(state, &copy, "bsize");
  _ogg_free(img);
}
#endif

static void od_split_superblocks(daala_enc_ctx *enc, int is_keyframe) {
  int nhsb;
  int nvsb;
  int i;
  int j;
  int k;
  int m;
  od_state *state;
  state = &enc->state;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  od_state_init_border(state);
  /* Allocate a blockSizeComp for scratch space and then calculate the block
     sizes eventually store them in bsize. */
  od_log_matrix_uchar(OD_LOG_GENERIC, OD_LOG_INFO, "bimg ",
   state->io_imgs[OD_FRAME_INPUT].planes[0].data -
   16*state->io_imgs[OD_FRAME_INPUT].planes[0].ystride - 16,
   state->io_imgs[OD_FRAME_INPUT].planes[0].ystride, (nvsb + 1)*32);
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_BLOCK_SIZE);
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_FRAME);
  for (i = 0; i < nvsb; i++) {
    unsigned char *bimg;
    unsigned char *rimg;
    int istride;
    int rstride;
    int bstride;
    bstride = state->bstride;
    istride = state->io_imgs[OD_FRAME_INPUT].planes[0].ystride;
    rstride = is_keyframe ? 0 :
     state->io_imgs[OD_FRAME_REC].planes[0].ystride;
    bimg = state->io_imgs[OD_FRAME_INPUT].planes[0].data + i*istride*32;
    rimg = state->io_imgs[OD_FRAME_REC].planes[0].data + i*rstride*32;
    for (j = 0; j < nhsb; j++) {
      int bsize[4][4];
      unsigned char *state_bsize;
      state_bsize = &state->bsize[i*4*state->bstride + j*4];
      od_split_superblock(enc->bs, bimg + j*32, istride,
       is_keyframe ? NULL : rimg + j*32, rstride, bsize, enc->quantizer[0]);
      /* Grab the 4x4 information returned from `od_split_superblock` in bsize
         and store it in the od_state bsize. */
      for (k = 0; k < 4; k++) {
        for (m = 0; m < 4; m++) {
          if (OD_LIMIT_BSIZE_MIN != OD_LIMIT_BSIZE_MAX) {
            state_bsize[k*bstride + m] =
             OD_MAXI(OD_MINI(bsize[k][m], OD_LIMIT_BSIZE_MAX),
             OD_LIMIT_BSIZE_MIN);
          }
          else {
            state_bsize[k*bstride + m] = OD_LIMIT_BSIZE_MIN;
          }
        }
      }
      if (OD_LIMIT_BSIZE_MIN != OD_LIMIT_BSIZE_MAX) {
        od_block_size_encode(&enc->ec, &enc->state.adapt, &state_bsize[0],
         bstride);
      }
    }
  }
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  od_log_matrix_uchar(OD_LOG_GENERIC, OD_LOG_INFO, "bsize ", state->bsize,
   state->bstride, (nvsb + 1)*4);
  for (i = 0; i < nvsb*4; i++) {
    for (j = 0; j < nhsb*4; j++) {
      OD_LOG_PARTIAL((OD_LOG_GENERIC, OD_LOG_INFO, "%d ",
       state->bsize[i*state->bstride + j]));
    }
    OD_LOG_PARTIAL((OD_LOG_GENERIC, OD_LOG_INFO, "\n"));
  }
#if defined(OD_DUMP_IMAGES)
  od_bsize_dump_img(&enc->state, nvsb, nhsb);
#endif
}

static void od_encode_mvs(daala_enc_ctx *enc) {
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  od_img *mvimg;
  int width;
  int height;
  int mv_res;
  od_mv_grid_pt *mvp;
  od_mv_grid_pt *other;
  od_mv_grid_pt **grid;
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_MOTION_VECTORS);
  nhmvbs = (enc->state.nhmbs + 1) << 2;
  nvmvbs = (enc->state.nvmbs + 1) << 2;
  mvimg = enc->state.io_imgs + OD_FRAME_REC;
  mv_res = enc->state.mv_res;
  OD_ASSERT(0 <= mv_res && mv_res < 3);
  od_ec_enc_uint(&enc->ec, mv_res, 3);
  width = (mvimg->width + 32) << (3 - mv_res);
  height = (mvimg->height + 32) << (3 - mv_res);
  grid = enc->state.mv_grid;
  /*Code the motion vectors and flags. At each level, the MVs are zero
    outside of the frame, so don't code them.*/
  /*Level 0.*/
  for (vy = 4; vy < nvmvbs; vy += 4) {
    for (vx = 4; vx < nhmvbs; vx += 4) {
      mvp = &grid[vy][vx];
      od_encode_mv(enc, mvp, vx, vy, 0, mv_res, width, height);
    }
  }
  /*od_ec_acct_add_label(&enc->ec.acct, "mvf-l1");
    od_ec_acct_add_label(&enc->ec.acct, "mvf-l2");
    od_ec_acct_add_label(&enc->ec.acct, "mvf-l3");
    od_ec_acct_add_label(&enc->ec.acct, "mvf-l4");*/
  /*Level 1.*/
  for (vy = 2; vy <= nvmvbs; vy += 4) {
    for (vx = 2; vx <= nhmvbs; vx += 4) {
      int p_invalid;
      p_invalid = od_mv_level1_probz(grid, vx, vy);
      mvp = &(grid[vy][vx]);
      /*od_ec_acct_record(&enc->ec.acct, "mvf-l1", mvp->valid, 2,
         od_mv_level1_ctx(grid, vx, vy));*/
      if (p_invalid >= 16384) {
        od_ec_encode_bool_q15(&enc->ec, mvp->valid, p_invalid);
      }
      else {
        od_ec_encode_bool_q15(&enc->ec, !mvp->valid, 32768 - p_invalid);
      }
      if (mvp->valid) {
        od_encode_mv(enc, mvp, vx, vy, 1, mv_res, width, height);
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
        /*od_ec_acct_record(&enc->ec.acct, "mvf-l2", mvp->valid, 2,
           od_mv_level2_ctx(grid, vx, vy));*/
        if (p_invalid >= 16384) {
          od_ec_encode_bool_q15(&enc->ec, mvp->valid, p_invalid);
        }
        else {
          od_ec_encode_bool_q15(&enc->ec, !mvp->valid, 32768 - p_invalid);
        }
        if (mvp->valid && vx >= 2 && vy >= 2 && vx <= nhmvbs - 2 &&
         vy <= nvmvbs - 2) {
          od_encode_mv(enc, mvp, vx, vy, 2, mv_res, width, height);
        }
      }
    }
  }
  /*Level 3.*/
  /*Level 3 motion vector flags outside the frame are specially coded
    since more information is known. On the grid edge, an L2 MV will only be
    valid if a L3 MV is needed outside of the frame. In the middle of the
    edge, this implies a tristate of the two possible child L3 MVs; they
    can't both be invalid. At the corner, one of the child L3 vectors will
    never appear, so an L2 MV directly implies the remaining L3 child.*/
  for (vy = 1; vy <= nvmvbs; vy += 2) {
    for (vx = 1; vx <= nhmvbs; vx += 2) {
      mvp = &grid[vy][vx];
      if (vy < 2 || vy > nvmvbs - 2) {
        if ((vx == 3 && grid[vy == 1 ? vy - 1 : vy + 1][vx - 1].valid)
         || (vx == nhmvbs - 3
         && grid[vy == 1 ? vy - 1 : vy + 1][vx + 1].valid)) {
          other = &grid[vy][vx == 3 ? vx - 2 : vx + 2];
          /*MVs are valid but will be zero.*/
          OD_ASSERT(mvp->valid && !mvp->mv[0] && !mvp->mv[1]
           && !other->valid);
        }
        else if (vx > 3 && vx < nhmvbs - 3) {
          other = &grid[vy][vx + 2];
          if (!(vx & 2) && grid[vy == 1 ? vy - 1 : vy + 1][vx + 1].valid) {
            /*0 = both valid, 1 = only this one, 2 = other one valid*/
            int s;
            s = mvp->valid && other->valid ? 0 : mvp->valid
             + (other->valid << 1);
            od_ec_encode_cdf_q15(&enc->ec, s, OD_UNIFORM_CDF_Q15(3), 3);
            /*MVs are valid but will be zero.*/
            OD_ASSERT((mvp->valid && !mvp->mv[0] && !mvp->mv[1])
             || (other->valid && !other->mv[0] && !other->mv[1]));
          }
          else if (!(vx & 2)) {
            OD_ASSERT(!mvp->valid && !other->valid);
          }
        }
        else {
          OD_ASSERT(!mvp->valid);
        }
      }
      else if (vx < 2 || vx > nhmvbs - 2) {
        od_mv_grid_pt *other;
        if ((vy == 3 && grid[vy - 1][vx == 1 ? vx - 1 : vx + 1].valid)
         || (vy == nvmvbs - 3
          && grid[vy + 1][vx == 1 ? vx - 1 : vx + 1].valid)) {
          other = &grid[vy == 3 ? vy - 2 : vy + 2][vx];
          /*MVs are valid but will be zero.*/
          OD_ASSERT(mvp->valid && !mvp->mv[0] && !mvp->mv[1]
           && !other->valid);
        }
        else if (vy > 3 && vy < nvmvbs - 3) {
          other = &grid[vy + 2][vx];
          if (!(vy & 2) && grid[vy + 1][vx == 1 ? vx - 1 : vx + 1].valid) {
            int s;
            s = mvp->valid && other->valid ? 0 : mvp->valid
             + (other->valid << 1);
            od_ec_encode_cdf_q15(&enc->ec, s, OD_UNIFORM_CDF_Q15(3), 3);
            /*MVs are valid but will be zero.*/
            OD_ASSERT((mvp->valid && !mvp->mv[0] && !mvp->mv[1])
             || (other->valid && !other->mv[0] && !other->mv[1]));
          }
          else if (!(vy & 2)) {
            OD_ASSERT(!mvp->valid && !other->valid);
          }
        }
        else {
          OD_ASSERT(!mvp->valid);
        }
      }
      else if (grid[vy - 1][vx - 1].valid && grid[vy - 1][vx + 1].valid
       && grid[vy + 1][vx + 1].valid && grid[vy + 1][vx - 1].valid) {
        int p_invalid;
        p_invalid = od_mv_level3_probz(grid, vx, vy);
        /*od_ec_acct_record(&enc->ec.acct, "mvf-l3", mvp->valid, 2,
           od_mv_level3_ctx(grid, vx, vy));*/
        if (p_invalid >= 16384) {
          od_ec_encode_bool_q15(&enc->ec, mvp->valid, p_invalid);
        }
        else {
          od_ec_encode_bool_q15(&enc->ec, !mvp->valid, 32768 - p_invalid);
        }
        if (mvp->valid) {
          od_encode_mv(enc, mvp, vx, vy, 3, mv_res, width, height);
        }
      }
      else {
        OD_ASSERT(!mvp->valid);
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
        /*od_ec_acct_record(&enc->ec.acct, "mvf-l4", mvp->valid, 2,
           od_mv_level4_ctx(grid, vx, vy));*/
        if (p_invalid >= 16384) {
          od_ec_encode_bool_q15(&enc->ec, mvp->valid, p_invalid);
        }
        else {
          od_ec_encode_bool_q15(&enc->ec, !mvp->valid, 32768 - p_invalid);
        }
        if (mvp->valid) {
          od_encode_mv(enc, mvp, vx, vy, 4, mv_res, width, height);
        }
        else {
          OD_ASSERT(!mvp->valid);
        }
      }
    }
  }
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
}

static void od_encode_residual(daala_enc_ctx *enc, od_mb_enc_ctx *mbctx) {
  int xdec;
  int ydec;
  int sby;
  int sbx;
  int h;
  int w;
  int y;
  int x;
  int pli;
  int nplanes;
  int frame_width;
  int frame_height;
  int nhsb;
  int nvsb;
  od_state *state = &enc->state;
  nplanes = state->info.nplanes;
  frame_width = state->frame_width;
  frame_height = state->frame_height;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  for (pli = 0; pli < nplanes; pli++) {
    xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_FRAME);
    OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_FRAME);
    /* TODO: We shouldn't be encoding the full, linear quantizer range. */
    od_ec_enc_uint(&enc->ec, enc->quantizer[pli], 512<<OD_COEFF_SHIFT);
    /*If the quantizer is zero (lossless), force scalar.*/
    OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
    OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_UNKNOWN);
  }
  for (pli = 0; pli < nplanes; pli++) {
    xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    /*Collect the image data needed for this plane.*/
    {
      unsigned char *data;
      unsigned char *mdata;
      int ystride;
      int coeff_shift;
      coeff_shift = enc->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
      data = state->io_imgs[OD_FRAME_INPUT].planes[pli].data;
      mdata = state->io_imgs[OD_FRAME_REC].planes[pli].data;
      ystride = state->io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
      for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
          state->ctmp[pli][y*w + x] = (data[ystride*y + x] - 128) <<
           coeff_shift;
          if (!mbctx->is_keyframe) {
            state->mctmp[pli][y*w + x] = (mdata[ystride*y + x] - 128)
             << coeff_shift;
          }
        }
      }
    }
    /*Apply the prefilter across the entire image.*/
    for (sby = 0; sby < nvsb; sby++) {
      for (sbx = 0; sbx < nhsb; sbx++) {
        od_apply_prefilter(state->ctmp[pli], w, sbx, sby, 3,
         state->bsize, state->bstride, xdec, ydec,
         (sbx > 0 ? OD_LEFT_EDGE : 0) |
         (sby < nvsb - 1 ? OD_BOTTOM_EDGE : 0));
        if (!mbctx->is_keyframe) {
          od_apply_prefilter(state->mctmp[pli], w, sbx, sby, 3, state->bsize,
           state->bstride, xdec, ydec, (sbx > 0 ? OD_LEFT_EDGE : 0) |
           (sby < nvsb - 1 ? OD_BOTTOM_EDGE : 0));
        }
      }
    }
  }
  for (sby = 0; sby < nvsb; sby++) {
    for (sbx = 0; sbx < nhsb; sbx++) {
      for (pli = 0; pli < nplanes; pli++) {
        OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_LUMA + pli);
        mbctx->c = state->ctmp[pli];
        mbctx->d = state->dtmp;
        mbctx->mc = state->mctmp[pli];
        mbctx->md = state->mdtmp[pli];
        mbctx->l = state->lbuf[pli];
        xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
        ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
        mbctx->nk = mbctx->k_total = mbctx->sum_ex_total_q8 = 0;
        mbctx->ncount = mbctx->count_total_q8 = mbctx->count_ex_total_q8 = 0;
        /*Need to update this to decay based on superblocks width.*/
        od_compute_dcts(enc, mbctx, pli, sbx, sby, 3, xdec, ydec);
        if (!OD_DISABLE_HAAR_DC && mbctx->is_keyframe) {
          od_quantize_haar_dc(enc, mbctx, pli, sbx, sby, 3, xdec, ydec, 0,
           0, sby > 0 && sbx < nhsb - 1);
        }
#if !defined(OD_DUMP_COEFFS)
        od_encode_recursive(enc, mbctx, pli, sbx, sby, 3, xdec, ydec);
#endif
      }
        OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_UNKNOWN);
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
         state->bstride, xdec, ydec,
         (sby > 0 ? OD_TOP_EDGE : 0) | (sbx < nhsb - 1 ? OD_RIGHT_EDGE : 0));
      }
    }
    {
      unsigned char *data;
      int ystride;
      int coeff_shift;
      coeff_shift = enc->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
      data = state->io_imgs[OD_FRAME_REC].planes[pli].data;
      ystride = state->io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
      for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
          data[ystride*y + x] = OD_CLAMP255(((state->ctmp[pli][y*w + x]
           + (1 << coeff_shift >> 1)) >> coeff_shift) + 128);
        }
      }
    }
  }
}

#if defined(OD_LOGGING_ENABLED)
static void od_dump_frame_metrics(od_state *state) {
  int pli;
  int nplanes;
  int frame_width;
  int frame_height;
  nplanes = state->info.nplanes;
  frame_width = state->frame_width;
  frame_height = state->frame_height;
  for (pli = 0; pli < nplanes; pli++) {
    unsigned char *data;
    ogg_int64_t enc_sqerr;
    ogg_uint32_t npixels;
    int ystride;
    int xdec;
    int ydec;
    int w;
    int h;
    int x;
    int y;
    enc_sqerr = 0;
    data = state->io_imgs[OD_FRAME_INPUT].planes[pli].data;
    ystride = state->io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
    xdec = state->io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = state->io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    npixels = w*h;
    for (y = 0; y < h; y++) {
      unsigned char *rec_row;
      unsigned char *inp_row;
      rec_row = state->io_imgs[OD_FRAME_REC].planes[pli].data +
       state->io_imgs[OD_FRAME_REC].planes[pli].ystride*y;
      inp_row = data + ystride*y;
      for (x = 0; x < w; x++) {
        int inp_val;
        int diff;
        inp_val = inp_row[x];
        diff = inp_val - rec_row[x];
        enc_sqerr += diff*diff;
      }
    }
    OD_LOG((OD_LOG_ENCODER, OD_LOG_DEBUG,
     "Encoded Plane %i, Squared Error: %12lli  Pixels: %6u  PSNR:  %5.4f",
     pli, (long long)enc_sqerr, npixels,
     10*log10(255*255.0*npixels/enc_sqerr)));
  }
}
#endif

static void od_interp_qm(unsigned char *out, int q, const od_qm_entry *entry1,
  const od_qm_entry *entry2) {
  int i;
  if (entry2 == NULL || entry2->qm_q4 == NULL
   || q < entry1->interp_q << OD_COEFF_SHIFT) {
    /* Use entry1. */
    for (i = 0; i < OD_QM_SIZE; i++) {
      out[i] = OD_MINI(255, entry1->qm_q4[i]*entry1->scale_q8 >> 8);
    }
  }
  else if (entry1 == NULL || entry1->qm_q4 == NULL
   || q > entry2->interp_q << OD_COEFF_SHIFT) {
    /* Use entry2. */
    for (i = 0; i < OD_QM_SIZE; i++) {
      out[i] = OD_MINI(255, entry2->qm_q4[i]*entry2->scale_q8 >> 8);
    }
  }
  else {
    /* Interpolate between entry1 and entry2. The interpolation is linear
       in terms of log(q) vs log(m*scale). Considering that we're ultimately
       multiplying the result it makes sense, but we haven't tried other
       interpolation methods. */
    double x;
    const unsigned char *m1;
    const unsigned char *m2;
    int q1;
    int q2;
    m1 = entry1->qm_q4;
    m2 = entry2->qm_q4;
    q1 = entry1->interp_q << OD_COEFF_SHIFT;
    q2 = entry2->interp_q << OD_COEFF_SHIFT;
    x = (log(q)-log(q1))/(log(q2)-log(q1));
    for (i = 0; i < OD_QM_SIZE; i++) {
      out[i] = OD_MINI(255, (int)floor(.5 + (1./256)*exp(
       x*log(m2[i]*entry2->scale_q8) + (1 - x)*log(m1[i]*entry1->scale_q8))));
    }
  }
}

int daala_encode_img_in(daala_enc_ctx *enc, od_img *img, int duration) {
  int refi;
  int nplanes;
  int pli;
  int frame_width;
  int frame_height;
  int pic_width;
  int pic_height;
  od_mb_enc_ctx mbctx;
#if defined(OD_ACCOUNTING)
  od_acct_reset(&enc->acct);
#endif
#if defined(OD_EC_ACCOUNTING)
  od_ec_acct_reset(&enc->ec.acct);
#endif
  if (enc == NULL || img == NULL) return OD_EFAULT;
  if (enc->packet_state == OD_PACKET_DONE) return OD_EINVAL;
  /*Check the input image dimensions to make sure they're compatible with the
     declared video size.*/
  nplanes = enc->state.info.nplanes;
  if (img->nplanes != nplanes) return OD_EINVAL;
  for (pli = 0; pli < nplanes; pli++) {
    if (img->planes[pli].xdec != enc->state.info.plane_info[pli].xdec
     || img->planes[pli].ydec != enc->state.info.plane_info[pli].ydec) {
      return OD_EINVAL;
    }
  }
  frame_width = enc->state.frame_width;
  frame_height = enc->state.frame_height;
  pic_width = enc->state.info.pic_width;
  pic_height = enc->state.info.pic_height;
  if (img->width != frame_width || img->height != frame_height) {
    /*The buffer does not match the frame size.
      Check to see if it matches the picture size.*/
    if (img->width != pic_width || img->height != pic_height) {
      /*It doesn't; we don't know how to handle it yet.*/
      return OD_EINVAL;
    }
  }
  od_img_copy_pad(&enc->state, img);

#if defined(OD_DUMP_IMAGES)
  if (od_logging_active(OD_LOG_GENERIC, OD_LOG_DEBUG)) {
    od_img_dump_padded(&enc->state);
  }
#endif
  /* Check if the frame should be a keyframe. */
  mbctx.is_keyframe = (enc->state.cur_time %
   (enc->state.info.keyframe_rate) == 0) ? 1 : 0;
  /*Update the buffer state.*/
  if (enc->state.ref_imgi[OD_FRAME_SELF] >= 0) {
    enc->state.ref_imgi[OD_FRAME_PREV] =
     enc->state.ref_imgi[OD_FRAME_SELF];
    /*TODO: Update golden frame.*/
    if (enc->state.ref_imgi[OD_FRAME_GOLD] < 0) {
      enc->state.ref_imgi[OD_FRAME_GOLD] =
       enc->state.ref_imgi[OD_FRAME_SELF];
      /*TODO: Mark keyframe timebase.*/
    }
  }
  /*Select a free buffer to use for this reference frame.*/
  for (refi = 0; refi == enc->state.ref_imgi[OD_FRAME_GOLD]
   || refi == enc->state.ref_imgi[OD_FRAME_PREV]
   || refi == enc->state.ref_imgi[OD_FRAME_NEXT]; refi++);
  enc->state.ref_imgi[OD_FRAME_SELF] = refi;
  /*We must be a keyframe if we don't have a reference.*/
  mbctx.is_keyframe |= !(enc->state.ref_imgi[OD_FRAME_PREV] >= 0);
  /*Initialize the entropy coder.*/
  od_ec_enc_reset(&enc->ec);
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_FRAME);
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_FRAME);
  /*Write a bit to mark this as a data packet.*/
  od_ec_encode_bool_q15(&enc->ec, 0, 16384);
  /*Code the keyframe bit.*/
  od_ec_encode_bool_q15(&enc->ec, mbctx.is_keyframe, 16384);
  for (pli = 0; pli < nplanes; pli++) {
    enc->quantizer[pli] = od_quantizer_from_quality(enc->quality[pli]);
  }
  if (mbctx.is_keyframe) {
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      int q;
      q = enc->quantizer[pli];
      if (q <= OD_DEFAULT_QMS[0][pli].interp_q << OD_COEFF_SHIFT) {
        od_interp_qm(&enc->state.pvq_qm_q4[pli][0], q, &OD_DEFAULT_QMS[0][pli],
         NULL);
      }
      else {
        i = 0;
        while (OD_DEFAULT_QMS[i + 1][pli].qm_q4 != NULL &&
         q > OD_DEFAULT_QMS[i + 1][pli].interp_q << OD_COEFF_SHIFT) {
          i++;
        }
        od_interp_qm(&enc->state.pvq_qm_q4[pli][0], q,
         &OD_DEFAULT_QMS[i][pli], &OD_DEFAULT_QMS[i + 1][pli]);
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      for (i = 0; i < OD_QM_SIZE; i++) {
        od_ec_enc_bits(&enc->ec, enc->state.pvq_qm_q4[pli][i], 8);
      }
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    /* At low rate, boost the keyframe quality by multiplying the quantizer
       by 29/32 (~0.9). */
    if (mbctx.is_keyframe && enc->quantizer[pli] > 20 << OD_COEFF_SHIFT) {
      enc->quantizer[pli] = (16+29*enc->quantizer[pli]) >> 5;
    }
  }
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "is_keyframe=%d", mbctx.is_keyframe));
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  OD_ENC_ACCT_UPDATE(enc, OD_ACCT_CAT_PLANE, OD_ACCT_TECH_UNKNOWN);
  /*TODO: Increment frame count.*/
  od_adapt_ctx_reset(&enc->state.adapt, mbctx.is_keyframe);
  if (!mbctx.is_keyframe) {
    od_predict_frame(enc);
    od_split_superblocks(enc, 0);
    od_encode_mvs(enc);
  }
  else {
    od_split_superblocks(enc, 1);
  }
  od_encode_residual(enc, &mbctx);
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  /*Dump YUV*/
  od_state_dump_yuv(&enc->state, enc->state.io_imgs + OD_FRAME_REC, "out");
#endif
#if defined(OD_LOGGING_ENABLED)
  od_dump_frame_metrics(&enc->state);
#endif
  enc->packet_state = OD_PACKET_READY;
  od_state_upsample8(&enc->state,
   enc->state.ref_imgs + enc->state.ref_imgi[OD_FRAME_SELF],
   enc->state.io_imgs + OD_FRAME_REC);
#if defined(OD_DUMP_IMAGES)
  /*Dump reference frame.*/
  /*od_state_dump_img(&enc->state,
   enc->state.ref_img + enc->state.ref_imigi[OD_FRAME_SELF], "ref");*/
#endif
#if defined(OD_ACCOUNTING)
  OD_ASSERT(enc->acct.last_frac_bits == od_ec_enc_tell_frac(&enc->ec));
  od_acct_write(&enc->acct, enc->state.cur_time);
#endif
#if defined(OD_EC_ACCOUNTING)
  od_ec_acct_write(&enc->ec.acct);
#endif
  if (enc->state.info.frame_duration == 0) enc->state.cur_time += duration;
  else enc->state.cur_time += enc->state.info.frame_duration;
  return 0;
}

#if defined(OD_ENCODER_CHECK)
static void daala_encoder_check(daala_enc_ctx *ctx, od_img *img,
 ogg_packet *op) {
  int pli;
  od_img dec_img;
  OD_ASSERT(ctx->dec);

  if (daala_decode_packet_in(ctx->dec, &dec_img, op) < 0) {
    fprintf(stderr,"decode failed!\n");
    return;
  }

  OD_ASSERT(img->nplanes == dec_img.nplanes);
  for (pli = 0; pli < img->nplanes; pli++) {
    int plane_width;
    int plane_height;
    int xdec;
    int ydec;
    int i;
    OD_ASSERT(img->planes[pli].xdec == dec_img.planes[pli].xdec);
    OD_ASSERT(img->planes[pli].ydec == dec_img.planes[pli].ydec);
    OD_ASSERT(img->planes[pli].ystride == dec_img.planes[pli].ystride);

    xdec = dec_img.planes[pli].xdec;
    ydec = dec_img.planes[pli].ydec;
    plane_width = ctx->dec->state.frame_width >> xdec;
    plane_height = ctx->dec->state.frame_height >> ydec;
    for (i = 0; i < plane_height; i++) {
      if (memcmp(img->planes[pli].data + img->planes[pli].ystride * i,
       dec_img.planes[pli].data + dec_img.planes[pli].ystride * i,
       plane_width)) {
        fprintf(stderr,"pixel mismatch in row %d\n", i);
      }
    }
  }
}
#endif

int daala_encode_packet_out(daala_enc_ctx *enc, int last, ogg_packet *op) {
  ogg_uint32_t nbytes;
  if (enc == NULL || op == NULL) return OD_EFAULT;
  else if (enc->packet_state <= 0 || enc->packet_state == OD_PACKET_DONE) {
    return 0;
  }
  op->packet = od_ec_enc_done(&enc->ec, &nbytes);
  op->bytes = nbytes;
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "Output Bytes: %ld (%ld Kbits)", op->bytes, op->bytes * 8 / 1024));
  op->b_o_s = 0;
  op->e_o_s = last;
  op->packetno = 0;
  op->granulepos = enc->state.cur_time;
  if (last) enc->packet_state = OD_PACKET_DONE;
  else enc->packet_state = OD_PACKET_EMPTY;

#if defined(OD_ENCODER_CHECK)
  /*Compare reconstructed frame against decoded frame.*/
  daala_encoder_check(enc, enc->state.io_imgs + OD_FRAME_REC, op);
#endif

  return 1;
}
