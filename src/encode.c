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

static double mode_bits = 0;
static double mode_count = 0;

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
    enc->quantizer[i] = od_quantizer_from_quality(10);
  }
  enc->mvest = od_mv_est_alloc(enc);
  enc->params.mv_level_min = 0;
  enc->params.mv_level_max = 4;
#if defined(OD_ACCOUNTING)
  od_acct_init(&enc->acct);
#endif
  enc->bs = (BlockSizeComp *)_ogg_malloc(sizeof(*enc->bs));
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
      OD_ASSERT(buf_sz == sizeof(*enc->quantizer));
      for (i = 0; i < OD_NPLANES_MAX; i++){
        enc->quantizer[i] = od_quantizer_from_quality(*(int *)buf);
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
  signed char *modes[OD_NPLANES_MAX];
  od_coeff *c;
  od_coeff **d;
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
typedef struct od_mb_enc_ctx od_mb_enc_ctx;

static void od_encode_compute_pred(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, od_coeff *pred,
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
  xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  ydec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
  frame_width = enc->state.frame_width;
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
        ogg_uint32_t mode_dist[OD_INTRA_NMODES];
        int m_l;
        int m_ul;
        int m_u;
        int mode;
        od_coeff *coeffs[4];
        int strides[4];
        /*Search predictors from the surrounding blocks.*/
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
        od_intra_pred_cdf(mode_cdf, enc->state.adapt.mode_probs[pli],
         OD_INTRA_NMODES, m_l, m_ul, m_u);
        (*OD_INTRA_DIST[ln])(mode_dist, d + (by << 2)*w + (bx << 2), w,
         coeffs, strides);
        /*Lambda = 1*/
#if OD_DISABLE_INTRA
        mode = 0;
#else
        /* Make lambda proportional to quantization step size, with exact
           factor based on quick experiments with subset1 (can be improved). */
        mode = od_intra_pred_search(mode_cdf, mode_dist, OD_INTRA_NMODES,
         OD_MINI(32767, enc->quantizer[pli] << 4));
#endif
        (*OD_INTRA_GET[ln])(pred, coeffs, strides, mode);
#if OD_DISABLE_INTRA
        OD_CLEAR(pred+1, n2-1);
#endif
        OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
         OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_INTRA_MODE);
#if !OD_DISABLE_INTRA
        od_ec_encode_cdf_unscaled(&enc->ec, mode, mode_cdf, OD_INTRA_NMODES);
#endif
        OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
         OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
        mode_bits -= M_LOG2E*log(
         (mode_cdf[mode] - (mode == 0 ? 0 : mode_cdf[mode - 1]))/
         (float)mode_cdf[OD_INTRA_NMODES - 1]);
        mode_count++;
        for (y = 0; y < (1 << ln); y++) {
          for (x = 0; x < (1 << ln); x++) {
            modes[(by + y)*(w >> 2) + bx + x] = mode;
          }
        }
        od_intra_pred_update(enc->state.adapt.mode_probs[pli], OD_INTRA_NMODES,
         mode, m_l, m_ul, m_u);
      }
      else {
        int mode;
        mode = modes[(by << ydec)*(frame_width >> 2) + (bx << xdec)];
        od_chroma_pred(pred, d, l, w, bx, by, ln, xdec, ydec,
         enc->state.bsize, enc->state.bstride,
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
        nsize = OD_BLOCK_SIZE4x4(enc->state.bsize, enc->state.bstride,
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
        nsize = OD_BLOCK_SIZE4x4(enc->state.bsize, enc->state.bstride,
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

static void od_single_band_scalar_quant(daala_enc_ctx *enc, int ln,
 od_coeff *scalar_out, const od_coeff *cblock, const od_coeff *predt,
 int q, int pli) {
  int *adapt;
  int vk;
  int zzi;
  int n2;
  ogg_int32_t adapt_curr[OD_NSB_ADAPT_CTXS];
  adapt = enc->state.adapt.pvq_adapt;
  vk = 0;
  n2 = 1 << (2*ln + 4);
  for (zzi = 1; zzi < n2; zzi++) {
    scalar_out[zzi] = OD_DIV_R0(cblock[zzi] - predt[zzi], q);
    vk += abs(scalar_out[zzi]);
  }
#if defined(OD_METRICS)
  pvq_frac_bits = od_ec_enc_tell_frac(&enc->ec);
#endif
  generic_encode(&enc->ec, &enc->state.adapt.model_g[pli], vk, -1,
   &enc->state.adapt.ex_g[pli][ln], 0);
  laplace_encode_vector(&enc->ec, scalar_out + 1, n2 - 1, vk, adapt_curr,
   adapt);
#if defined(OD_METRICS)
  enc->state.bit_metrics[OD_METRIC_PVQ] += od_ec_enc_tell_frac(&enc->ec) -
   pvq_frac_bits;
#endif
  for (zzi = 1; zzi < n2; zzi++) {
    scalar_out[zzi] = scalar_out[zzi]*q + predt[zzi];
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

void od_block_encode(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int ln,
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
  od_coeff cblock[16*16];
  od_coeff scalar_out[16*16];
  int run_pvq;
  int quant;
  int dc_quant;
#ifndef USE_BAND_PARTITIONS
  unsigned char const *zig;
#endif
#if defined(OD_OUTPUT_PRED)
  od_coeff preds[16*16];
  int zzi;
#endif
  OD_ASSERT(ln >= 0 && ln <= 2);
  n = 1 << (ln + 2);
  run_pvq = ctx->run_pvq[pli];
  bx <<= ln;
  by <<= ln;
#ifndef USE_BAND_PARTITIONS
  zig = OD_DCT_ZIGS[ln];
#endif
  xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
  frame_width = enc->state.frame_width;
  w = frame_width >> xdec;
  c = ctx->c;
  d = ctx->d[pli];
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  tf = ctx->tf[pli];
  md = ctx->md;
  mc = ctx->mc;
  /* Apply forward transform. */
  if (!ctx->is_keyframe) {
    (*OD_FDCT_2D[ln])(md + (by << 2)*w + (bx << 2), w,
     mc + (by << 2)*w + (bx << 2), w);
  }
  od_encode_compute_pred(enc, ctx, pred, ln, pli, bx, by, has_ur);
#if defined(OD_OUTPUT_PRED)
  for (zzi = 0; zzi < (n*n); zzi++) preds[zzi] = pred[zzi];
#endif
  /* Change ordering for encoding. */
#ifdef USE_BAND_PARTITIONS
  od_raster_to_coding_order(cblock,  n, &d[((by << 2))*w + (bx << 2)], w,
   !run_pvq);
  od_raster_to_coding_order(predt,  n, &pred[0], n, !run_pvq);
#else
  /*Zig-zag*/
  {
    int y;
    for (y = 0; y < n; y++) {
      int x;
      for (x = 0; x < n; x++) {
        cblock[zig[y*n + x]] = d[((by << 2) + y)*w + (bx << 2) + x];
        predt[zig[y*n + x]] = pred[y*n + x];
      }
    }
  }
#endif
  /* Lossless encoding uses an actual quantizer of 1, but is signalled
     with a 'quantizer' of 0. */
  quant = OD_MAXI(1, enc->quantizer[pli]);
  if (run_pvq)
    dc_quant = OD_MAXI(1, quant*OD_PVQ_QM_Q4[pli][ln][0] >> 4);
  else
    dc_quant = (pli==0 || enc->quantizer[pli]==0) ? quant : (quant + 1) >> 1;
  if (OD_DISABLE_HAAR_DC || !ctx->is_keyframe) {
    if (abs(cblock[0] - predt[0]) < dc_quant * 141 / 256) { /* 0.55 */
      scalar_out[0] = 0;
    } else {
      scalar_out[0] = OD_DIV_R0(cblock[0] - predt[0], dc_quant);
    }
  }
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
    OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_AC_COEFFS);
  if (run_pvq) {
    pvq_encode(enc, predt, cblock, scalar_out, quant, pli, ln,
     OD_PVQ_QM_Q4[pli][ln], OD_PVQ_BETA[pli][ln],
     OD_PVQ_INTER_BAND_MASKING[ln], ctx->is_keyframe);
  }
  else {
    od_single_band_scalar_quant(enc, ln, scalar_out, cblock, predt, quant,
     pli);
  }
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  if (OD_DISABLE_HAAR_DC || !ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && run_pvq;
    OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
     OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_DC_COEFF);
    if (!has_dc_skip || scalar_out[0]) {
      generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli],
       abs(scalar_out[0]) - has_dc_skip, -1, &enc->state.adapt.ex_dc[pli][ln][0], 2);
    }
    if (scalar_out[0]) od_ec_enc_bits(&enc->ec, scalar_out[0] < 0, 1);
    OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
     OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
    scalar_out[0] = scalar_out[0]*dc_quant;
    scalar_out[0] += predt[0];
  }
  else {
    scalar_out[0] = cblock[0];
  }
#ifdef USE_BAND_PARTITIONS
  od_coding_order_to_raster(&d[((by << 2))*w + (bx << 2)], w, scalar_out, n,
   !run_pvq);
#else
  /*De-zigzag*/
  {
    int y;
    for (y = 0; y < n; y++) {
      int x;
      for (x = 0; x < n; x++) {
        d[((by << 2) + y)*w + (bx << 2) + x] = scalar_out[zig[y*n + x]];
      }
    }
  }
#endif
  /*Update the TF'd luma plane with CfL, or all the planes without CfL.*/
  if (ctx->is_keyframe && (pli == 0 || OD_DISABLE_CFL)) {
    od_convert_block_down(tf + (by << 2)*w + (bx << 2), w,
     d + (by << 2)*w + (bx << 2), w, ln, 0, 0);
  }
  /*Apply the inverse transform.*/
#if !defined(OD_OUTPUT_PRED)
  (*OD_IDCT_2D[ln])(c + (by << 2)*w + (bx << 2), w, d + (by << 2)*w
   + (bx << 2), w);
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
  (*OD_IDCT_2D[ln])(c + (by << 2)*w + (bx << 2), w, preds, n);
#endif
}

static void od_32x32_encode(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int ln,
 int pli, int bx, int by, int has_ur) {
  bx <<= 1;
  by <<= 1;
  od_block_encode(enc, ctx, ln - 1, pli, bx + 0, by + 0, 1);
  od_block_encode(enc, ctx, ln - 1, pli, bx + 1, by + 0, has_ur);
  od_block_encode(enc, ctx, ln - 1, pli, bx + 0, by + 1, 1);
  od_block_encode(enc, ctx, ln - 1, pli, bx + 1, by + 1, 0);
}

typedef void (*od_enc_func)(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int ln,
 int pli, int bx, int by, int has_ur);

const od_enc_func OD_ENCODE_BLOCK[OD_NBSIZES + 2] = {
  od_block_encode,
  od_block_encode,
  od_block_encode,
  od_32x32_encode
};

static void od_compute_dcts(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int pli,
  int bx, int by, int l, int xdec, int ydec) {
  int od;
  int d;
  int w;
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
    (*OD_FDCT_2D[d])(c + (by << (2 + d))*w + (bx << (2 + d)), w,
      ctx->c + (by << (2 + d))*w + (bx << (2 + d)), w);
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
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
     OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_DC_COEFF);
  if (l == 3) {
    int nhsb;
    int quant;
    int dc0;
    int l2;
    od_coeff sb_dc_pred;
    od_coeff sb_dc_curr;
    od_coeff *sb_dc_mem;
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
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
}
#endif

static void od_encode_block(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int pli,
 int bx, int by, int l, int xdec, int ydec, int has_ur) {
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
    (*OD_ENCODE_BLOCK[d])(enc, ctx, d, pli, bx, by, has_ur);
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_encode_block(enc, ctx, pli, bx + 0, by + 0, l, xdec, ydec, 1);
    od_encode_block(enc, ctx, pli, bx + 1, by + 0, l, xdec, ydec, has_ur);
    od_encode_block(enc, ctx, pli, bx + 0, by + 1, l, xdec, ydec, 1);
    od_encode_block(enc, ctx, pli, bx + 1, by + 1, l, xdec, ydec, 0);
  }
}

static void od_encode_mv(daala_enc_ctx *enc, od_mv_grid_pt *mvg, int vx,
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
  od_state_get_predictor(&enc->state, pred, vx, vy, level, mv_res);
  ox = (mvg->mv[0] >> mv_res) - pred[0];
  oy = (mvg->mv[1] >> mv_res) - pred[1];
  /*Interleave positive and negative values.*/
  mv_ex = enc->state.adapt.mv_ex;
  mv_ey = enc->state.adapt.mv_ey;
  model = &enc->state.adapt.mv_model;
  ex = mv_ex[level] >> mv_res;
  ey = mv_ex[level] >> mv_res;
  id = OD_MINI(abs(oy), 3)*4 + OD_MINI(abs(ox), 3);
  od_encode_cdf_adapt(&enc->ec, id, enc->state.adapt.mv_small_cdf, 16,
   enc->state.adapt.mv_small_increment);
  if (abs(ox) >= 3) generic_encode(&enc->ec, model, abs(ox) - 3, width << (3 - mv_res), &ex, 2);
  if (abs(oy) >= 3) generic_encode(&enc->ec, model, abs(oy) - 3, height << (3 - mv_res), &ey, 2);
  if (abs(ox)) od_ec_enc_bits(&enc->ec, ox < 0, 1);
  if (abs(oy)) od_ec_enc_bits(&enc->ec, oy < 0, 1);
  mv_ex[level] -= (mv_ex[level] - (abs(ox) << mv_res << 16)) >> 6;
  mv_ey[level] -= (mv_ey[level] - (abs(oy) << mv_res << 16)) >> 6;
}

int daala_encode_img_in(daala_enc_ctx *enc, od_img *img, int duration) {
  int refi;
  int nplanes;
  int pli;
  int frame_width;
  int frame_height;
  int pic_width;
  int pic_height;
  int i;
  int j;
  int k;
  int m;
  int nhsb;
  int nvsb;
  od_mb_enc_ctx mbctx;
#if defined(OD_ACCOUNTING)
  od_acct_reset(&enc->acct);
#endif
  if (enc == NULL || img == NULL) return OD_EFAULT;
  if (enc->packet_state == OD_PACKET_DONE) return OD_EINVAL;
  /*Check the input image dimensions to make sure they're compatible with the
     declared video size.*/
  nplanes = enc->state.info.nplanes;
  for (pli = 0; pli < OD_NPLANES_MAX; pli++) {
    mbctx.tf[pli] = enc->state.tf[pli];
    mbctx.modes[pli] = enc->state.modes[pli];
  }
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
  nhsb = enc->state.nhsb;
  nvsb = enc->state.nvsb;
  if (img->width != frame_width || img->height != frame_height) {
    /*The buffer does not match the frame size.
      Check to see if it matches the picture size.*/
    if (img->width != pic_width || img->height != pic_height) {
      /*It doesn't; we don't know how to handle it yet.*/
      return OD_EINVAL;
    }
  }
  /* Copy and pad the image. */
  for (pli = 0; pli < nplanes; pli++) {
    od_img_plane plane;
    int plane_width;
    int plane_height;
    *&plane = *(img->planes + pli);
    plane_width = ((pic_width + (1 << plane.xdec) - 1) >> plane.xdec);
    plane_height = ((pic_height + (1 << plane.ydec) - 1) >>
     plane.ydec);
    od_img_plane_copy_pad8(&enc->state.io_imgs[OD_FRAME_INPUT].planes[pli],
     frame_width >> plane.xdec, frame_height >> plane.ydec,
     &plane, plane_width, plane_height);
    od_img_plane_edge_ext8(&enc->state.io_imgs[OD_FRAME_INPUT].planes[pli],
     frame_width >> plane.xdec, frame_height >> plane.ydec,
     OD_UMV_PADDING >> plane.xdec, OD_UMV_PADDING >> plane.ydec);
  }
#if defined(OD_DUMP_IMAGES)
  if (od_logging_active(OD_LOG_GENERIC, OD_LOG_DEBUG)) {
    daala_info *info;
    od_img img;
    info = &enc->state.info;
    /*Modify the image offsets to include the padding.*/
    *&img = *(enc->state.io_imgs+OD_FRAME_INPUT);
    for (pli = 0; pli < nplanes; pli++) {
      img.planes[pli].data -= (OD_UMV_PADDING>>info->plane_info[pli].xdec)
       +img.planes[pli].ystride*(OD_UMV_PADDING>>info->plane_info[pli].ydec);
    }
    img.width += OD_UMV_PADDING<<1;
    img.height += OD_UMV_PADDING<<1;
    od_state_dump_img(&enc->state, &img, "pad");
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
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_FRAME);
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_FRAME);
  /*Write a bit to mark this as a data packet.*/
  od_ec_encode_bool_q15(&enc->ec, 0, 16384);
  /*Code the keyframe bit.*/
  od_ec_encode_bool_q15(&enc->ec, mbctx.is_keyframe, 16384);
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "is_keyframe=%d", mbctx.is_keyframe));
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_PLANE, OD_ACCT_TECH_UNKNOWN);
  /*TODO: Incrment frame count.*/
  /*Motion estimation and compensation.*/
  if (!mbctx.is_keyframe) {
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
  od_adapt_ctx_reset(&enc->state.adapt, mbctx.is_keyframe);
  /*Block size switching.*/
  od_state_init_border(&enc->state);
  /* Allocate a blockSizeComp for scratch space and then calculate the block sizes
     eventually store them in bsize. */
  od_log_matrix_uchar(OD_LOG_GENERIC, OD_LOG_INFO, "bimg ",
   enc->state.io_imgs[OD_FRAME_INPUT].planes[0].data -
   16*enc->state.io_imgs[OD_FRAME_INPUT].planes[0].ystride - 16,
   enc->state.io_imgs[OD_FRAME_INPUT].planes[0].ystride, (nvsb + 1)*32);
   OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
    OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_BLOCK_SIZE);
   OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
    OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_FRAME);
  for (i = 0; i < nvsb; i++) {
    unsigned char *bimg;
    unsigned char *rimg;
    int istride;
    int rstride;
    int bstride;
    int kf;
    kf = mbctx.is_keyframe;
    bstride = enc->state.bstride;
    istride = enc->state.io_imgs[OD_FRAME_INPUT].planes[0].ystride;
    rstride = kf ? 0 : enc->state.io_imgs[OD_FRAME_REC].planes[0].ystride;
    bimg = enc->state.io_imgs[OD_FRAME_INPUT].planes[0].data + i*istride*32;
    rimg = enc->state.io_imgs[OD_FRAME_REC].planes[0].data + i*rstride*32;
    for (j = 0; j < nhsb; j++) {
      int bsize[4][4];
      unsigned char *state_bsize;
      state_bsize = &enc->state.bsize[i*4*enc->state.bstride + j*4];
      process_block_size32(enc->bs, bimg + j*32, istride,
       kf ? NULL : rimg + j*32, rstride, bsize, enc->quantizer[0]);
      /* Grab the 4x4 information returned from process_block_size32 in bsize
         and store it in the od_state bsize. */
      for (k = 0; k < 4; k++) {
        for (m = 0; m < 4; m++) {
          if (OD_LIMIT_LOG_BSIZE_MIN != OD_LIMIT_LOG_BSIZE_MAX) {
            state_bsize[k*bstride + m] =
             OD_MAXI(OD_MINI(bsize[k][m], OD_LIMIT_LOG_BSIZE_MAX
             - OD_LOG_BSIZE0), OD_LIMIT_LOG_BSIZE_MIN - OD_LOG_BSIZE0);
          }
          else {
            state_bsize[k*bstride + m] =
             OD_LIMIT_LOG_BSIZE_MIN - OD_LOG_BSIZE0;
          }
        }
      }
      if (OD_LIMIT_LOG_BSIZE_MIN != OD_LIMIT_LOG_BSIZE_MAX) {
        od_block_size_encode(&enc->ec, &enc->state.adapt, &state_bsize[0], bstride);
      }
    }
  }
  OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
   OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  od_log_matrix_uchar(OD_LOG_GENERIC, OD_LOG_INFO, "bsize ", enc->state.bsize,
   enc->state.bstride, (nvsb + 1)*4);
  for (i = 0; i < nvsb*4; i++) {
    for (j = 0; j < nhsb*4; j++) {
      OD_LOG_PARTIAL((OD_LOG_GENERIC, OD_LOG_INFO, "%d ",
       enc->state.bsize[i*enc->state.bstride + j]));
    }
    OD_LOG_PARTIAL((OD_LOG_GENERIC, OD_LOG_INFO, "\n"));
  }
  /*Code the motion vectors.*/
  if (!mbctx.is_keyframe) {
    int nhmvbs;
    int nvmvbs;
    int vx;
    int vy;
    od_img *mvimg;
    int width;
    int height;
    int mv_res;
    od_mv_grid_pt *mvp;
    od_mv_grid_pt **grid;
    OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
     OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_MOTION_VECTORS);
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
    /*Level 1.*/
    for (vy = 2; vy <= nvmvbs; vy += 4) {
      for (vx = 2; vx <= nhmvbs; vx += 4) {
        int p_invalid;
        p_invalid = od_mv_level1_prob(grid, vx, vy);
        mvp = &(grid[vy][vx]);
        od_ec_encode_bool_q15(&enc->ec, mvp->valid, p_invalid);
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
          od_ec_encode_bool_q15(&enc->ec, mvp->valid, 13684);
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
            /*MVs are valid but will be zero.*/
          }
          else if (vx > 3 && vx < nhmvbs - 3) {
            if (!(vx & 2) && grid[vx == 3 ? vy - 1 : vy + 1][vx + 1].valid) {
              /*0 = both valid, 1 = only this one, 2 = other one valid*/
              int s;
              s = mvp->valid && grid[vy][vx + 2].valid ? 0 : mvp->valid
               + (grid[vy][vx + 2].valid << 1);
              od_ec_encode_cdf_q15(&enc->ec, s, OD_UNIFORM_CDF_Q15(3), 3);
              /*MVs are valid but will be zero.*/
            }
            else if (!(vx & 2)) {
              OD_ASSERT(!mvp->valid && !grid[vy][vx + 2].valid);
            }
          }
          else {
            OD_ASSERT(!mvp->valid);
          }
        }
        else if (vx < 2 || vx > nhmvbs - 2) {
          if ((vy == 3 && grid[vy - 1][vx == 1 ? vx - 1 : vx + 1].valid)
           || (vy == nvmvbs - 3
           && grid[vy - 1][vx == 1 ? vx - 1 : vx + 1].valid)) {
            /*MVs are valid but will be zero.*/
          }
          else if (!(vy & 2) && grid[vy + 1][vx == 1 ? vx - 1 : vx + 1].valid) {
            int s;
            s = mvp->valid && grid[vy + 2][vx].valid ? 0 : mvp->valid
             + (grid[vy + 2][vx].valid << 1);
            od_ec_encode_cdf_q15(&enc->ec, s, OD_UNIFORM_CDF_Q15(3), 3);
            /*MVs are valid but will be zero.*/
          }
          else if (!(vy & 2)) {
            OD_ASSERT(!mvp->valid && !grid[vy + 2][vx].valid);
          }
        }
        else if (grid[vy - 1][vx - 1].valid && grid[vy - 1][vx + 1].valid
         && grid[vy + 1][vx + 1].valid && grid[vy + 1][vx - 1].valid) {
          od_ec_encode_bool_q15(&enc->ec, mvp->valid, 16384);
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
          od_ec_encode_bool_q15(&enc->ec, mvp->valid, 16384);
          if (mvp->valid) {
            od_encode_mv(enc, mvp, vx, vy, 4, mv_res, width, height);
          }
          else {
            OD_ASSERT(!mvp->valid);
          }
        }
      }
    }
    OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
     OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
  }
  {
    int xdec;
    int ydec;
    int sby;
    int sbx;
    int h;
    int w;
    int y;
    int x;
    for (pli = 0; pli < nplanes; pli++) {
      xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /* Set this to 1 to enable the new (experimental, encode-only) PVQ
         implementation */
      mbctx.run_pvq[pli] = !OD_DISABLE_PVQ;
      OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
       OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_FRAME);
      OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
       OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_FRAME);
      /* TODO: We shouldn't be encoding the full, linear quantizer range. */
      od_ec_enc_uint(&enc->ec, enc->quantizer[pli], 512<<OD_COEFF_SHIFT);
      /*If the quantizer is zero (lossless), force scalar.*/
      if (!enc->quantizer[pli]) mbctx.run_pvq[pli] = 0;
      else od_ec_encode_bool_q15(&enc->ec, mbctx.run_pvq[pli], 16384);
      OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
       OD_ACCT_CAT_TECHNIQUE, OD_ACCT_TECH_UNKNOWN);
      OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
       OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_UNKNOWN);
    }
    for (pli = 0; pli < nplanes; pli++) {
      xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /*Collect the image data needed for this plane.*/
      {
        unsigned char *data;
        unsigned char *mdata;
        int ystride;
        int coeff_shift;
        coeff_shift = enc->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
        data = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].data;
        mdata = enc->state.io_imgs[OD_FRAME_REC].planes[pli].data;
        ystride = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
        for (y = 0; y < h; y++) {
          for (x = 0; x < w; x++) {
            enc->state.ctmp[pli][y*w + x] = (data[ystride*y + x] - 128) <<
             coeff_shift;
            if (!mbctx.is_keyframe) {
              enc->state.mctmp[pli][y*w + x] = (mdata[ystride*y + x] - 128)
               << coeff_shift;
            }
          }
        }
      }
      /*Apply the prefilter across the entire image.*/
      for (sby = 0; sby < nvsb; sby++) {
        for (sbx = 0; sbx < nhsb; sbx++) {
          od_apply_prefilter(enc->state.ctmp[pli], w, sbx, sby, 3,
           enc->state.bsize, enc->state.bstride, xdec, ydec,
           (sbx > 0 ? OD_LEFT_EDGE : 0) |
           (sby < nvsb - 1 ? OD_BOTTOM_EDGE : 0));
          if (!mbctx.is_keyframe) {
            od_apply_prefilter(enc->state.mctmp[pli], w, sbx, sby, 3, enc->state.bsize,
             enc->state.bstride, xdec, ydec, (sbx > 0 ? OD_LEFT_EDGE : 0) |
             (sby < nvsb - 1 ? OD_BOTTOM_EDGE : 0));
          }
        }
      }
    }
    for (sby = 0; sby < nvsb; sby++) {
      for (sbx = 0; sbx < nhsb; sbx++) {
        for (pli = 0; pli < nplanes; pli++) {
          OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
           OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_LUMA + pli);
          mbctx.c = enc->state.ctmp[pli];
          mbctx.d = enc->state.dtmp;
          mbctx.mc = enc->state.mctmp[pli];
          mbctx.md = enc->state.mdtmp[pli];
          mbctx.l = enc->state.lbuf[pli];
          xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
          ydec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
          mbctx.nk = mbctx.k_total = mbctx.sum_ex_total_q8 = 0;
          mbctx.ncount = mbctx.count_total_q8 = mbctx.count_ex_total_q8 = 0;
          /*Need to update this to decay based on superblocks width.*/
          od_compute_dcts(enc, &mbctx, pli, sbx, sby, 3, xdec, ydec);
          if (!OD_DISABLE_HAAR_DC && mbctx.is_keyframe) {
            od_quantize_haar_dc(enc, &mbctx, pli, sbx, sby, 3, xdec, ydec, 0,
             0, sby > 0 && sbx < nhsb - 1);
          }
          od_encode_block(enc, &mbctx, pli, sbx, sby, 3, xdec, ydec,
           sby > 0 && sbx < nhsb - 1);
        }
          OD_ACCT_UPDATE(&enc->acct, od_ec_enc_tell_frac(&enc->ec),
           OD_ACCT_CAT_PLANE, OD_ACCT_PLANE_UNKNOWN);
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /*Apply the postfilter across the entire image.*/
      for (sby = 0; sby < nvsb; sby++) {
        for (sbx = 0; sbx < nhsb; sbx++) {
          od_apply_postfilter(enc->state.ctmp[pli], w, sbx, sby, 3, enc->state.bsize,
           enc->state.bstride, xdec, ydec, (sby > 0 ? OD_TOP_EDGE : 0) |
           (sbx < nhsb - 1 ? OD_RIGHT_EDGE : 0));
        }
      }
      {
        unsigned char *data;
        int ystride;
        int coeff_shift;
        coeff_shift = enc->quantizer[pli] == 0 ? 0 : OD_COEFF_SHIFT;
        data = enc->state.io_imgs[OD_FRAME_REC].planes[pli].data;
        ystride = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
        for (y = 0; y < h; y++) {
          for (x = 0; x < w; x++) {
            data[ystride*y + x] = OD_CLAMP255(((enc->state.ctmp[pli][y*w + x]
             + (1 << coeff_shift >> 1)) >> coeff_shift) + 128);
          }
        }
      }
    }
  }
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  /*Dump YUV*/
  od_state_dump_yuv(&enc->state, enc->state.io_imgs + OD_FRAME_REC, "out");
#endif
#if defined(OD_LOGGING_ENABLED)
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
    data = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].data;
    ystride = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
    xdec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec = enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    npixels = w*h;
    for (y = 0; y < h; y++) {
      unsigned char *rec_row;
      unsigned char *inp_row;
      rec_row = enc->state.io_imgs[OD_FRAME_REC].planes[pli].data +
       enc->state.io_imgs[OD_FRAME_REC].planes[pli].ystride*y;
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
#endif
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO,
   "mode bits: %f/%f=%f", mode_bits, mode_count, mode_bits/mode_count));
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
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "Output Bytes: %ld", op->bytes));
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
