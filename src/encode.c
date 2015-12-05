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
#include "pvq_encoder.h"
#include "block_size.h"
#include "logging.h"
#include "tf.h"
#include "state.h"
#include "mcenc.h"
#include "quantizer.h"
#if defined(OD_X86ASM)
# include "x86/x86int.h"
#endif

/* These are the PVQ equivalent of quantization matrices, except that
   the values are per-band. */
#define OD_MASKING_DISABLED 0
#define OD_MASKING_ENABLED 1

#define OD_GOLDEN_FRAME_INTERVAL 10

static const unsigned char OD_LUMA_QM_Q4[2][OD_QM_SIZE] = {
/* Flat quantization for PSNR. The DC component isn't 16 because the DC
   magnitude compensation is done here for inter (Haar DC doesn't need it).
   Masking disabled: */
 {
  21, 16,
  18, 16, 16, 16,
  17, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
 },
/* The non-flat AC coefficients compensate for the non-linear scaling caused
   by activity masking. The values are currently hand-tuned so that the rate
   of each band remains roughly constant when enabling activity masking
   on intra.
   Masking enabled: */
 {
  21, 16,
  18, 18, 28, 32,
  17, 14, 20, 20, 28, 32,
  16, 11, 16, 14, 16, 16, 23, 28
 }
};

static const unsigned char OD_CHROMA_QM_Q4[2][OD_QM_SIZE] = {
/* Chroma quantization is different because of the reduced lapping.
   FIXME: Use the same matrix as luma for 4:4:4.
   Masking disabled: */
 {
  21, 16,
  18, 16, 16, 16,
  17, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
 },
/* The AC part is flat for chroma because it has no activity masking.
   Masking enabled: */
 {
  21, 16,
  18, 16, 16, 16,
  17, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
 }
};

typedef struct od_qm_entry {
  int interp_q;
  int scale_q8;
  const unsigned char *qm_q4;
} od_qm_entry;

/* No interpolation, always use od_flat_qm_q4, but use a different scale for
   each plane.
   FIXME: Add interpolation and properly tune chroma. */
static const od_qm_entry OD_DEFAULT_QMS[2][3][OD_NPLANES_MAX] = {
 /* Masking disabled */
 {{{4, 256, OD_LUMA_QM_Q4[OD_MASKING_DISABLED]},
   {4, 448, OD_CHROMA_QM_Q4[OD_MASKING_DISABLED]},
   {4, 320, OD_CHROMA_QM_Q4[OD_MASKING_DISABLED]}},
  {{318, 256, OD_LUMA_QM_Q4[OD_MASKING_DISABLED]},
   {318, 140, OD_CHROMA_QM_Q4[OD_MASKING_DISABLED]},
   {318, 100, OD_CHROMA_QM_Q4[OD_MASKING_DISABLED]}},
  {{0, 0, NULL},
   {0, 0, NULL},
   {0, 0, NULL}}},
 /* Masking enabled */
 {{{4, 256, OD_LUMA_QM_Q4[OD_MASKING_ENABLED]},
   {4, 448, OD_CHROMA_QM_Q4[OD_MASKING_ENABLED]},
   {4, 320, OD_CHROMA_QM_Q4[OD_MASKING_ENABLED]}},
  {{318, 256, OD_LUMA_QM_Q4[OD_MASKING_ENABLED]},
   {318, 140, OD_CHROMA_QM_Q4[OD_MASKING_ENABLED]},
   {318, 100, OD_CHROMA_QM_Q4[OD_MASKING_ENABLED]}},
  {{0, 0, NULL},
   {0, 0, NULL},
   {0, 0, NULL}}}
};

/*TODO: This makes little sense with the coded quantizer mapping
   changes, but that's a problem for later.
  Maintain current quality setting handling both here and in the
   encode_ctl.*/
static int od_quantizer_from_quality(int quality) {
  return quality == 0 ? 0 :
   (quality << OD_COEFF_SHIFT >> OD_QUALITY_SHIFT) +
   (1 << OD_COEFF_SHIFT >> 1);
}

void od_enc_opt_vtbl_init_c(od_enc_ctx *enc) {
  if (enc->state.full_precision_references) {
    enc->opt_vtbl.mc_compute_sad_4x4 =
      od_mc_compute_sad16_4x4_c;
    enc->opt_vtbl.mc_compute_sad_8x8 =
      od_mc_compute_sad16_8x8_c;
    enc->opt_vtbl.mc_compute_sad_16x16 =
      od_mc_compute_sad16_16x16_c;
    enc->opt_vtbl.mc_compute_sad_32x32 =
      od_mc_compute_sad16_32x32_c;
    enc->opt_vtbl.mc_compute_sad_64x64 =
      od_mc_compute_sad16_64x64_c;
    enc->opt_vtbl.mc_compute_satd_4x4 =
      od_mc_compute_satd16_4x4_c;
    enc->opt_vtbl.mc_compute_satd_8x8 =
      od_mc_compute_satd16_8x8_c;
    enc->opt_vtbl.mc_compute_satd_16x16 =
      od_mc_compute_satd16_16x16_c;
    enc->opt_vtbl.mc_compute_satd_32x32 =
      od_mc_compute_satd16_32x32_c;
    enc->opt_vtbl.mc_compute_satd_64x64 =
      od_mc_compute_satd16_64x64_c;
  }
  else {
    enc->opt_vtbl.mc_compute_sad_4x4 =
      od_mc_compute_sad8_4x4_c;
    enc->opt_vtbl.mc_compute_sad_8x8 =
      od_mc_compute_sad8_8x8_c;
    enc->opt_vtbl.mc_compute_sad_16x16 =
      od_mc_compute_sad8_16x16_c;
    enc->opt_vtbl.mc_compute_sad_32x32 =
      od_mc_compute_sad8_32x32_c;
    enc->opt_vtbl.mc_compute_sad_64x64 =
      od_mc_compute_sad8_64x64_c;
    enc->opt_vtbl.mc_compute_satd_4x4 =
      od_mc_compute_satd8_4x4_c;
    enc->opt_vtbl.mc_compute_satd_8x8 =
      od_mc_compute_satd8_8x8_c;
    enc->opt_vtbl.mc_compute_satd_16x16 =
      od_mc_compute_satd8_16x16_c;
    enc->opt_vtbl.mc_compute_satd_32x32 =
      od_mc_compute_satd8_32x32_c;
    enc->opt_vtbl.mc_compute_satd_64x64 =
      od_mc_compute_satd8_64x64_c;
  }
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
  int pli;
  od_img *img;
  od_img_plane *iplane;
  size_t data_sz;
  unsigned char *input_img_data;
  int frame_buf_width;
  int frame_buf_height;
  int plane_buf_width;
  int plane_buf_height;
  int reference_bytes;
  int reference_bits;
  int ret;
#if defined(OD_DUMP_BSIZE_DIST)
  char dist_fname[1024];
  const char *suf;
#endif
  ret = od_state_init(&enc->state, info);
  if (ret < 0) return ret;
  enc->use_satd = 0;
  od_enc_opt_vtbl_init(enc);
  oggbyte_writeinit(&enc->obb);
  od_ec_enc_init(&enc->ec, 65025);
  enc->packet_state = OD_PACKET_INFO_HDR;
  for (i = 0; i < OD_NPLANES_MAX; i++){
    enc->quality[i] = 10;
  }
  enc->complexity = 7;
  enc->use_activity_masking = 1;
  enc->qm = OD_HVS_QM;
  enc->use_haar_wavelet = OD_USE_HAAR_WAVELET;
  enc->mvest = od_mv_est_alloc(enc);
  if (OD_UNLIKELY(!enc->mvest)) {
    return OD_EFAULT;
  }
  enc->params.mv_level_min = 0;
  enc->params.mv_level_max = 4;
  enc->bs = (od_block_size_comp *)malloc(sizeof(*enc->bs));
  data_sz = 0;
  reference_bytes = enc->state.full_precision_references ? 2 : 1;
  reference_bits =
   enc->state.full_precision_references ? 8 + OD_COEFF_SHIFT : 8;
  /*TODO: Check for overflow before allocating.*/
  frame_buf_width = enc->state.frame_width + (OD_BUFFER_PADDING << 1);
  frame_buf_height = enc->state.frame_height + (OD_BUFFER_PADDING << 1);
  for (pli = 0; pli < info->nplanes; pli++) {
    plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
    plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
    /*Reserve space for input plane buffer */
    data_sz += plane_buf_width*plane_buf_height*reference_bytes;
#if defined(OD_DUMP_IMAGES)
    /*Reserve space for this plane in 1 visualization image.*/
    data_sz += plane_buf_width*plane_buf_height << 2;
    /*Reserve space for this plane in 1 temporary image used to obtain
       the visualization image.*/
    data_sz += plane_buf_width*plane_buf_height << 2;
    /*Reserve space for the line buffer in the up-sampler.*/
    /*Upsampler is vis-only and always runs at 8 bits.*/
    data_sz += (frame_buf_width << 1) * 8;
#endif
  }
  enc->input_img_data = input_img_data =
    (unsigned char *)od_aligned_malloc(data_sz, 32);
  if (OD_UNLIKELY(!input_img_data)) {
    return OD_EFAULT;
  }
  /*Fill in the input img structure.*/
  img = &enc->input_img;
  img->nplanes = info->nplanes;
  img->width = enc->state.frame_width;
  img->height = enc->state.frame_height;
  for (pli = 0; pli < img->nplanes; pli++) {
    plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
    plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
    iplane = img->planes + pli;
    iplane->xdec = info->plane_info[pli].xdec;
    iplane->ydec = info->plane_info[pli].ydec;
    iplane->bitdepth = reference_bits;
    /*Internal buffers are always planar.*/
    iplane->xstride = reference_bytes;
    iplane->ystride = plane_buf_width*iplane->xstride;
    iplane->data = input_img_data
     + iplane->xstride*(OD_BUFFER_PADDING >> info->plane_info[pli].xdec)
     + iplane->ystride*(OD_BUFFER_PADDING >> info->plane_info[pli].ydec);
    input_img_data += plane_buf_height*iplane->ystride;
  }
#if defined(OD_DUMP_IMAGES)
  /*Fill in the line buffers.*/
  {
    int y;
    for (y = 0; y < 8; y++) {
      enc->upsample_line_buf[y] = input_img_data + (OD_BUFFER_PADDING << 1);
      input_img_data += frame_buf_width << 1;
    }
  }
  /*Fill in the visualization image structure.*/
  img = &enc->vis_img;
  img->nplanes = info->nplanes;
  img->width = frame_buf_width << 1;
  img->height = frame_buf_height << 1;
  for (pli = 0; pli < img->nplanes; pli++) {
    iplane = img->planes + pli;
    plane_buf_width = img->width >> info->plane_info[pli].xdec;
    plane_buf_height = img->height >> info->plane_info[pli].ydec;
    iplane->data = input_img_data;
    input_img_data += plane_buf_width*plane_buf_height;
    iplane->xdec = info->plane_info[pli].xdec;
    iplane->ydec = info->plane_info[pli].ydec;
    iplane->bitdepth = 8;
    iplane->xstride = 1;
    iplane->ystride = plane_buf_width;
  }
  /*Fill in the temporary image structure.*/
  img = &enc->tmp_vis_img;
  img->nplanes = info->nplanes;
  img->width = frame_buf_width << 1;
  img->height = frame_buf_height << 1;
  for (pli = 0; pli < img->nplanes; pli++) {
    iplane = img->planes + pli;
    plane_buf_width = img->width >> info->plane_info[pli].xdec;
    plane_buf_height = img->height >> info->plane_info[pli].ydec;
    iplane->data = input_img_data;
    input_img_data += plane_buf_width*plane_buf_height;
    iplane->xdec = info->plane_info[pli].xdec;
    iplane->ydec = info->plane_info[pli].ydec;
    iplane->bitdepth = 8;
    iplane->xstride = 1;
    iplane->ystride = plane_buf_width;
  }
#endif
#if defined(OD_ENCODER_CHECK)
  enc->dec = daala_decode_create(info, NULL);
#endif
#if defined(OD_DUMP_BSIZE_DIST)
  suf = getenv("OD_DUMP_BSIZE_DIST_SUFFIX");
  if (!suf) {
    suf = "dist";
  }
  sprintf(dist_fname, "%08i-%s.out",
   (int)daala_granule_basetime(&enc->state, enc->state.cur_time), suf);
  for (i = 0; i < OD_NPLANES_MAX; i++){
    enc->bsize_dist[i] = 0.0;
    enc->bsize_dist_total[i] = 0.0;
  }
  /* FIXME: something unique-ish maybe? */
  enc->bsize_dist_file = fopen(dist_fname, "wb");
  if (OD_UNLIKELY(enc->bsize_dist_file == NULL)) {
    return OD_EFAULT;
  }
#endif
  return 0;
}

static void od_enc_clear(od_enc_ctx *enc) {
  od_mv_est_free(enc->mvest);
  od_ec_enc_clear(&enc->ec);
  oggbyte_writeclear(&enc->obb);
  od_aligned_free(enc->input_img_data);
  od_state_clear(&enc->state);
}

daala_enc_ctx *daala_encode_create(const daala_info *info) {
  od_enc_ctx *enc;
  if (info == NULL) return NULL;
  enc = (od_enc_ctx *)malloc(sizeof(*enc));
  if (od_enc_init(enc, info) < 0) {
    free(enc);
    return NULL;
  }
  return enc;
}

void daala_encode_free(daala_enc_ctx *enc) {
#if defined(OD_DUMP_BSIZE_DIST)
    int i;
    fprintf(enc->bsize_dist_file, "Total: ");
    for (i = 0; i < enc->state.info.nplanes; i++){
      fprintf(enc->bsize_dist_file, "%-7G\t",
       10*log10(enc->bsize_dist_total[i]));
    }
    fprintf(enc->bsize_dist_file, "\n");
    fclose(enc->bsize_dist_file);
#endif
  if (enc != NULL) {
#if defined(OD_ENCODER_CHECK)
    if (enc->dec != NULL) {
      daala_decode_free(enc->dec);
    }
#endif
    free(enc->bs);
    od_enc_clear(enc);
    free(enc);
  }
}

int daala_encode_ctl(daala_enc_ctx *enc, int req, void *buf, size_t buf_sz) {
  (void)buf;
  (void)buf_sz;
  switch (req) {
    case OD_SET_QUANT:
    {
      int i;
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(*enc->quality), OD_EINVAL);
      for (i = 0; i < OD_NPLANES_MAX; i++){
        int tmp = *(int *)buf;
        enc->quality[i] = tmp > 0 ? (tmp << OD_QUALITY_SHIFT) - 8 : 0;
      }
      return OD_SUCCESS;
    }
    case OD_SET_COMPLEXITY: {
      int complexity;
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(enc->complexity), OD_EINVAL);
      complexity = *(const int *)buf;
      if (complexity < 0 || complexity > 10) return OD_EINVAL;
      enc->complexity = complexity;
      return OD_SUCCESS;
    }
    case OD_GET_COMPLEXITY: {
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(enc->complexity), OD_EINVAL);
      *(int *)buf = enc->complexity;
      return OD_SUCCESS;
    }
    case OD_SET_MC_CHROMA:
    {
      int mc_use_chroma;
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(mc_use_chroma), OD_EINVAL);
      mc_use_chroma = *(int *)buf;
      if (mc_use_chroma) {
        enc->mvest->flags |= OD_MC_USE_CHROMA;
      }
      else {
        enc->mvest->flags &= ~OD_MC_USE_CHROMA;
      }
      return OD_SUCCESS;
    }
    case OD_SET_MC_SATD: {
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(enc->use_satd), OD_EINVAL);
      enc->use_satd = !!*(const int *)buf;
      return OD_SUCCESS;
    }
    case OD_SET_ACTIVITY_MASKING: {
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(enc->use_activity_masking), OD_EINVAL);
      enc->use_activity_masking = !!*(const int *)buf;
      return OD_SUCCESS;
    }
    case OD_SET_DERING: {
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(enc->use_dering), OD_EINVAL);
      enc->use_dering = !!*(const int *)buf;
      return OD_SUCCESS;
    }
    case OD_SET_QM: {
      int qm;
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(qm), OD_EINVAL);
      qm = *(const int *)buf;
      if (qm < OD_FLAT_QM || qm > OD_HVS_QM) {
          return OD_EINVAL;
      }
      enc->qm = qm;
      return OD_SUCCESS;
    }
    case OD_SET_MV_RES_MIN:
    {
      int mv_res_min;
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(mv_res_min), OD_EINVAL);
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
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(mv_level_min), OD_EINVAL);
      mv_level_min = *(int *)buf;
      if (mv_level_min < 0 || mv_level_min > OD_MC_LEVEL_MAX) {
        return OD_EINVAL;
      }
      enc->params.mv_level_min = mv_level_min;
      return OD_SUCCESS;
    }
    case OD_SET_MV_LEVEL_MAX:
    {
      int mv_level_max;
      OD_RETURN_CHECK(enc, OD_EFAULT);
      OD_RETURN_CHECK(buf, OD_EFAULT);
      OD_RETURN_CHECK(buf_sz == sizeof(mv_level_max), OD_EINVAL);
      mv_level_max = *(int *)buf;
      if (mv_level_max < 0 || mv_level_max > OD_MC_LEVEL_MAX) {
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

static void od_img_plane_copy_pad(od_img *dst,
 int plane_width, int plane_height, od_img *src,
 int pic_width, int pic_height, int pli) {
  od_img_plane *dst_p;
  unsigned char *dst_data;
  int dst_xstride;
  int dst_ystride;
  int x;
  int y;
  dst_p = dst->planes + pli;
  dst_xstride = dst_p->xstride;
  dst_ystride = dst_p->ystride;

  /*If we have _no_ data, just encode a dull green.*/
  if (pic_width == 0 || pic_height == 0) {
    dst_data = dst_p->data;
    for (y = 0; y < plane_height; y++) {
      OD_CLEAR(dst_data, plane_width*dst_xstride);
      dst_data += dst_ystride;
    }
  }
  else {
    /*Otherwise, Step 1: Copy the data we do have.*/
    od_img_plane_copy(dst, src, pli);
    /*Step 2: Perform a low-pass extension into the padding region.*/
    /*Right side.*/
    for (x = pic_width; x < plane_width; x++) {
      dst_data = dst_p->data + (x - 1)*dst_xstride;
      if (dst_xstride == 1) {
        for (y = 0; y < pic_height; y++) {
          unsigned char uppercase_u;
          unsigned char uppercase_c;
          unsigned char uppercase_d;
          uppercase_c = *dst_data;
          uppercase_u = *(dst_data - (dst_ystride & -(y > 0)));
          uppercase_d = *(dst_data + (dst_ystride & -(y + 1 < pic_height)));
          dst_data[1] = (2*uppercase_c + uppercase_u + uppercase_d + 2) >> 2;
          dst_data += dst_ystride;
        }
      }
      else {
        for (y = 0; y < pic_height; y++) {
          uint16_t uppercase_u;
          uint16_t uppercase_c;
          uint16_t uppercase_d;
          uppercase_c = *(uint16_t *)dst_data;
          uppercase_u = *(uint16_t *)(dst_data - (dst_ystride & -(y > 0)));
          uppercase_d = *(uint16_t *)(dst_data +
           (dst_ystride & -(y + 1 < pic_height)));
          ((uint16_t *)dst_data)[1] =
           (2*uppercase_c + uppercase_u + uppercase_d + 2) >> 2;
          dst_data += dst_ystride;
        }
      }
    }
    /*Bottom.*/
    dst_data = dst_p->data + dst_ystride*pic_height;
    for (y = pic_height; y < plane_height; y++) {
      if (dst_xstride == 1) {
        for (x = 0; x < plane_width; x++) {
          unsigned char uppercase_l;
          unsigned char uppercase_c;
          unsigned char uppercase_r;
          uppercase_c = (dst_data - dst_ystride)[x];
          uppercase_l = (dst_data - dst_ystride)[x - (x > 0)];
          uppercase_r = (dst_data - dst_ystride)[x + (x + 1 < plane_width)];
          dst_data[x] = (2*uppercase_c + uppercase_l + uppercase_r + 2) >> 2;
        }
      }
      else{
        for (x = 0; x < plane_width; x++) {
          uint16_t uppercase_l;
          uint16_t uppercase_c;
          uint16_t uppercase_r;
          uppercase_c = ((uint16_t *)(dst_data - dst_ystride))[x];
          uppercase_l = ((uint16_t *)(dst_data - dst_ystride))[x - (x > 0)];
          uppercase_r =
           ((uint16_t *)(dst_data - dst_ystride))[x + (x + 1 < plane_width)];
          ((uint16_t *)dst_data)[x] =
           (2*uppercase_c + uppercase_l + uppercase_r + 2) >> 2;
        }
      }
      dst_data += dst_ystride;
    }
  }
}

/*Block-level encoder context information.
  Global encoder context information is in od_enc_ctx.*/
struct od_mb_enc_ctx {
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
  int q_scaling;
};
typedef struct od_mb_enc_ctx od_mb_enc_ctx;

static void od_encode_compute_pred(daala_enc_ctx *enc, od_mb_enc_ctx *ctx,
 od_coeff *pred, const od_coeff *d, int bs, int pli, int bx, int by) {
  int n;
  int xdec;
  int w;
  int bo;
  int y;
  int x;
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  n = 1 << (bs + OD_LOG_BSIZE0);
  xdec = enc->input_img.planes[pli].xdec;
  w = enc->state.frame_width >> xdec;
  bo = (by << OD_LOG_BSIZE0)*w + (bx << OD_LOG_BSIZE0);
  /*We never use tf on the chroma planes, but if we do it will blow up, which
    is better than always using luma's tf.*/
  if (ctx->is_keyframe) {
    if (pli == 0 || OD_DISABLE_CFL || ctx->use_haar_wavelet) {
      OD_CLEAR(pred, n*n);
      if (pli == 0 && !ctx->use_haar_wavelet) {
        od_hv_intra_pred(pred, d, w, bx, by, enc->state.bsize,
         enc->state.bstride, bs);
      }
    }
    else {
      for (y = 0; y < n; y++) {
        for (x = 0; x < n; x++) {
          pred[n*y + x] = ctx->l[n*y + x];
        }
      }
    }
  }
  else {
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        pred[n*y + x] = ctx->md[bo + y*w + x];
      }
    }
  }
}

/* Compute the sum of the tree (parent and descendents) at each level. */
static int od_compute_max_tree(od_coeff tree_sum[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int x, int y, const od_coeff *c, int ln) {
  int n;
  int maxval;
  n = 1 << ln;
  maxval = 0;
  if (2*x < n && 2*y < n) {
    int tmp;
    tmp = od_compute_max_tree(tree_sum, 2*x, 2*y, c, ln);
    maxval = maxval + tmp;
    tmp = od_compute_max_tree(tree_sum, 2*x + 1, 2*y, c, ln);
    maxval = maxval + tmp;
    tmp = od_compute_max_tree(tree_sum, 2*x, 2*y + 1, c, ln);
    maxval = maxval + tmp;
    tmp = od_compute_max_tree(tree_sum, 2*x + 1, 2*y + 1, c, ln);
    maxval = maxval + tmp;
  }
  maxval = maxval + abs(c[y*n + x]);
  tree_sum[y][x] = maxval;
  return maxval;
}

/* Encode with unary (Rice) code.
   FIXME: This should go away. */
static void od_ec_enc_unary(od_ec_enc *ec, int x) {
  if (x) od_ec_enc_bits(ec, 0, x);
  od_ec_enc_bits(ec, 1, 1);
}

/* Encodes the magnitude of the coefficient at the root of a tree based
   on the sum for the entire tree (distribution can be highly biased). */
static void od_encode_coeff_split(daala_enc_ctx *enc, int a, int sum,
 int ctx) {
  int shift;
  if (sum == 0) return;
  shift = OD_MAXI(0, OD_ILOG(sum) - 4);
  if (shift) {
    od_ec_enc_bits(&enc->ec, a & ((1 << shift) - 1), shift);
    a >>= shift;
    sum >>= shift;
  }
  od_encode_cdf_adapt(&enc->ec, a, enc->state.adapt.haar_coeff_cdf[15*ctx + sum
   - 1], sum + 1, enc->state.adapt.haar_coeff_increment);
}

/* Encodes the magnitude of one side of a tree split based on the sum of the
   two sides. The distribution should be roughly symmetric. */
static void od_encode_tree_split(daala_enc_ctx *enc, int a, int sum, int ctx) {
  int shift;
  if (sum == 0) return;
  shift = OD_MAXI(0, OD_ILOG(sum) - 4);
  if (shift) {
    od_ec_enc_bits(&enc->ec, a & ((1 << shift) - 1), shift);
    a >>= shift;
    sum >>= shift;
  }
  od_encode_cdf_adapt(&enc->ec, a, enc->state.adapt.haar_split_cdf[15*
   (2*ctx + OD_MINI(shift, 1)) + sum - 1], sum + 1,
   enc->state.adapt.haar_split_increment);
}

static void od_encode_sum_tree(daala_enc_ctx *enc, const od_coeff *c, int ln,
 od_coeff tree_sum[OD_BSIZE_MAX][OD_BSIZE_MAX], int x, int y, int dir,
 int pli) {
  int n;
  int coeff_mag;
  int children_sum;
  n = 1 << ln;
  if (tree_sum[y][x] == 0) return;
  coeff_mag = abs(c[y*n + x]);
  od_encode_coeff_split(enc, coeff_mag, tree_sum[y][x], dir
   + 3*(OD_ILOG(OD_MAXI(x,y)) - 1));
  /* Encode max of each four children relative to tree. */
  children_sum = tree_sum[2*y][2*x] + tree_sum[2*y][2*x + 1]
   + tree_sum[2*y + 1][2*x] + tree_sum[2*y + 1][2*x + 1];
  if (children_sum) {
    OD_ASSERT(coeff_mag + children_sum == tree_sum[y][x]);
    /* Split in a different order depending on direction. */
    if (dir == 0) {
      od_encode_tree_split(enc, tree_sum[2*y][2*x] + tree_sum[2*y][2*x + 1],
       children_sum, 0);
      od_encode_tree_split(enc, tree_sum[2*y][2*x], tree_sum[2*y][2*x]
       + tree_sum[2*y][2*x + 1], 2);
      od_encode_tree_split(enc, tree_sum[2*y + 1][2*x], tree_sum[2*y + 1][2*x]
       + tree_sum[2*y + 1][2*x + 1], 2);
    }
    else {
      od_encode_tree_split(enc, tree_sum[2*y][2*x] + tree_sum[2*y + 1][2*x],
       children_sum, 1);
      od_encode_tree_split(enc, tree_sum[2*y][2*x], tree_sum[2*y][2*x]
       + tree_sum[2*y + 1][2*x], 2);
      od_encode_tree_split(enc, tree_sum[2*y][2*x + 1], tree_sum[2*y][2*x + 1]
       + tree_sum[2*y + 1][2*x + 1], 2);
    }
  }
  if (4*x < n && 4*y < n) {
    /* Recursive calls. */
    od_encode_sum_tree(enc, c, ln, tree_sum, 2*x, 2*y, dir, pli);
    od_encode_sum_tree(enc, c, ln, tree_sum, 2*x + 1, 2*y, dir, pli);
    od_encode_sum_tree(enc, c, ln, tree_sum, 2*x, 2*y + 1, dir, pli);
    od_encode_sum_tree(enc, c, ln, tree_sum, 2*x + 1, 2*y + 1, dir, pli);
  }
}

static int od_wavelet_quantize(daala_enc_ctx *enc, int ln,
 od_coeff *out, const od_coeff *dblock, const od_coeff *predt,
 int quant, int pli) {
  int n;
  int i, j;
  int dir;
  od_coeff tree_sum[OD_BSIZE_MAX][OD_BSIZE_MAX];
  n = 1 << ln;
  /* Quantize everything but DC. */
  for (dir = 0; dir < 3; dir++) {
    int level;
    for (level = 0; level < ln; level++) {
      int bo;
      int q;
      bo = (((dir + 1) >> 1) << level)*n + (((dir + 1) & 1) << level);
      if (quant == 0) q = 1;
      else q = quant*OD_HAAR_QM[dir == 2][level] >> 4;
      for (i = 0; i < 1 << level; i++) {
        for (j = 0; j < 1 << level; j++) {
          out[bo + i*n + j] = OD_DIV_R0(dblock[bo + i*n + j]
           - predt[bo + i*n + j], q);
        }
      }
    }
  }
  /* Compute magnitude at each level of each tree. */
  od_compute_max_tree(tree_sum, 1, 0, out, ln);
  od_compute_max_tree(tree_sum, 0, 1, out, ln);
  od_compute_max_tree(tree_sum, 1, 1, out, ln);
  /* Encode magnitude for the top of each tree */
  tree_sum[0][0] = tree_sum[0][1] + tree_sum[1][0] + tree_sum[1][1];
  {
    int bits;
    bits = OD_ILOG(tree_sum[0][0]);
    /* Encode the sum of quantized coefficients for the entire block as log2(x)
       followed by the LSBs. */
    od_encode_cdf_adapt(&enc->ec, OD_MINI(bits, 15),
     enc->state.adapt.haar_bits_cdf[pli], 16,
     enc->state.adapt.haar_bits_increment);
    if (bits >= 15) od_ec_enc_unary(&enc->ec, bits - 15);
    if (bits > 1) {
      od_ec_enc_bits(&enc->ec, tree_sum[0][0] & ((1 << (bits - 1)) - 1),
       bits - 1);
    }
    /* Encode diagonal tree sum. */
    od_encode_tree_split(enc, tree_sum[1][1], tree_sum[0][0], 3);
    /* Horizontal vs vertical. */
    od_encode_tree_split(enc, tree_sum[0][1], tree_sum[0][0] - tree_sum[1][1],
     4);
  }
  /* Encode all 3 trees. */
  od_encode_sum_tree(enc, out, ln, tree_sum, 1, 0, 0, pli);
  od_encode_sum_tree(enc, out, ln, tree_sum, 0, 1, 1, pli);
  od_encode_sum_tree(enc, out, ln, tree_sum, 1, 1, 2, pli);
  /* For all significant coeffs, encode sign. */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) if (i + j) {
      od_coeff in;
      in = out[i*n + j];
      if (in) od_ec_enc_bits(&enc->ec, in < 0, 1);
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
        for (j = 0; j < 1 << level; j++)
          out[bo + i*n + j] = q*out[bo + i*n + j] + predt[bo + i*n + j];
      }
    }
  }
  return 0;
}

static int od_compute_var_4x4(od_coeff *x, int stride) {
  int sum;
  int s2;
  int i;
  sum = 0;
  s2 = 0;
  for (i = 0; i < 4; i++) {
    int j;
    for (j = 0; j < 4; j++) {
      int t;
      /* Avoids overflow in the sum^2 below because the pre-filtered input
         can be much larger than +/-128 << OD_COEFF_SHIFT. Shifting the sum
         itself is a bad idea because it leads to large error on low
         variance. */
      t = x[i*stride + j] >> 2;
      sum += t;
      s2 += t*t;
    }
  }
  return (s2 - (sum*sum >> 4));
}

static double od_compute_dist_8x8(daala_enc_ctx *enc, od_coeff *x, od_coeff *y,
 int stride, int bs) {
  od_coeff e[8*8];
  od_coeff et[8*8];
  double sum;
  int min_var;
  double mean_var;
  double var_stat;
  double activity;
  double calibration;
  int i;
  int j;
  double vardist;
  vardist = 0;
  OD_ASSERT(enc->qm != OD_FLAT_QM);
#if 1
  min_var = INT_MAX;
  mean_var = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      int varx;
      int vary;
      double diff;
      varx = od_compute_var_4x4(x + 2*i*stride + 2*j, stride);
      vary = od_compute_var_4x4(y + 2*i*stride + 2*j, stride);
      min_var = OD_MINI(min_var, varx);
      mean_var += 1./(1 + varx);
      diff = sqrt(varx) - sqrt(vary);
      vardist += diff*diff;
    }
  }
  /* We use a different variance statistic depending on whether activity
     masking is used, since the harmonic mean appeared slghtly worse with
     masking off. The calibration constant just ensures that we preserve the
     rate compared to activity=1. */
  if (enc->use_activity_masking) {
    calibration = 1.95;
    var_stat = 9./mean_var;
  }
  else {
    calibration = 1.62;
    var_stat = min_var;
  }
  /* 1.62 is a calibration constant, 0.25 is a noise floor and 1/6 is the
     activity masking constant. */
  activity = calibration*pow(.25 + var_stat/(1 << 2*OD_COEFF_SHIFT), -1./6);
#else
  activity = 1;
#endif
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) e[8*i + j] = x[i*stride + j] - y[i*stride + j];
  }
  (*enc->state.opt_vtbl.fdct_2d[OD_BLOCK_8X8])(&et[0], 8, &e[0], 8);
  sum = 0;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      double mag;
      mag = 16./OD_QM8_Q4_HVS[i*8 + j];
      /* We attempt to consider the basis magnitudes here, though that's not
         perfect for block size 16x16 and above since only some edges are
         filtered then. */
      mag *= OD_BASIS_MAG[0][bs][i << (bs - 1)]*
       OD_BASIS_MAG[0][bs][j << (bs - 1)];
      mag *= mag;
      sum += et[8*i + j]*(double)et[8*i + j]*mag;
    }
  }
  return activity*activity*(sum + vardist);
}

static double od_compute_dist(daala_enc_ctx *enc, od_coeff *x, od_coeff *y,
 int n, int bs) {
  int i;
  double sum;
  sum = 0;
  if (enc->qm == OD_FLAT_QM) {
    for (i = 0; i < n*n; i++) {
      double tmp;
      tmp = x[i] - y[i];
      sum += tmp*tmp;
    }
  }
  else {
    for (i = 0; i < n; i += 8) {
      int j;
      for (j = 0; j < n; j += 8) {
        sum += od_compute_dist_8x8(enc, &x[i*n + j], &y[i*n + j], n, bs);
      }
    }
    /* Compensate for the fact that the quantization matrix lowers the
       distortion value. We tried a half-dozen values and picked the one where
       we liked the ntt-short1 curves best. The tuning is approximate since
       the different metrics go in different directions. */
    sum *= 1.7;
  }
  return sum;
}

/* Computes block size RDO lambda (for 1/8 bits) from the quantizer. */
static double od_bs_rdo_lambda(int q) {
  return OD_PVQ_LAMBDA*(1./(1 << OD_BITRES))*q*q;
}

/* Returns 1 if the block is skipped, zero otherwise. */
static int od_block_encode(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int bs,
 int pli, int bx, int by, int rdo_only) {
  int n;
  int xdec;
  int w;
  int bo;
  int frame_width;
  int use_masking;
  od_coeff *c;
  od_coeff *d;
  od_coeff *md;
  od_coeff *mc;
  od_coeff pred[OD_BSIZE_MAX*OD_BSIZE_MAX];
  od_coeff predt[OD_BSIZE_MAX*OD_BSIZE_MAX];
  od_coeff dblock[OD_BSIZE_MAX*OD_BSIZE_MAX];
  od_coeff scalar_out[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int quant;
  int dc_quant;
  int lossless;
  int skip;
  const int *qm;
  double dist_noskip;
  int tell;
  int i;
  int j;
  int has_late_skip_rdo;
  od_rollback_buffer pre_encode_buf;
  od_coeff *c_orig;
  od_coeff *mc_orig;
  qm = ctx->qm == OD_HVS_QM ? OD_QM8_Q4_HVS : OD_QM8_Q4_FLAT;
#if defined(OD_OUTPUT_PRED)
  od_coeff preds[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int zzi;
#endif
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  tell = 0;
  n = 1 << (bs + 2);
  bx <<= bs;
  by <<= bs;
  xdec = enc->input_img.planes[pli].xdec;
  frame_width = enc->state.frame_width;
  use_masking = enc->use_activity_masking;
  w = frame_width >> xdec;
  bo = (by << 2)*w + (bx << 2);
  c = ctx->c;
  d = ctx->d[pli];
  md = ctx->md;
  mc = ctx->mc;
  lossless = OD_LOSSLESS(enc, pli);
  c_orig = enc->block_c_orig;
  mc_orig = enc->block_mc_orig;
  has_late_skip_rdo = !ctx->is_keyframe && !ctx->use_haar_wavelet && bs > 0;
  if (has_late_skip_rdo) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) c_orig[n*i + j] = c[bo + i*w + j];
    }
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) mc_orig[n*i + j] = mc[bo + i*w + j];
    }
    tell = od_ec_enc_tell_frac(&enc->ec);
    od_encode_checkpoint(enc, &pre_encode_buf);
  }
  /* Apply forward transform. */
  if (ctx->use_haar_wavelet) {
    if (rdo_only || !ctx->is_keyframe) {
      od_haar(d + bo, w, c + bo, w, bs + 2);
    }
    if (!ctx->is_keyframe) {
      od_haar(md + bo, w, mc + bo, w, bs + 2);
    }
  }
  else {
    if (rdo_only || !ctx->is_keyframe) {
      int quantized_dc;
      quantized_dc = d[bo];
      (*enc->state.opt_vtbl.fdct_2d[bs])(d + bo, w, c + bo, w);
      if (ctx->is_keyframe) d[bo] = quantized_dc;
      od_apply_qm(d + bo, w, d + bo, w, bs, xdec, 0, qm);
    }
    if (!ctx->is_keyframe) {
      (*enc->state.opt_vtbl.fdct_2d[bs])(md + bo, w, mc + bo, w);
      od_apply_qm(md + bo, w, md + bo, w, bs, xdec, 0, qm);
    }
  }
  od_encode_compute_pred(enc, ctx, pred, d, bs, pli, bx, by);
#if defined(OD_OUTPUT_PRED)
  for (zzi = 0; zzi < (n*n); zzi++) preds[zzi] = pred[zzi];
#endif
  if (ctx->use_haar_wavelet) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        dblock[i*n + j] = d[bo + i*w + j];
        predt[i*n + j] = pred[i*n + j];
      }
    }
  }
  else {
    /* Change ordering for encoding. */
    od_raster_to_coding_order(dblock,  n, &d[bo], w);
    od_raster_to_coding_order(predt,  n, &pred[0], n);
  }
  /* Lossless encoding uses an actual quantizer of 1, but is signalled
     with a 'quantizer' of 0. */
  quant = OD_MAXI(1, enc->state.quantizer[pli]);
  if (lossless) dc_quant = quant;
  else {
    dc_quant = OD_MAXI(1, quant*
     enc->state.pvq_qm_q4[pli][od_qm_get_index(bs, 0)] >> 4);
  }
  /* This quantization may be overridden in the PVQ code for full RDO. */
  if (!ctx->is_keyframe) {
    if (abs(dblock[0] - predt[0]) < dc_quant*141/256) { /* 0.55 */
      scalar_out[0] = 0;
    }
    else {
      scalar_out[0] = OD_DIV_R0(dblock[0] - predt[0], dc_quant);
    }
  }
  if (ctx->use_haar_wavelet) {
    skip = od_wavelet_quantize(enc, bs + 2, scalar_out, dblock, predt,
     enc->state.quantizer[pli], pli);
  }
  else {
    skip = od_pvq_encode(enc, predt, dblock, scalar_out, quant, pli, bs,
     OD_PVQ_BETA[use_masking][pli][bs], OD_ROBUST_STREAM, ctx->is_keyframe,
     ctx->q_scaling, bx, by);
  }
  if (!ctx->is_keyframe) {
    int has_dc_skip;
    has_dc_skip = !ctx->is_keyframe && !ctx->use_haar_wavelet;
    if (!has_dc_skip || scalar_out[0]) {
      generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli],
       abs(scalar_out[0]) - has_dc_skip, -1,
       &enc->state.adapt.ex_dc[pli][bs][0], 2);
    }
    if (scalar_out[0]) {
      od_ec_enc_bits(&enc->ec, scalar_out[0] < 0, 1);
      skip = 0;
    }
    scalar_out[0] = scalar_out[0]*dc_quant;
    scalar_out[0] += predt[0];
  }
  else {
    scalar_out[0] = dblock[0];
  }
  if (ctx->use_haar_wavelet) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        d[bo + i*w + j] = scalar_out[i*n + j];
      }
    }
  }
  else {
    /*Safely initialize d since some coeffs are skipped by PVQ.*/
    od_init_skipped_coeffs(d, pred, ctx->is_keyframe, bo, n, w);
    od_coding_order_to_raster(&d[bo], w, scalar_out, n);
  }
  /*Apply the inverse transform.*/
#if !defined(OD_OUTPUT_PRED)
  if (ctx->use_haar_wavelet) {
    od_haar_inv(c + bo, w, d + bo, w, bs + 2);
  }
  else {
    od_apply_qm(d + bo, w, d + bo, w, bs, xdec, 1, qm);
    (*enc->state.opt_vtbl.idct_2d[bs])(c + bo, w, d + bo, w);
  }
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
  (*enc->state.opt_vtbl.idct_2d[bs])(c + bo, w, preds, n);
#endif
  /* Allow skipping if it helps the RDO metric, even if the PVQ metric didn't
     skip. */
  if (!skip && has_late_skip_rdo) {
    double lambda;
    double dist_skip;
    double rate_skip;
    int rate_noskip;
    od_coeff *c_noskip;
    c_noskip = enc->block_c_noskip;
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) c_noskip[n*i + j] = c[bo + i*w + j];
    }
    dist_noskip = od_compute_dist(enc, c_orig, c_noskip, n, bs);
    lambda = od_bs_rdo_lambda(enc->state.quantizer[pli]);
    rate_noskip = od_ec_enc_tell_frac(&enc->ec) - tell;
    dist_skip = od_compute_dist(enc, c_orig, mc_orig, n, bs);
    rate_skip = (1 << OD_BITRES)*od_encode_cdf_cost(2,
     enc->state.adapt.skip_cdf[2*bs + (pli != 0)],
     4 + (pli == 0 && bs > 0));
    if (dist_skip + lambda*rate_skip < dist_noskip + lambda*rate_noskip) {
      od_encode_rollback(enc, &pre_encode_buf);
      /* Code the "skip this block" symbol (2). */
      od_encode_cdf_adapt(&enc->ec, 2,
       enc->state.adapt.skip_cdf[2*bs + (pli != 0)], 4 + (pli == 0 && bs > 0),
       enc->state.adapt.skip_increment);
#if OD_SIGNAL_Q_SCALING
      if (bs == (OD_NBSIZES - 1) && pli == 0) {
        od_encode_quantizer_scaling(enc, 0,
         bx >> (OD_NBSIZES - 1), by >> (OD_NBSIZES - 1), 1);
      }
#endif
      skip = 1;
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
          d[bo + i*w + j] = md[bo + i*w + j];
        }
      }
      od_apply_qm(d + bo, w, d + bo, w, bs, xdec, 1, qm);
      (*enc->state.opt_vtbl.idct_2d[bs])(c + bo, w, d + bo, w);
    }
  }
  return skip;
}

static void od_compute_dcts(daala_enc_ctx *enc, od_mb_enc_ctx *ctx, int pli,
  int bx, int by, int bsi, int xdec, int ydec, int use_haar) {
  int obs;
  int bs;
  int w;
  int bo;
  int lossless;
  od_coeff *d;
  const int *qm;
  qm = ctx->qm == OD_HVS_QM ? OD_QM8_Q4_HVS : OD_QM8_Q4_FLAT;
  lossless = OD_LOSSLESS(enc, pli);
  d = ctx->d[pli];
  w = enc->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  obs = OD_BLOCK_SIZE4x4(enc->state.bsize,
   enc->state.bstride, bx << bsi, by << bsi);
  bs = OD_MAXI(obs, xdec);
  OD_ASSERT(bs <= bsi);
  if (bs == bsi) {
    bs -= xdec;
    bo = (by << (OD_LOG_BSIZE0 + bs))*w + (bx << (OD_LOG_BSIZE0 + bs));
    if (use_haar) {
      od_haar(d + bo, w, ctx->c + bo, w, bs + 2);
    }
    else {
      (*enc->state.opt_vtbl.fdct_2d[bs])(d + bo, w, ctx->c + bo, w);
      if (!lossless) od_apply_qm(d + bo, w, d + bo, w, bs, xdec, 0, qm);
    }
  }
  else {
    int f;
    bs = bsi - xdec;
    f = OD_FILT_SIZE(bs - 1, xdec);
    bo = (by << (OD_LOG_BSIZE0 + bs))*w + (bx << (OD_LOG_BSIZE0 + bs));
    od_prefilter_split(ctx->c + bo, w, bs, f);
    bsi--;
    bx <<= 1;
    by <<= 1;
    od_compute_dcts(enc, ctx, pli, bx + 0, by + 0, bsi, xdec, ydec, use_haar);
    od_compute_dcts(enc, ctx, pli, bx + 1, by + 0, bsi, xdec, ydec, use_haar);
    od_compute_dcts(enc, ctx, pli, bx + 0, by + 1, bsi, xdec, ydec, use_haar);
    od_compute_dcts(enc, ctx, pli, bx + 1, by + 1, bsi, xdec, ydec, use_haar);
    if (ctx->is_keyframe) {
      od_coeff x[4];
      int ln;
      ln = bsi - xdec + 2;
      x[0] = d[(by << ln)*w + (bx << ln)];
      x[1] = d[(by << ln)*w + ((bx + 1) << ln)];
      x[2] = d[((by + 1) << ln)*w + (bx << ln)];
      x[3] = d[((by + 1) << ln)*w + ((bx + 1) << ln)];
      OD_HAAR_KERNEL(x[0], x[2], x[1], x[3]);
      d[(by << ln)*w + (bx << ln)] = x[0];
      d[(by << ln)*w + ((bx + 1) << ln)] = x[1];
      d[((by + 1) << ln)*w + (bx << ln)] = x[2];
      d[((by + 1) << ln)*w + ((bx + 1) << ln)] = x[3];
    }
  }
}

#if !OD_DISABLE_HAAR_DC
/** Encode all DCs from a superblock by applying a Haar transform at each level.
 *  This code also computes the distortion and rate to allow block size RDO
 *  on intra. This is made possible because the combination of 4 DC basis
 *  functions of size NxN is exactly equal (outside of rounding error) to the
 *  DC basis function of size 2Nx2N. We also save the quantized DC (dc) and
 *  rates (dc_rate) in recursive raster order.
 * @param [in,out] enc    Daala encoder context
 * @param [in,out] ctx    encoder scratch context
 * @param [in]     pli    plane index
 * @param [in]     bx     horizontal block index (superblock at toplevel)
 * @param [in]     by     vertical block index (superblock at toplevel)
 * @param [in]     l      transform level (starts at 3, going down)
 * @param [in]     xdec   horizontal subsampling
 * @param [in]     ydec   vertical subsampling
 * @param [in]     hgrad  horizontal gradient (ignored for l=3)
 * @param [in]     vgrad  vertical gradient (ignored for l=3)
 * @param [in]     has_ur whether the up-right superblock is available
 * @param [in,out] dc     buffer for recursively saving the quantized DC for
 * RDO purposes
 * @param [in,out] dc_rate buffer for recursively saving the rate of the DC for
 * RDO purposes
 */
static void od_quantize_haar_dc_sb(daala_enc_ctx *enc, od_mb_enc_ctx *ctx,
 int pli, int bx, int by, int xdec, int ydec, int has_ur,
 od_coeff *ohgrad, od_coeff *ovgrad) {
  int w;
  int dc_quant;
  od_coeff *d;
  int nhsb;
  int quant;
  int dc0;
  int ln;
  od_coeff sb_dc_pred;
  od_coeff sb_dc_curr;
  od_coeff *sb_dc_mem;
  (void)ydec;
  d = ctx->d[pli];
  w = enc->state.frame_width >> xdec;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  if (OD_LOSSLESS(enc, pli)) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, enc->state.quantizer[pli]*
     enc->state.pvq_qm_q4[pli][od_qm_get_index(OD_NBSIZES - 1, 0)] >> 4);
  }
  nhsb = enc->state.nhsb;
  sb_dc_mem = enc->state.sb_dc_mem[pli];
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
  dc0 = d[(by << ln)*w + (bx << ln)] - sb_dc_pred;
  quant = OD_DIV_R0(dc0, dc_quant);
  generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli], abs(quant), -1,
   &enc->state.adapt.ex_sb_dc[pli], 2);
  if (quant) od_ec_enc_bits(&enc->ec, quant < 0, 1);
  sb_dc_curr = quant*dc_quant + sb_dc_pred;
  d[(by << ln)*w + (bx << ln)] = sb_dc_curr;
  sb_dc_mem[by*nhsb + bx] = sb_dc_curr;
  if (by > 0) *ovgrad = sb_dc_mem[(by - 1)*nhsb + bx] - sb_dc_curr;
  if (bx > 0) *ohgrad = sb_dc_mem[by*nhsb + bx - 1]- sb_dc_curr;
}
#endif

static void od_quantize_haar_dc_level(daala_enc_ctx *enc, od_mb_enc_ctx *ctx,
  int pli, int bx, int by, int bsi, int xdec, od_coeff *hgrad,
  od_coeff *vgrad) {
  od_coeff x[4];
  int ln;
  int ac_quant[2];
  int i;
  int dc_quant;
  int w;
  w = enc->state.frame_width >> xdec;
  if (OD_LOSSLESS(enc, pli)) dc_quant = 1;
  else {
    dc_quant = OD_MAXI(1, enc->state.quantizer[pli]*
     enc->state.pvq_qm_q4[pli][od_qm_get_index(OD_NBSIZES - 1, 0)] >> 4);
  }
  if (OD_LOSSLESS(enc, pli)) ac_quant[0] = ac_quant[1] = 1;
  else {
    ac_quant[0] = (dc_quant*OD_DC_QM[bsi - xdec][0] + 8) >> 4;
    ac_quant[1] = (dc_quant*OD_DC_QM[bsi - xdec][1] + 8) >> 4;
  }
  ln = bsi - xdec + 2;
  x[0] = ctx->d[pli][(by << ln)*w + (bx << ln)];
  x[1] = ctx->d[pli][(by << ln)*w + ((bx + 1) << ln)];
  x[2] = ctx->d[pli][((by + 1) << ln)*w + (bx << ln)];
  x[3] = ctx->d[pli][((by + 1) << ln)*w + ((bx + 1) << ln)];
  x[1] -= *hgrad/5;
  x[2] -= *vgrad/5;
  for (i = 1; i < 4; i++) {
    int quant;
    int sign;
    double cost;
    int q;
    q = ac_quant[i == 3];
    sign = x[i] < 0;
    x[i] = abs(x[i]);
#if 1 /* Set to zero to disable RDO. */
    quant = x[i]/q;
    cost = generic_encode_cost(&enc->state.adapt.model_dc[pli], quant + 1,
     -1, &enc->state.adapt.ex_dc[pli][bsi][i-1]);
    cost -= generic_encode_cost(&enc->state.adapt.model_dc[pli], quant,
     -1, &enc->state.adapt.ex_dc[pli][bsi][i-1]);
    /* Count cost of sign bit. */
    if (quant == 0) cost += 1;
    if (q*q - 2*q*(x[i] - quant*q) + q*q*OD_PVQ_LAMBDA*cost < 0) quant++;
#else
    quant = OD_DIV_R0(x[i], q);
#endif
    generic_encode(&enc->ec, &enc->state.adapt.model_dc[pli], quant, -1,
     &enc->state.adapt.ex_dc[pli][bsi][i-1], 2);
    if (quant) od_ec_enc_bits(&enc->ec, sign, 1);
    x[i] = quant*ac_quant[i == 3];
    if (sign) x[i] = -x[i];
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

/* Returns 1 if the block is skipped, zero otherwise. */
static int od_encode_recursive(daala_enc_ctx *enc, od_mb_enc_ctx *ctx,
 int pli, int bx, int by, int bsi, int xdec, int ydec, int rdo_only,
 od_coeff hgrad, od_coeff vgrad) {
  int obs;
  int bs;
  int frame_width;
  int w;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  obs = OD_BLOCK_SIZE4x4(enc->state.bsize,
   enc->state.bstride, bx << bsi, by << bsi);
  frame_width = enc->state.frame_width;
  w = frame_width >> xdec;
  bs = OD_MAXI(obs, xdec);
  OD_ASSERT(bs <= bsi);
  if (bs == bsi) {
    int skip;
    int i;
    int j;
    bs -= xdec;
    /*Construct the luma predictors for chroma planes.*/
    if (ctx->l != NULL) {
      OD_ASSERT(pli > 0);
      od_resample_luma_coeffs(ctx->l, 1 << (bs + OD_LOG_BSIZE0),
       ctx->d[0] + (by << (2 + bsi))*frame_width + (bx << (2 + bsi)),
       frame_width, xdec, ydec, bs, obs);
    }
    skip = od_block_encode(enc, ctx, bs, pli, bx, by, rdo_only);
    for (i = 0; i < 1 << bs; i++) {
      for (j = 0; j < 1 << bs; j++) {
        enc->state.bskip[pli][((by << bs) + i)*enc->state.skip_stride
         + (bx << bs) + j] = skip && !ctx->is_keyframe;
      }
    }
    return skip;
  }
  else {
    int f;
    int bo;
    int n;
    int tell;
    int skip_split;
    int skip_nosplit;
    int skip_block;
    od_rollback_buffer pre_encode_buf;
    od_rollback_buffer post_nosplit_buf;
    od_coeff *mc_orig;
    od_coeff *c_orig;
    od_coeff *nosplit;
    od_coeff *split;
    int rate_nosplit;
    int rate_split;
    c_orig = enc->c_orig[bsi - 1];
    mc_orig = enc->mc_orig[bsi - 1];
    nosplit = enc->nosplit[bsi - 1];
    split = enc->split[bsi - 1];
    /* Silence gcc -Wmaybe-uninitialized */
    rate_nosplit = skip_nosplit = 0;
    bs = bsi - xdec;
    bo = (by << (OD_LOG_BSIZE0 + bs))*w + (bx << (OD_LOG_BSIZE0 + bs));
    n = 4 << bs;
    if (rdo_only && bsi <= OD_LIMIT_BSIZE_MAX) {
      int i;
      int j;
      od_coeff dc_orig[(OD_BSIZE_MAX/4)*(OD_BSIZE_MAX/4)];
      tell = od_ec_enc_tell_frac(&enc->ec);
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) c_orig[n*i + j] = ctx->c[bo + i*w + j];
      }
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) mc_orig[n*i + j] = ctx->mc[bo + i*w + j];
      }
      /* Save only the DCs from the transform coeffs. */
      for (i = 0; i < n/4; i++) {
        for (j = 0; j < n/4; j++) {
          dc_orig[n/4*i + j] = ctx->d[pli][bo + 4*i*w + 4*j];
        }
      }
      od_encode_checkpoint(enc, &pre_encode_buf);
      skip_nosplit = od_block_encode(enc, ctx, bs, pli, bx, by, rdo_only);
      rate_nosplit = od_ec_enc_tell_frac(&enc->ec) - tell;
      od_encode_checkpoint(enc, &post_nosplit_buf);
      od_encode_rollback(enc, &pre_encode_buf);
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) nosplit[n*i + j] = ctx->c[bo + i*w + j];
      }
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) ctx->c[bo + i*w + j] = c_orig[n*i + j];
      }
      for (i = 0; i < n/4; i++) {
        for (j = 0; j < n/4; j++) {
          ctx->d[pli][bo + 4*i*w + 4*j] = dc_orig[n/4*i + j];
        }
      }
    }
    f = OD_FILT_SIZE(bs - 1, xdec);
    od_prefilter_split(ctx->c + bo, w, bs, f);
    if (!ctx->is_keyframe) od_prefilter_split(ctx->mc + bo, w, bs, f);
    skip_split = 1;
    if (pli == 0) {
      /* Code the "split this block" symbol (4). */
      od_encode_cdf_adapt(&enc->ec, 4,
       enc->state.adapt.skip_cdf[2*bs + (pli != 0)], 5,
       enc->state.adapt.skip_increment);
#if OD_SIGNAL_Q_SCALING
      if (bs == (OD_NBSIZES - 1)) {
        od_encode_quantizer_scaling(enc, ctx->q_scaling, bx, by, 0);
      }
#endif
    }
    if (ctx->is_keyframe) {
      od_quantize_haar_dc_level(enc, ctx, pli, 2*bx, 2*by, bsi - 1, xdec,
       &hgrad, &vgrad);
    }
    skip_split &= od_encode_recursive(enc, ctx, pli, 2*bx + 0, 2*by + 0,
     bsi - 1, xdec, ydec, rdo_only, hgrad, vgrad);
    skip_split &= od_encode_recursive(enc, ctx, pli, 2*bx + 1, 2*by + 0,
     bsi - 1, xdec, ydec, rdo_only, hgrad, vgrad);
    skip_split &= od_encode_recursive(enc, ctx, pli, 2*bx + 0, 2*by + 1,
     bsi - 1, xdec, ydec, rdo_only, hgrad, vgrad);
    skip_split &= od_encode_recursive(enc, ctx, pli, 2*bx + 1, 2*by + 1,
     bsi - 1, xdec, ydec, rdo_only, hgrad, vgrad);
    skip_block = skip_split;
    od_postfilter_split(ctx->c + bo, w, bs, f, enc->state.coded_quantizer[pli],
     &enc->state.bskip[pli][(by << bs)*enc->state.skip_stride + (bx << bs)],
     enc->state.skip_stride);
    if (rdo_only && bsi <= OD_LIMIT_BSIZE_MAX) {
      int i;
      int j;
      double lambda;
      double dist_split;
      double dist_nosplit;
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) split[n*i + j] = ctx->c[bo + i*w + j];
      }
      rate_split = od_ec_enc_tell_frac(&enc->ec) - tell;
      dist_split = od_compute_dist(enc, c_orig, split, n, bs);
      dist_nosplit = od_compute_dist(enc, c_orig, nosplit, n, bs);
      lambda = od_bs_rdo_lambda(enc->state.quantizer[pli]);
      if (skip_split || dist_nosplit + lambda*rate_nosplit < dist_split
       + lambda*rate_split) {
        /* This rollback call leaves the entropy coder in an inconsistent state
           because the bytes in the buffer are not being copied back. This is
           not a problem here because we are only tracking the rate and we will
           rollback everything at the end of the RDO stage anyway. */
        od_encode_rollback(enc, &post_nosplit_buf);
        for (i = 0; i < n; i++) {
          for (j = 0; j < n; j++) ctx->c[bo + i*w + j] = nosplit[n*i + j];
        }
        for (i = 0; i < 1 << (bs - 1); i++) {
          for (j = 0; j < 1 << (bs - 1); j++) {
            enc->state.bsize[((by << bsi >> 1) + i)*enc->state.bstride
             + (bx << bsi >> 1) + j] = bs;
          }
        }
        for (i = 0; i < 1 << bs; i++) {
          for (j = 0; j < 1 << bs; j++) {
            enc->state.bskip[pli][((by << bs) + i)*enc->state.skip_stride
             + (bx << bs) + j] = skip_nosplit && !ctx->is_keyframe;
          }
        }
        skip_block = skip_nosplit;
#if defined(OD_DUMP_BSIZE_DIST)
        if (bsi == OD_NBSIZES - 2) {
          enc->bsize_dist[pli] += dist_nosplit;
        }
#endif
      }
#if defined(OD_DUMP_BSIZE_DIST)
      else if (bsi == OD_NBSIZES - 2) {
        enc->bsize_dist[pli] += dist_split;
      }
#endif
      for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) ctx->mc[bo + i*w + j] = mc_orig[n*i + j];
      }
    }
    return skip_block && rdo_only;
  }
}

static void od_encode_mv(daala_enc_ctx *enc, int num_refs, od_mv_grid_pt *mvg,
 int vx, int vy, int level, int mv_res, int mv_range_x, int mv_range_y) {
  generic_encoder *model;
  int pred[2];
  int ox;
  int oy;
  int id;
  int equal_mvs;
  int ref_pred;
  if (num_refs > 1) {
    /* Code reference index. */
    ref_pred = od_mc_get_ref_predictor(&enc->state, vx, vy, level);
    OD_ASSERT(ref_pred >= 0);
    OD_ASSERT(ref_pred < num_refs);
    OD_ASSERT(mvg->ref < num_refs);
    od_encode_cdf_adapt(&enc->ec, mvg->ref,
     enc->state.adapt.mv_ref_cdf[ref_pred], num_refs, 256);
  }
  equal_mvs = od_state_get_predictor(&enc->state, pred, vx, vy, level,
   mv_res, mvg->ref);
  ox = (mvg->mv[0] >> mv_res) - pred[0];
  oy = (mvg->mv[1] >> mv_res) - pred[1];
  /*Interleave positive and negative values.*/
  model = &enc->state.adapt.mv_model;
  id = OD_MINI(abs(oy), 3)*4 + OD_MINI(abs(ox), 3);
  od_encode_cdf_adapt(&enc->ec, id, enc->state.adapt.mv_small_cdf[equal_mvs],
   16, enc->state.adapt.mv_small_increment);
  if (abs(ox) >= 3) {
    generic_encode(&enc->ec, model, abs(ox) - 3, mv_range_x,
     &enc->state.adapt.mv_ex[level], 6);
  }
  if (abs(oy) >= 3) {
    generic_encode(&enc->ec, model, abs(oy) - 3, mv_range_y,
     &enc->state.adapt.mv_ey[level], 6);
  }
  if (abs(ox)) od_ec_enc_bits(&enc->ec, ox < 0, 1);
  if (abs(oy)) od_ec_enc_bits(&enc->ec, oy < 0, 1);
}

static void od_img_copy_pad(daala_enc_ctx *enc, od_img *img) {
  od_state *state;
  int pli;
  int nplanes;
  state = &enc->state;
  nplanes = img->nplanes;
  /* Copy and pad the image. */
  for (pli = 0; pli < nplanes; pli++) {
    int plane_width;
    int plane_height;
    int xdec;
    int ydec;
    xdec = img->planes[pli].xdec;
    ydec = img->planes[pli].ydec;
    plane_width = ((state->info.pic_width + (1 << xdec) - 1) >> xdec);
    plane_height = ((state->info.pic_height + (1 << ydec) - 1) >> ydec);
    od_img_plane_copy_pad(&enc->input_img,
     state->frame_width >> xdec, state->frame_height >> ydec,
     img, plane_width, plane_height, pli);
  }
  od_img_edge_ext(&enc->input_img);
}

#if defined(OD_DUMP_IMAGES)
/*State visualization.
  None of this is particularly fast.*/

/*static const unsigned char OD_YCbCr_BORDER[3] = {49, 109, 184};*/
static const unsigned char OD_YCbCr_BORDER[3] = {113, 72, 137};
static const unsigned char OD_YCbCr_EDGE[3] = {41, 240, 110};
static const unsigned char OD_YCbCr_MV[3] = {81, 90, 240};

#define OD_SRCVAL_TO_8(x) (siplane->bitdepth == 8 ? src[(x)] : \
 (((int16_t *)src)[(x)] + (1 << siplane->bitdepth - 9) \
 >> siplane->bitdepth - 8))

/*Upsamples the reconstructed image (any depth) to an 8-bit reference image.
  Output img is 2x larger in both horizontal and vertical dimensions.
  TODO: Pipeline with reconstruction.*/
static void od_img_upsample8(daala_enc_ctx *enc, od_img *dimg,
 const od_img *simg) {
  int pli;
  for (pli = 0; pli < simg->nplanes; pli++) {
    const od_img_plane *siplane;
    od_img_plane *diplane;
    const unsigned char *src;
    unsigned char *dst;
    int xpad;
    int ypad;
    int w;
    int h;
    int x;
    int y;
    siplane = simg->planes + pli;
    diplane = dimg->planes + pli;
    xpad = OD_BUFFER_PADDING >> siplane->xdec;
    ypad = OD_BUFFER_PADDING >> siplane->ydec;
    w = simg->width >> siplane->xdec;
    h = simg->height >> siplane->ydec;
    src = siplane->data;
    dst = diplane->data - (diplane->ystride << 1)*ypad;
    for (y = -ypad; y < h + ypad + 3; y++) {
      /*Horizontal filtering:*/
      if (y < h + ypad) {
        unsigned char *buf;
        buf = enc->upsample_line_buf[y & 7];
        memset(buf - (xpad << 1), OD_SRCVAL_TO_8(0), (xpad - 2) << 1);
        /*for (x = -xpad; x < -2; x++) {
          *(buf + (x << 1)) = OD_SRCVAL_TO_8(0);
          *(buf + (x << 1 | 1)) = OD_SRCVAL_TO_8(0);
        }*/
        *(buf - 4) = OD_SRCVAL_TO_8(0);
        *(buf - 3) = OD_CLAMP255((31*OD_SRCVAL_TO_8(0) + OD_SRCVAL_TO_8(1)
         + 16) >> 5);
        *(buf - 2) = OD_SRCVAL_TO_8(0);
        *(buf - 1) = OD_CLAMP255((36*OD_SRCVAL_TO_8(0) - 5*OD_SRCVAL_TO_8(1)
         + OD_SRCVAL_TO_8(1) + 16) >> 5);
        buf[0] = OD_SRCVAL_TO_8(0);
        buf[1] = OD_CLAMP255((20*(OD_SRCVAL_TO_8(0) + OD_SRCVAL_TO_8(1))
         - 5*(OD_SRCVAL_TO_8(0) + OD_SRCVAL_TO_8(2)) + OD_SRCVAL_TO_8(0)
         + OD_SRCVAL_TO_8(3) + 16) >> 5);
        buf[2] = OD_SRCVAL_TO_8(1);
        buf[3] = OD_CLAMP255((20*(OD_SRCVAL_TO_8(1) + OD_SRCVAL_TO_8(2))
         - 5*(OD_SRCVAL_TO_8(0) + OD_SRCVAL_TO_8(3)) + OD_SRCVAL_TO_8(0)
         + OD_SRCVAL_TO_8(4) + 16) >> 5);
        for (x = 2; x < w - 3; x++) {
          buf[x << 1] = OD_SRCVAL_TO_8(x);
          buf[x << 1 | 1] = OD_CLAMP255((20*(OD_SRCVAL_TO_8(x)
           + OD_SRCVAL_TO_8(x + 1)) - 5*(OD_SRCVAL_TO_8(x - 1)
           + OD_SRCVAL_TO_8(x + 2)) + OD_SRCVAL_TO_8(x - 2)
           + OD_SRCVAL_TO_8(x + 3) + 16) >> 5);
        }
        buf[x << 1] = OD_SRCVAL_TO_8(x);
        buf[x << 1 | 1] = OD_CLAMP255((20*(OD_SRCVAL_TO_8(x)
         + OD_SRCVAL_TO_8(x + 1)) - 5*(OD_SRCVAL_TO_8(x - 1)
         + OD_SRCVAL_TO_8(x + 2)) + OD_SRCVAL_TO_8(x - 2)
         + OD_SRCVAL_TO_8(x + 2) + 16) >> 5);
        x++;
        buf[x << 1] = OD_SRCVAL_TO_8(x);
        buf[x << 1 | 1] = OD_CLAMP255((20*(OD_SRCVAL_TO_8(x)
         + OD_SRCVAL_TO_8(x + 1)) - 5*(OD_SRCVAL_TO_8(x - 1)
         + OD_SRCVAL_TO_8(x + 1)) + OD_SRCVAL_TO_8(x - 2)
         + OD_SRCVAL_TO_8(x + 1) + 16) >> 5);
        x++;
        buf[x << 1] = OD_SRCVAL_TO_8(x);
        buf[x << 1 | 1] =
         OD_CLAMP255((36*OD_SRCVAL_TO_8(x) - 5*OD_SRCVAL_TO_8(x - 1)
          + OD_SRCVAL_TO_8(x - 2) + 16) >> 5);
        x++;
        buf[x << 1] = OD_SRCVAL_TO_8(w - 1);
        buf[x << 1 | 1] = OD_CLAMP255((31*OD_SRCVAL_TO_8(w - 1)
         + OD_SRCVAL_TO_8(w - 2) + 16) >> 5);
        memset(buf + (++x << 1), OD_SRCVAL_TO_8(w - 1), (xpad - 1) << 1);
        /*for (x++; x < w + xpad; x++) {
          buf[x << 1] = OD_SRCVAL_TO_8(w - 1);
          buf[x << 1 | 1]=OD_SRCVAL_TO_8(w - 1);
        }*/
        if (y >= 0 && y + 1 < h) src += siplane->ystride;
      }
      /*Vertical filtering:*/
      if (y >= -ypad + 3) {
        if (y < 1 || y > h + 3) {
          OD_COPY(dst - (xpad << 1),
           enc->upsample_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          OD_COPY(dst - (xpad << 1),
           enc->upsample_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
        else {
          unsigned char *buf[6];
          buf[0] = enc->upsample_line_buf[(y - 5) & 7];
          buf[1] = enc->upsample_line_buf[(y - 4) & 7];
          buf[2] = enc->upsample_line_buf[(y - 3) & 7];
          buf[3] = enc->upsample_line_buf[(y - 2) & 7];
          buf[4] = enc->upsample_line_buf[(y - 1) & 7];
          buf[5] = enc->upsample_line_buf[(y - 0) & 7];
          OD_COPY(dst - (xpad << 1),
           enc->upsample_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            *(dst + x) = OD_CLAMP255((20*(*(buf[2] + x) + *(buf[3] + x))
             - 5*(*(buf[1] + x) + *(buf[4] + x))
             + *(buf[0] + x) + *(buf[5] + x) + 16) >> 5);
          }
          /*fprintf(stderr, "%3i: ", (y - 3) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
      }
    }
  }
}

static void od_img_draw_point(od_img *img, int x, int y,
 const unsigned char ycbcr[3]) {
  if (x < 0 || y < 0 || x >= img->width || y >= img->height) return;
  *(img->planes[0].data + img->planes[0].ystride*(y >> img->planes[0].ydec)
   + (x >> img->planes[0].xdec)) = ycbcr[0];
  if (img->nplanes >= 3) {
    *(img->planes[1].data + img->planes[1].ystride*(y >> img->planes[1].ydec)
     + (x >> img->planes[1].xdec)) = ycbcr[1];
    *(img->planes[2].data + img->planes[2].ystride*(y >> img->planes[2].ydec)
     + (x >> img->planes[2].xdec)) = ycbcr[2];
  }
}

void od_img_draw_line(od_img *img,
 int x0, int y0, int x1, int y1, const unsigned char ycbcr[3]) {
  int p0[2];
  int p1[2];
  int dx[2];
  int step[2];
  int steep;
  int err;
  int derr;
  steep = abs(y1 - y0) > abs(x1 - x0);
  p0[0] = x0;
  p0[1] = y0;
  p1[0] = x1;
  p1[1] = y1;
  dx[0] = abs(p1[0] - p0[0]);
  dx[1] = abs(p1[1] - p0[1]);
  err = 0;
  derr = dx[1 - steep];
  step[0] = ((p0[0] < p1[0]) << 1) - 1;
  step[1] = ((p0[1] < p1[1]) << 1) - 1;
  od_img_draw_point(img, p0[0], p0[1], ycbcr);
  while (p0[steep] != p1[steep]) {
    p0[steep] += step[steep];
    err += derr;
    if (err << 1 > dx[steep]) {
      p0[1 - steep] += step[1 - steep];
      err -= dx[steep];
    }
    od_img_draw_point(img, p0[0], p0[1], ycbcr);
  }
}

static void od_state_draw_mv_grid_block(daala_enc_ctx *enc,
 int vx, int vy, int log_mvb_sz) {
  od_state *state;
  int half_mvb_sz;
  state = &enc->state;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_state_draw_mv_grid_block(enc, vx, vy, log_mvb_sz - 1);
    od_state_draw_mv_grid_block(enc, vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_state_draw_mv_grid_block(enc, vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_state_draw_mv_grid_block(enc, vx + half_mvb_sz, vy + half_mvb_sz,
     log_mvb_sz - 1);
  }
  else {
    int mvb_sz;
    int x0;
    int y0;
    mvb_sz = 1 << log_mvb_sz;
    x0 = (vx << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
    y0 = (vy << (OD_LOG_MVBSIZE_MIN + 1)) + (OD_BUFFER_PADDING << 1);
    od_img_draw_line(&enc->vis_img, x0, y0,
     x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), y0, OD_YCbCr_EDGE);
    od_img_draw_line(&enc->vis_img,
     x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), y0,
     x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
     y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), OD_YCbCr_EDGE);
    od_img_draw_line(&enc->vis_img,
     x0, y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
     x0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)),
     y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), OD_YCbCr_EDGE);
    od_img_draw_line(&enc->vis_img, x0, y0,
     x0, y0 + (mvb_sz << (OD_LOG_MVBSIZE_MIN + 1)), OD_YCbCr_EDGE);
  }
}

static void od_state_draw_mv_grid(daala_enc_ctx *enc) {
  od_state *state;
  int vx;
  int vy;
  int nhmvbs;
  int nvmvbs;
  state = &enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy < nvmvbs; vy += OD_MVB_DELTA0) {
    for (vx = 0; vx < nhmvbs; vx += OD_MVB_DELTA0) {
      od_state_draw_mv_grid_block(enc, vx, vy, OD_LOG_MVB_DELTA0);
    }
  }
}

static void od_state_draw_mvs_block(daala_enc_ctx *enc,
 int vx, int vy, int log_mvb_sz) {
  od_state *state;
  int half_mvb_sz;
  state = &enc->state;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_state_draw_mvs_block(enc, vx, vy, log_mvb_sz - 1);
    od_state_draw_mvs_block(enc, vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_state_draw_mvs_block(enc, vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_state_draw_mvs_block(enc, vx + half_mvb_sz, vy + half_mvb_sz,
     log_mvb_sz - 1);
  }
  else {
    od_mv_grid_pt *grid[4];
    const int *dxp;
    const int *dyp;
    int x0;
    int y0;
    int k;
    int oc;
    int s;
    if (log_mvb_sz < OD_LOG_MVB_DELTA0) {
      int mask;
      int s1vx;
      int s1vy;
      int s3vx;
      int s3vy;
      mask = (1 << (log_mvb_sz + 1)) - 1;
      oc = !!(vx & mask);
      if (vy & mask) oc = 3 - oc;
      s1vx = vx + (OD_VERT_DX[(oc + 1) & 3] << log_mvb_sz);
      s1vy = vy + (OD_VERT_DY[(oc + 1) & 3] << log_mvb_sz);
      s3vx = vx + (OD_VERT_DX[(oc + 3) & 3] << log_mvb_sz);
      s3vy = vy + (OD_VERT_DY[(oc + 3) & 3] << log_mvb_sz);
      s = state->mv_grid[s1vy][s1vx].valid |
       state->mv_grid[s3vy][s3vx].valid << 1;
    }
    else {
      oc = 0;
      s = 3;
    }
    dxp = OD_VERT_SETUP_DX[oc][s];
    dyp = OD_VERT_SETUP_DY[oc][s];
    for (k = 0; k < 4; k++) {
      grid[k] = state->mv_grid[vy + (dyp[k] << log_mvb_sz)]
       + vx + (dxp[k] << log_mvb_sz);
    }
    for (k = 0; k < 4; k++) {
      x0 = ((vx + (dxp[k] << log_mvb_sz)) << (OD_LOG_MVBSIZE_MIN + 1))
       + (OD_BUFFER_PADDING << 1);
      y0 = ((vy + (dyp[k] << log_mvb_sz)) << (OD_LOG_MVBSIZE_MIN + 1))
       + (OD_BUFFER_PADDING << 1);
      /*od_img_draw_point(&enc->vis_img, x0, y0, OD_YCbCr_MV);*/
      od_img_draw_line(&enc->vis_img, x0, y0,
       x0 + OD_DIV_ROUND_POW2(grid[k]->mv[0], 2, 2),
       y0 + OD_DIV_ROUND_POW2(grid[k]->mv[1], 2, 2), OD_YCbCr_MV);
    }
  }
}

void od_state_draw_mvs(daala_enc_ctx *enc) {
  od_state *state;
  int vx;
  int vy;
  int nhmvbs;
  int nvmvbs;
  state = &enc->state;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy < nvmvbs; vy += OD_MVB_DELTA0) {
    for (vx = 0; vx < nhmvbs; vx += OD_MVB_DELTA0) {
      od_state_draw_mvs_block(enc, vx, vy, OD_LOG_MVB_DELTA0);
    }
  }
}

void od_encode_fill_vis(od_enc_ctx *enc) {
  od_state *state;
  od_img *img;
  od_img *ref_img;
  int pli;
  int xdec;
  int ydec;
  int border;
  int x;
  int y;
  state = &enc->state;
  img = &enc->vis_img;
  /*Upsample the reconstructed image for better quality.*/
  /*Adjust the data pointers so that the padding works like the reference
     images.*/
  border = OD_BUFFER_PADDING << 1;
  for (pli = 0; pli < img->nplanes; pli++) {
    img->planes[pli].data += (border >> img->planes[pli].xdec)
     + img->planes[pli].ystride*(border >> img->planes[pli].ydec);
  }
  od_img_upsample8(enc, img,
   state->ref_imgs + state->ref_imgi[OD_FRAME_SELF]);
  /*Upsample the input image, as well, and subtract it to get a difference
     image.*/
  ref_img = &enc->tmp_vis_img;
  od_img_upsample8(enc, ref_img, &enc->input_img);
  xdec = state->info.plane_info[0].xdec;
  ydec = state->info.plane_info[0].ydec;
  for (y = 0; y < ref_img->height; y++) {
    for (x = 0; x < ref_img->width; x++) {
      int diff;
      int px;
      int py;
      px = x >> xdec;
      py = y >> ydec;
      diff = *(img->planes[0].data + img->planes[0].ystride*py + px)
       - *(ref_img->planes[0].data + ref_img->planes[0].ystride*py + px);
      /*Scale the differences by 2 to make them visible.*/
      diff = OD_CLAMP255((diff << 1) + 128);
      *(img->planes[0].data + img->planes[0].ystride*py + px) =
       (unsigned char)diff;
    }
  }
  /*Undo the adjustment.*/
  for (pli = 0; pli < img->nplanes; pli++) {
    img->planes[pli].data -= (border >> img->planes[pli].xdec) +
     img->planes[pli].ystride*(border >> img->planes[pli].ydec);
  }
  /*Clear the border region.*/
  for (y = 0; y < (border >> ydec); y++) {
    OD_CLEAR(img->planes[0].data + (img->planes[0].ystride)*y,
     img->width >> xdec);
  }
  for (; y < (img->height - border) >> ydec; y++) {
    OD_CLEAR(img->planes[0].data + img->planes[0].ystride*y, border >> xdec);
    OD_CLEAR(img->planes[0].data + img->planes[0].ystride*y
     + ((img->width - border) >> xdec), border >> xdec);
  }
  for (; y < img->height >> ydec; y++) {
    OD_CLEAR(img->planes[0].data + (img->planes[0].ystride)*y,
     img->width >> xdec);
  }
  /*Clear the chroma planes.*/
  for (pli = 1; pli < img->nplanes; pli++) {
    memset(img->planes[pli].data, 128, (img->height >> img->planes[pli].ydec)*
     (img->width >> img->planes[pli].xdec));
  }
  od_img_draw_line(img, border - 1, border - 1,
   img->width - border, border - 1, OD_YCbCr_BORDER);
  od_img_draw_line(img, border - 2, border - 2,
   img->width - border + 1, border - 2, OD_YCbCr_BORDER);
  od_img_draw_line(img, img->width - border, border - 1,
   img->width - border, img->height - border, OD_YCbCr_BORDER);
  od_img_draw_line(img, img->width - border + 1, border - 2,
   img->width - border + 1, img->height - border + 1, OD_YCbCr_BORDER);
  od_img_draw_line(img, border - 1, img->height - border,
   img->width - border, img->height - border, OD_YCbCr_BORDER);
  od_img_draw_line(img, border - 2, img->height - border + 1,
   img->width - border + 1, img->height - border + 1, OD_YCbCr_BORDER);
  od_img_draw_line(img, border - 1, border - 1,
   border - 1, img->height - border, OD_YCbCr_BORDER);
  od_img_draw_line(img, border - 2, border - 2,
   border - 2, img->height - border + 1, OD_YCbCr_BORDER);
  od_state_draw_mv_grid(enc);
  od_state_draw_mvs(enc);
}

static void od_img_dump_padded(daala_enc_ctx *enc) {
  od_state *state;
  daala_info *info;
  od_img img;
  int nplanes;
  int pli;
  state = &enc->state;
  info = &state->info;
  nplanes = info->nplanes;
  /*Modify the image offsets to include the padding.*/
  *&img = enc->input_img;
  for (pli = 0; pli < nplanes; pli++) {
    img.planes[pli].data -=
     img.planes[pli].ystride*(OD_BUFFER_PADDING>>info->plane_info[pli].ydec)
     + img.planes[pli].xstride*(OD_BUFFER_PADDING>>info->plane_info[pli].xdec);
  }
  img.width += OD_BUFFER_PADDING<<1;
  img.height += OD_BUFFER_PADDING<<1;
  od_state_dump_img(state, &img, "pad");
}
#endif

static void od_predict_frame(daala_enc_ctx *enc) {
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  enc->ani_iter = 0;
#endif
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "Predicting frame %i:",
   (int)daala_granule_basetime(enc, enc->state.cur_time)));
  /*4640000 ~= 0.55313 (or sqrt(0.30595)) in Q23.
   The lower bound of 40 is there because we do not yet consider PVQ noref
    flags during the motion search, so we waste far too many bits trying to
    predict unpredictable areas when lambda is too small.
   Hopefully when we fix that, we can remove the limit.*/
  od_mv_est(enc->mvest,
   OD_MAXI(((4640000 + (((1 << OD_COEFF_SHIFT) - 1) >> 1)) >> OD_COEFF_SHIFT)*
   enc->state.quantizer[0] >> (23 - OD_LAMBDA_SCALE), 40));
  od_state_mc_predict(&enc->state,
   enc->state.ref_imgs + enc->state.ref_imgi[OD_FRAME_SELF]);
  /*Do edge extension here because the block-size analysis needs to read
    outside the frame, but otherwise isn't read from.*/
  od_img_edge_ext(enc->state.ref_imgs + enc->state.ref_imgi[OD_FRAME_SELF]);
#if defined(OD_DUMP_IMAGES)
  od_encode_fill_vis(enc);
  od_state_dump_img(&enc->state, &enc->vis_img, "vis");
#endif
}

static void od_split_superblocks(daala_enc_ctx *enc, int is_keyframe) {
  int nhsb;
  int nvsb;
  int i;
  int j;
  int k;
  int m;
  od_state *state;
  state = &enc->state;
  OD_ASSERT(enc->input_img.planes[0].xstride == 1);
  OD_ASSERT(state->ref_imgs[state->ref_imgi[OD_FRAME_SELF]].
   planes[0].xstride == 1);
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  od_state_init_border(state);
  /* Allocate a blockSizeComp for scratch space and then calculate the block
     sizes eventually store them in bsize. */
  od_log_matrix_uchar(OD_LOG_GENERIC, OD_LOG_INFO, "bimg ",
   enc->input_img.planes[0].data - 16*enc->input_img.planes[0].ystride - 16,
   enc->input_img.planes[0].ystride, (nvsb + 1)*32);
  for (i = 0; i < nvsb; i++) {
    unsigned char *bimg;
    unsigned char *rimg;
    int istride;
    int rstride;
    int bstride;
    bstride = state->bstride;
    istride = enc->input_img.planes[0].ystride;
    rstride = is_keyframe ? 0 :
     state->ref_imgs[state->ref_imgi[OD_FRAME_SELF]].planes[0].ystride;
    bimg = enc->input_img.planes[0].data
     + i*istride*OD_BSIZE_MAX;
    rimg = state->ref_imgs[state->ref_imgi[OD_FRAME_SELF]].planes[0].data
     + i*rstride*OD_BSIZE_MAX;
    for (j = 0; j < nhsb; j++) {
      int bsize[4][4];
      unsigned char *state_bsize;
      state_bsize = &state->bsize[i*4*state->bstride + j*4];
      od_split_superblock(enc->bs, bimg + j*OD_BSIZE_MAX, istride,
       is_keyframe ? NULL : rimg + j*OD_BSIZE_MAX, rstride, bsize,
       state->quantizer[0]);
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
    }
  }
}

static void od_encode_mvs(daala_enc_ctx *enc, int num_refs) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  od_img *mvimg;
  int width;
  int height;
  int mv_res;
  int log_mvb_sz;
  int level;
  od_mv_grid_pt *mvp;
  od_mv_grid_pt **grid;
  uint16_t *cdf;
  state = &enc->state;
  nhmvbs = enc->state.nhmvbs;
  nvmvbs = enc->state.nvmvbs;
  mvimg = state->ref_imgs + state->ref_imgi[OD_FRAME_SELF];
  mv_res = enc->state.mv_res;
  OD_ASSERT(0 <= mv_res && mv_res < 3);
  od_ec_enc_uint(&enc->ec, mv_res, 3);
  width = (mvimg->width + 32) << ((3 - mv_res) + 1); /* delta mvx range */
  height = (mvimg->height + 32) << ((3 - mv_res) + 1);/* delta mvy range */
  grid = enc->state.mv_grid;
  /*Code the motion vectors and flags. At each level, the MVs are zero
    outside of the frame, so don't code them.*/
  /*Level 0.*/
  for (vy = 0; vy <= nvmvbs; vy += OD_MVB_DELTA0) {
    for (vx = 0; vx <= nhmvbs; vx += OD_MVB_DELTA0) {
      mvp = grid[vy] + vx;
      od_encode_mv(enc, num_refs, mvp, vx, vy, 0, mv_res, width, height);
    }
  }
  /*od_ec_acct_add_label(&enc->ec.acct, "mvf-l1");
    od_ec_acct_add_label(&enc->ec.acct, "mvf-l2");
    od_ec_acct_add_label(&enc->ec.acct, "mvf-l3");
    od_ec_acct_add_label(&enc->ec.acct, "mvf-l4");*/
  for (log_mvb_sz = OD_LOG_MVB_DELTA0, level = 1; log_mvb_sz-- > 0; level++) {
    int mvb_sz;
    mvb_sz = 1 << log_mvb_sz;
    /*Odd levels.*/
    for (vy = mvb_sz; vy <= nvmvbs; vy += 2*mvb_sz) {
      for (vx = mvb_sz; vx <= nhmvbs; vx += 2*mvb_sz) {
        mvp = grid[vy] + vx;
        if (grid[vy - mvb_sz][vx - mvb_sz].valid
         && grid[vy - mvb_sz][vx + mvb_sz].valid
         && grid[vy + mvb_sz][vx + mvb_sz].valid
         && grid[vy + mvb_sz][vx - mvb_sz].valid) {
          cdf = od_mv_split_flag_cdf(&enc->state, vx, vy, level);
          od_encode_cdf_adapt(&enc->ec, mvp->valid,
           cdf, 2, enc->state.adapt.split_flag_increment);
          if (mvp->valid) {
            od_encode_mv(enc, num_refs, mvp, vx, vy, level, mv_res,
             width, height);
          }
        }
        else {
          OD_ASSERT(!mvp->valid);
        }
      }
    }
    level++;
    /*Even levels.*/
    for (vy = 0; vy <= nvmvbs; vy += mvb_sz) {
      for (vx = mvb_sz*!(vy & mvb_sz); vx <= nhmvbs; vx += 2*mvb_sz) {
        mvp = grid[vy] + vx;
        if ((vy - mvb_sz < 0 || grid[vy - mvb_sz][vx].valid)
         && (vx - mvb_sz < 0 || grid[vy][vx - mvb_sz].valid)
         && (vy + mvb_sz > nvmvbs || grid[vy + mvb_sz][vx].valid)
         && (vx + mvb_sz > nhmvbs || grid[vy][vx + mvb_sz].valid)) {
          cdf = od_mv_split_flag_cdf(&enc->state, vx, vy, level);
          od_encode_cdf_adapt(&enc->ec, mvp->valid,
           cdf, 2, enc->state.adapt.split_flag_increment);
          if (mvp->valid) {
            od_encode_mv(enc, num_refs, mvp, vx, vy, level, mv_res,
             width, height);
          }
        }
        else {
          OD_ASSERT(!mvp->valid);
        }
      }
    }
  }
}

/* FIXME: add a real scaling calculation */
static int od_compute_superblock_q_scaling(daala_enc_ctx *enc, od_coeff *x,
 int stride) {
  (void)enc;
  (void)x;
  (void)stride;
  return 0;
}

#define OD_ENCODE_REAL (0)
#define OD_ENCODE_RDO (1)
static void od_encode_coefficients(daala_enc_ctx *enc, od_mb_enc_ctx *mbctx,
 int rdo_only) {
  int xdec;
  int ydec;
  int sby;
  int sbx;
  int w;
  int y;
  int x;
  int pli;
  int nplanes;
  int frame_width;
  int nhsb;
  int nvsb;
  od_state *state;
  od_img *rec;
  state = &enc->state;
  nplanes = state->info.nplanes;
  if (rdo_only) nplanes = 1;
  frame_width = state->frame_width;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  rec = state->ref_imgs + state->ref_imgi[OD_FRAME_SELF];
  for (pli = 0; pli < nplanes; pli++) {
    od_ec_enc_uint(&enc->ec, state->coded_quantizer[pli],
     OD_N_CODED_QUANTIZERS);
  }
  for (pli = 0; pli < nplanes; pli++) {
    /*Collect the image data needed for this plane.*/
    xdec = enc->input_img.planes[pli].xdec;
    ydec = enc->input_img.planes[pli].ydec;
    w = frame_width >> xdec;
    od_ref_plane_to_coeff(state, state->ctmp[pli],
     OD_LOSSLESS(enc, pli), &enc->input_img, pli);
    if (!mbctx->use_haar_wavelet) {
      od_apply_prefilter_frame_sbs(state->ctmp[pli],
       w, nhsb, nvsb, xdec, ydec);
    }
    if (!mbctx->is_keyframe) {
      od_ref_plane_to_coeff(state,
       state->mctmp[pli], OD_LOSSLESS(enc, pli), rec, pli);
      if (!mbctx->use_haar_wavelet) {
        od_apply_prefilter_frame_sbs(state->mctmp[pli], w, nhsb, nvsb, xdec,
         ydec);
      }
    }
  }
  for (sby = 0; sby < nvsb; sby++) {
    for (sbx = 0; sbx < nhsb; sbx++) {
      for (pli = 0; pli < nplanes; pli++) {
        od_coeff *c_orig;
        int i;
        int j;
        int width;
        od_rollback_buffer buf;
        od_coeff hgrad;
        od_coeff vgrad;
        width = enc->state.frame_width;
        hgrad = vgrad = 0;
        c_orig = enc->c_orig[0];
        mbctx->c = state->ctmp[pli];
        mbctx->d = state->dtmp;
        mbctx->mc = state->mctmp[pli];
        mbctx->md = state->mdtmp[pli];
        mbctx->l = state->lbuf[pli];
        xdec = enc->input_img.planes[pli].xdec;
        ydec = enc->input_img.planes[pli].ydec;
        if (pli == 0 || (rdo_only && mbctx->is_keyframe)) {
          for (i = 0; i < OD_BSIZE_MAX; i++) {
            for (j = 0; j < OD_BSIZE_MAX; j++) {
              c_orig[i*OD_BSIZE_MAX + j] =
               mbctx->c[(OD_BSIZE_MAX*sby + i)*width + OD_BSIZE_MAX*sbx + j];
            }
          }
        }
        if (mbctx->is_keyframe) {
          if (rdo_only) {
            od_encode_checkpoint(enc, &buf);
          }
          od_compute_dcts(enc, mbctx, pli, sbx, sby, OD_NBSIZES - 1, xdec,
           ydec, mbctx->use_haar_wavelet && !rdo_only);
          od_quantize_haar_dc_sb(enc, mbctx, pli, sbx, sby, xdec, ydec,
           sby > 0 && sbx < nhsb - 1, &hgrad, &vgrad);
          if (rdo_only) {
            od_encode_rollback(enc, &buf);
            for (i = 0; i < OD_BSIZE_MAX; i++) {
              for (j = 0; j < OD_BSIZE_MAX; j++) {
                mbctx->c[(OD_BSIZE_MAX*sby + i)*width + OD_BSIZE_MAX*sbx + j] =
                 c_orig[i*OD_BSIZE_MAX + j];
              }
            }
          }
        }
        if (pli == 0 && state->quantizer[pli] != 0) {
          mbctx->q_scaling =
           od_compute_superblock_q_scaling(enc, c_orig, OD_BSIZE_MAX);
        }
        od_encode_recursive(enc, mbctx, pli, sbx, sby, OD_NBSIZES - 1, xdec,
         ydec, rdo_only, hgrad, vgrad);
      }
    }
  }
#if defined(OD_DUMP_IMAGES)
  if (!rdo_only) {
    /*Dump the lapped frame (before the postfilter has been applied)*/
    for (pli = 0; pli < nplanes; pli++) {
      od_coeff_to_ref_plane(state, rec, pli,
       state->ctmp[pli], OD_LOSSLESS(enc, pli));
    }
    od_state_dump_img(&enc->state, rec, "lapped");
  }
#endif
  for (pli = 0; pli < nplanes; pli++) {
    xdec = enc->input_img.planes[pli].xdec;
    ydec = enc->input_img.planes[pli].ydec;
    w = frame_width >> xdec;
    if (!mbctx->use_haar_wavelet) {
      od_apply_postfilter_frame_sbs(state->ctmp[pli], w, nhsb, nvsb, xdec,
       ydec, state->coded_quantizer[pli], &enc->state.bskip[pli][0],
       enc->state.skip_stride);
    }
  }
  if (!rdo_only && state->quantizer[0] > 0) {
    int nhdr;
    int nvdr;
    /* We copy ctmp to dtmp so we can use it as an unmodified input
       and avoid filtering some pixels twice. */
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      int size;
      xdec = enc->input_img.planes[pli].xdec;
      ydec = enc->input_img.planes[pli].ydec;
      size = nvsb*nhsb*OD_BSIZE_MAX*OD_BSIZE_MAX >> xdec >> ydec;
      for (i = 0; i < size; i++) {
        state->etmp[pli][i] = state->ctmp[pli][i];
      }
    }
    nhdr = state->frame_width >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
    nvdr = state->frame_height >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
    for (sby = 0; sby < nvdr; sby++) {
      for (sbx = 0; sbx < nhdr; sbx++) {
        int ln;
        int n;
        int16_t buf[OD_BSIZE_MAX*OD_BSIZE_MAX];
        od_coeff orig[OD_BSIZE_MAX*OD_BSIZE_MAX];
        double unfiltered_error;
        double filtered_error;
        int ystride;
        int xstride;
        unsigned char *input;
        od_coeff *output;
        int filtered;
        int up;
        int left;
        int c;
        int q2;
        double filtered_rate;
        double unfiltered_rate;
        int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
        int i;
        int j;
        unsigned char *bskip;
        state->dering_flags[sby*nhdr + sbx] = 0;
        bskip = enc->state.bskip[0] +
         (sby << OD_LOG_DERING_GRID)*enc->state.skip_stride +
         (sbx << OD_LOG_DERING_GRID);
        for (j = 0; j < 1 << OD_LOG_DERING_GRID; j++) {
          for (i = 0; i < 1 << OD_LOG_DERING_GRID; i++) {
            if (!bskip[j*enc->state.skip_stride + i]) {
              state->dering_flags[sby*nhdr + sbx] = 1;
            }
          }
        }
        if (!state->dering_flags[sby*nhdr + sbx]) {
          continue;
        }
        pli = 0;
        xdec = enc->input_img.planes[pli].xdec;
        ydec = enc->input_img.planes[pli].ydec;
        w = frame_width >> xdec;
        OD_ASSERT(xdec == ydec);
        ln = OD_LOG_DERING_GRID + OD_LOG_BSIZE0 - xdec;
        n = 1 << ln;
        od_dering(state, buf, n, &state->etmp[pli][(sby << ln)*w +
         (sbx << ln)], w, ln, sbx, sby, nhdr, nvdr, state->quantizer[0],
         xdec, dir, pli, &enc->state.bskip[pli]
         [(sby << (OD_LOG_DERING_GRID - ydec))*enc->state.skip_stride
         + (sbx << (OD_LOG_DERING_GRID - xdec))], enc->state.skip_stride);
        xstride = enc->input_img.planes[pli].xstride;
        ystride = enc->input_img.planes[pli].ystride;
        input = (unsigned char *)&enc->input_img.planes[pli].
         data[(sby << ln)*ystride + (sbx << ln)*xstride];
        output = &state->ctmp[pli][(sby << ln)*w + (sbx << ln)];
        unfiltered_error = 0;
        filtered_error = 0;
        od_ref_buf_to_coeff(state, orig, n, 0, input, xstride, ystride, n, n);
#if 0
        /* Optimize deringing for PSNR. */
        for (y = 0; y < n; y++) {
          for (x = 0; x < n; x++) {
            od_coeff r;
            od_coeff p;
            od_coeff o;
            r = orig[y*OD_BSIZE_MAX + x];
            p = buf[y*OD_BSIZE_MAX + x];
            filtered_error += (r - p)*(double)(r - p);
            o = output[y*w + x];
            unfiltered_error += (r - o)*(double)(r - o);
          }
        }
#else
        /* Optimize deringing for the block size decision metric. */
        {
          od_coeff out[OD_BSIZE_MAX*OD_BSIZE_MAX];
          od_coeff buf32[OD_BSIZE_MAX*OD_BSIZE_MAX];
          for (y = 0; y < n; y++) {
            for (x = 0; x < n; x++) {
              out[y*n + x] = output[y*w + x];
              buf32[y*n + x] = buf[y*n + x];
            }
          }
          unfiltered_error = od_compute_dist(enc, orig, out, n, 3);
          filtered_error = od_compute_dist(enc, orig, buf32, n, 3);
        }
#endif
        up = 0;
        if (sby > 0) {
          up = state->dering_flags[(sby - 1)*nhdr + sbx];
        }
        left = 0;
        if (sbx > 0) {
          left = state->dering_flags[sby*nhdr + (sbx - 1)];
        }
        c = (up << 1) + left;
        filtered_rate = od_encode_cdf_cost(1, state->adapt.clpf_cdf[c], 2);
        unfiltered_rate = od_encode_cdf_cost(0, state->adapt.clpf_cdf[c], 2);
        q2 = state->quantizer[0] * state->quantizer[0];
        filtered = (filtered_error + OD_PVQ_LAMBDA*q2*filtered_rate) <
         (unfiltered_error + OD_PVQ_LAMBDA*q2*unfiltered_rate);
        /*When use_dering is 0, force the deringing filter off.*/
        if (!enc->use_dering) {
          filtered = 0;
        }
        state->dering_flags[sby*nhdr + sbx] = filtered;
        od_encode_cdf_adapt(&enc->ec, filtered, state->adapt.clpf_cdf[c], 2,
         state->adapt.clpf_increment);
        if (filtered) {
          for (y = 0; y < n; y++) {
            for (x = 0; x < n; x++) {
              output[y*w + x] = buf[y*n + x];
            }
          }
          for (pli = 1; pli < nplanes; pli++) {
            xdec = enc->input_img.planes[pli].xdec;
            ydec = enc->input_img.planes[pli].ydec;
            w = frame_width >> xdec;
            ln = OD_LOG_DERING_GRID + OD_LOG_BSIZE0 - xdec;
            n = 1 << ln;
            od_dering(state, buf, n, &state->etmp[pli][(sby << ln)*w +
             (sbx << ln)], w, ln, sbx, sby, nhdr, nvdr,
             state->quantizer[pli], xdec, dir, pli, &enc->state.bskip[pli]
             [(sby << (OD_LOG_DERING_GRID - ydec))*enc->state.skip_stride
             + (sbx << (OD_LOG_DERING_GRID - xdec))], enc->state.skip_stride);
            output = &state->ctmp[pli][(sby << ln)*w + (sbx << ln)];
            for (y = 0; y < n; y++) {
              for (x = 0; x < n; x++) {
                output[y*w + x] = buf[y*n + x];
              }
            }
          }
        }
      }
    }
  }
  for (pli = 0; pli < nplanes; pli++) {
    xdec = enc->input_img.planes[pli].xdec;
    ydec = enc->input_img.planes[pli].ydec;
    w = frame_width >> xdec;
    if (!rdo_only && state->quantizer[0] > 0) {
      for (sby = 0; sby < nvsb; sby++) {
        for (sbx = 0; sbx < nhsb; sbx++) {
          if (mbctx->is_keyframe) {
            od_smooth_recursive(state->ctmp[pli], enc->state.bsize,
             enc->state.bstride, sbx, sby, OD_NBSIZES - 1, w, xdec, ydec,
             OD_BLOCK_32X32, state->quantizer[pli], pli);
          }
        }
      }
    }
    if (!rdo_only) {
      od_coeff_to_ref_plane(state, rec, pli,
       state->ctmp[pli], OD_LOSSLESS(enc, pli));
    }
  }
}

#if defined(OD_LOGGING_ENABLED)
static void od_dump_frame_metrics(daala_enc_ctx *enc) {
  od_state *state;
  int pli;
  int nplanes;
  int frame_width;
  int frame_height;
  state = &enc->state;
  nplanes = state->info.nplanes;
  frame_width = state->frame_width;
  frame_height = state->frame_height;
  for (pli = 0; pli < nplanes; pli++) {
    unsigned char *data;
    int64_t enc_sqerr;
    uint32_t npixels;
    int ystride;
    int xdec;
    int ydec;
    int w;
    int h;
    int x;
    int y;
    enc_sqerr = 0;
    data = enc->input_img.planes[pli].data;
    ystride = enc->input_img.planes[pli].ystride;
    xdec = enc->input_img.planes[pli].xdec;
    ydec = enc->input_img.planes[pli].ydec;
    w = frame_width >> xdec;
    h = frame_height >> ydec;
    npixels = w*h;
    for (y = 0; y < h; y++) {
      unsigned char *rec_row;
      unsigned char *inp_row;
      rec_row = state->ref_imgs[state->ref_imgi[OD_FRAME_SELF]].planes[pli].data
       + state->ref_imgs[state->ref_imgi[OD_FRAME_SELF]].planes[pli].ystride*y;
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

static void od_split_superblocks_rdo(daala_enc_ctx *enc,
 od_mb_enc_ctx *mbctx) {
  od_rollback_buffer rbuf;
  od_encode_checkpoint(enc, &rbuf);
  od_encode_coefficients(enc, mbctx, OD_ENCODE_RDO);
  od_encode_rollback(enc, &rbuf);
}

int daala_encode_img_in(daala_enc_ctx *enc, od_img *img, int duration) {
  int refi;
  int nplanes;
  int pli;
  int frame_width;
  int frame_height;
  int pic_width;
  int pic_height;
  int use_masking;
  od_mb_enc_ctx mbctx;
  od_img *ref_img;
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
  use_masking = enc->use_activity_masking;
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
  od_img_copy_pad(enc, img);

#if defined(OD_DUMP_IMAGES)
  if (od_logging_active(OD_LOG_GENERIC, OD_LOG_DEBUG)) {
    od_img_dump_padded(enc);
  }
#endif
  /* Check if the frame should be a keyframe. */
  mbctx.is_keyframe = (enc->state.cur_time %
   (enc->state.info.keyframe_rate) == 0) ? 1 : 0;
  mbctx.is_golden_frame = (enc->state.cur_time %
   (OD_GOLDEN_FRAME_INTERVAL) == 0) ? 1 : 0;
  if (enc->state.ref_imgi[OD_FRAME_GOLD] < 0) {
    mbctx.is_golden_frame = 1;
  }
  /*Update the buffer state.*/
  if (enc->state.ref_imgi[OD_FRAME_SELF] >= 0) {
    enc->state.ref_imgi[OD_FRAME_PREV] =
     enc->state.ref_imgi[OD_FRAME_SELF];
  }
  /*Select a free buffer to use for this reference frame.*/
  for (refi = 0; refi == enc->state.ref_imgi[OD_FRAME_GOLD]
   || refi == enc->state.ref_imgi[OD_FRAME_PREV]
   || refi == enc->state.ref_imgi[OD_FRAME_NEXT]; refi++);
  enc->state.ref_imgi[OD_FRAME_SELF] = refi;
  /*We must be a keyframe if we don't have a reference.*/
  mbctx.is_keyframe |= !(enc->state.ref_imgi[OD_FRAME_PREV] >= 0);
  mbctx.num_refs = mbctx.is_keyframe ? 0 : OD_MAX_CODED_REFS;
  /* FIXME: This should be dynamic */
  mbctx.use_activity_masking = enc->use_activity_masking;
  mbctx.qm = enc->qm;
  /* Use Haar for lossless since 1) it's more efficient than the DCT and 2)
     PVQ isn't lossless. We only look at luma quality based on the assumption
     that it's silly to have just some planes be lossless. */
  mbctx.use_haar_wavelet = enc->use_haar_wavelet || enc->quality[0] == 0;
  /*Initialize the entropy coder.*/
  od_ec_enc_reset(&enc->ec);
  /*Write a bit to mark this as a data packet.*/
  od_ec_encode_bool_q15(&enc->ec, 0, 16384);
  /*Code the keyframe bit.*/
  od_ec_encode_bool_q15(&enc->ec, mbctx.is_keyframe, 16384);
  /* Code the number of references. */
  if (!mbctx.is_keyframe) {
    od_ec_enc_uint(&enc->ec, mbctx.num_refs - 1, OD_MAX_CODED_REFS);
  }
  /*Code whether or not activity masking is being used.*/
  od_ec_encode_bool_q15(&enc->ec, mbctx.use_activity_masking, 16384);
  /*Code whether flat or hvs quantization matrices are being used.
   * FIXME: will need to be a wider type if other QMs get added */
  od_ec_encode_bool_q15(&enc->ec, mbctx.qm, 16384);
  od_ec_encode_bool_q15(&enc->ec, mbctx.use_haar_wavelet, 16384);
  od_ec_encode_bool_q15(&enc->ec, mbctx.is_golden_frame, 16384);
  for (pli = 0; pli < nplanes; pli++) {
    enc->state.coded_quantizer[pli] =
     od_quantizer_to_codedquantizer(
      od_quantizer_from_quality(enc->quality[pli]));
    enc->state.quantizer[pli] =
     od_codedquantizer_to_quantizer(enc->state.coded_quantizer[pli]);
  }
  if (mbctx.is_keyframe) {
    for (pli = 0; pli < nplanes; pli++) {
      int i;
      int q;
      q = enc->state.quantizer[pli];
      if (q <= OD_DEFAULT_QMS[use_masking][0][pli].interp_q << OD_COEFF_SHIFT) {
        od_interp_qm(&enc->state.pvq_qm_q4[pli][0], q,
         &OD_DEFAULT_QMS[use_masking][0][pli], NULL);
      }
      else {
        i = 0;
        while (OD_DEFAULT_QMS[use_masking][i + 1][pli].qm_q4 != NULL &&
         q > OD_DEFAULT_QMS[use_masking][i + 1][pli].interp_q
         << OD_COEFF_SHIFT) {
          i++;
        }
        od_interp_qm(&enc->state.pvq_qm_q4[pli][0], q,
         &OD_DEFAULT_QMS[use_masking][i][pli],
         &OD_DEFAULT_QMS[use_masking][i + 1][pli]);
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
    /*Boost the keyframe quality slightly (one coded quantizer
      step is the minimum possible).*/
    if ((mbctx.is_keyframe || mbctx.is_golden_frame)
     && enc->state.coded_quantizer[pli] != 0) {
      enc->state.coded_quantizer[pli] =
       OD_MAXI(1, enc->state.coded_quantizer[pli] - 3);
      enc->state.quantizer[pli] =
       od_codedquantizer_to_quantizer(enc->state.coded_quantizer[pli]);
    }
  }
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "is_keyframe=%d", mbctx.is_keyframe));
  /*TODO: Increment frame count.*/
  od_adapt_ctx_reset(&enc->state.adapt, mbctx.is_keyframe);
  if (!mbctx.is_keyframe) {
    od_predict_frame(enc);
    od_encode_mvs(enc, mbctx.num_refs);
  }
  if (mbctx.use_haar_wavelet) {
    od_state_init_superblock_split(&enc->state, OD_BLOCK_32X32);
  } else {
    /* Enable block size RDO for all but complexity 0 and 1. We might want to
       revise that choice if we get a better open-loop block size algorithm. */
    od_state_init_superblock_split(&enc->state, OD_LIMIT_BSIZE_MIN);
    if (enc->complexity >= 2) od_split_superblocks_rdo(enc, &mbctx);
    else od_split_superblocks(enc, mbctx.is_keyframe);
  }
  od_encode_coefficients(enc, &mbctx, OD_ENCODE_REAL);
  ref_img = enc->state.ref_imgs + enc->state.ref_imgi[OD_FRAME_SELF];
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  /*Dump YUV*/
  od_state_dump_yuv(&enc->state, ref_img, "out");
#endif
#if defined(OD_LOGGING_ENABLED)
  od_dump_frame_metrics(enc);
#endif
  enc->packet_state = OD_PACKET_READY;
  od_img_edge_ext(ref_img);
  if (mbctx.is_golden_frame) {
    enc->state.ref_imgi[OD_FRAME_GOLD] =
     enc->state.ref_imgi[OD_FRAME_SELF];
  }
#if defined(OD_DUMP_IMAGES)
  /*Dump reference frame.*/
  /*od_state_dump_img(&enc->state,
   enc->state.ref_img + enc->state.ref_imigi[OD_FRAME_SELF], "ref");*/
#endif
  if (enc->state.info.frame_duration == 0) enc->state.cur_time += duration;
  else enc->state.cur_time += enc->state.info.frame_duration;
#if defined(OD_DUMP_BSIZE_DIST)
  for (pli = 0; pli < nplanes; pli++){
    /* Write value for this frame and reset it */
    fprintf(enc->bsize_dist_file, "%-7G\t",
     10*log10(enc->bsize_dist[pli]));
    enc->bsize_dist_total[pli] += enc->bsize_dist[pli];
    enc->bsize_dist[pli] = 0.0;
  }
  fprintf(enc->bsize_dist_file, "\n");
#endif
  return 0;
}

#if defined(OD_ENCODER_CHECK)
static void daala_encoder_check(daala_enc_ctx *ctx, od_img *img,
 daala_packet *op) {
  int pli;
  od_img out_img;
  od_img dec_img;
  OD_ASSERT(ctx->dec);

  if (daala_decode_packet_in(ctx->dec, &out_img, op) < 0) {
    fprintf(stderr,"decode failed!\n");
    return;
  }
  dec_img = ctx->dec->state.ref_imgs[ctx->dec->state.ref_imgi[OD_FRAME_SELF]];

  OD_ASSERT(img->nplanes == dec_img.nplanes);
  for (pli = 0; pli < img->nplanes; pli++) {
    int plane_width;
    int plane_height;
    int xdec;
    int ydec;
    int i;
    OD_ASSERT(img->planes[pli].xdec == dec_img.planes[pli].xdec);
    OD_ASSERT(img->planes[pli].ydec == dec_img.planes[pli].ydec);
    OD_ASSERT(img->planes[pli].xstride == dec_img.planes[pli].xstride);
    OD_ASSERT(img->planes[pli].ystride == dec_img.planes[pli].ystride);

    xdec = dec_img.planes[pli].xdec;
    ydec = dec_img.planes[pli].ydec;
    plane_width = ctx->dec->state.frame_width >> xdec;
    plane_height = ctx->dec->state.frame_height >> ydec;
    for (i = 0; i < plane_height; i++) {
      if (memcmp(img->planes[pli].data + img->planes[pli].ystride * i,
       dec_img.planes[pli].data + dec_img.planes[pli].ystride * i,
       plane_width)) {
        OD_LOG((OD_LOG_ENCODER, OD_LOG_ERR,
         "Pixel mismatch in row:%d, plane:%d\n", i, pli));
      }
    }
  }
}
#endif

int daala_encode_packet_out(daala_enc_ctx *enc, int last, daala_packet *op) {
  uint32_t nbytes;
  if (enc == NULL || op == NULL) return OD_EFAULT;
  else if (enc->packet_state <= 0 || enc->packet_state == OD_PACKET_DONE) {
    return 0;
  }
  op->packet = od_ec_enc_done(&enc->ec, &nbytes);
  op->bytes = nbytes;
  OD_LOG((OD_LOG_ENCODER, OD_LOG_INFO, "Output Bytes: %ld (%ld Kbits)",
   op->bytes, op->bytes*8/1024));
  op->b_o_s = 0;
  op->e_o_s = last;
  op->packetno = 0;
  op->granulepos = enc->state.cur_time;
  if (last) enc->packet_state = OD_PACKET_DONE;
  else enc->packet_state = OD_PACKET_EMPTY;

#if defined(OD_ENCODER_CHECK)
  /*Compare reconstructed frame against decoded frame.*/
  daala_encoder_check(enc,
   enc->state.ref_imgs + enc->state.ref_imgi[OD_FRAME_SELF], op);
#endif

  return 1;
}
