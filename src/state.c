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
#include <string.h>
#include "state.h"
#include "util.h"
#if defined(OD_X86ASM)
# include "x86/x86int.h"
#endif
#include "block_size.h"

/* Scaling compensation for the Haar equivalent basis function. Left is
   for horizontal/vertical. Right is for diagonal. */
#if OD_DISABLE_FILTER || OD_DEBLOCKING
const od_coeff OD_DC_QM[OD_NBSIZES - 1][2] = {
  {16, 16}, {16, 16}, {16, 16}
};
#else
const od_coeff OD_DC_QM[OD_NBSIZES - 1][2] = {
  {21, 25}, {18, 20}, {17, 18}
};
#endif

/* Haar "quantization matrix" for each decomposition level (starting from LF).
   */
const int OD_HAAR_QM[2][OD_LOG_BSIZE_MAX] = {
  /* horizontal/vertical direction. */
  {16, 16, 16, 24, 32},
  /* "diagonal" direction. */
  {16, 16, 24, 32, 48},
};

void *od_aligned_malloc(size_t _sz,size_t _align) {
  unsigned char *p;
  if (_align - 1 > UCHAR_MAX || (_align & (_align - 1))
   || _sz > ~(size_t)0 - _align) {
    return NULL;
  }
  p = (unsigned char *)malloc(_sz + _align);
  if (p != NULL) {
    int offs;
    offs = ((p - (unsigned char *)0) - 1) & (_align - 1);
    p[offs] = offs;
    p += offs + 1;
  }
  return p;
}

void od_aligned_free(void *_ptr) {
  unsigned char *p;
  p = (unsigned char *)_ptr;
  if (p != NULL) {
    int offs;
    offs = *--p;
    free(p - offs);
  }
}

/*This is a smart copy that copies the intersection of the two img planes
   and performs any needed bitdepth or packing conversion.
  Note that this copy performs clamping as it's only used for copying into
   and out of the codec, not for internal reference->reference copies.
  Does not touch any padding/border/unintersected area.*/
void od_img_plane_copy(od_img* dst, od_img* src, int pli) {
  od_img_plane *dst_plane;
  od_img_plane *src_plane;
  unsigned char *dst_data;
  unsigned char *src_data;
  int dst_xstride;
  int dst_ystride;
  int src_xstride;
  int src_ystride;
  int dst_xdec;
  int dst_ydec;
  int src_xdec;
  int src_ydec;
  int dst_plane_width;
  int dst_plane_height;
  int src_plane_width;
  int src_plane_height;
  int w;
  int h;
  int x;
  int y;
  if (pli >= dst->nplanes || pli >= src->nplanes) {
    return;
  }
  dst_plane = dst->planes+pli;
  src_plane = src->planes+pli;
  dst_xstride = dst_plane->xstride;
  dst_ystride = dst_plane->ystride;
  src_xstride = src_plane->xstride;
  src_ystride = src_plane->ystride;

  OD_ASSERT(dst_xstride == 1 || dst_xstride == 2);
  OD_ASSERT(src_xstride == 1 || src_xstride == 2);
  dst_xdec = dst_plane->xdec;
  dst_ydec = dst_plane->ydec;
  src_xdec = src_plane->xdec;
  src_ydec = src_plane->ydec;
  dst_data = dst_plane->data;
  src_data = src_plane->data;
  dst_plane_width = ((dst->width + (1 << dst_xdec) - 1) >> dst_xdec);
  dst_plane_height = ((dst->height + (1 << dst_ydec) - 1) >> dst_ydec);
  src_plane_width = ((src->width + (1 << src_xdec) - 1) >> src_xdec);
  src_plane_height = ((src->height + (1 << src_ydec) - 1) >> src_ydec);
  w = OD_MINI(dst_plane_width, src_plane_width);
  h = OD_MINI(dst_plane_height, src_plane_height);
  for (y = 0; y < h; y++) {
    unsigned char *src_ptr;
    unsigned char *dst_ptr;
    src_ptr = src_data;
    dst_ptr = dst_data;
    if (src_plane->bitdepth > 8) {
      if (dst_plane->bitdepth > 8) {
        /*Both source and destination are multibyte.*/
        if (dst_plane->bitdepth >= src_plane->bitdepth) {
          /*Shift source up into destination.*/
          int upshift;
          upshift = dst_plane->bitdepth - src_plane->bitdepth;
          for (x = 0; x < w; x++) {
            *(int16_t *)dst_ptr = OD_CLAMPI(0,
             *(int16_t *)src_ptr << upshift,
             (1 << dst_plane->bitdepth) - 1);
            src_ptr += src_xstride;
            dst_ptr += dst_xstride;
          }
        }
        else {
          /*Round source down into destination.*/
          int dnshift;
          dnshift = src_plane->bitdepth - dst_plane->bitdepth;
          for (x = 0; x < w; x++) {
            *(int16_t *)dst_ptr = OD_CLAMPI(0,
             (*(int16_t *)src_ptr + (1 << dnshift >> 1)) >> dnshift,
             (1 << dst_plane->bitdepth) - 1);
            src_ptr += src_xstride;
            dst_ptr += dst_xstride;
          }
        }
      }
      else {
        /*Round multibyte source down into 8-bit destination.*/
        int dnshift;
        dnshift = src_plane->bitdepth - 8;
        OD_ASSERT(dst_plane->bitdepth == 8);
        for (x = 0; x < w; x++) {
          *dst_ptr = OD_CLAMP255((*(int16_t *)src_ptr
           + (1 << dnshift >> 1)) >> dnshift);
          src_ptr += src_xstride;
          dst_ptr += dst_xstride;
        }
      }
    }
    else {
      OD_ASSERT(src_plane->bitdepth == 8);
      if (dst_plane->bitdepth > 8) {
        /*Shift 8-bit source up into multibyte destination.*/
        int upshift;
        upshift = dst_plane->bitdepth - 8;
        for (x = 0; x < w; x++) {
          *(int16_t *)dst_ptr = OD_CLAMPI(0,
           *src_ptr << upshift, (1 << dst_plane->bitdepth) - 1);
          src_ptr += src_xstride;
          dst_ptr += dst_xstride;
        }
      }
      else {
        /*Source and destination are both 8-bit buffers.*/
        OD_ASSERT(dst_plane->bitdepth == 8);
        if (src_xstride == 1 && dst_xstride == 1) {
          /*Neither source nor dest are packed; this can be a simple copy.*/
          OD_COPY(dst_data, src_data, w);
        }
        else {
          /*General purpose case: One or both of source and destination have
            an xstride > 1 despite being 8-bit buffers.*/
          for (x = 0; x < w; x++) {
            *dst_ptr = *src_ptr;
            src_ptr += src_xstride;
            dst_ptr += dst_xstride;
          }
        }
      }
    }
    dst_data += dst_ystride;
    src_data += src_ystride;
  }
}

/*This is a smart copy that copies the intersection of the two od_img
   and performs any needed bitdepth conversion.
  Does not touch any padding/border/unintersected area, nor unintersected
   planes.*/
void od_img_copy(od_img* dst, od_img* src) {
  int pli;
  for (pli = 0; pli < OD_MINI(dst->nplanes, src->nplanes); pli++) {
    od_img_plane_copy(dst, src, pli);
  }
}

/*Initializes the buffers used for reference frames.
  These buffers are padded with 16 extra pixels on each side, to allow
   (relatively) unrestricted motion vectors without special casing reading
   outside the image boundary.
  If chroma is decimated in either direction, the padding is reduced by an
   appropriate factor on the appropriate sides.*/
static int od_state_ref_imgs_init(od_state *state, int nrefs) {
  daala_info *info;
  od_img *img;
  od_img_plane *iplane;
  unsigned char *ref_img_data;
  size_t data_sz;
  int frame_buf_width;
  int frame_buf_height;
  int plane_buf_width;
  int plane_buf_height;
  int reference_bytes;
  int reference_bits;
  int imgi;
  int pli;
  OD_ASSERT(nrefs == 4);
  info = &state->info;
  data_sz = 0;
  /*Reference bit depth is constant over the lifetime of od_state, and by
    extension, an encode or decode ctx.*/
  reference_bytes = state->full_precision_references ? 2 : 1;
  reference_bits = state->full_precision_references ? 8 + OD_COEFF_SHIFT : 8;
  /*TODO: Check for overflow before allocating.*/
  frame_buf_width = state->frame_width + (OD_BUFFER_PADDING << 1);
  frame_buf_height = state->frame_height + (OD_BUFFER_PADDING << 1);
  /*Reserve space for the motion comp buffers.*/
  data_sz += OD_MVBSIZE_MAX*OD_MVBSIZE_MAX*reference_bytes*5;
  for (pli = 0; pli < info->nplanes; pli++) {
    /*Reserve space for this plane in nrefs reference images.*/
    plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
    plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
    data_sz += plane_buf_width*plane_buf_height*reference_bytes*nrefs;
  }
  state->ref_img_data = ref_img_data =
    (unsigned char *)od_aligned_malloc(data_sz, 32);
  if (OD_UNLIKELY(!ref_img_data)) {
    return OD_EFAULT;
  }
  /*Fill in the motion comp buffers.*/
  for (imgi = 0; imgi < 5; imgi++) {
    state->mc_buf[imgi] = ref_img_data;
    ref_img_data += OD_MVBSIZE_MAX*OD_MVBSIZE_MAX*reference_bytes;
  }
  /*Fill in the reference image structures.*/
  for (imgi = 0; imgi < nrefs; imgi++) {
    img = state->ref_imgs + imgi;
    img->nplanes = info->nplanes;
    img->width = state->frame_width;
    img->height = state->frame_height;
    for (pli = 0; pli < img->nplanes; pli++) {
      iplane = img->planes + pli;
      iplane->xdec = info->plane_info[pli].xdec;
      iplane->ydec = info->plane_info[pli].ydec;
      plane_buf_width = frame_buf_width >> iplane->xdec;
      plane_buf_height = frame_buf_height >> iplane->ydec;
      iplane->bitdepth = reference_bits;
      iplane->xstride = reference_bytes;
      iplane->ystride = plane_buf_width*reference_bytes;
      iplane->data = ref_img_data
       + (OD_BUFFER_PADDING >> iplane->xdec)*iplane->xstride
       + (OD_BUFFER_PADDING >> iplane->ydec)*iplane->ystride;
      ref_img_data += plane_buf_height*iplane->ystride;
    }
  }
  /*Mark all of the reference image buffers available.*/
  for (imgi = 0; imgi < nrefs; imgi++) state->ref_imgi[imgi] = -1;
  return OD_SUCCESS;
}

static int od_state_mvs_init(od_state *state) {
  state->mv_grid = (od_mv_grid_pt **)od_calloc_2d(state->nvmvbs + 1,
   state->nhmvbs + 1, sizeof(**state->mv_grid));
  if (OD_UNLIKELY(!state->mv_grid)) {
    return OD_EFAULT;
  }
  return OD_SUCCESS;
}

static void od_restore_fpu_c(void) {}

void od_restore_fpu(od_state *state) {
  (*state->opt_vtbl.restore_fpu)();
}

void od_state_opt_vtbl_init_c(od_state *state) {
  if (state->full_precision_references) {
    state->opt_vtbl.mc_predict1fmv = od_mc_predict1fmv16_c;
    state->opt_vtbl.mc_blend_full = od_mc_blend_full16_c;
    state->opt_vtbl.mc_blend_full_split = od_mc_blend_full_split16_c;
    state->opt_vtbl.mc_blend_multi = od_mc_blend_multi16_c;
    state->opt_vtbl.mc_blend_multi_split = od_mc_blend_multi_split16_c;
    OD_COPY(state->opt_vtbl.od_copy_nxn,
     OD_COPY_NXN_16_C, OD_LOG_COPYBSIZE_MAX + 1);
  }
  else {
    state->opt_vtbl.mc_predict1fmv = od_mc_predict1fmv8_c;
    state->opt_vtbl.mc_blend_full = od_mc_blend_full8_c;
    state->opt_vtbl.mc_blend_full_split = od_mc_blend_full_split8_c;
    state->opt_vtbl.mc_blend_multi = od_mc_blend_multi8_c;
    state->opt_vtbl.mc_blend_multi_split = od_mc_blend_multi_split8_c;
    OD_COPY(state->opt_vtbl.od_copy_nxn,
     OD_COPY_NXN_8_C, OD_LOG_COPYBSIZE_MAX + 1);
  }
  OD_COPY(state->opt_vtbl.filter_dering_direction, OD_DERING_DIRECTION_C,
   OD_DERINGSIZES);
  OD_COPY(state->opt_vtbl.filter_dering_orthogonal, OD_DERING_ORTHOGONAL_C,
   OD_DERINGSIZES);
  state->opt_vtbl.restore_fpu = od_restore_fpu_c;
  OD_COPY(state->opt_vtbl.fdct_2d, OD_FDCT_2D_C, OD_NBSIZES + 1);
  OD_COPY(state->opt_vtbl.idct_2d, OD_IDCT_2D_C, OD_NBSIZES + 1);
}

static void od_state_opt_vtbl_init(od_state *state) {
#if defined(OD_X86ASM)
  od_state_opt_vtbl_init_x86(state);
#else
  od_state_opt_vtbl_init_c(state);
#endif
}

static int od_state_init_impl(od_state *state, const daala_info *info) {
  int nplanes;
  int pli;
  /*First validate the parameters.*/
  if (info == NULL) return OD_EFAULT;
  nplanes = info->nplanes;
  if (nplanes <= 0 || nplanes > OD_NPLANES_MAX) return OD_EINVAL;
  /*The first plane (the luma plane) must not be subsampled.*/
  if (info->plane_info[0].xdec || info->plane_info[0].ydec) return OD_EINVAL;
  /*The bitdepth is restricted to a few allowed values.*/
  if (info->bitdepth_mode < OD_BITDEPTH_MODE_8
   || info->bitdepth_mode > OD_BITDEPTH_MODE_12) {
    return OD_EINVAL;
  }
  OD_CLEAR(state, 1);
  OD_COPY(&state->info, info, 1);
  /*Frame size is a multiple of a super block.*/
  state->frame_width = (info->pic_width + (2*OD_BSIZE_MAX - 1)) &
   ~(2*OD_BSIZE_MAX - 1);
  state->frame_height = (info->pic_height + (2*OD_BSIZE_MAX - 1)) &
   ~(2*OD_BSIZE_MAX - 1);
  state->nhmvbs = state->frame_width >> OD_LOG_MVBSIZE_MIN;
  state->nvmvbs = state->frame_height >> OD_LOG_MVBSIZE_MIN;
  /*The master switch; once FPR is ready to go, this can be set to always-on,
     or on only when encoding high-depth.
    FPR must be on for high-depth, including lossless high-depth.
    When FPR is on for 8-bit or 10-bit content, lossless frames are still
     stored with 8 + OD_COEFF_SHIFT bit depth to allow streams with mixed lossy
     and lossless frames.
    FIXME: Switch on when FPR SIMD lands.*/
  state->full_precision_references = 0;
  od_state_opt_vtbl_init(state);
  if (OD_UNLIKELY(od_state_ref_imgs_init(state, 4))) {
    return OD_EFAULT;
  }
  if (OD_UNLIKELY(od_state_mvs_init(state))) {
    return OD_EFAULT;
  }
  state->nhsb = state->frame_width >> OD_LOG_BSIZE_MAX;
  state->nvsb = state->frame_height >> OD_LOG_BSIZE_MAX;
  for (pli = 0; pli < nplanes; pli++) {
    int xdec;
    int ydec;
    int w;
    int h;
    state->sb_dc_mem[pli] = (od_coeff*)malloc(
     sizeof(state->sb_dc_mem[pli][0])*state->nhsb*state->nvsb);
    if (OD_UNLIKELY(!state->sb_dc_mem[pli])) {
      return OD_EFAULT;
    }
    xdec = info->plane_info[pli].xdec;
    ydec = info->plane_info[pli].ydec;
    w = state->frame_width >> xdec;
    h = state->frame_height >> ydec;
    state->ctmp[pli] = (od_coeff *)malloc(w*h*sizeof(*state->ctmp[pli]));
    if (OD_UNLIKELY(!state->ctmp[pli])) {
      return OD_EFAULT;
    }
    state->dtmp[pli] = (od_coeff *)malloc(w*h*sizeof(*state->dtmp[pli]));
    if (OD_UNLIKELY(!state->dtmp[pli])) {
      return OD_EFAULT;
    }
    state->etmp[pli] = (int16_t *)malloc(w*h*sizeof(*state->etmp[pli]));
    if (OD_UNLIKELY(!state->etmp[pli])) {
      return OD_EFAULT;
    }
    state->mctmp[pli] = (od_coeff *)malloc(w*h*sizeof(*state->mctmp[pli]));
    if (OD_UNLIKELY(!state->mctmp[pli])) {
      return OD_EFAULT;
    }
    state->mdtmp[pli] = (od_coeff *)malloc(w*h*sizeof(*state->mdtmp[pli]));
    if (OD_UNLIKELY(!state->mdtmp[pli])) {
      return OD_EFAULT;
    }
    /*We predict chroma planes from the luma plane.  Since chroma can be
      subsampled, we cache subsampled versions of the luma plane in the
      frequency domain.  We can share buffers with the same subsampling.*/
    if (pli > 0) {
      int plj;
      for (plj = 1; plj < pli; plj++) {
        if (xdec == info->plane_info[plj].xdec
          && ydec == info->plane_info[plj].ydec) {
          state->ltmp[pli] = NULL;
          state->lbuf[pli] = state->ltmp[plj];
        }
      }
      if (plj >= pli) {
        state->lbuf[pli] = state->ltmp[pli] = (od_coeff *)malloc(OD_BSIZE_MAX*
         OD_BSIZE_MAX*sizeof(*state->ltmp[pli]));
        if (OD_UNLIKELY(!state->lbuf[pli])) {
          return OD_EFAULT;
        }
      }
    }
    else state->lbuf[pli] = state->ltmp[pli] = NULL;
    state->bskip[pli] = (unsigned char *)malloc(sizeof(*state->bskip)*
     state->nhsb*state->nvsb<<(2*(OD_NBSIZES-1) - xdec - ydec));
    if (OD_UNLIKELY(!state->bskip[pli])) {
      return OD_EFAULT;
    }
  }
  state->bsize = (unsigned char *)malloc(
   sizeof(*state->bsize)*(state->nhsb + 2)*4*(state->nvsb + 2)*4);
  if (OD_UNLIKELY(!state->bsize)) {
    return OD_EFAULT;
  }
  state->bstride = (state->nhsb + 2)*4;
  state->bsize += 4*state->bstride + 4;
  state->skip_stride = state->nhsb << (OD_NBSIZES - 1);
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  state->dump_tags = 0;
  state->dump_files = 0;
#endif
  {
    int nhdr;
    int nvdr;
    nhdr = state->frame_width >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
    nvdr = state->frame_height >> (OD_LOG_DERING_GRID + OD_LOG_BSIZE0);
    state->dering_flags = (unsigned char *)malloc(nhdr * nvdr);
    if (OD_UNLIKELY(!state->dering_flags)) {
      return OD_EFAULT;
    }
  }
  state->sb_q_scaling = (unsigned char *)malloc(state->nhsb * state->nvsb);
  if (OD_UNLIKELY(!state->sb_q_scaling)) {
    return OD_EFAULT;
  }
  return OD_SUCCESS;
}

int od_state_init(od_state *state, const daala_info *info) {
  int ret;
  ret = od_state_init_impl(state, info);
  if (OD_UNLIKELY(ret < 0)) {
    od_state_clear(state);
  }
  return ret;
}

void od_state_clear(od_state *state) {
  int pli;
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  int i;
  if (state->dump_tags > 0) {
    for (i = 0; i < state->dump_tags; i++) fclose(state->dump_files[i].fd);
    free(state->dump_files);
    state->dump_files = 0;
    state->dump_tags = 0;
  }
#endif
  od_free_2d(state->mv_grid);
  od_aligned_free(state->ref_img_data);
  state->bsize -= 4*state->bstride + 4;
  for (pli = 0; pli < state->info.nplanes; pli++) {
    free(state->sb_dc_mem[pli]);
    free(state->ltmp[pli]);
    free(state->dtmp[pli]);
    free(state->etmp[pli]);
    free(state->ctmp[pli]);
    free(state->mctmp[pli]);
    free(state->mdtmp[pli]);
  }
  free(state->bsize);
  for (pli = 0; pli < 3; pli++) free(state->bskip[pli]);
  free(state->dering_flags);
  free(state->sb_q_scaling);
}

/*Probabilities that a motion vector is not coded given two neighbors and the
  consistency of the nearby motion field. Which MVs are used varies by
  level due to the grid geometry, but critically, we never look at MVs in
  blocks to our right or below.

  This data was compiled from video-subset1-short using
  tools/collect_mvf_ec.sh and tools/mv_ec_stats.jl:

  LEVEL 1:
   Probabilties:
    [30512 31715 32546
     19755 22768 25170
     8822 11180 13710]
   Totals:
    [2865820 1304318 1970291
     844820 196641 65383
     311051 50215 11931]
  LEVEL 2:
   Probabilties:
    [15025 11377 11630
     11771 13799 17357
     9106 12384 14943]
   Totals:
    [429074 69682 11586
     329998 44923 6228
     104242 9983 978]
  LEVEL 3:
   Probabilties:
    [20517 21744 24679
     12351 12900 16429
     8029 9085 12245]
   Totals:
    [188607 29632 5623
     145680 24551 3265
     44420 6687 867]
  LEVEL 4:
   Probabilties:
    [9803 8953 10887
     11962 12496 18801
     11424 17400 24094]
   Totals:
    [107908 10753 1833
     51862 4492 671
     7260 371 68]

  These statistics should be regenerated if the number of levels or the size
   of the levels change.*/

static const uint16_t OD_MV_SPLIT_FLAG_PROBZ_Q15[OD_MC_LEVEL_MAX][9] = {
  { 30512, 31715, 32546, 19755, 22768, 25170, 8822, 11180, 13710 },
  { 15025, 11377, 11630, 11771, 13799, 17357, 9106, 12384, 14943 },
  { 20517, 21744, 24679, 12351, 12900, 16429, 8029, 9085, 12245 },
  { 9803, 8953, 10887, 11962, 12496, 18801, 11424, 17400, 24094 },
  { 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384 },
  { 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384, 16384 }
};

void od_adapt_ctx_reset(od_adapt_ctx *state, int is_keyframe) {
  int i;
  int level;
  int pli;
  od_adapt_pvq_ctx_reset(&state->pvq, is_keyframe);
  generic_model_init(&state->mv_model);
  OD_CDFS_INIT(state->mv_ref_cdf, 128);
  state->skip_increment = 128;
  OD_CDFS_INIT(state->skip_cdf, state->skip_increment >> 2);
  state->mv_small_increment = 128;
  OD_CDFS_INIT_FIRST(state->mv_small_cdf, state->mv_small_increment,
   10*state->mv_small_increment);
  state->split_flag_increment = 128;
  for (level = 0; level < OD_MC_LEVEL_MAX; level++) {
    for (i = 0; i < 9; i++) {
      state->split_flag_cdf[level][i][0] = (uint16_t)(
       (uint32_t)OD_MV_SPLIT_FLAG_PROBZ_Q15[level][i]*
       (state->split_flag_increment >> 1) >> 15);
      state->split_flag_cdf[level][i][1] = state->split_flag_increment >> 1;
    }
  }
  state->haar_coeff_increment = 128;
  OD_CDFS_INIT(state->haar_coeff_cdf, state->haar_coeff_increment >> 2);
  state->haar_split_increment = 128;
  OD_CDFS_INIT(state->haar_split_cdf, state->haar_split_increment >> 2);
  state->haar_bits_increment = 128;
  OD_CDFS_INIT(state->haar_bits_cdf, state->haar_bits_increment >> 2);
  for (pli = 0; pli < OD_NPLANES_MAX; pli++) {
    generic_model_init(&state->model_dc[pli]);
    for (i = 0; i < OD_NBSIZES; i++) {
      state->ex_g[pli][i] = 8;
    }
    state->ex_sb_dc[pli] = pli > 0 ? 8 : 32768;
    for (i = 0; i < 4; i++) {
      int j;
      for (j = 0; j < 3; j++) {
        state->ex_dc[pli][i][j] = pli > 0 ? 8 : 32768;
      }
    }
  }
  state->clpf_increment = 128;
  OD_CDFS_INIT(state->clpf_cdf, state->clpf_increment >> 2);
  state->q_increment = 128;
  OD_CDFS_INIT(state->q_cdf, state->q_increment >> 2);
}

void od_state_set_mv_res(od_state *state, int mv_res) {
  int i;
  state->mv_res = mv_res;
  for (i = 0; i < OD_MC_NLEVELS; i++) {
    state->adapt.mv_ex[i] = state->adapt.mv_ey[i] = (24 << 16) >> mv_res;
  }
}

/*The data used to build the following two arrays.*/
const int OD_VERT_D[22] = {
/*0  1        4  5        8  9        12  13          17  18*/
  0, 0, 1, 1, 0, 0, 1, 2, 0, 0, 2, 1, 0, -1, 1, 1, 0, -1, 0, 1, 1, -1
};

/*The vector offsets in the X direction for each vertex from the upper-left,
   indexed by [exterior corner][split state][vertex].*/
const int *const OD_VERT_SETUP_DX[4][4] = {
  {
    OD_VERT_D + 9, OD_VERT_D + 1, OD_VERT_D + 9, OD_VERT_D + 1,
  },
  {
    OD_VERT_D + 13, OD_VERT_D + 13, OD_VERT_D + 1, OD_VERT_D + 1,
  },
  {
    OD_VERT_D + 18, OD_VERT_D + 1, OD_VERT_D + 18, OD_VERT_D + 1,
  },
  {
    OD_VERT_D + 5, OD_VERT_D + 5, OD_VERT_D + 1, OD_VERT_D + 1,
  }
};

/*The vector offsets in the Y direction for each vertex from the upper-left,
   indexed by [exterior corner][split state][vertex].*/
const int *const OD_VERT_SETUP_DY[4][4] = {
  {
    OD_VERT_DY + 4, OD_VERT_DY + 4, OD_VERT_DY + 0, OD_VERT_DY + 0,
  },
  {
    OD_VERT_DY + 8, OD_VERT_DY + 0, OD_VERT_DY + 8, OD_VERT_DY + 0,
  },
  {
    OD_VERT_DY + 12, OD_VERT_DY + 12, OD_VERT_DY + 0, OD_VERT_DY + 0,
  },
  {
    OD_VERT_DY + 17, OD_VERT_DY + 0, OD_VERT_DY + 17, OD_VERT_DY + 0,
  }
};

void od_state_pred_block_from_setup(od_state *state,
 unsigned char *buf, int ystride, int pli,
 int vx, int vy, int oc, int s, int log_mvb_sz) {
  od_img_plane *iplane;
  od_mv_grid_pt *grid[4];
  int32_t mvx[4];
  int32_t mvy[4];
  const unsigned char *src[4];
  const int *dxp;
  const int *dyp;
  int x;
  int y;
  int k;
  int xdec;
  int ydec;
  /* Assumes that xdec and ydec are the same on all references. */
  xdec = state->ref_imgs[state->ref_imgi[OD_FRAME_PREV]].planes[pli].xdec;
  ydec = state->ref_imgs[state->ref_imgi[OD_FRAME_PREV]].planes[pli].ydec;
  dxp = OD_VERT_SETUP_DX[oc][s];
  dyp = OD_VERT_SETUP_DY[oc][s];
  for (k = 0; k < 4; k++) {
    grid[k] = state->mv_grid[vy + (dyp[k]*(1 << log_mvb_sz))]
     + vx + (dxp[k]*(1 << log_mvb_sz));
    mvx[k] = (int32_t)OD_DIV_POW2_RE(grid[k]->mv[0], xdec);
    mvy[k] = (int32_t)OD_DIV_POW2_RE(grid[k]->mv[1], ydec);
    iplane = state->ref_imgs[state->ref_imgi[grid[k]->ref]].planes+pli;
    x = vx << (OD_LOG_MVBSIZE_MIN - iplane->xdec);
    y = vy << (OD_LOG_MVBSIZE_MIN - iplane->ydec);
    src[k] = iplane->data + y*iplane->ystride + x*iplane->xstride;
  }
  od_mc_predict(state, buf, ystride, src,
   iplane->ystride, mvx, mvy, oc, s,
   log_mvb_sz + OD_LOG_MVBSIZE_MIN - iplane->xdec,
   log_mvb_sz + OD_LOG_MVBSIZE_MIN - iplane->ydec);
}

void od_state_pred_block(od_state *state,
 unsigned char *buf, int ystride, int xstride,
 int pli, int vx, int vy, int log_mvb_sz) {
  int half_mvb_sz;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_img_plane *iplane;
    int half_xblk_sz;
    int half_yblk_sz;
    iplane = state->ref_imgs[state->ref_imgi[OD_FRAME_PREV]].planes + pli;
    half_xblk_sz = 1 << (log_mvb_sz + OD_LOG_MVBSIZE_MIN - 1 - iplane->xdec);
    half_yblk_sz = 1 << (log_mvb_sz + OD_LOG_MVBSIZE_MIN - 1 - iplane->ydec);
    od_state_pred_block(state, buf,
     ystride, xstride, pli, vx, vy, log_mvb_sz - 1);
    od_state_pred_block(state, buf + half_xblk_sz*xstride,
     ystride, xstride, pli, vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_state_pred_block(state, buf + half_yblk_sz*ystride,
     ystride, xstride, pli, vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_state_pred_block(state,
     buf + half_yblk_sz*ystride + half_xblk_sz*xstride, ystride, xstride,
     pli, vx + half_mvb_sz, vy + half_mvb_sz, log_mvb_sz - 1);
  }
  else {
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
    od_state_pred_block_from_setup(state,
     buf, ystride, pli, vx, vy, oc, s, log_mvb_sz);
  }
}

int od_state_dump_yuv(od_state *state, od_img *img, const char *tag) {
  static const char *CHROMA_TAGS[4] = {
    " C420jpeg", "", " C422jpeg", " C444"
  };
  char fname[1024];
  FILE *fp;
  int pic_width;
  int pic_height;
  int x;
  int y;
  int pli;
  int needs_header;
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  int i;
  needs_header = 0;
  for (i = 0; i < state->dump_tags &&
    strcmp(tag,state->dump_files[i].tag) != 0; i++);
  if(i>=state->dump_tags) {
    const char *suf;
    OD_ASSERT(strlen(tag)<16);
    state->dump_tags++;
    state->dump_files = realloc(state->dump_files,
     state->dump_tags*sizeof(od_yuv_dumpfile));
    OD_ASSERT(state->dump_files);
    strncpy(state->dump_files[i].tag,tag,16);
#else
  {
    const char *suf;
#endif
    needs_header = 1;
    suf = getenv("OD_DUMP_IMAGES_SUFFIX");
    if (!suf) {
      suf="";
    }
    sprintf(fname, "%08i%s-%s.y4m",
     (int)daala_granule_basetime(state, state->cur_time), tag, suf);
#if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
    state->dump_files[i].fd = fopen(fname, "wb");
  }
  fp = state->dump_files[i].fd;
#else
    fp = fopen(fname, "wb");
  }
#endif
  pic_width = state->info.pic_width;
  pic_height = state->info.pic_height;
  OD_ASSERT(img->nplanes != 2);
  if (needs_header) {
    int fps_num;
    int fps_denom;
    const char *chroma;
    fps_num = state->info.timebase_numerator;
    fps_denom = state->info.timebase_denominator*state->info.frame_duration;
    chroma = img->nplanes == 1 ? " Cmono" :
     CHROMA_TAGS[(img->planes[1].xdec == 0) + (img->planes[1].ydec == 0)*2];
    fprintf(fp, "YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
     pic_width, pic_height, fps_num, fps_denom,
     state->info.pixel_aspect_numerator, state->info.pixel_aspect_denominator,
     chroma);
  }
  fprintf(fp, "FRAME\n");
  for (pli = 0; pli < OD_MINI(img->nplanes, 3); pli++) {
    int xdec;
    int ydec;
    int xstride;
    int ystride;
    xdec = img->planes[pli].xdec;
    ydec = img->planes[pli].ydec;
    xstride = img->planes[pli].xstride;
    ystride = img->planes[pli].ystride;
    for (y = 0; y < (pic_height + ydec) >> ydec; y++) {
      if (xstride > 1) {
        for (x = 0; x < (pic_width + xdec) >> xdec; x++) {
          int value;
          value = (*((int16_t *)(img->planes[pli].data + ystride*y + xstride*x))
           + (1 << (img->planes[pli].bitdepth - 9)))
           >> (img->planes[pli].bitdepth - 8);
          if (fputc(OD_CLAMP255(value), fp) == EOF) {
            fprintf(stderr, "Error writing to \"%s\".\n", fname);
            return OD_EFAULT;
          }
        }
      }
      else {
        if (fwrite(img->planes[pli].data + ystride*y,
                   (pic_width + xdec) >> xdec, 1, fp) < 1) {
          fprintf(stderr, "Error writing to \"%s\".\n", fname);
          return OD_EFAULT;
        }
      }
    }
  }
  return 0;
}

#if defined(OD_DUMP_IMAGES)
# include <png.h>
# include <zlib.h>

#define OD_SRCVAL_TO_8FP(pli,x) (img->planes[pli].bitdepth == 8 ? \
  ((float)p[(pli)][(x)]) : \
  (((int16_t *)p[(pli)])[(x)]*(1.0F/(1 << img->planes[pli].bitdepth - 8))))

/*Dump a PNG of the reconstructed image, or a reference frame.*/
int od_state_dump_img(od_state *state, od_img *img, const char *tag) {
  png_structp png;
  png_infop info;
  png_bytep *data;
  FILE *fp;
  char fname[1024];
  unsigned char *p_rows[3];
  unsigned char *p[3];
  int nplanes;
  int pli;
  int x;
  int y;
  char *suf;
  suf = getenv("OD_DUMP_IMAGES_SUFFIX");
  if (!suf) {
    suf="";
  }
  sprintf(fname, "%08i%s%s.png",
   (int)daala_granule_basetime(state, state->cur_time), tag, suf);
  fp = fopen(fname, "wb");
  png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png == NULL) {
    fclose(fp);
    return OD_EFAULT;
  }
  info = png_create_info_struct(png);
  if (info == NULL) {
    png_destroy_write_struct(&png, NULL);
    fclose(fp);
    return OD_EFAULT;
  }
  if (setjmp(png_jmpbuf(png))) {
    png_destroy_write_struct(&png, &info);
    fclose(fp);
    return OD_EFAULT;
  }
  data = (png_bytep *)od_malloc_2d(img->height, 6*img->width, sizeof(**data));
  if (img->nplanes < 3) nplanes = 1;
  else nplanes = 3;
  for (pli = 0; pli < nplanes; pli++) p_rows[pli] = img->planes[pli].data;
  /*Chroma up-sampling is just done with a box filter.
    This is very likely what will actually be used in practice on a real
     display, and also removes one more layer to search in for the source of
     artifacts.
    As an added bonus, it's dead simple.*/
  for (y = 0; y < img->height; y++) {
    int mask;
    /*LOOP VECTORIZES.*/
    for (pli = 0; pli < nplanes; pli++) p[pli] = p_rows[pli];
    for (x = 0; x < img->width; x++) {
      float yval;
      float cbval;
      float crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      /*This is intentionally slow and very accurate.*/
      yval = (OD_SRCVAL_TO_8FP(0, 0) - 16)*(1.0F/219);
      if (nplanes >= 3) {
        cbval = (OD_SRCVAL_TO_8FP(1, 0) - 128)*(2*(1 - 0.114F)/224);
        crval = (OD_SRCVAL_TO_8FP(2, 0) - 128)*(2*(1 - 0.299F)/224);
      }
      else cbval = crval = 0;
      rval = OD_CLAMPI(0, (int)(65535*(yval + crval) + 0.5F), 65535);
      gval = OD_CLAMPI(0, (int)(65535*
       (yval - cbval*(0.114F/0.587F) - crval*(0.299F/0.587F)) + 0.5F), 65535);
      bval = OD_CLAMPI(0, (int)(65535*(yval + cbval) + 0.5F), 65535);
      data[y][6*x + 0] = (unsigned char)(rval >> 8);
      data[y][6*x + 1] = (unsigned char)(rval & 0xFF);
      data[y][6*x + 2] = (unsigned char)(gval >> 8);
      data[y][6*x + 3] = (unsigned char)(gval & 0xFF);
      data[y][6*x + 4] = (unsigned char)(bval >> 8);
      data[y][6*x + 5] = (unsigned char)(bval & 0xFF);
      for (pli = 0; pli < nplanes; pli++) {
        mask = (1 << img->planes[pli].xdec) - 1;
        p[pli] += ((x & mask) == mask)*img->planes[pli].xstride;
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      mask = (1 << img->planes[pli].ydec) - 1;
      p_rows[pli] += ((y & mask) == mask)*img->planes[pli].ystride;
    }
  }
  png_init_io(png, fp);
  png_set_compression_level(png, Z_DEFAULT_COMPRESSION);
  png_set_IHDR(png, info, img->width, img->height, 16, PNG_COLOR_TYPE_RGB,
   PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  /*TODO: Define real colorspace.*/
  /*png_set_gAMA(png, info, 2.2);
  png_set_cHRM_fixed(png, info, 31006, 31616, 67000, 32000,
   21000, 71000, 14000, 8000);*/
  png_set_pHYs(png, info, state->info.pixel_aspect_numerator,
   state->info.pixel_aspect_denominator, 0);
  png_set_rows(png, info, data);
  png_write_png(png, info, PNG_TRANSFORM_IDENTITY, NULL);
  png_write_end(png, info);
  png_destroy_write_struct(&png, &info);
  od_free_2d(data);
  fclose(fp);
  return 0;
}
#endif

void od_state_mc_predict(od_state *state, od_img *img_dst) {
  int nhmvbs;
  int nvmvbs;
  int pli;
  int vx;
  int vy;
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy < nvmvbs; vy += OD_MVB_DELTA0) {
    for (vx = 0; vx < nhmvbs; vx += OD_MVB_DELTA0) {
      for (pli = 0; pli < img_dst->nplanes; pli++) {
        od_img_plane *iplane_dst;
        int blk_x;
        int blk_y;
        int xstride;
        int ystride;
        iplane_dst = img_dst->planes + pli;
        blk_x = vx << OD_LOG_MVBSIZE_MIN >> iplane_dst->xdec;
        blk_y = vy << OD_LOG_MVBSIZE_MIN >> iplane_dst->ydec;
        xstride = iplane_dst->xstride;
        ystride = iplane_dst->ystride;
        od_state_pred_block(state,
         iplane_dst->data + blk_y*ystride + blk_x*xstride,
         ystride, xstride, pli, vx, vy, OD_LOG_MVB_DELTA0);
      }
    }
  }
}

/* Initialize superblock split decisions. */
void od_state_init_superblock_split(od_state *state, unsigned char bsize) {
  int nhsb;
  int nvsb;
  int i;
  int j;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  for (i = 0; i < 4*nvsb; i++) {
    for (j = 0; j < 4*nhsb; j++) {
      state->bsize[i*state->bstride + j] = bsize;
    }
  }
}

/*To avoiding having to special-case superblocks on the edges of the image,
   one superblock of padding is maintained on each side of the image.
  These "dummy" superblocks are notionally not subdivided.
  See the comment for the `bsize` member of `od_state` for more information
   about the data layout and meaning.*/
void od_state_init_border(od_state *state) {
  int i;
  int j;
  int nhsb;
  int nvsb;
  unsigned char *bsize;
  int bstride;
  nhsb = state->nhsb;
  nvsb = state->nvsb;
  bsize = state->bsize;
  bstride = state->bstride;
  for (i = -4; i < (nhsb+1)*4; i++) {
    for (j = -4; j < 0; j++) {
      bsize[(j*bstride) + i] = OD_LIMIT_BSIZE_MAX;
    }
    for (j = nvsb*4; j < (nvsb+1)*4; j++) {
      bsize[(j*bstride) + i] = OD_LIMIT_BSIZE_MAX;
    }
  }
  for (j = -4; j < (nvsb+1)*4; j++) {
    for (i = -4; i < 0; i++) {
      bsize[(j*bstride) + i] = OD_LIMIT_BSIZE_MAX;
    }
    for (i = nhsb*4; i < (nhsb+1)*4; i++) {
      bsize[(j*bstride) + i] = OD_LIMIT_BSIZE_MAX;
    }
  }
}

int64_t daala_granule_basetime(void *encdec, int64_t granpos) {
  od_state *state;
  state = (od_state *)encdec;
  if (granpos >= 0) {
    int64_t key_time;
    int64_t delta_time;
    key_time = granpos >> state->info.keyframe_granule_shift;
    delta_time = granpos - (key_time << state->info.keyframe_granule_shift);
    return key_time + delta_time;
  }
  return -1;
}

double daala_granule_time(void *encdec, int64_t granpos) {
  od_state *state;
  int64_t base_time;
  state = (od_state *)encdec;
  base_time = daala_granule_basetime(encdec, granpos);
  if (base_time >= 0) {
    return base_time*(double)state->info.timebase_denominator/
     state->info.timebase_numerator;
  }
  return -1;
}

/*Extend the edge into the padding.*/
/*This is used on internal and thus planar buffers only.
  We can assume depth == 8 implies xstride == 1 and depth > 8 implies
   xstride == 2.*/
static void od_img_plane_edge_ext(od_img_plane *dst_p,
 int plane_width, int plane_height, int horz_padding, int vert_padding) {
  ptrdiff_t xstride;
  ptrdiff_t ystride;
  unsigned char *dst_data;
  unsigned char *dst;
  int x;
  int y;
  xstride = dst_p->xstride;
  ystride = dst_p->ystride;
  dst_data = dst_p->data;

  OD_ASSERT((horz_padding&1) == 0);
  OD_ASSERT((xstride == 1 && dst_p->bitdepth == 8)
   || (xstride == 2 && dst_p->bitdepth > 8));
  /*Left side.*/
  for (y = 0; y < plane_height; y++) {
    dst = dst_data + ystride*y;
    if(xstride == 1){
      for (x = 1; x <= horz_padding; x++) {
        (dst - x)[0] = dst[0];
      }
    }else{
      for (x = 1; x <= horz_padding; x++) {
        (dst - (x << 1))[0] = dst[0];
        (dst - (x << 1))[1] = dst[1];
      }
    }
  }
  /*Right side.*/
  for (y = 0; y < plane_height; y++) {
    dst = dst_data + ystride*y + xstride*(plane_width - 1);
    if(xstride == 1){
      for (x = 1; x <= horz_padding; x++) {
        dst[x] = dst[0];
      }
    }else{
      for (x = 1; x <= horz_padding; x++) {
        (dst + (x << 1))[0] = dst[0];
        (dst + (x << 1))[1] = dst[1];
      }
    }
  }
  /*Top.*/
  dst = dst_data - horz_padding*xstride;
  for (y = 0; y < vert_padding; y++) {
    for (x = 0; x < (plane_width + 2*horz_padding)*xstride; x++) {
      (dst - ystride)[x] = dst[x];
    }
    dst -= ystride;
  }
  /*Bottom.*/
  dst = dst_data - horz_padding*xstride + plane_height*ystride;
  for (y = 0; y < vert_padding; y++) {
    for (x = 0; x < (plane_width + 2*horz_padding)*xstride; x++) {
      dst[x] = (dst - ystride)[x];
    }
    dst += ystride;
  }
}

void od_img_edge_ext(od_img* src) {
  int pli;
  for (pli = 0; pli < src->nplanes; pli++) {
    int xdec;
    int ydec;
    xdec = (src->planes + pli)->xdec;
    ydec = (src->planes + pli)->ydec;
    od_img_plane_edge_ext(&src->planes[pli],
     src->width >> xdec, src->height >> ydec,
     OD_BUFFER_PADDING >> xdec, OD_BUFFER_PADDING >> ydec);
  }
}

/*General purpose reference od_img block to coefficient block
  conversion routine.*/
void od_ref_buf_to_coeff(od_state *state,
 od_coeff *dst, int dst_ystride, int lossless_p,
 unsigned char *src, int src_xstride, int src_ystride,
 int w, int h) {
  int x;
  int y;
  int coeff_shift;
  OD_ASSERT(src_xstride == 1 || src_xstride == 2);
  if (src_xstride == 1) {
    /*The references are running at 8 bits.
      The transforms and coefficients may be operating at 8 bits (during
       lossless coding) or 8 + OD_COEFF_SHIFT bits (during lossy coding).*/
    coeff_shift = lossless_p ?
     (state->info.bitdepth_mode - OD_BITDEPTH_MODE_8)*2 :
     OD_COEFF_SHIFT;
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
        dst[x] = (src[x] - 128)*(1 << coeff_shift);
      }
      dst += dst_ystride;
      src += src_ystride;
    }
  }
  else {
    /*The references are running at greater than 8 bits, implying FPR.
      The transforms and coefficients may be operating at any supported
       bit depth (8, 10 or 12) */
    coeff_shift = lossless_p ?
     OD_COEFF_SHIFT - (state->info.bitdepth_mode - OD_BITDEPTH_MODE_8)*2 :
     0;
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
        dst[x] = (((int16_t *)src)[x] - (1 << (8 + OD_COEFF_SHIFT) >> 1)
         + (1 << coeff_shift >> 1)) >> coeff_shift;
      }
      dst += dst_ystride;
      src += src_ystride;
    }
  }
}

/*Convenience version of the general reference->coeff copying routine
   for a whole-visible-buffer copy.*/
void od_ref_plane_to_coeff(od_state *state, od_coeff *dst, int lossless_p,
 od_img *src, int pli) {
  od_img_plane *iplane;
  int xdec;
  int ydec;
  int h;
  int w;
  int src_xstride;
  int src_ystride;
  iplane = src->planes + pli;
  xdec = iplane->xdec;
  ydec = iplane->ydec;
  w = src->width >> xdec;
  h = src->height >> ydec;
  src_xstride = iplane->xstride;
  src_ystride = iplane->ystride;
  od_ref_buf_to_coeff(state, dst, w, lossless_p,
   iplane->data, src_xstride, src_ystride, w, h);
}

/*General purpose coefficient block to reference od_img block
  conversion routine.*/
void od_coeff_to_ref_buf(od_state *state,
 unsigned char *dst, int dst_xstride, int dst_ystride,
 od_coeff *src, int src_ystride, int lossless_p,
 int w, int h) {
  int x;
  int y;
  int coeff_shift;
  OD_ASSERT(dst_xstride == 1 || dst_xstride == 2);
  if (dst_xstride == 1) {
    /*The references are running at 8 bits.
      The transforms and coefficients may be operating at 8 bits (during
       lossless coding) or 8 + OD_COEFF_SHIFT bits (during lossy coding).*/
    coeff_shift = lossless_p ?
     (state->info.bitdepth_mode - OD_BITDEPTH_MODE_8)*2 :
     OD_COEFF_SHIFT;
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
        *(dst + x) =
         OD_CLAMP255(((src[x] + (1 << coeff_shift >> 1)) >> coeff_shift) +
         128);
      }
      dst += dst_ystride;
      src += src_ystride;
    }
  }
  else {
    /*The references are running at greater than 8 bits, implying FPR.
      An FPR reference must run at full depth (8 + OD_COEFF_SHIFT).
      The transforms and coefficients may be operating at any supported
       bit depth (8, 10 or 12) */
    coeff_shift = lossless_p
     ? OD_COEFF_SHIFT - (state->info.bitdepth_mode - OD_BITDEPTH_MODE_8)*2
     : 0;
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x++) {
        ((int16_t *)dst)[x] =
          OD_CLAMPFPR((src[x]*(1 << coeff_shift)) + (128 << OD_COEFF_SHIFT));
      }
      dst += dst_ystride;
      src += src_ystride;
    }
  }
}

/*Convenience version of the general coeff->reference copying routine
   for a whole-visible-buffer copy.*/
void od_coeff_to_ref_plane(od_state *state, od_img *dst, int pli,
 od_coeff *src, int lossless_p) {
  od_img_plane *iplane;
  int xdec;
  int ydec;
  int h;
  int w;
  int dst_xstride;
  int dst_ystride;
  iplane = dst->planes + pli;
  xdec = iplane->xdec;
  ydec = iplane->ydec;
  w = dst->width >> xdec;
  h = dst->height >> ydec;
  dst_xstride = iplane->xstride;
  dst_ystride = iplane->ystride;
  od_coeff_to_ref_buf(state, iplane->data, dst_xstride, dst_ystride,
   src, w, lossless_p, w, h);
}

void od_init_skipped_coeffs(od_coeff *d, od_coeff *pred, int is_keyframe,
 int bo, int n, int w) {
  int i;
  int j;
  if (is_keyframe) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        /* skip DC */
        if (i || j) d[bo + i*w + j] = 0;
      }
    }
  } else {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        d[bo + i*w + j] = pred[i*n + j];
      }
    }
  }
}
