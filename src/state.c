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
#include <stdio.h>
#include "state.h"
#if defined(OD_X86ASM)
# include "x86/x86int.h"
#endif

/*Initializes the buffers used for reference frames.
  These buffers are padded with 16 extra pixels on each side, to allow
   (relatively) unrestricted motion vectors without special casing reading
   outside the image boundary.
  If chroma is decimated in either direction, the padding is reduced by an
   appropriate factor on the appropriate sides.*/
static void od_state_ref_imgs_init(od_state *state, int nrefs, int nio) {
  daala_info *info;
  od_img *img;
  od_img_plane *iplane;
  unsigned char *ref_img_data;
  size_t data_sz;
  int frame_buf_width;
  int frame_buf_height;
  int plane_buf_width;
  int plane_buf_height;
  int imgi;
  int pli;
  int y;
  OD_ASSERT(nrefs >= 3);
  OD_ASSERT(nrefs <= 4);
  OD_ASSERT(nio == 2);
  info = &state->info;
  data_sz = 0;
  /*TODO: Check for overflow before allocating.*/
  frame_buf_width = state->frame_width + (OD_UMV_PADDING << 1);
  frame_buf_height = state->frame_height + (OD_UMV_PADDING << 1);
  for (pli = 0; pli < info->nplanes; pli++) {
    /*Reserve space for this plane in nrefs reference images.*/
    plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
    plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
    data_sz += plane_buf_width*plane_buf_height*nrefs << 2;
#if defined(OD_DUMP_IMAGES)
    /*Reserve space for this plane in 1 visualization image.*/
    data_sz += plane_buf_width*plane_buf_height << 2;
#endif
    /*Reserve space for this plane in nio input/output images.*/
    data_sz += plane_buf_width*plane_buf_height*nio;
  }
  /*Reserve space for the line buffer in the up-sampler.*/
  data_sz += (frame_buf_width << 1)*8;
  state->ref_img_data = ref_img_data = (unsigned char *)_ogg_malloc(data_sz);
  /*Fill in the reference image structures.*/
  for (imgi = 0; imgi < nrefs; imgi++) {
    img = state->ref_imgs + imgi;
    img->nplanes = info->nplanes;
    img->width = state->frame_width << 1;
    img->height = state->frame_height << 1;
    for (pli = 0; pli < img->nplanes; pli++) {
      plane_buf_width = (frame_buf_width << 1) >> info->plane_info[pli].xdec;
      plane_buf_height = (frame_buf_height << 1) >> info->plane_info[pli].ydec;
      iplane = img->planes + pli;
      iplane->data = ref_img_data
       + ((OD_UMV_PADDING << 1) >> info->plane_info[pli].xdec)
       + plane_buf_width*((OD_UMV_PADDING << 1) >> info->plane_info[pli].ydec);
      ref_img_data += plane_buf_width*plane_buf_height;
      iplane->xdec = info->plane_info[pli].xdec;
      iplane->ydec = info->plane_info[pli].ydec;
      iplane->xstride = 1;
      iplane->ystride = plane_buf_width;
    }
  }
  /*Fill in the reconstruction image structure.*/
  for (imgi = 0; imgi < nio; imgi++) {
    img = state->io_imgs + imgi;
    img->nplanes = info->nplanes;
    img->width = state->frame_width;
    img->height = state->frame_height;
    for (pli = 0; pli < img->nplanes; pli++) {
      plane_buf_width = frame_buf_width >> info->plane_info[pli].xdec;
      plane_buf_height = frame_buf_height >> info->plane_info[pli].ydec;
      iplane = img->planes + pli;
      iplane->data = ref_img_data
       + (OD_UMV_PADDING >> info->plane_info[pli].xdec)
       + plane_buf_width*(OD_UMV_PADDING >> info->plane_info[pli].ydec);
      ref_img_data += plane_buf_width*plane_buf_height;
      iplane->xdec = info->plane_info[pli].xdec;
      iplane->ydec = info->plane_info[pli].ydec;
      iplane->xstride = 1;
      iplane->ystride = plane_buf_width;
    }
  }
  /*Fill in the line buffers.*/
  for (y = 0; y < 8; y++) {
    state->ref_line_buf[y] = ref_img_data + (OD_UMV_PADDING << 1);
    ref_img_data += frame_buf_width << 1;
  }
  /*Mark all of the reference image buffers available.*/
  for (imgi = 0; imgi < nrefs; imgi++) state->ref_imgi[imgi] = -1;
#if defined(OD_DUMP_IMAGES)
  /*Fill in the visualization image structure.*/
  img = &state->vis_img;
  img->nplanes = info->nplanes;
  img->width = frame_buf_width << 1;
  img->height = frame_buf_height << 1;
  for (pli = 0; pli < img->nplanes; pli++) {
    iplane = img->planes + pli;
    plane_buf_width = img->width >> info->plane_info[pli].xdec;
    plane_buf_height = img->height >> info->plane_info[pli].ydec;
    iplane->data = ref_img_data;
    ref_img_data += plane_buf_width*plane_buf_height;
    iplane->xdec = info->plane_info[pli].xdec;
    iplane->ydec = info->plane_info[pli].ydec;
    iplane->xstride = 1;
    iplane->ystride = plane_buf_width;
  }
#endif
}

static void od_state_mvs_init(od_state *state) {
  int nhmvbs;
  int nvmvbs;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  state->mv_grid = (od_mv_grid_pt **)od_calloc_2d(nvmvbs + 1, nhmvbs + 1,
   sizeof(**state->mv_grid));
}

void od_state_opt_vtbl_init_c(od_state *state) {
  state->opt_vtbl.mc_predict1imv8 = od_mc_predict1imv8_c;
  state->opt_vtbl.mc_predict1fmv8 = od_mc_predict1fmv8_c;
  state->opt_vtbl.mc_blend_full8 = od_mc_blend_full8_c;
  state->opt_vtbl.mc_blend_full_split8 = od_mc_blend_full_split8_c;
}

static void od_state_opt_vtbl_init(od_state *state) {
#if defined(OD_X86ASM)
  od_state_opt_vtbl_init_x86(state);
#else
  od_state_opt_vtbl_init_c(state);
#endif
}

int od_state_init(od_state *state, const daala_info *info) {
  int nplanes;
  int pli;
  /*First validate the parameters.*/
  if (info == NULL) return OD_EFAULT;
  nplanes = info->nplanes;
  if (nplanes <= 0 || nplanes > OD_NPLANES_MAX) return OD_EINVAL;
  /*The first plane (the luma plane) must not be subsampled.*/
  if (info->plane_info[0].xdec || info->plane_info[0].ydec) return OD_EINVAL;
  memset(state, 0, sizeof(*state));
  memcpy(&state->info, info, sizeof(*info));
  /*Frame size is a multiple of a super block.*/
  state->frame_width = (info->pic_width + (OD_SUPERBLOCK_SIZE - 1)) &
   ~(OD_SUPERBLOCK_SIZE - 1);
  state->frame_height = (info->pic_height + (OD_SUPERBLOCK_SIZE - 1)) &
   ~(OD_SUPERBLOCK_SIZE - 1);
  state->nhmbs = state->frame_width >> 4;
  state->nvmbs = state->frame_height >> 4;
  od_state_opt_vtbl_init(state);
  od_state_ref_imgs_init(state, 4, 2);
  od_state_mvs_init(state);
  state->nhsb = state->frame_width >> 5;
  state->nvsb = state->frame_height >> 5;
  for (pli = 0; pli < nplanes; pli++) {
    od_adapt_init(&state->adapt_sb[pli], state->nhsb,
     OD_NSB_ADAPT_CTXS, OD_SB_ADAPT_PARAMS);
  }
  state->bsize = (unsigned char *)_ogg_malloc(
   sizeof(*state->bsize)*(state->nhsb + 2)*4*(state->nvsb + 2)*4);
  state->bstride = (state->nhsb + 2)*4;
  state->bsize += 4*state->bstride + 4;
  generic_model_init(&state->pvq_gain_model);
  state->pvq_adapt[OD_ADAPT_K_Q8] = 384;
  state->pvq_adapt[OD_ADAPT_SUM_EX_Q8] = 256;
  state->pvq_adapt[OD_ADAPT_COUNT_Q8] = 104;
  state->pvq_adapt[OD_ADAPT_COUNT_EX_Q8] = 128;
  state->pvq_exg = 2<<16;
  state->pvq_ext = 2<<16;
  return 0;
}

void od_state_clear(od_state *state) {
  int nplanes;
  int pli;
  nplanes = state->info.nplanes;
  for (pli = nplanes; pli-- > 0;) od_adapt_clear(&state->adapt_sb[pli]);
  od_free_2d(state->mv_grid);
  _ogg_free(state->ref_img_data);
  state->bsize -= 4*state->bstride + 4;
  _ogg_free(state->bsize);
}

#if 0
/*Upsamples the reconstructed image to a reference image.
  TODO: Pipeline with reconstruction.*/
void od_state_upsample8(od_state *state, int refi) {
  int pli;
  for (pli = 0; pli < state->io_imgs[OD_FRAME_REC].nplanes; pli++) {
    od_img_plane *siplane;
    od_img_plane *diplane;
    unsigned char *src;
    unsigned char *dst;
    int xpad;
    int ypad;
    int w;
    int h;
    int x;
    int y;
    siplane = state->io_imgs[OD_FRAME_REC].planes + pli;
    diplane = state->ref_imgs[refi].planes + pli;
    xpad = OD_UMV_PADDING >> siplane->xdec;
    ypad = OD_UMV_PADDING >> siplane->ydec;
    w = state->io_imgs[OD_FRAME_REC].width >> siplane->xdec;
    h = state->io_imgs[OD_FRAME_REC].height >> siplane->ydec;
    src = siplane->data;
    dst = diplane->data - (diplane->ystride << 1)*ypad;
    for (y = -ypad; y < h + ypad + 2; y++) {
      /*Horizontal filtering:*/
      if (y < h + ypad) {
        unsigned char *buf;
        buf = state->ref_line_buf[y & 3];
        for (x = -xpad; x < -1; x++) {
          *(buf + (x << 1)) = src[0];
          *(buf + (x << 1| 1)) = src[0];
        }
        *(buf - 2) = src[0];
        *(buf - 1) = OD_CLAMP255((132*src[0] - 4*src[1] + 64) >> 7);
        buf[0] = OD_CLAMP255((121*src[0] + 7*src[1] + 64) >> 7);
        /*buf[0] = src[0];*/
        buf[1] =
         OD_CLAMP255((68*(src[0] + src[1]) - 4*(src[0] + src[2]) + 64) >> 7);
        for (x = 1; x < w - 2; x++) {
          buf[x << 1] =
           OD_CLAMP255((7*src[x - 1] + 114*src[x] + 7*src[x + 1] + 64) >> 7);
          /*buf[x << 1] = src[x];*/
          buf[x << 1 | 1] = OD_CLAMP255((68*(src[x] + src[x + 1])
           - 4*(src[x - 1] + src[x + 2]) + 64) >> 7);
        }
        buf[x << 1] =
         OD_CLAMP255((7*src[x - 1] + 114*src[x] + 7*src[x + 1] + 64) >> 7);
        /*buf[x << 1] = src[x];*/
        buf[x << 1 | 1] = OD_CLAMP255((68*(src[x] + src[x + 1])
         - 4*(src[x - 1] + src[x + 1]) + 64) >> 7);
        x++;
        buf[x << 1] =
         OD_CLAMP255((7*src[x - 1] + 114*src[x] + 7*src[x] + 64) >> 7);
        /*buf[x << 1] = src[x];*/
        buf[x << 1 | 1] = OD_CLAMP255((132*src[x] - 4*src[x - 1] + 64) >> 7);
        for (x++; x < w + xpad; x++) {
          buf[x << 1] = src[w - 1];
          buf[x << 1 | 1] = src[w - 1];
        }
        if (y >= 0 && y + 1 < h) src += siplane->ystride;
      }
      /*Vertical filtering:*/
      if (y >= -ypad + 2) {
        if (y < 1 || y > h + 1) {
          memcpy(dst - (xpad << 1),
           state->ref_line_buf[(y - 2) & 3] - (xpad << 1),
           w + (xpad << 1) << 1);
          /*fprintf(stderr, "%3i: ", (y - 2) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          memcpy(dst - (xpad << 1),
           state->ref_line_buf[(y - 2) & 3] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 2) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
        else {
          unsigned char *buf[4];
          buf[0] = state->ref_line_buf[(y - 3) & 3];
          buf[1] = state->ref_line_buf[(y - 2) & 3];
          buf[2] = state->ref_line_buf[(y - 1) & 3];
          buf[3] = state->ref_line_buf[(y - 0) & 3];
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            *(dst + x) = OD_CLAMP255((7**(buf[0] + x) + 114**(buf[1] + x)
             + 7**(buf[2] + x) + 64) >> 7);
            /**(dst + x) = *(buf[1] + x);*/
          }
          /*fprintf(stderr, "%3i: ", (y - 2) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            *(dst + x) = OD_CLAMP255((68*(*(buf[1] + x) + *(buf[2] + x))
             - 4*(*(buf[0] + x) + *(buf[3] + x)) + 64) >> 7);
          }
          /*fprintf(stderr, "%3i: ", (y - 2) << 1 | 1);
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
#else
/*Upsamples the reconstructed image to a reference image.
  TODO: Pipeline with reconstruction.*/
void od_state_upsample8(od_state *state, od_img *dimg, const od_img *simg) {
  int pli;
  for (pli = 0; pli < state->io_imgs[OD_FRAME_REC].nplanes; pli++) {
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
    xpad = OD_UMV_PADDING >> siplane->xdec;
    ypad = OD_UMV_PADDING >> siplane->ydec;
    w = simg->width >> siplane->xdec;
    h = simg->height >> siplane->ydec;
    src = siplane->data;
    dst = diplane->data - (diplane->ystride << 1)*ypad;
    for (y = -ypad; y < h + ypad + 3; y++) {
      /*Horizontal filtering:*/
      if (y < h + ypad) {
        unsigned char *buf;
        buf = state->ref_line_buf[y & 7];
        memset(buf - (xpad << 1), src[0], (xpad - 2) << 1);
        /*for (x = -xpad; x < -2; x++) {
          *(buf + (x << 1)) = src[0];
          *(buf + (x << 1 | 1)) = src[0];
        }*/
        *(buf - 4) = src[0];
        *(buf - 3) = OD_CLAMP255((31*src[0] + src[1] + 16) >> 5);
        *(buf - 2) = src[0];
        *(buf - 1) = OD_CLAMP255((36*src[0] - 5*src[1] + src[1] + 16) >> 5);
        buf[0] = src[0];
        buf[1] = OD_CLAMP255((20*(src[0] + src[1])
         - 5*(src[0] + src[2]) + src[0] + src[3] + 16) >> 5);
        buf[2] = src[1];
        buf[3] = OD_CLAMP255((20*(src[1] + src[2])
         - 5*(src[0] + src[3]) + src[0] + src[4] + 16) >> 5);
        for (x = 2; x < w - 3; x++) {
          buf[x << 1] = src[x];
          buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
           - 5*(src[x - 1] + src[x + 2]) + src[x - 2] + src[x + 3] + 16) >> 5);
        }
        buf[x << 1] = src[x];
        buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
         - 5*(src[x - 1] + src[x + 2]) + src[x - 2] + src[x + 2] + 16) >> 5);
        x++;
        buf[x << 1] = src[x];
        buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
         - 5*(src[x - 1] + src[x + 1]) + src[x - 2] + src[x + 1] + 16) >> 5);
        x++;
        buf[x << 1] = src[x];
        buf[x << 1 | 1] =
         OD_CLAMP255((36*src[x] - 5*src[x - 1] + src[x - 2] + 16) >> 5);
        x++;
        buf[x << 1] = src[w - 1];
        buf[x << 1 | 1] = OD_CLAMP255((31*src[w - 1] + src[w - 2] + 16) >> 5);
        memset(buf + (++x << 1), src[w - 1], (xpad - 1) << 1);
        /*for (x++; x < w + xpad; x++) {
          buf[x << 1] = src[w - 1];
          buf[x << 1 | 1]=src[w - 1];
        }*/
        if (y >= 0 && y + 1 < h) src += siplane->ystride;
      }
      /*Vertical filtering:*/
      if (y >= -ypad + 3) {
        if (y < 1 || y > h + 3) {
          memcpy(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          memcpy(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
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
          buf[0] = state->ref_line_buf[(y - 5) & 7];
          buf[1] = state->ref_line_buf[(y - 4) & 7];
          buf[2] = state->ref_line_buf[(y - 3) & 7];
          buf[3] = state->ref_line_buf[(y - 2) & 7];
          buf[4] = state->ref_line_buf[(y - 1) & 7];
          buf[5] = state->ref_line_buf[(y - 0) & 7];
          memcpy(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
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
#endif

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
 unsigned char *buf, int ystride, int ref, int pli,
 int vx, int vy, int oc, int s, int log_mvb_sz) {
  od_img_plane *iplane;
  od_mv_grid_pt *grid[4];
  ogg_int32_t mvx[4];
  ogg_int32_t mvy[4];
  const int *dxp;
  const int *dyp;
  int etype;
  int x;
  int y;
  int k;
  iplane = state->ref_imgs[state->ref_imgi[ref]].planes+pli;
  dxp = OD_VERT_SETUP_DX[oc][s];
  dyp = OD_VERT_SETUP_DY[oc][s];
  for (k = 0; k < 4; k++) {
    grid[k] = state->mv_grid[vy + (dyp[k] << log_mvb_sz)]
     + vx + (dxp[k] << log_mvb_sz);
    mvx[k] = (ogg_int32_t)grid[k]->mv[0] << (14 - iplane->xdec);
    mvy[k] = (ogg_int32_t)grid[k]->mv[1] << (14 - iplane->ydec);
  }
  /*Note: We pull the edge types from the grid points _after_ they've been
    offset to account for unsplit edges.
   This means that we do not rely on all the flags that would be implicitly set
    to the edge type of the unsplit edge to be up to date.
   This might also cause some other flags to be set incorrectly, but those are
    precisely the ones we override below.*/
  etype = grid[0]->right | grid[1]->down << 1 |
   grid[3]->right << 2 | grid[0]->down << 3;
  if (!(s & 1)) {
    etype &= ~(1 << ((oc + 1) & 3));
    /*if (oc == 3) etype |= (etype & 8) >> 3;
    else etype |= (etype & 1 << oc) << 1;*/
    etype |= (etype | etype << 4) >> 3 & 1 << ((oc + 1) & 3);
    if (etype >> oc & 1) {
      mvx[(oc + 1) & 3] = mvx[oc] + mvx[(oc + 1) & 3] >> 1;
      mvy[(oc + 1) & 3] = mvy[oc] + mvy[(oc + 1) & 3] >> 1;
    }
  }
  if (!(s & 2)) {
    etype &= ~(1 << ((oc + 2) & 3));
    /*if (oc == 1) etype |= (etype & 1) << 3;
    else etype |= (etype & 1 << ((oc + 3) & 3)) >> 1;*/
    etype |= (etype | etype << 4) >> 1 & 1 << ((oc + 2) & 3);
    if (etype << (-oc & 3) & 8) {
      mvx[(oc + 3) & 3] = mvx[oc] + mvx[(oc + 3) & 3] >> 1;
      mvy[(oc + 3) & 3] = mvy[oc] + mvy[(oc + 3) & 3] >> 1;
    }
  }
  /*fprintf(stderr, "interpolation type: %c%c%c%c (0x%X) ",
   etype & 1 ? 'V' : 'B', etype & 2 ? 'V' : 'B',
   etype & 4 ? 'V' : 'B', etype & 8 ? 'V' : 'B', etype);*/
  x = (vx - 2) << (3 - iplane->xdec);
  y = (vy - 2) << (3 - iplane->ydec);
  od_mc_predict8(state, buf, ystride, iplane->data + y*iplane->ystride + x,
   iplane->ystride, mvx, mvy, etype, oc, s,
   log_mvb_sz + 2 - iplane->xdec, log_mvb_sz + 2 - iplane->ydec);
}

void od_state_pred_block(od_state *state, unsigned char *buf, int ystride,
 int ref, int pli, int vx, int vy, int log_mvb_sz) {
  int half_mvb_sz;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_img_plane *iplane;
    int half_xblk_sz;
    int half_yblk_sz;
    iplane = state->ref_imgs[state->ref_imgi[ref]].planes + pli;
    half_xblk_sz = 1 << (log_mvb_sz + 1 - iplane->xdec);
    half_yblk_sz = 1 << (log_mvb_sz + 1 - iplane->ydec);
    od_state_pred_block(state, buf,
     ystride, ref, pli, vx, vy, log_mvb_sz - 1);
    od_state_pred_block(state, buf + half_xblk_sz,
     ystride, ref, pli, vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_state_pred_block(state, buf + half_yblk_sz*ystride,
     ystride, ref, pli, vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_state_pred_block(state, buf + half_yblk_sz*ystride + half_xblk_sz,
     ystride, ref, pli, vx + half_mvb_sz, vy + half_mvb_sz, log_mvb_sz - 1);
  }
  else {
    int oc;
    int s;
    if (log_mvb_sz < 2) {
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
     buf, ystride, ref, pli, vx, vy, oc, s, log_mvb_sz);
  }
}

int od_state_dump_yuv(od_state *state, od_img *img, const char *suf) {
  static const char *CHROMA_TAGS[4] = {
    " C420jpeg", "", " C422jpeg", " C444"
  };
  char fname[128];
  FILE *fp;
  int pic_width;
  int pic_height;
  int y;
  int pli;
  sprintf(fname, "%08i%s.y4m",
   (int)daala_granule_basetime(state, state->cur_time), suf);
  fp = fopen(fname, "wb");
  pic_width = state->info.pic_width;
  pic_height = state->info.pic_height;
  OD_ASSERT(img->nplanes >= 3);
  fprintf(fp, "YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   pic_width, pic_height, 1, 1, 0, 0,
   CHROMA_TAGS[(img->planes[1].xdec == 0) + (img->planes[1].ydec == 0)*2]);
  fprintf(fp, "FRAME\n");
  for (pli = 0; pli < 3; pli++) {
    int xdec;
    int ydec;
    int ystride;
    xdec = img->planes[pli].xdec;
    ydec = img->planes[pli].ydec;
    ystride = img->planes[pli].ystride;
    for (y = 0; y < (pic_height + ydec) >> ydec; y++) {
      if (fwrite(img->planes[pli].data + ystride*y,
       (pic_width + xdec) >> xdec, 1, fp) < 1) {
        fprintf(stderr, "Error writing to \"%s\".\n", fname);
        return OD_EFAULT;
      }
    }
  }
  fclose(fp);
  return 0;
}

#if defined(OD_DUMP_IMAGES)
# include <png.h>
# include <zlib.h>

/*State visualization.
  None of this is particularly fast.*/

/*static const unsigned char OD_YCbCr_BORDER[3] = {49, 109, 184};*/
static const unsigned char OD_YCbCr_BORDER[3] = {113, 72, 137};
static const unsigned char OD_YCbCr_BEDGE[3] = {41, 240, 110};
static const unsigned char OD_YCbCr_VEDGE[3] = {145, 54, 34};
static const unsigned char OD_YCbCr_MV[3] = {81, 90, 240};

void od_img_draw_point(od_img *img, int x, int y,
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

static void od_state_draw_mv_grid_block(od_state *state,
 int vx, int vy, int log_mvb_sz) {
  int half_mvb_sz;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_state_draw_mv_grid_block(state, vx, vy, log_mvb_sz - 1);
    od_state_draw_mv_grid_block(state, vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_state_draw_mv_grid_block(state, vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_state_draw_mv_grid_block(state, vx + half_mvb_sz, vy + half_mvb_sz,
     log_mvb_sz - 1);
  }
  else {
    od_mv_grid_pt *grid[4];
    ogg_int32_t mvx[4];
    ogg_int32_t mvy[4];
    const int *dxp;
    const int *dyp;
    int mvb_sz;
    int etype;
    int x0;
    int y0;
    int k;
    int oc;
    int s;
    mvb_sz = 1 << log_mvb_sz;
    if (log_mvb_sz < 2) {
      int mask;
      int s1vx;
      int s1vy;
      int s3vx;
      int s3vy;
      mask = (1 << log_mvb_sz + 1) - 1;
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
      mvx[k] = grid[k]->mv[0];
      mvy[k] = grid[k]->mv[1];
    }
    etype = grid[0]->right | grid[1]->down << 1 |
     grid[3]->right << 2 | grid[0]->down << 3;
    if (!(s & 1)) {
      etype &= ~(1 << ((oc + 1) & 3));
      etype |= (etype | etype << 4) >> 3 & 1 << ((oc + 1) & 3);
      if (etype >> oc & 1) {
        mvx[(oc + 1) & 3] = mvx[oc] + mvx[(oc + 1) & 3] >> 1;
        mvy[(oc + 1) & 3] = mvy[oc] + mvy[(oc + 1) & 3] >> 1;
      }
    }
    if (!(s & 2)) {
      etype &= ~(1 << ((oc + 2) & 3));
      etype |= (etype | etype << 4) >> 1 & 1 << ((oc + 2) & 3);
      if (etype << (-oc & 3) & 8) {
        mvx[(oc + 3) & 3] = mvx[oc] + mvx[(oc + 3) & 3] >> 1;
        mvy[(oc + 3) & 3] = mvy[oc] + mvy[(oc + 3) & 3] >> 1;
      }
    }
    x0 = ((vx - 2) << 3) + (OD_UMV_PADDING << 1);
    y0 = ((vy - 2) << 3) + (OD_UMV_PADDING << 1);
    od_img_draw_line(&state->vis_img, x0, y0, x0 + (mvb_sz << 3), y0,
     etype & 1 ? OD_YCbCr_VEDGE : OD_YCbCr_BEDGE);
    od_img_draw_line(&state->vis_img, x0 + (mvb_sz << 3), y0,
     x0 + (mvb_sz << 3), y0 + (mvb_sz << 3),
     etype & 2 ? OD_YCbCr_VEDGE : OD_YCbCr_BEDGE);
    od_img_draw_line(&state->vis_img, x0, y0 + (mvb_sz << 3),
     x0 + (mvb_sz << 3), y0 + (mvb_sz << 3),
     etype & 4 ? OD_YCbCr_VEDGE : OD_YCbCr_BEDGE);
    od_img_draw_line(&state->vis_img, x0, y0, x0, y0 + (mvb_sz << 3),
     etype & 8 ? OD_YCbCr_VEDGE : OD_YCbCr_BEDGE);
  }
}

void od_state_draw_mv_grid(od_state *state) {
  int vx;
  int vy;
  int nhmvbs;
  int nvmvbs;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  for (vy = 0; vy < nvmvbs; vy += 4) {
    for (vx = 0; vx < nhmvbs; vx += 4) {
      od_state_draw_mv_grid_block(state, vx, vy, 2);
    }
  }
}

static void od_state_draw_mvs_block(od_state *state,
 int vx, int vy, int log_mvb_sz) {
  int half_mvb_sz;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_state_draw_mvs_block(state, vx, vy, log_mvb_sz - 1);
    od_state_draw_mvs_block(state, vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_state_draw_mvs_block(state, vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_state_draw_mvs_block(state, vx + half_mvb_sz, vy + half_mvb_sz,
     log_mvb_sz - 1);
  }
  else {
    od_mv_grid_pt *grid[4];
    ogg_int32_t mvx[4];
    ogg_int32_t mvy[4];
    const int *dxp;
    const int *dyp;
    int etype;
    int x0;
    int y0;
    int k;
    int oc;
    int s;
    if (log_mvb_sz < 2) {
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
      mvx[k] = grid[k]->mv[0];
      mvy[k] = grid[k]->mv[1];
    }
    etype = grid[0]->right | grid[1]->down << 1 |
     grid[3]->right << 2 | grid[0]->down << 3;
    if (!(s & 1)) {
      etype &= ~(1 << ((oc + 1) & 3));
      etype |= (etype | etype << 4) >> 3 & 1 << ((oc + 1) & 3);
      if (etype >> oc & 1) {
        mvx[(oc + 1) & 3] = (mvx[oc] + mvx[(oc + 1) & 3]) >> 1;
        mvy[(oc + 1) & 3] = (mvy[oc] + mvy[(oc + 1) & 3]) >> 1;
      }
    }
    if (!(s & 2)) {
      etype &= ~(1 << ((oc + 2) & 3));
      etype |= (etype | etype << 4) >> 1 & 1 << ((oc + 2) & 3);
      if (etype << (-oc & 3) & 8) {
        mvx[(oc + 3) & 3] = (mvx[oc] + mvx[(oc + 3) & 3]) >> 1;
        mvy[(oc + 3) & 3] = (mvy[oc] + mvy[(oc + 3) & 3]) >> 1;
      }
    }
    x0 = ((vx - 2) << 3) + (OD_UMV_PADDING << 1);
    y0 = ((vy - 2) << 3) + (OD_UMV_PADDING << 1);
    for (k = 0; k < 4; k++) {
      x0 = (vx - 2 + (dxp[k] << log_mvb_sz) << 3) + (OD_UMV_PADDING << 1);
      y0 = (vy - 2 + (dyp[k] << log_mvb_sz) << 3) + (OD_UMV_PADDING << 1);
      /*od_img_draw_point(&state->vis_img, x0, y0, OD_YCbCr_MV);*/
      od_img_draw_line(&state->vis_img, x0, y0,
       x0 + OD_DIV_ROUND_POW2(grid[k]->mv[0], 2, 2),
       y0 + OD_DIV_ROUND_POW2(grid[k]->mv[1], 2, 2), OD_YCbCr_MV);
    }
  }
}

void od_state_draw_mvs(od_state *state) {
  int vx;
  int vy;
  int nhmvbs;
  int nvmvbs;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  for (vy = 0; vy < nvmvbs; vy += 4) {
    for (vx = 0; vx < nhmvbs; vx += 4) {
      od_state_draw_mvs_block(state, vx, vy, 2);
    }
  }
}

void od_state_fill_vis(od_state *state) {
  od_img *img;
  od_img *ref_img;
  int pli;
  int xdec;
  int ydec;
  int border;
  int x;
  int y;
  img = &state->vis_img;
  /*Upsample the reconstructed image for better quality.*/
  /*Adjust the data pointers so that the padding works like the reference
     images.*/
  border = OD_UMV_PADDING << 1;
  for (pli = 0; pli < img->nplanes; pli++) {
    img->planes[pli].data += (border >> img->planes[pli].xdec)
     + img->planes[pli].ystride*(border >> img->planes[pli].ydec);
  }
  od_state_upsample8(state, img, state->io_imgs + OD_FRAME_REC);
  /*Upsample the input image, as well, and subtract it to get a difference
     image.*/
  ref_img = state->ref_imgs + state->ref_imgi[OD_FRAME_SELF];
  od_state_upsample8(state, ref_img, &state->input);
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
    memset(img->planes[0].data + (img->planes[0].ystride)*y,
     0, img->width >> xdec);
  }
  for (; y < (img->height - border) >> ydec; y++) {
    memset(img->planes[0].data + img->planes[0].ystride*y, 0, border >> xdec);
    memset(img->planes[0].data + img->planes[0].ystride*y
     + ((img->width - border) >> xdec), 0, border >> xdec);
  }
  for (; y < img->height >> ydec; y++) {
    memset(img->planes[0].data + (img->planes[0].ystride)*y, 0,
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
  od_state_draw_mv_grid(state);
  od_state_draw_mvs(state);
}

/*Dump a PNG of the reconstructed image, or a reference frame.*/
int od_state_dump_img(od_state *state, od_img *img, const char *suf) {
  png_structp png;
  png_infop info;
  png_bytep *data;
  FILE *fp;
  char fname[128];
  unsigned char *p_rows[3];
  unsigned char *p[3];
  int nplanes;
  int pli;
  int x;
  int y;
  sprintf(fname, "%08i%s.png",
   (int)daala_granule_basetime(state, state->cur_time), suf);
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
      yval = (p[0][0] - 16)*(1.0F/219);
      if (nplanes >= 3) {
        cbval = (p[1][0] - 128)*(2*(1 - 0.114F)/224);
        crval = (p[2][0] - 128)*(2*(1 - 0.299F)/224);
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
        p[pli] += (x & mask) == mask;
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      mask = (1 << img->planes[pli].ydec) - 1;
      p_rows[pli] += ((y & mask) == mask)*img->planes[pli].ystride;
    }
  }
  png_init_io(png, fp);
  png_set_compression_level(png, Z_BEST_COMPRESSION);
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

void od_state_mc_predict(od_state *state, int ref) {
  unsigned char __attribute__((aligned(16))) buf[16][16];
  od_img *img;
  int nhmvbs;
  int nvmvbs;
  int pli;
  int vx;
  int vy;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  img = state->io_imgs + OD_FRAME_REC;
  for (vy = 0; vy < nvmvbs; vy += 4) {
    for (vx = 0; vx < nhmvbs; vx += 4) {
      for (pli = 0; pli < img->nplanes; pli++) {
        od_img_plane *iplane;
        unsigned char *p;
        int blk_w;
        int blk_h;
        int blk_x;
        int blk_y;
        int y;
        od_state_pred_block(state,
         buf[0], sizeof(buf[0]), ref, pli, vx, vy, 2);
        /*Copy the predictor into the image, with clipping.*/
        iplane = img->planes + pli;
        blk_w = 16 >> iplane->xdec;
        blk_h = 16 >> iplane->ydec;
        blk_x = (vx - 2) << (2 - iplane->xdec);
        blk_y = (vy - 2) << (2 - iplane->ydec);
        p = buf[0];
        if (blk_x < 0) {
          blk_w += blk_x;
          p -= blk_x;
          blk_x = 0;
        }
        if (blk_y < 0) {
          blk_h += blk_y;
          p -= blk_y*sizeof(buf[0]);
          blk_y = 0;
        }
        if (blk_x + blk_w > img->width >> iplane->xdec) {
          blk_w = (img->width >> iplane->xdec) - blk_x;
        }
        if (blk_y + blk_h > img->height >> iplane->ydec) {
          blk_h = (img->height >> iplane->ydec) - blk_y;
        }
        for (y = blk_y; y < blk_y + blk_h; y++) {
          memcpy(iplane->data + y*iplane->ystride + blk_x, p, blk_w);
          p += sizeof(buf[0]);
        }
      }
    }
  }
}

/*To avoiding having to special-case superblocks on the edges of the image,
   one superblock of padding is maintained on each side of the image.
  These "dummy" superblocks are notionally not subdivided.
  See the comment for the `bsize` member of `od_state` for more information
   about the data layout and meaning.*/
void od_state_init_border_as_32x32(od_state *state) {
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
      bsize[(j*bstride) + i] = 3;
    }
    for (j = nvsb*4; j < (nvsb+1)*4; j++) {
      bsize[(j*bstride) + i] = 3;
    }
  }
  for (j = -4; j < (nvsb+1)*4; j++) {
    for (i = -4; i < 0; i++) {
      bsize[(j*bstride) + i] = 3;
    }
    for (i = nhsb*4; i < (nhsb+1)*4; i++) {
      bsize[(j*bstride) + i] = 3;
    }
  }
}

ogg_int64_t daala_granule_basetime(void *encdec, ogg_int64_t granpos) {
  od_state *state;
  state = (od_state *)encdec;
  if (granpos >= 0) {
    ogg_int64_t key_time;
    ogg_int64_t delta_time;
    key_time = granpos >> state->info.keyframe_granule_shift;
    delta_time = granpos - (key_time << state->info.keyframe_granule_shift);
    return key_time + delta_time;
  }
  return -1;
}

double daala_granule_time(void *encdec, ogg_int64_t granpos) {
  od_state *state;
  ogg_int64_t base_time;
  state = (od_state *)encdec;
  base_time = daala_granule_basetime(encdec, granpos);
  if (base_time >= 0) {
    return base_time*(double)state->info.timebase_denominator/
     state->info.timebase_numerator;
  }
  return -1;
}
