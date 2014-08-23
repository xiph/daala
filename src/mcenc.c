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
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>

/*Encoder-only motion compensation routines.*/

/*TODO:
 - Develop a real encoding and measure real bits.
 - Thresholds for DP.
   + How do we calculate them?
   + How do they propagate between frames (block sizes change)*/

/*The frame number to animate.*/
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
# define ANI_FRAME (69)
#endif

typedef struct od_mv_err_node od_mv_err_node;

#include "logging.h"
#include "mcenc.h"

typedef int od_offset[2];
typedef int od_pattern[8];

/*The number of bits to reduce chroma SADs by, if used.*/
#define OD_MC_CHROMA_SCALE (2)

/*The subdivision level of a MV in the mesh, given its position (mod 4).*/
static const int OD_MC_LEVEL[4][4] = {
  { 0, 4, 2, 4 },
  { 4, 3, 4, 3 },
  { 2, 4, 1, 4 },
  { 4, 3, 4, 3 }
};

/*Ancestor lists for a vertex.
  These are stored as lists of offsets to the vertices in the domain.
  Level 0 ancestors are not included, as they cannot be decimated.*/
/*Lists for level 2 vertices.*/
static const od_offset OD_ANCESTORS2[2][2] = {
  { { 0, -2}, { 0,  2} },
  { {-2,  0}, { 2,  0} },
};
/*Lists for level 3 vertices.*/
static const od_offset OD_ANCESTORS3[4][5] = {
  { { 1, -1}, {-1, 1}, { 1, -3}, {-3,  1}, { 1,  1} },
  { {-1, -1}, { 1, 1}, {-1, -3}, {-1,  1}, { 3,  1} },
  { {-1, -1}, { 1, 1}, {-3, -1}, { 1, -1}, { 1,  3} },
  { { 1, -1}, {-1, 1}, {-1, -1}, { 3, -1}, {-1,  3} },
};
/*Lists for level 4 vertices.*/
static const od_offset OD_ANCESTORS4[8][9] = {
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2},
    {-3, -2}, { 1, -2}, {-3,  2}, { 1,  2}
  },
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2},
    {-1, -2}, { 3, -2}, {-1,  2}, { 3,  2}
  },
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1},
    {-2, -3}, { 2, -3}, {-2,  1}, { 2,  1}
  },
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1},
    { 0, -3}, {-4,  1}, { 0,  1}, { 4,  1}
  },
  {
    { 0, -1}, { 0,  1}, { 1, -2}, {-1,  0}, { 1,  2},
    { 1, -4}, {-3,  0}, { 1,  0}, { 1,  4}
  },
  {
    { 0, -1}, { 0,  1}, {-1, -2}, { 1,  0}, {-1,  2},
    {-1, -4}, {-1,  0}, { 3,  0}, {-1,  4}
  },
  {
    {-1,  0}, { 1,  0}, { 0, -1}, {-2,  1}, { 2,  1},
    {-2, -1}, { 2, -1}, {-2,  3}, { 2,  3}
  },
  {
    {-1,  0}, { 1,  0}, {-2, -1}, { 2, -1}, { 0,  1},
    {-4, -1}, { 0, -1}, { 4, -1}, { 0,  3}
  },
};

/*The number of ancestors in each list in the grid pattern.*/
static const int OD_NANCESTORS[4][4] = {
  { 0, 9, 2, 9 },
  { 9, 5, 9, 5 },
  { 2, 9, 0, 9 },
  { 9, 5, 9, 5 }
};
/*The lists for each vertex in the grid pattern.*/
static const od_offset *OD_ANCESTORS[4][4] = {
  {             NULL, OD_ANCESTORS4[0], OD_ANCESTORS2[0], OD_ANCESTORS4[1] },
  { OD_ANCESTORS4[2], OD_ANCESTORS3[0], OD_ANCESTORS4[3], OD_ANCESTORS3[1] },
  { OD_ANCESTORS2[1], OD_ANCESTORS4[4],             NULL, OD_ANCESTORS4[5] },
  { OD_ANCESTORS4[6], OD_ANCESTORS3[2], OD_ANCESTORS4[7], OD_ANCESTORS3[3] }
};

/*Computes the Sum of Absolute Differences: slow path.*/
int od_mc_compute_sad_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride, int dxstride, int w, int h) {
  const unsigned char *ref0;
  int i;
  int j;
  int ret;
  ret = 0;
  ref0 = ref;
  for (j = 0; j < h; j++) {
    ref = ref0;
    for (i = 0; i < w; i++) {
      ret += abs(ref[0] - src[i]);
      ref += dxstride;
    }
    src += systride;
    ref0 += dystride;
  }
  return ret;
}

int od_mc_compute_sad_4x4_xstride_1_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad_c(src, systride, ref, dystride, 1, 4, 4);
}

int od_mc_compute_sad_8x8_xstride_1_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad_c(src, systride, ref, dystride, 1, 8, 8);
}

int od_mc_compute_sad_16x16_xstride_1_c(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride) {
  return od_mc_compute_sad_c(src, systride, ref, dystride, 1, 16, 16);
}

/*Computes the SAD of the input image against the given predictor.*/
static ogg_int32_t od_enc_sad8(od_enc_ctx *enc, const unsigned char *p,
 int pystride, int pxstride, int pli, int x, int y, int log_blk_sz) {
  od_state *state;
  od_img_plane *iplane;
  unsigned char *src;
  int clipx;
  int clipy;
  int clipw;
  int cliph;
  int w;
  int h;
  ogg_int32_t ret;
  state = &enc->state;
  iplane = state->io_imgs[OD_FRAME_INPUT].planes + pli;
  /*Compute the block dimensions in the target image plane.*/
  x >>= iplane->xdec;
  y >>= iplane->ydec;
  w = 1 << (log_blk_sz - iplane->xdec);
  h = 1 << (log_blk_sz - iplane->ydec);
  /*Clip the block against the active picture region.*/
  clipx = -x;
  if (clipx > 0) {
    w -= clipx;
    p += clipx*pxstride;
    x += clipx;
  }
  clipy = -y;
  if (clipy > 0) {
    h -= clipy;
    p += clipy*pystride;
    y += clipy;
  }
  clipw = ((state->info.pic_width + (1 << iplane->xdec) - 1) >> iplane->xdec)
   - x;
  w = OD_MINI(w, clipw);
  cliph = ((state->info.pic_height + (1 << iplane->ydec) - 1) >> iplane->ydec)
   - y;
  h = OD_MINI(h, cliph);
  /*OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "[%i, %i]x[%i, %i]", x, y, w, h));*/
  /*Compute the SAD.*/
  src = iplane->data + y*iplane->ystride + x*iplane->xstride;
  if (pxstride != 1) {
    /*Default C implementation.*/
    ret = od_mc_compute_sad_c(src, iplane->ystride,
     p, pystride, pxstride, w, h);
  }
  else if (w == 4 && h == 4) {
    ret = (*enc->opt_vtbl.mc_compute_sad_4x4_xstride_1)(src, iplane->ystride,
     p, pystride);
  }
  else if (w == 8 && h == 8) {
    ret = (*enc->opt_vtbl.mc_compute_sad_8x8_xstride_1)(src, iplane->ystride,
     p, pystride);
  }
  else if (w == 16 && h == 16) {
    ret = (*enc->opt_vtbl.mc_compute_sad_16x16_xstride_1)(src, iplane->ystride,
     p, pystride);
  }
  else {
    /*Default C implementation.*/
    ret = od_mc_compute_sad_c(src, iplane->ystride,
     p, pystride, pxstride, w, h);
  }
  return ret;
}

static void od_mv_est_init(od_mv_est_ctx *est, od_enc_ctx *enc) {
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  est->enc = enc;
  nhmvbs = (enc->state.nhmbs + 1) << 2;
  nvmvbs = (enc->state.nvmbs + 1) << 2;
  est->sad_cache[1] = (od_sad4 **)od_malloc_2d(nvmvbs >> 1, nhmvbs >> 1,
   sizeof(est->sad_cache[1][0][0]));
  est->sad_cache[0] = (od_sad4 **)od_malloc_2d(nvmvbs, nhmvbs,
   sizeof(est->sad_cache[1][0][0]));
  est->mvs = (od_mv_node **)od_calloc_2d(nvmvbs + 1, nhmvbs + 1,
   sizeof(est->mvs[0][0]));
  est->refine_grid = (od_mv_grid_pt **)od_malloc_2d(nvmvbs + 1, nhmvbs + 1,
   sizeof(est->refine_grid[0][0]));
  est->dp_nodes = (od_mv_dp_node *)_ogg_malloc(
   sizeof(od_mv_dp_node)*(OD_MAXI(nhmvbs, nvmvbs) + 1));
  est->row_counts =
   (unsigned *)_ogg_malloc(sizeof(*est->row_counts)*(nvmvbs + 1));
  est->col_counts =
   (unsigned *)_ogg_malloc(sizeof(*est->col_counts)*(nhmvbs + 1));
  for (vy = 0; vy <= nvmvbs; vy++) {
    for (vx = 0; vx <= nhmvbs; vx++) {
      est->mvs[vy][vx].vx = vx;
      est->mvs[vy][vx].vy = vy;
      est->mvs[vy][vx].heapi = -1;
      enc->state.mv_grid[vy][vx].valid = 1;
    }
  }
  est->dec_heap = (od_mv_node **)_ogg_malloc(
   sizeof(*est->dec_heap)*(nvmvbs + 1)*(nhmvbs + 1));
  est->hit_bit = 0;
  est->mv_res_min = 0;
  est->flags = OD_MC_USE_CHROMA;
}

static void od_mv_est_clear(od_mv_est_ctx *est) {
  _ogg_free(est->dec_heap);
  _ogg_free(est->col_counts);
  _ogg_free(est->row_counts);
  _ogg_free(est->dp_nodes);
  od_free_2d(est->refine_grid);
  od_free_2d(est->mvs);
  od_free_2d(est->sad_cache[0]);
  od_free_2d(est->sad_cache[1]);
}

/*STAGE 1: INITIAL MV ESTIMATES (via EPZS^2).*/

/*The amount to right shift the minimum error by when inflating it for
   computing the second maximum SAD threshold.*/
#define OD_MC_THRESH2_SCALE_BITS (3)

/*The vector offsets in the X direction for each search site in the various
   patterns.*/
static const int OD_SITE_DX[13] = {
  -1, 0, 1, -1, 0, 1, -1, 0, 1, -2, 0, 2, 0
};
/*The vector offsets in the Y direction for each search site in the various
   patterns.*/
static const int OD_SITE_DY[13] = {
  -1, -1, -1, 0, 0, 0, 1, 1, 1, 0, -2, 0, 2
};

/*Set to do the initial search with a square search pattern instead of the 3-D
   predict hexagon state machine.*/
#undef OD_USE_SQUARE_SEARCH
/*Set to use logarithmicly scaled search radius instead of a fixed search
   radius during motion vector refinement.*/
#undef OD_USE_LOGARITHMIC_REFINEMENT

#if defined(OD_USE_SQUARE_SEARCH) || defined(OD_USE_LOGARITHMIC_REFINEMENT)

/*The number of sites to search of each boundary condition in the square
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_SQUARE_NSITES[11] = { 8, 5, 5, 0, 5, 3, 3, 0, 5, 3, 3 };
/*The list of sites to search for each boundary condition in the square
   pattern.*/
static const od_pattern OD_SQUARE_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 0, 1, 2, 3, 5, 6, 7, 8 },
  /*-32 == dx,        -32 < dy < 31*/
  { 1, 2, 5, 7, 8 },
  /*       dx == 31,  -32 < dy < 31*/
  { 0, 1, 3, 6, 7 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 3, 5, 6, 7, 8 },
  /*-32 == dx,       -32 == dy*/
  { 5, 7, 8 },
  /*       dx == 31, -32 == dy*/
  { 3, 6, 7 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 0, 1, 2, 3, 5 },
  /*-32 == dx,              dy == 31*/
  { 1, 2, 5 },
  /*       dx == 31,        dy == 31*/
  { 0, 1, 3 }
};

#endif

/*The number of sites to search of each boundary condition in the diamond
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_DIAMOND_NSITES[11] = { 4, 3, 3, 0, 3, 2, 2, 0, 3, 2, 2 };
/*The list of sites to search for each boundary condition in the square
   pattern.*/
static const od_pattern OD_DIAMOND_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 1, 3, 5, 7 },
  /*-32 == dx,        -32 < dy < 31*/
  { 1, 5, 7 },
  /*       dx == 31,  -32 < dy < 31*/
  { 1, 3, 7 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 3, 5, 7 },
  /*-32 == dx,       -32 == dy*/
  { 5, 7 },
  /*       dx == 31, -32 == dy*/
  { 3, 7 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 1, 3, 5 },
  /*-32 == dx,              dy == 31*/
  { 1, 5 },
  /*       dx == 31,        dy == 31*/
  { 1, 3 }
};

#if !defined(OD_USE_SQUARE_SEARCH)

/*The number of sites to search of each boundary condition in the horizontal
   hex pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_HHEX_NSITES[11] = { 6, 3, 3, 0, 4, 2, 2, 0, 4, 2, 2 };
/*The list of sites to search for each boundary condition in the horizontal
   hex pattern.*/
static const od_pattern OD_HHEX_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 0, 2, 6, 8, 9, 11 },
  /*-32 == dx,        -32 < dy < 31*/
  { 2, 8, 11 },
  /*       dx == 31,  -32 < dy < 31*/
  { 0, 6, 9 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 6, 8, 9, 11 },
  /*-32 == dx,       -32 == dy*/
  { 8, 11 },
  /*       dx == 31, -32 == dy*/
  { 6, 9 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 0, 2, 9, 11 },
  /*-32 == dx,              dy == 31*/
  { 2, 11 },
  /*     dx == 31,          dy == 31*/
  { 0, 9 }
};

/*The number of sites to search of each boundary condition in the vertical hex
   pattern.
  Bit flags for the boundary conditions are as follows:
  1: -32 == dx
  2: dx == 31
  4: -32 == dy
  8: dy == 31*/
static const int OD_VHEX_NSITES[11] = { 6, 4, 4, 0, 3, 2, 2, 0, 3, 2, 2 };
/*The list of sites to search for each boundary condition in the vertical hex
   pattern.*/
static const od_pattern OD_VHEX_SITES[11] = {
  /* -32 < dx < 31,   -32 < dy < 31*/
  { 0, 2, 6, 8, 10, 12 },
  /*-32 == dx,        -32 < dy < 31*/
  { 2, 8, 10, 12 },
  /*       dx == 31,  -32 < dy < 31*/
  { 0, 6, 10, 12 },
  /*-32 == dx == 31,  -32 < dy < 31*/
  { -1 },
  /* -32 < dx < 31,  -32 == dy*/
  { 6, 8, 12 },
  /*-32 == dx,       -32 == dy*/
  { 8, 12 },
  /*       dx == 31, -32 == dy*/
  { 6, 12 },
  /*-32 == dx == 31, -32 == dy*/
  { -1 },
  /* -32 < dx < 31,         dy == 31*/
  { 0, 2, 10 },
  /*-32 == dx,              dy == 31*/
  { 2, 10 },
  /*       dx == 31,        dy == 31*/
  { 0, 10 }
};

#endif

#if defined(OD_USE_SQUARE_SEARCH)
/*The search state indicating we found a local minimum.*/
# define OD_SEARCH_STATE_DONE (1)

/*The number of sites in the pattern to use for each state.*/
static const int *const OD_SEARCH_NSITES[1] = {
  OD_SQUARE_NSITES
};

/*The sites in the pattern to use for each state.*/
static const od_pattern *const OD_SEARCH_SITES[1] = {
  OD_SQUARE_SITES
};

/*The successor state given the current state and the terminating site.*/
static const int OD_SEARCH_STATES[1][13] = {
  /*Just use a square pattern for the whole search.*/
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, -1, -1, -1 }
};
#else
/*The state machine used to select which patterns to use in the gradient
   descent BMA search.
  These are based on the procedure described in
  @ARTICLE{TP06,
    author="Tsung-Han Tsai and Yu-Nan Pan",
    title="A Novel 3-D Predict Hexagon Search Algorithm for Fast Block Motion
     Estimation on H.264 Video Coding",
    journal="{IEEE} Transactions on Circuits and Systems for Video Technology",
    volume=16,
    number=12,
    pages="1542--1549",
    month=Dec,
    year=2006
  }*/

/*The search state indicating we found a local minimum.*/
# define OD_SEARCH_STATE_DONE (6)

/*The number of sites in the pattern to use for each state.*/
static const int *const OD_SEARCH_NSITES[6] = {
  OD_DIAMOND_NSITES,
  OD_DIAMOND_NSITES,
  OD_DIAMOND_NSITES,
  OD_HHEX_NSITES,
  OD_VHEX_NSITES,
  OD_DIAMOND_NSITES
};

/*The sites in the pattern to use for each state.*/
static const od_pattern *const OD_SEARCH_SITES[6] = {
  OD_DIAMOND_SITES,
  OD_DIAMOND_SITES,
  OD_DIAMOND_SITES,
  OD_HHEX_SITES,
  OD_VHEX_SITES,
  OD_DIAMOND_SITES
};

/*The successor state given the current state and the terminating site.*/
static const int OD_SEARCH_STATES[6][13] = {
  /*Start with a small diamond in the first step.*/
  { -1, 2, -1, 1, 6, 1, -1, 2, -1, -1, -1, -1, -1 },
  /*Use a small diamond for the second step, too, but remember if we took a
     horizontal step in the first step...*/
  { -1, 3, -1, 3, 6, 3, -1, 3, -1, -1, -1, -1, -1 },
  /*...or a vertical one.*/
  { -1, 4, -1, 4, 6, 4, -1, 4, -1, -1, -1, -1, -1 },
  /*Then switch to a horizontal hex pattern for all remaining steps...*/
  { 3, -1, 3, -1, 5, -1, 3, -1, 3, 3, -1, 3, -1 },
  /*...or a vertical hex pattern, depending.*/
  { 4, -1, 4, -1, 5, -1, 4, -1, 4, -1, 4, -1, 4 },
  /*And revert back to a small diamond for the last step.*/
  { -1, 6, -1, 6, 6, 6, -1, 6, -1, -1, -1, -1, -1 }
};
#endif

/*Clear the cache of motion vectors we've examined.*/
static void od_mv_est_clear_hit_cache(od_mv_est_ctx *est) {
  if (est->hit_bit++ == 0) memset(est->hit_cache, 0, sizeof(est->hit_cache));
  else est->hit_bit &= UCHAR_MAX;
}

/*Test if a motion vector has been examined.*/
static int od_mv_est_is_hit(od_mv_est_ctx *est, int mvx, int mvy) {
  return est->hit_cache[mvy + 32][mvx + 32] == est->hit_bit;
}

/*Mark a motion vector examined.*/
static void od_mv_est_set_hit(od_mv_est_ctx *est, int mvx, int mvy) {
  est->hit_cache[mvy + 32][mvx + 32] = (unsigned char)est->hit_bit;
}

/*Estimated rate (in units of OD_BITRES) of a MV component of a given length.
  TODO: This should be replaced with more realistic estimates.*/
static const int OD_MV_EST_RATE[256] = {
    8,  32,  32,  48,  48,  48,  48,  64,
   64,  64,  64,  64,  64,  64,  64,  80,
   80,  80,  80,  80,  80,  80,  80,  80,
   80,  80,  80,  80,  80,  80,  80,  96,
   96,  96,  96,  96,  96,  96,  96,  96,
   96,  96,  96,  96,  96,  96,  96,  96,
   96,  96,  96,  96,  96,  96,  96,  96,
   96,  96,  96,  96,  96,  96,  96, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 112,
  112, 112, 112, 112, 112, 112, 112, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 128,
  128, 128, 128, 128, 128, 128, 128, 144
};

/*Estimate the number of bits that will be used to encode the given MV and its
   predictor.*/
static int od_mv_est_bits(int dx, int dy, int predx, int predy) {
  int pdx;
  int pdy;
  pdx = OD_MINI(abs(dx - predx), 255);
  pdy = OD_MINI(abs(dy - predy), 255);
  dx = OD_MINI(abs(dx), 255);
  dy = OD_MINI(abs(dy), 255);
  return (1 << OD_BITRES) + OD_MINI(OD_MV_EST_RATE[dx] + OD_MV_EST_RATE[dy],
   OD_MV_EST_RATE[pdx] + OD_MV_EST_RATE[pdy]);
}

/*Computes the SAD of a whole-pel BMA block with the given parameters.*/
static ogg_int32_t od_mv_est_bma_sad8(od_mv_est_ctx *est,
 int ref, int bx, int by, int mvx, int mvy, int log_mvb_sz) {
  od_state *state;
  od_img_plane *iplane;
  ogg_int32_t ret;
  int refi;
  int pmvx;
  int pmvy;
  int pbx;
  int pby;
  int dx;
  int dy;
  state = &est->enc->state;
  refi = state->ref_imgi[ref];
  iplane = state->ref_imgs[refi].planes + 0;
  pmvx = OD_DIV_POW2_RE(mvx << 1, iplane->xdec);
  pmvy = OD_DIV_POW2_RE(mvy << 1, iplane->ydec);
  pbx = (bx + (1 << iplane->xdec) - 1) & ~((1 << iplane->xdec) - 1);
  pby = (by + (1 << iplane->ydec) - 1) & ~((1 << iplane->ydec) - 1);
  dx = (pbx << 1 >> iplane->xdec) + pmvx;
  dy = (pby << 1 >> iplane->ydec) + pmvy;
  ret = od_enc_sad8(est->enc, iplane->data + dy*iplane->ystride + dx,
   iplane->ystride << 1, 2, 0, pbx, pby, log_mvb_sz + 2);
  if (est->flags & OD_MC_USE_CHROMA) {
    int pli;
    for (pli = 1; pli < state->io_imgs[OD_FRAME_INPUT].nplanes; pli++) {
      iplane = state->ref_imgs[refi].planes + pli;
      pmvx = OD_DIV_POW2_RE(mvx << 1, iplane->xdec);
      pmvy = OD_DIV_POW2_RE(mvy << 1, iplane->ydec);
      pbx = (bx + (1 << iplane->xdec) - 1) & ~((1 << iplane->xdec) - 1);
      pby = (by + (1 << iplane->ydec) - 1) & ~((1 << iplane->ydec) - 1);
      dx = (pbx << 1 >> iplane->xdec) + pmvx;
      dy = (pby << 1 >> iplane->ydec) + pmvy;
      ret += od_enc_sad8(est->enc, iplane->data + dy*iplane->ystride + dx,
       iplane->ystride << 1, 2, pli, pbx, pby, log_mvb_sz + 2) >>
       OD_MC_CHROMA_SCALE;
    }
  }
  return ret;
}

/*Computes the SAD of a block with the given parameters.*/
static ogg_int32_t od_mv_est_sad8(od_mv_est_ctx *est,
 int ref, int vx, int vy, int oc, int s, int log_mvb_sz) {
  unsigned char __attribute((aligned(16))) pred[16][16];
  od_state *state;
  ogg_int32_t ret;
  state = &est->enc->state;
  od_state_pred_block_from_setup(state, pred[0], sizeof(pred[0]), ref, 0,
   vx, vy, oc, s, log_mvb_sz);
  ret = od_enc_sad8(est->enc, pred[0], sizeof(pred[0]), 1, 0,
   (vx - 2) << 2, (vy - 2) << 2, log_mvb_sz + 2);
  if (est->flags & OD_MC_USE_CHROMA) {
    int pli;
    for (pli = 1; pli < state->io_imgs[OD_FRAME_INPUT].nplanes; pli++) {
      od_state_pred_block_from_setup(state, pred[0], sizeof(pred[0]), ref, pli,
       vx, vy, oc, s, log_mvb_sz);
      ret += od_enc_sad8(est->enc, pred[0], sizeof(pred[0]), 1, pli,
       (vx - 2) << 2, (vy - 2) << 2, log_mvb_sz + 2) >> OD_MC_CHROMA_SCALE;
    }
  }
  return ret;
}

#if defined(OD_LOGGING_ENABLED)
/*Checks to make sure our current mv_rate and sad values are correct.
  This is used for debugging only.*/
void od_mv_est_check_rd_block_state(od_mv_est_ctx *est,
 int ref, int vx, int vy, int log_mvb_sz) {
  od_state *state;
  int half_mvb_sz;
  state = &est->enc->state;
  half_mvb_sz = 1 << log_mvb_sz >> 1;
  if (log_mvb_sz > 0
   && state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
    od_mv_est_check_rd_block_state(est, ref, vx, vy, log_mvb_sz - 1);
    od_mv_est_check_rd_block_state(est, ref,
     vx + half_mvb_sz, vy, log_mvb_sz - 1);
    od_mv_est_check_rd_block_state(est, ref,
     vx, vy + half_mvb_sz, log_mvb_sz - 1);
    od_mv_est_check_rd_block_state(est, ref,
     vx + half_mvb_sz, vy + half_mvb_sz, log_mvb_sz - 1);
  }
  else {
    od_mv_node *block;
    ogg_int32_t sad;
    int oc;
    int s;
    block = est->mvs[vy] + vx;
    if (block->log_mvb_sz != log_mvb_sz) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
       "Failure at node (%i, %i): log_mvb_sz should be %i (is %i)",
       vx, vy, log_mvb_sz, block->log_mvb_sz));
    }
    if (log_mvb_sz < 2) {
      int mask;
      int s1vx;
      int s1vy;
      int s3vx;
      int s3vy;
      mask = (1 << (log_mvb_sz + 1)) - 1;
      oc = !!(vx & mask);
      if (vy & mask) oc = 3 - oc;
      if (block->oc != oc) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
         "Failure at node (%i, %i): oc should be %i (is %i)",
         vx, vy, oc, block->oc));
      }
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
    if (block->s != s) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
       "Failure at node (%i, %i): s should be %i (is %i)",
       vx, vy, s, block->s));
    }
    sad = od_mv_est_sad8(est, ref, vx, vy, oc, s, log_mvb_sz);
    if (block->sad != sad) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
       "Failure at node (%i, %i): sad should be %i (is %i)",
       vx, vy, sad, block->sad));
    }
  }
}

/*Checks to make sure our current mv_rate and sad values are correct.
  This is used for debugging only.*/
void od_mv_est_check_rd_state(od_mv_est_ctx *est, int ref, int mv_res) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  for (vy = 0; vy < nvmvbs; vy += 4) {
    for (vx = 0; vx < nhmvbs; vx += 4) {
      od_mv_est_check_rd_block_state(est, ref, vx, vy, 2);
    }
  }
  for (vy = 0; vy < nvmvbs; vy++) {
    for (vx = 0; vx < nhmvbs; vx++) {
      od_mv_grid_pt *mvg;
      od_mv_node *mv;
      int pred[2];
      int mv_rate;
      mvg = state->mv_grid[vy] + vx;
      if (!mvg->valid) continue;
      mv = est->mvs[vy] + vx;
      if (vx >= 2 && vx <= nhmvbs - 2 && vy >= 2 && vy <= nvmvbs - 2) {
        od_state_get_predictor(state, pred, vx, vy,
         OD_MC_LEVEL[vy & 3][vx & 3], mv_res);
        mv_rate = od_mv_est_bits(mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res,
         pred[0], pred[1]);
      }
      else pred[0] = pred[1] = mv_rate = 0;
      if (mv_rate != mv->mv_rate) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
         "Failure at node (%i, %i): mv_rate should be %i (is %i)",
         vx, vy, mv_rate, mv->mv_rate));
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN,
         "Predictor was: (%i, %i)   MV was: (%i, %i)",
         pred[0], pred[1], mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res));
      }
    }
  }
}
#endif

#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
static const unsigned char OD_YCbCr_MVCAND[3] = { 210, 16, 214 };
#endif

static void od_mv_est_init_mv(od_mv_est_ctx *est, int ref, int vx, int vy) {
  od_state *state;
  od_mv_grid_pt *mvg;
  od_mv_node *mv;
  od_mv_node *cneighbors[4];
  od_mv_node *pneighbors[4];
  ogg_int32_t t2;
  ogg_int32_t best_sad;
  ogg_int32_t best_cost;
  int best_rate;
  int cands[6][2];
  int best_vec[2];
  int a[4][2];
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int bx;
  int by;
  int ncns;
  int mvxmin;
  int mvxmax;
  int mvymin;
  int mvymax;
  int candx;
  int candy;
  int predx;
  int predy;
  int ci;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  int animating;
  int x0;
  int y0;
#endif
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Initial search for MV (%i, %i):", vx, vy));
  state = &est->enc->state;
  mv = est->mvs[vy] + vx;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  mvg = state->mv_grid[vy] + vx;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  animating = daala_granule_basetime(state, state->cur_time) == ANI_FRAME;
  if (animating) {
    od_state_mc_predict(state, ref);
    od_state_fill_vis(state);
    x0 = ((vx - 2) << 3) + (OD_UMV_PADDING << 1);
    y0 = ((vy - 2) << 3) + (OD_UMV_PADDING << 1);
  }
#endif
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Level %i (%ix%i block)", level, mvb_sz << 2, mvb_sz << 2));
  bx = (vx - 2) << 2;
  by = (vy - 2) << 2;
  mvxmin = OD_MAXI(bx - (mvb_sz << 2) -32, -16) - (bx - (mvb_sz << 2));
  mvxmax = OD_MINI(bx + (mvb_sz << 2) + 32, state->frame_width + 16)
   - (bx + (mvb_sz << 2)) - 1;
  mvymin = OD_MAXI(by - (mvb_sz << 2) - 32, -16) - (by - (mvb_sz << 2));
  mvymax = OD_MINI(by + (mvb_sz << 2) + 32, state->frame_height + 16)
   - (by + (mvb_sz << 2)) - 1;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "(%i, %i): Search range: [%i, %i]x[%i, %i]",
   bx, by, mvxmin, mvymin, mvxmax, mvymax));
  bx -= mvb_sz << 1;
  by -= mvb_sz << 1;
  ncns = 4;
  if (level == 0) {
    cneighbors[0] = est->mvs[vy - 4] + vx - 4;
    cneighbors[1] = est->mvs[vy - 4] + vx;
    cneighbors[2] = est->mvs[vy - 4] + vx + 4;
    cneighbors[3] = est->mvs[vy] + vx - 4;
    pneighbors[0] = est->mvs[vy - 4] + vx;
    pneighbors[1] = est->mvs[vy] + vx - 4;
    pneighbors[2] = est->mvs[vy] + vx + 4;
    pneighbors[3] = est->mvs[vy + 4] + vx;
  }
  else {
    if (level & 1) {
      pneighbors[0] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      pneighbors[1] = est->mvs[vy - mvb_sz] + vx + mvb_sz;
      pneighbors[2] = est->mvs[vy + mvb_sz] + vx - mvb_sz;
      pneighbors[3] = est->mvs[vy + mvb_sz] + vx + mvb_sz;
      memcpy(cneighbors, pneighbors, sizeof(cneighbors));
    }
    else {
      pneighbors[0] = est->mvs[vy - mvb_sz] + vx;
      pneighbors[1] = est->mvs[vy] + vx - mvb_sz;
      pneighbors[2] = est->mvs[vy] + vx + mvb_sz;
      pneighbors[3] = est->mvs[vy + mvb_sz] + vx;
      cneighbors[0] = pneighbors[0];
      cneighbors[1] = pneighbors[1];
      /*NOTE: Only one of these candidatss can be excluded at a time, so
         there will always be at least 3.*/
      if (vx + mvb_sz > ((vx + 3) & ~3)) ncns--;
      else cneighbors[2] = pneighbors[2];
      if (vy + mvb_sz > ((vy + 3) & ~3)) ncns--;
      else cneighbors[ncns - 1] = pneighbors[3];
    }
  }
  /*Spatially correlated predictors (from the current frame):*/
  for (ci = 0; ci < ncns; ci++) {
    a[ci][0] = cneighbors[ci]->bma_mvs[0][ref][0];
    a[ci][1] = cneighbors[ci]->bma_mvs[0][ref][1];
    cands[ci][0] = OD_CLAMPI(mvxmin, a[ci][0], mvxmax);
    cands[ci][1] = OD_CLAMPI(mvymin, a[ci][1], mvymax);
  }
  /*Compute the median predictor:*/
  if (ncns > 3) {
    /*Median-of-4.*/
    OD_SORT2I(a[0][0], a[1][0]);
    OD_SORT2I(a[0][1], a[1][1]);
    OD_SORT2I(a[2][0], a[3][0]);
    OD_SORT2I(a[2][1], a[3][1]);
    OD_SORT2I(a[0][0], a[2][0]);
    OD_SORT2I(a[0][1], a[2][1]);
    OD_SORT2I(a[1][0], a[3][0]);
    OD_SORT2I(a[1][1], a[3][1]);
    predx = a[1][0] + a[2][0];
    predy = a[1][1] + a[2][1];
    candx = OD_CLAMPI(mvxmin, OD_DIV2_RE(predx), mvxmax);
    candy = OD_CLAMPI(mvymin, OD_DIV2_RE(predy), mvymax);
  }
  else {
    /*Median-of-3.*/
    OD_SORT2I(a[0][0], a[1][0]);
    OD_SORT2I(a[0][1], a[1][1]);
    OD_SORT2I(a[1][0], a[2][0]);
    OD_SORT2I(a[1][1], a[2][1]);
    OD_SORT2I(a[0][0], a[1][0]);
    OD_SORT2I(a[0][1], a[1][1]);
    predx = a[1][0] << 1;
    predy = a[1][1] << 1;
    candx = OD_CLAMPI(mvxmin, a[1][0], mvxmax);
    candy = OD_CLAMPI(mvymin, a[1][1], mvymax);
  }
  od_mv_est_clear_hit_cache(est);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  if (animating) {
    od_img_draw_line(&state->vis_img, x0, y0,
     x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
  }
#endif
  best_sad = od_mv_est_bma_sad8(est, ref, bx, by, candx, candy, log_mvb_sz);
  best_rate = od_mv_est_bits(candx << 1, candy << 1, predx, predy);
  best_cost = (best_sad << OD_ERROR_SCALE) + best_rate*est->lambda;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Median predictor: (%i, %i)   Cost: %i", candx, candy, best_cost));
  od_mv_est_set_hit(est, candx, candy);
  best_vec[0] = candx;
  best_vec[1] = candy;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Threshold: %i", est->thresh1[log_mvb_sz]));
  if (best_sad > est->thresh1[log_mvb_sz]) {
    ogg_int32_t sad;
    ogg_int32_t cost;
    int rate;
    /*Compute the early termination threshold for set B.*/
    t2 = mv->bma_sad;
    for (ci = 0; ci < ncns; ci++) {
      int log_cnb_sz;
      int cvx;
      int cvy;
      cvx = cneighbors[ci]->vy;
      cvy = cneighbors[ci]->vx;
      log_cnb_sz = (4 - OD_MC_LEVEL[cvy & 3][cvx & 3]) >> 1;
      t2 = OD_MINI(t2,
       cneighbors[ci]->bma_sad >> ((log_cnb_sz - log_mvb_sz) << 1));
    }
    t2 = t2 + (t2 >> OD_MC_THRESH2_SCALE_BITS) + est->thresh2_offs[log_mvb_sz];
    /*Constant velocity predictor:*/
    cands[ncns][0] =
     OD_CLAMPI(mvxmin, mv->bma_mvs[1][ref][0], mvxmax);
    cands[ncns][1] =
     OD_CLAMPI(mvymin, mv->bma_mvs[1][ref][1], mvymax);
    ncns++;
    /*Zero predictor.*/
    cands[ncns][0] = 0;
    cands[ncns][1] = 0;
    ncns++;
    /*Examine the candidates in Set B.*/
    for (ci = 0; ci < ncns; ci++) {
      candx = cands[ci][0];
      candy = cands[ci][1];
      if (od_mv_est_is_hit(est, candx, candy)) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Set B predictor %i: (%i, %i) ...Skipping.", ci, candx, candy));
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "...Skipping."));
        continue;
      }
      od_mv_est_set_hit(est, candx, candy);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (animating) {
        od_img_draw_line(&state->vis_img, x0, y0,
         x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
      }
#endif
      sad = od_mv_est_bma_sad8(est, ref, bx, by, candx, candy, log_mvb_sz);
      rate = od_mv_est_bits(candx << 1, candy << 1, predx, predy);
      cost = (sad << OD_ERROR_SCALE) + rate*est->lambda;
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Set B predictor %i: (%i, %i)    Cost: %i", ci, candx, candy, cost));
      if (cost < best_cost) {
        best_sad = sad;
        best_rate = rate;
        best_cost = cost;
        best_vec[0] = candx;
        best_vec[1] = candy;
      }
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Threshold: %i", t2));
    if (best_sad > t2) {
      /*Constant velocity predictors from the previous frame:*/
      for (ci = 0; ci < 4; ci++) {
        cands[ci][0] =
         OD_CLAMPI(mvxmin, pneighbors[ci]->bma_mvs[1][ref][0], mvxmax);
        cands[ci][1] =
         OD_CLAMPI(mvymin, pneighbors[ci]->bma_mvs[1][ref][1], mvymax);
      }
      /*The constant acceleration predictor:*/
      cands[4][0] = OD_CLAMPI(mvxmin,
       OD_DIV_ROUND_POW2(mv->bma_mvs[1][ref][0]*est->mvapw[ref][0]
       - mv->bma_mvs[2][ref][0]*est->mvapw[ref][1], 16, 0x8000), mvxmax);
      cands[4][1] = OD_CLAMPI(mvymin,
       OD_DIV_ROUND_POW2(mv->bma_mvs[1][ref][1]*est->mvapw[ref][0]
       - mv->bma_mvs[2][ref][1]*est->mvapw[ref][1], 16, 0x8000), mvymax);
      /*Examine the candidates in Set C.*/
      for (ci = 0; ci < 5; ci++) {
        candx = cands[ci][0];
        candy = cands[ci][1];
        if (od_mv_est_is_hit(est, candx, candy)) {
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Set C predictor %i: (%i, %i) ...Skipping.", ci, candx, candy));
          continue;
        }
        /*if (od_mv_est_is_hit(est, candx, candy)) continue;*/
        od_mv_est_set_hit(est, candx, candy);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
        if (animating) {
          od_img_draw_line(&state->vis_img, x0, y0,
           x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
        }
#endif
        sad = od_mv_est_bma_sad8(est, ref, bx, by, candx, candy, log_mvb_sz);
        rate = od_mv_est_bits(candx << 1, candy << 1, predx, predy);
        cost = (sad << OD_ERROR_SCALE) + rate*est->lambda;
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Set C predictor %i: (%i, %i)    Cost: %i", ci, candx, candy, cost));
        if (cost < best_cost) {
          best_sad = sad;
          best_rate = rate;
          best_cost = cost;
          best_vec[0] = candx;
          best_vec[1] = candy;
        }
      }
      /*Use the same threshold for Set C as in Set B.*/
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Threshold: %i", t2));
      if (best_sad > t2) {
        const int *pattern;
        int mvstate;
        int best_site;
        int nsites;
        int sitei;
        int site;
        int b;
        /*Gradient descent pattern search.*/
        mvstate = 0;
        for (;;) {
          best_site = 4;
          b = (best_vec[0] <= mvxmin) | (best_vec[0] >= mvxmax) << 1 |
           (best_vec[1] <= mvymin) << 2 | (best_vec[1] >= mvymax) << 3;
          pattern = OD_SEARCH_SITES[mvstate][b];
          nsites = OD_SEARCH_NSITES[mvstate][b];
          for (sitei = 0; sitei < nsites; sitei++) {
            site = pattern[sitei];
            candx = best_vec[0] + OD_SITE_DX[site];
            candy = best_vec[1] + OD_SITE_DY[site];
            /*For the large search patterns, our simple mechanism to move
               bounds checking out of the inner loop doesn't work (it would
               need 2 more bits, or 4 times as much table storage, and require
               4 extra compares, when there are often fewer than 4 sites).
              If the displacement is larger than +/-1 in any direction (which
               happens when site > 8), check the bounds explicitly.*/
            if (site > 8 && (candx < mvxmin || candx > mvxmax
             || candy < mvymin || candy > mvymax)) {
              continue;
            }
            if (od_mv_est_is_hit(est, candx, candy)) {
              OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
               "Pattern search %i: (%i, %i) ...Skipping.",
               site, candx, candy));
              continue;
            }
            od_mv_est_set_hit(est, candx, candy);
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
            if (animating) {
              od_img_draw_line(&state->vis_img, x0, y0,
               x0 + (candx << 1), y0 + (candy << 1), OD_YCbCr_MVCAND);
            }
#endif
            sad = od_mv_est_bma_sad8(est,
             ref, bx, by, candx, candy, log_mvb_sz);
            rate = od_mv_est_bits(candx << 1, candy << 1, predx, predy);
            cost = (sad << OD_ERROR_SCALE) + rate*est->lambda;
            OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
             "Pattern search %i: (%i, %i)    Cost: %i",
             site, candx, candy, cost));
            if (cost < best_cost) {
              best_sad = sad;
              best_rate = rate;
              best_cost = cost;
              best_site = site;
            }
          }
          mvstate = OD_SEARCH_STATES[mvstate][best_site];
          best_vec[0] += OD_SITE_DX[best_site];
          best_vec[1] += OD_SITE_DY[best_site];
          if (mvstate == OD_SEARCH_STATE_DONE) break;
        }
      }
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Finished. Best vector: (%i, %i)  Best cost %i",
   best_vec[0], best_vec[1], best_cost));
  mv->bma_mvs[0][ref][0] = best_vec[0];
  mv->bma_mvs[0][ref][1] = best_vec[1];
  mvg->mv[0] = best_vec[0] << 3;
  mvg->mv[1] = best_vec[1] << 3;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  if (animating) {
    char iter_label[16];
    const od_offset *anc;
    od_mv_grid_pt *amvg;
    int nanc;
    int ai;
    int ax;
    int ay;
    mvg->valid = 1;
    nanc = OD_NANCESTORS[vy & 3][vx & 3];
    anc = OD_ANCESTORS[vy & 3][vx & 3];
    for (ai = 0; ai < nanc; ai++) {
      ax = vx + anc[ai][0];
      if (ax < 0 || ax > ((state->nhmbs + 1) << 2)) continue;
      ay = vy + anc[ai][1];
      if (ay < 0 || ay > ((state->nvmbs + 1) << 2)) continue;
      amvg = state->mv_grid[ay] + ax;
      amvg->valid = 1;
    }
    sprintf(iter_label, "ani%08i", state->ani_iter++);
    od_state_dump_img(state, &state->vis_img, iter_label);
  }
#endif
  mv->bma_sad = best_sad;
  mv->mv_rate = best_rate;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Initialized MV (%2i, %2i): (%3i, %3i), SAD: %i",
   vx, vy, best_vec[0], best_vec[1], best_sad));
  /*od_state_get_predictor(state, a[0], vx, vy, level, 2);
  if (a[0][0] != predx || a[0][1] != predy) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Failure in MV predictor init: (%i, %i) != (%i, %i)",
     a[0][0], a[0][1], predx, predy));
  }
  mv->mv_rate = od_mv_est_bits(mvg->mv[0] >> 2, mvg->mv[1] >> 2,
   a[0][0], a[0][1]);
  if (mv->mv_rate != best_rate) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Failure in MV rate init: %i != %i",
     mv->mv_rate, best_rate));
  }*/
}

static void od_mv_est_init_mvs(od_mv_est_ctx *est, int ref) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  /*Move the motion vector predictors back a frame.*/
  for (vy = 2; vy <= nvmvbs - 2; vy++) {
    for (vx = 2; vx <= nhmvbs - 2; vx++) {
      od_mv_node *mv;
      mv = est->mvs[vy] + vx;
      OD_MOVE(mv->bma_mvs + 1, mv->bma_mvs + 0, 2);
    }
  }
  /*We initialize MVs a MVB at a time for cache coherency.
    Proceeding level-by-level would involve less branching and less complex
     code, but the SADs dominate.*/
  for (vy = 4; vy <= nvmvbs; vy += 4) {
    for (vx = 4; vx <= nhmvbs; vx += 4) {
      int b;
      /*Keep track of what edges we're on.
        We only need to know about the bottom/right to start with, since the
         top/left are already initialized (or do not need to be).*/
      b = (vx < nhmvbs) << 1 | (vy < nvmvbs) << 3;
      /*Initialization MUST proceed in order by level to ensure the necessary
         predictors are available.
        Order within a level does not matter, except for level 0.
        Level 0 is the only level to use predictors outside the current MVB,
         and must proceed in raster order.*/
      /*Level 0 vertex.*/
      if ((b & 0xA) == 0xA) od_mv_est_init_mv(est, ref, vx, vy);
      if (est->level_max < 1) continue;
      /*Level 1 vertex.*/
      od_mv_est_init_mv(est, ref, vx - 2, vy - 2);
      if (est->level_max < 2) continue;
      /*Level 2 vertices.*/
      if (b & 2) od_mv_est_init_mv(est, ref, vx, vy - 2);
      if (b & 8) od_mv_est_init_mv(est, ref, vx - 2, vy);
      if (est->level_max < 3) continue;
      /*Level 3 vertices.*/
      /*Add in flags for the top/left edges.*/
      b |= (vx > 4) | (vy > 4) << 2;
      if (b & 4) {
        if (b & 1) od_mv_est_init_mv(est, ref, vx - 3, vy - 3);
        if (b & 2) od_mv_est_init_mv(est, ref, vx - 1, vy - 3);
      }
      if (b & 8) {
        if (b & 1) od_mv_est_init_mv(est, ref, vx - 3, vy - 1);
        if (b & 2) od_mv_est_init_mv(est, ref, vx - 1, vy - 1);
      }
      if (est->level_max < 4) continue;
      /*Level 4 vertices.*/
      if (b & 1) od_mv_est_init_mv(est, ref, vx - 3, vy - 2);
      if (b & 2) od_mv_est_init_mv(est, ref, vx - 1, vy - 2);
      if (b & 4) {
        od_mv_est_init_mv(est, ref, vx - 2, vy - 3);
        if (b & 2) od_mv_est_init_mv(est, ref, vx, vy - 3);
      }
      if (b & 8) {
        od_mv_est_init_mv(est, ref, vx - 2, vy - 1);
        if (b & 1) od_mv_est_init_mv(est, ref, vx - 3, vy);
        if (b & 2) {
          od_mv_est_init_mv(est, ref, vx, vy - 1);
          od_mv_est_init_mv(est, ref, vx - 1, vy);
        }
      }
    }
  }
}

/*STAGE 2: DECIMATION.*/

/*Merging domains.
  These are stored as lists of offsets to the vertices in the domain.
  Note that vertices in the merging domain must appear in order from finest
   scale (largest level) to coarsest (smallest level).
  Each list ends with the vertex (0, 0), the actual vertex be decimated.*/
/*Level 4 vertex:
            4*/
static const od_offset OD_MERGEDOM4[1] = {
  { 0, 0 },
};

/*Level 3 vertex:
            4
          4-3-4
            4*/
static const od_offset OD_MERGEDOM3[5] = {
  {  0, -1 }, { -1,  0 }, {  1,  0 }, {  0,  1 }, {  0,  0 }
};

/*Level 2 vertex:
          4   4
          |   |
        4-3-4-3-4
          | | |
          4-2-4
          | | |
        4-3-4-3-4
          |   |
          4   4*/
static const od_offset OD_MERGEDOM2[17] = {
  { -1, -2 }, {  1, -2 }, { -2, -1 }, {  0, -1 }, {  2, -1 }, { -1,  0 },
  {  1,  0 }, { -2,  1 }, {  0,  1 }, {  2,  1 }, { -1,  2 }, {  1,  2 },
  { -1, -1 }, {  1, -1 }, { -1,  1 }, {  1,  1 }, {  0,  0 }
};

/*Level 1 vertex:
          4   4
          |   |
        4-3-4-3-4
          | | |
      4   4-2-4   4
      |   | | |   |
    4-3-4-3-4-3-4-3-4
      | | | | | | |
      4-2-4-1-4-2-4
      | | | | | | |
    4-3-4-3-4-3-4-3-4
      |   | | |   |
      4   4-2-4   4
          | | |
        4-3-4-3-4
          |   |
          4   4
*/
static const od_offset OD_MERGEDOM1[49] = {
  { -1, -4 }, {  1, -4 }, { -2, -3 }, {  0, -3 }, {  2, -3 }, { -3, -2 },
  { -1, -2 }, {  1, -2 }, {  3, -2 }, { -4, -1 }, { -2, -1 }, {  0, -1 },
  {  2, -1 }, {  4, -1 }, { -3,  0 }, { -1,  0 }, {  1,  0 }, {  3,  0 },
  { -4,  1 }, { -2,  1 }, {  0,  1 }, {  2,  1 }, {  4,  1 }, { -3,  2 },
  { -1,  2 }, {  1,  2 }, {  3,  2 }, { -2,  3 }, {  0,  3 }, {  2,  3 },
  { -1,  4 }, {  1,  4 }, { -1, -3 }, {  1, -3 }, { -3, -1 }, { -1, -1 },
  {  1, -1 }, {  3, -1 }, { -3,  1 }, { -1,  1 }, {  1,  1 }, {  3,  1 },
  { -1,  3 }, {  1,  3 }, {  0, -2 }, { -2,  0 }, {  2,  0 }, {  0,  2 },
  {  0,  0 }
};

/*The merging domain for a vertex, indexed by level-1.*/
static const od_offset *OD_MERGEDOM[4] = {
  OD_MERGEDOM1,
  OD_MERGEDOM2,
  OD_MERGEDOM3,
  OD_MERGEDOM4
};

/*Error support regions.
  These are the blocks whose SAD will change after decimating a vertex at a
   given level, assuming no other vertices in the mesh have been decimated.
  Vertices in the figures at a higher level than the one removed illustrate one
   possible configuration; there may be others.*/
struct od_mv_err_node {
  int dx;
  int dy;
  int log_mvb_sz;
};

/*Level 4 support:
          4-3-4
          |/|\|
          2-.-1
          |\|/|
          4-3-4*/
static const od_mv_err_node OD_ERRDOM4[4] = {
  { -1, -1, 0 }, {  0, -1, 0 }, { -1,  0, 0 }, {  0,  0, 0 }
};

/*Level 3 support:
          4-3-4
          |/|\|
        4-0-.-2-4
        |/|   |\|
        3-.   .-3
        |\|   |/|
        4-2-.-1-4
          |\|/|
          4-3-4*/
static const od_mv_err_node OD_ERRDOM3[9] = {
  { -1, -2, 0 }, {  0, -2, 0 }, { -2, -1, 0 }, {  1, -1, 0 },
  { -2,  0, 0 }, {  1,  0, 0 }, { -1,  1, 0 }, {  0,  1, 0 },
  { -1, -1, 1 }
};

/*Level 2 support:
        4-3-4-3-4
        |/|\|/|\|
      4-2-.-1-.-2-4
      |/|  /|\  |\|
      3-. / | \ .-3
      |\|/  |  \|/|
      4-0---.---0-4
      |/|\  |  /|\|
      3-. \ | / .-3
      |\|  \|/  |/|
      4-2-.-1-.-2-4
        |\|/|\|/|
        4-3-4-3-4*/
static const od_mv_err_node OD_ERRDOM2[20] = {
  { -2, -3, 0 }, { -1, -3, 0 }, {  0, -3, 0 }, {  1, -3, 0 },
  { -3, -2, 0 }, {  2, -2, 0 }, { -3, -1, 0 }, {  2, -1, 0 },
  { -3,  0, 0 }, {  2,  0, 0 }, { -3,  1, 0 }, {  2,  1, 0 },
  { -2,  2, 0 }, { -1,  2, 0 }, {  0,  2, 0 }, {  1,  2, 0 },
  { -2, -2, 1 }, {  0, -2, 1 }, { -2,  0, 1 }, {  0,  0, 1 }
};

/*Level 1 support:
        4-3-4-3-4
        |/|\|/|\|
      4-2-.-1-.-2-4
      |/|  /|\  |\|
    4-3-. / | \ .-3-4
    |/| |/  |  \| |\|
  4-2-.-0---.---0-.-2-4
  |/|  /|       |\  |\|
  3-. / |       | \ .-3
  |\|/  |       |  \|/|
  4-1---.       .---1-4
  |/|\  |       |  /|\|
  3-. \ |       | / .-3
  |\|  \|       |/  |/|
  4-2-.-0---.---0-.-2-4
    |\| |\  |  /| |/|
    4-3-. \ | / .-3-4
      |\|  \|/  |/|
      4-2-.-1-.-2-4
        |\|/|\|/|
        4-3-4-3-4*/
static const od_mv_err_node OD_ERRDOM1[37] = {
  { -2, -5, 0 }, { -1, -5, 0 }, {  0, -5, 0 }, {  1, -5, 0 },
  { -3, -4, 0 }, {  2, -4, 0 }, { -4, -3, 0 }, { -3, -3, 0 },
  {  2, -3, 0 }, {  3, -3, 0 }, { -5, -2, 0 }, {  4, -2, 0 },
  { -5, -1, 0 }, {  4, -1, 0 }, { -5,  0, 0 }, {  4,  0, 0 },
  { -5,  1, 0 }, {  4,  1, 0 }, { -4,  2, 0 }, { -3,  2, 0 },
  {  2,  2, 0 }, {  3,  2, 0 }, { -3,  3, 0 }, {  2,  3, 0 },
  { -2,  4, 0 }, { -1,  4, 0 }, {  0,  4, 0 }, {  1,  4, 0 },
  { -2, -4, 1 }, {  0, -4, 1 }, { -4, -2, 1 }, {  2, -2, 1 },
  { -4,  0, 1 }, {  2,  0, 1 }, { -2,  2, 1 }, {  0,  2, 1 },
  { -2, -2, 2 }
};

/*The number of blocks in each decimated error domain.*/
static const int OD_NERRDOM[4] = { 37, 20, 9, 4 };
/*The error domain for a vertex, indexed by level-1.*/
static const od_mv_err_node *OD_ERRDOM[4] = {
  OD_ERRDOM1,
  OD_ERRDOM2,
  OD_ERRDOM3,
  OD_ERRDOM4
};

/*Returns a negative value, 0, or a positive value, depending on whether
  -dd1/dr1 is less, equal or greater than -dd2/dr2.*/
static int od_mv_dddr_cmp(ogg_int32_t dd1, int dr1,
 ogg_int32_t dd2, int dr2) {
  ogg_int64_t diff;
  /*dr == 0 and dd != 0 should not be possible, but we check for it anyway just
     in case, to prevent a bug from trashing the whole optimization process.*/
  if (dr1 == 0) {
    return dr2 == 0 ? OD_SIGNI(dd2 - dd1) : (OD_SIGNI(dd1) << 1) - 1;
  }
  else if (dr2 == 0) return (OD_SIGNI(-dd2) << 1) + 1;
  diff = dd2*(ogg_int64_t)dr1 - dd1*(ogg_int64_t)dr2;
  return OD_SIGNI(diff);
}

/*Compare two nodes on the decimation heap.*/
static int od_mv_dec_cmp(od_mv_node *n1, od_mv_node *n2) {
  return od_mv_dddr_cmp(n1->dd, n1->dr, n2->dd, n2->dr);
}

/*Swap the two nodes on the decimation heap at indices p and q.*/
static void od_mv_dec_heap_swap(od_mv_node **heap, int p, int q) {
  od_mv_node *t;
  heap[p]->heapi = q;
  heap[q]->heapi = p;
  t = heap[p];
  heap[p] = heap[q];
  heap[q] = t;
}

/*Convert the list of nodes to be decimated to a heap.*/
static void od_mv_dec_heapify(od_mv_est_ctx *est) {
  od_mv_node **heap;
  int l;
  int r;
  int i;
  heap = est->dec_heap;
  l = est->dec_nheap >> 1;
  r = est->dec_nheap - 1;
  for (i = l; i-- > 0;) {
    int p;
    p = i;
    do {
      int q;
      q = (p << 1) + 1;
      if (q < r && od_mv_dec_cmp(heap[q], heap[q + 1]) >= 0) q++;
      if (od_mv_dec_cmp(heap[p], heap[q]) <= 0) break;
      od_mv_dec_heap_swap(heap, p, q);
      p = q;
    }
    while (p < l);
  }
}

/*Restore the heap structure at the given index by moving it down the heap.*/
static void od_mv_dec_heap_down(od_mv_est_ctx *est, int heapi) {
  od_mv_node **heap;
  int l;
  int r;
  int p;
  heap = est->dec_heap;
  l = est->dec_nheap >> 1;
  r = est->dec_nheap - 1;
  p = heapi;
  while (p < l) {
    int q;
    q = (p << 1) + 1;
    if (q < r && od_mv_dec_cmp(heap[q], heap[q + 1]) >= 0) q++;
    if (od_mv_dec_cmp(heap[p], heap[q]) <= 0) break;
    od_mv_dec_heap_swap(heap, p, q);
    p = q;
  }
}

/*Restore the heap structure at the given index by moving it up the heap.*/
static void od_mv_dec_heap_up(od_mv_est_ctx *est, int heapi) {
  od_mv_node **heap;
  int p;
  heap = est->dec_heap;
  p = heapi;
  while (p > 0) {
    int q;
    q = p;
    p = ((q + 1) >> 1) - 1;
    if (od_mv_dec_cmp(heap[p], heap[q]) <= 0) break;
    od_mv_dec_heap_swap(heap, p, q);
  }
}

/*Retrieve the item at the top of the heap.
  Returns NULL if there are no more nodes to decimate.*/
static od_mv_node *od_mv_dec_heap_delhead(od_mv_est_ctx *est) {
  od_mv_node *ret;
  if (est->dec_nheap <= 0) return NULL;
  ret = est->dec_heap[0];
  ret->heapi = -1;
  if (--est->dec_nheap > 0) {
    est->dec_heap[0] = est->dec_heap[est->dec_nheap];
    est->dec_heap[0]->heapi = 0;
    od_mv_dec_heap_down(est, 0);
  }
  return ret;
}

static void od_mv_dec_heap_del(od_mv_est_ctx *est, od_mv_node *node) {
  int heapi;
  heapi = node->heapi;
  if (heapi >= 0) {
    node->heapi = -1;
    est->dec_nheap--;
    if (est->dec_nheap > heapi) {
      est->dec_heap[heapi] = est->dec_heap[est->dec_nheap];
      est->dec_heap[heapi]->heapi = heapi;
      if (od_mv_dec_cmp(node, est->dec_heap[heapi]) >= 0) {
        od_mv_dec_heap_up(est, heapi);
      }
      else od_mv_dec_heap_down(est, heapi);
    }
    else est->dec_heap[est->dec_nheap] = NULL;
  }
}

/*Sets the dd and dr values of the given node, restoring the heap structure
   afterwards.*/
static void od_mv_dec_update(od_mv_est_ctx *est, od_mv_node *node,
 int dd, int dr) {
  int diff;
  diff = od_mv_dddr_cmp(dd, dr, node->dd, node->dr);
  node->dd = dd;
  node->dr = dr;
  if (node->heapi >= 0) {
    if (diff <= 0) od_mv_dec_heap_up(est, node->heapi);
    else od_mv_dec_heap_down(est, node->heapi);
  }
}

static void od_mv_est_init_nodes(od_mv_est_ctx *est) {
  od_state *state;
  od_mv_node *mv_row;
  od_mv_grid_pt *grid;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  for (vy = 0; vy <= nvmvbs; vy++) {
    mv_row = est->mvs[vy];
    grid = state->mv_grid[vy];
    for (vx = 0; vx <= nhmvbs; vx++) {
      int level;
      level = OD_MC_LEVEL[vy & 3][vx & 3];
      if (level <= est->level_max) {
        int flag_rate;
        /*While we're here, reset the MV state.*/
        grid[vx].valid = 1;
        est->row_counts[vy]++;
        est->col_counts[vx]++;
        /*Inbetween the level limits, vertices require on average 2 bits to
           indicate the presence of children.
          TODO: Fix a more exact representation.
          TODO: Fix-up child flags for blocks outside the frame border.*/
        flag_rate = (est->level_min <= level && level < est->level_max) <<
         (1 + OD_BITRES);
        mv_row[vx].dr = -mv_row[vx].mv_rate - flag_rate;
      }
      else grid[vx].valid = 0;
    }
  }
}

/*Computes the SAD of all blocks at all scales with all possible edge
   splittings, using OBMC.
  These are what will drive the error of the adaptive subdivision process.*/
static void od_mv_est_calc_sads(od_mv_est_ctx *est, int ref) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  int oc;
  int s;
  state = &est->enc->state;
  /*TODO: Interleaved evaluation would probably provide better cache
     coherency.*/
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  if (est->level_max >= 3) {
    for (vy = 0; vy < nvmvbs; vy++) {
      od_mv_node *mv_row;
      mv_row = est->mvs[vy];
      for (vx = 0; vx < nhmvbs; vx++) {
        oc = (vx & 1) ^ ((vy & 1) << 1 | (vy & 1));
        /*While we're here, fill in the block's setup state.*/
        mv_row[vx].oc = oc;
        mv_row[vx].log_mvb_sz = 0;
        if (est->level_max >= 4) {
          for (s = 0; s < 4; s++) {
            est->sad_cache[0][vy][vx][s] =
             (ogg_uint16_t)od_mv_est_sad8(est, ref, vx, vy, oc, s, 0);
          }
          mv_row[vx].s = 3;
          mv_row[vx].sad = est->sad_cache[0][vy][vx][3];
        }
        else {
          mv_row[vx].s = 0;
          mv_row[vx].sad = od_mv_est_sad8(est, ref, vx, vy, oc, 0, 0);
        }
      }
    }
  }
  nhmvbs >>= 1;
  nvmvbs >>= 1;
  if (est->level_max >= 1) {
    if (est->level_min < 3) {
      for (vy = 0; vy < nvmvbs; vy++) {
        od_mv_node *mv_row;
        mv_row = est->mvs[vy << 1];
        for (vx = 0; vx < nhmvbs; vx++) {
          oc = (vx & 1) ^ ((vy & 1) << 1 | (vy & 1));
          if (est->level_max >= 2) {
            for (s = 0; s < 4; s++) {
              est->sad_cache[1][vy][vx][s] =
               (ogg_uint16_t)od_mv_est_sad8(est,
               ref, vx << 1, vy << 1, oc, s, 1);
            }
            if (est->level_max <= 2) {
              mv_row[vx << 1].oc = oc;
              mv_row[vx << 1].s = 3;
              mv_row[vx << 1].log_mvb_sz = 1;
              mv_row[vx << 1].sad = est->sad_cache[1][vy][vx][3];
            }
          }
          else {
            mv_row[vx << 1].oc = oc;
            mv_row[vx << 1].s = 0;
            mv_row[vx << 1].log_mvb_sz = 1;
            mv_row[vx << 1].sad = od_mv_est_sad8(est,
             ref, vx << 1, vy << 1, oc, 0, 1);
          }
        }
      }
    }
  }
  else {
    nhmvbs >>= 1;
    nvmvbs >>= 1;
    for (vy = 0; vy < nvmvbs; vy++) {
      od_mv_node *mv_row;
      mv_row = est->mvs[vy << 2];
      for (vx = 0; vx < nhmvbs; vx++) {
        mv_row[vx << 2].oc = 0;
        mv_row[vx << 2].s = 3;
        mv_row[vx << 2].log_mvb_sz = 2;
        mv_row[vx << 2].sad = od_mv_est_sad8(est,
         ref, vx << 2, vy << 2, 0, 3, 2);
      }
    }
  }
}

static void od_mv_est_init_du(od_mv_est_ctx *est, int ref, int vx, int vy) {
  od_state *state;
  od_mv_node *dec;
  od_mv_node *merge;
  const od_mv_err_node *errdom;
  int nerrdom;
  const od_offset *mergedom;
  int nhmvbs;
  int nvmvbs;
  int level;
  int dlev;
  int log_mvb_sz_min;
  int log_mvb_sz;
  int di;
  int dvx;
  int dvy;
  int dx;
  int dy;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Computing du's for (%i, %i)", vx, vy));
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  dec = est->mvs[vy] + vx;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  dlev = est->level_max <= 2;
  log_mvb_sz_min = (5 - est->level_max) >> 1;
  errdom = OD_ERRDOM[level - 1 + (dlev << 1)];
  nerrdom = OD_NERRDOM[level - 1 + (dlev << 1)];
  mergedom = OD_MERGEDOM[level - 1 + (dlev << 1)];
  dec->dd = 0;
  /*Subtract off the error before decimation.*/
  for (di = 0; di < nerrdom; di++) {
    dvx = vx + (errdom[di].dx << dlev);
    dvy = vy + (errdom[di].dy << dlev);
    if (dvx >= 0 && dvy >= 0 && dvx < nhmvbs && dvy < nvmvbs) {
      int mvb_sz;
      log_mvb_sz = errdom[di].log_mvb_sz + dlev;
      if (log_mvb_sz < log_mvb_sz_min) continue;
      mvb_sz = 1 << (log_mvb_sz - dlev);
      for (dy = 0; dy < mvb_sz; dy++) {
        for (dx = 0; dx < mvb_sz; dx++) {
          dec->dd -= est->mvs[dvy + (dy << dlev)][dvx + (dx << dlev)].sad;
          /*dec->dd -= est->sad_cache[dlev][(dvy >> dlev) + dy]
           [(dvx >> dlev) + dx][undecs];*/
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Added error (%i, %i) [%ix%i]: %i", dvx + (dx << dlev),
           dvy + (dy << dlev), 4 << dlev, 4 << dlev, dec->dd));
        }
      }
    }
    else {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "(%i, %i) outside [%i, %i]x[%i, %i]",
       dvx, dvy, 0, 0, nhmvbs, nvmvbs));
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Subtracted initial error: %i", dec->dd));
  /*Decimate the vertices in the merging domain.
    Also sum up the rate changes while we do it.*/
  for (di = 0;; di++) {
    dvx = vx + (mergedom[di][0] << dlev);
    if (dvx < 0 || dvx > nhmvbs) continue;
    dvy = vy + (mergedom[di][1] << dlev);
    if (dvy < 0 || dvy > nvmvbs) continue;
    if (OD_MC_LEVEL[dvy & 3][dvx & 3] > est->level_max) continue;
    state->mv_grid[dvy][dvx].valid = 0;
    merge = est->mvs[dvy] + dvx;
    if (merge == dec) break;
    dec->dr += merge->dr;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Merged vertex (%2i, %2i), dr: %i", dvx, dvy, dec->dr));
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Merged vertex (%2i, %2i)", dvx, dvy));
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Decimated vertices in merging domain."));
  /*Add in the error after decimation.*/
  for (di = 0; di < nerrdom; di++) {
    dvx = vx + (errdom[di].dx << dlev);
    dvy = vy + (errdom[di].dy << dlev);
    if (dvx >= 0 && dvy >= 0 && dvx < nhmvbs && dvy < nvmvbs) {
      log_mvb_sz = errdom[di].log_mvb_sz + dlev;
      if (log_mvb_sz < log_mvb_sz_min) continue;
      else if (log_mvb_sz < 2) {
        int mask;
        int s1vx;
        int s1vy;
        int s3vx;
        int s3vy;
        int oc;
        int s;
        mask = (1 << (log_mvb_sz + 1)) - 1;
        oc = !!(dvx & mask);
        if (dvy & mask) oc = 3 - oc;
        s1vx = dvx + (OD_VERT_DX[(oc + 1) & 3] << log_mvb_sz);
        s1vy = dvy + (OD_VERT_DY[(oc + 1) & 3] << log_mvb_sz);
        s3vx = dvx + (OD_VERT_DX[(oc + 3) & 3] << log_mvb_sz);
        s3vy = dvy + (OD_VERT_DY[(oc + 3) & 3] << log_mvb_sz);
        s = state->mv_grid[s1vy][s1vx].valid |
         state->mv_grid[s3vy][s3vx].valid << 1;
        dec->dd +=
         est->sad_cache[log_mvb_sz][dvy >> log_mvb_sz][dvx >> log_mvb_sz][s];
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Added error (%i, %i) [%ix%i] {%i, %i}: %i", dvx, dvy,
         1 << (log_mvb_sz + 2), 1 << (log_mvb_sz + 2), oc, s, dec->dd));
      }
      else {
        /*Cache the SAD for top-level blocks in the dd field, which is
           otherwise unused (since they cannot be decimated).*/
        est->mvs[dvy][dvx].dd = od_mv_est_sad8(est, ref, dvx, dvy, 0, 3, 2);
        dec->dd += est->mvs[dvy][dvx].dd;
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Added error (%i, %i) [%ix%i]: %i", dvx, dvy,
         1 << (log_mvb_sz + 2), 1 << (log_mvb_sz + 2), dec->dd));
      }
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Total merging error: %i", dec->dd));
  /*Restore the vertices in the merging domain.*/
  for (di = 0;; di++) {
    dvx = vx + (mergedom[di][0] << dlev);
    if (dvx < 0 || dvx > nhmvbs) continue;
    dvy = vy + (mergedom[di][1] << dlev);
    if (dvy < 0 || dvy > nvmvbs) continue;
    if (OD_MC_LEVEL[dvy & 3][dvx & 3] > est->level_max) continue;
    state->mv_grid[dvy][dvx].valid = 1;
    if (dvx == vx && dvy == vy) break;
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Restored vertices in merging domain."));
  /*Add this node to the heap.*/
  dec->heapi = est->dec_nheap;
  est->dec_heap[est->dec_nheap++] = dec;
}

static void od_mv_est_init_dus(od_mv_est_ctx *est, int ref) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  OD_CLEAR(est->row_counts, nvmvbs + 1);
  OD_CLEAR(est->col_counts, nhmvbs + 1);
  od_mv_est_init_nodes(est);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Finished MV bits."));
  od_mv_est_calc_sads(est, ref);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Finished SADs."));
  /*Clear the merge heap.*/
  est->dec_nheap = 0;
  est->dec_heap[0] = NULL;
  /*The initialization is destructive to dr, and so must proceed by level from
     top to bottom.*/
  if (est->level_max >= 1) {
    /*Level 1 vertices.*/
    if (est->level_min < 1) {
      for (vy = 2; vy <= nvmvbs; vy += 4) {
        for (vx = 2; vx <= nhmvbs; vx += 4) {
          od_mv_est_init_du(est, ref, vx, vy);
        }
      }
    }
    if (est->level_max >= 2) {
      /*Level 2 vertices.*/
      if (est->level_min < 2) {
        for (vy = 0;; vy += 2) {
          for (vx = 2; vx <= nhmvbs; vx += 4) {
            od_mv_est_init_du(est, ref, vx, vy);
          }
          vy += 2;
          if (vy > nvmvbs) break;
          for (vx = 0; vx <= nhmvbs; vx += 4) {
            od_mv_est_init_du(est, ref, vx, vy);
          }
        }
      }
      if (est->level_max >= 3) {
        if (est->level_min < 3) {
          /*Level 3 vertices.*/
          for (vy = 1; vy <= nvmvbs; vy += 2) {
            for (vx = 1; vx <= nhmvbs; vx += 2) {
              od_mv_est_init_du(est, ref, vx, vy);
            }
          }
        }
        if (est->level_max >= 4) {
          /*Level 4 vertices.*/
          if (est->level_min < 4) {
            for (vy = 0;; vy++) {
              for (vx = 1; vx <= nhmvbs; vx += 2) {
                od_mv_est_init_du(est, ref, vx, vy);
              }
              vy++;
              if (vy > nvmvbs) break;
              for (vx = 0; vx <= nhmvbs; vx += 2) {
                od_mv_est_init_du(est, ref, vx, vy);
              }
            }
          }
        }
      }
    }
  }
  /*Make the node list into a proper heap.*/
  od_mv_dec_heapify(est);
}

static void od_mv_est_decimate(od_mv_est_ctx *est, int ref) {
  od_mv_node *dec;
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int dlev;
  int vx;
  int vy;
  od_mv_est_init_dus(est, ref);
#if defined(OD_LOGGING_ENABLED)
  if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN)) {
    od_mv_est_check_rd_state(est, ref, 2);
  }
#endif
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  dlev = est->level_max <= 2;
  for (;;) {
    const od_offset *mergedom;
    int level;
    int di;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
    if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
      char iter_label[16];
      od_state_mc_predict(state, ref);
      od_state_fill_vis(state);
      sprintf(iter_label, "ani%08i", state->ani_iter++);
      od_state_dump_img(state, &state->vis_img, iter_label);
    }
#endif
    dec = od_mv_dec_heap_delhead(est);
    /*Stop if we've fully decimated the mesh, or if this decimation would not
       improve R-D performance at the current lambda.*/
    if (dec == NULL
     || dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE) > 0) {
      break;
    }
    level = OD_MC_LEVEL[dec->vy & 3][dec->vx & 3];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Merging node (%2i, %2i), level %i, dd %5i, dr %5i, dopt %5i:",
     dec->vx, dec->vy, level, dec->dd, dec->dr,
     dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE)));
    mergedom = OD_MERGEDOM[level - 1 + (dlev << 1)];
    for (di = 0;; di++) {
      od_mv_node *merge;
      od_mv_node *ancestor;
      od_mv_node *block;
      const od_offset *anc;
      int nanc;
      int ai;
      int ax;
      int ay;
      int bx;
      int by;
      int log_mvb_sz;
      int mask;
      /*Don't decimate vertices outside of the mesh.*/
      vx = dec->vx + (mergedom[di][0] << dlev);
      if (vx < 0 || vx > nhmvbs) continue;
      vy = dec->vy + (mergedom[di][1] << dlev);
      if (vy < 0 || vy > nvmvbs) continue;
      merge = est->mvs[vy] + vx;
      /*Don't decimate vertices that have already been decimated.*/
      if (!state->mv_grid[vy][vx].valid) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Skipping node (%i, %i) (already merged).", vx, vy));
        continue;
      }
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Merging node (%2i, %2i), dd %5i, dr %5i:",
       vx, vy, merge->dd, merge->dr));
      /*Update the deltas for this vertex in the merging domain.
        The simple rule applied below handles overlapped domains with an
         inclusion-exclusion approach.
        See Balmelli 2001 for details.*/
      nanc = OD_NANCESTORS[vy & 3][vx & 3];
      anc = OD_ANCESTORS[vy & 3][vx & 3];
      for (ai = 0; ai < nanc; ai++) {
        ax = vx + anc[ai][0];
        if (ax < 0 || ax > nhmvbs) continue;
        ay = vy + anc[ai][1];
        if (ay < 0 || ay > nvmvbs) continue;
        ancestor = est->mvs[ay] + ax;
        od_mv_dec_update(est, ancestor,
         ancestor->dd - merge->dd, ancestor->dr - merge->dr);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Updated ancestor (%2i, %2i) of (%2i, %2i): dd %5i, dr %5i",
         ax, ay, vx, vy, ancestor->dd, ancestor->dr));
      }
      state->mv_grid[vy][vx].valid = 0;
      od_mv_dec_heap_del(est, merge);
      est->row_counts[vy]--;
      est->col_counts[vx]--;
      level = OD_MC_LEVEL[vy & 3][vx & 3];
      log_mvb_sz = (4 - level) >> 1;
      /*Account for quadrilaterals which may have only partially belonged to
         the merging domain (e.g., that would not have belonged were we using
         triangles).*/
      if (!(level & 1)) {
        static const int OD_CDX[4] = { -1, 1, -1, 1 };
        static const int OD_CDY[4] = { -1, -1, 1, 1 };
        int k;
        mask = (1 << (log_mvb_sz + 1)) - 1;
        for (k = 0; k < 4; k++) {
          int cx;
          int cy;
          int ddd;
          int s;
          cx = vx + (OD_CDX[k] << log_mvb_sz);
          if (cx < 0 || cx > nhmvbs) continue;
          cy = vy + (OD_CDY[k] << log_mvb_sz);
          if (cy < 0 || cy > nvmvbs) continue;
          bx = vx + (OD_ERRDOM4[k].dx << log_mvb_sz);
          by = vy + (OD_ERRDOM4[k].dy << log_mvb_sz);
          block = est->mvs[by] + bx;
          by >>= log_mvb_sz;
          bx >>= log_mvb_sz;
          if (!state->mv_grid[cy][cx].valid) {
            block->s = 0;
            block->sad = est->sad_cache[log_mvb_sz][by][bx][0];
            /*If the opposing corner has already been decimated, the remaining
               adjustments have already been made.*/
            continue;
          }
          /*s is the split state of the error block with (vx, vy) decimated,
             and (cx, cy) undecimated.*/
          s = 1 << ((((k + 3) & 3) >> 1) ^ !!(vx & mask));
          block->s = s;
          block->sad = est->sad_cache[log_mvb_sz][by][bx][s];
          /*Replace the old decimation error change with the new one.*/
          ddd = est->sad_cache[log_mvb_sz][by][bx][0]
           - est->sad_cache[log_mvb_sz][by][bx][s ^ 3]
           + est->sad_cache[log_mvb_sz][by][bx][3]
           - est->sad_cache[log_mvb_sz][by][bx][s];
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Checking opposing corner (%2i, %2i): ddd %i", cx, cy, ddd));
          /*This happens in regions of constant motion.*/
          if (ddd == 0) continue;
          ancestor = est->mvs[cy] + cx;
          od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Updated corner (%2i, %2i): dd %5i, dr %5i",
           cx, cy, ancestor->dd, ancestor->dr));
          /*Update the opposing corner's ancestors, which also, of
             necessity, must contain the affected quadrilateral, and must
             not have been decimated yet.*/
          nanc = OD_NANCESTORS[cy & 3][cx & 3];
          anc = OD_ANCESTORS[cy & 3][cx & 3];
          for (ai = 0; ai < nanc; ai++) {
            ax = cx + anc[ai][0];
            if (ax < 0 || ax > nhmvbs) continue;
            ay = cy + anc[ai][1];
            if (ay < 0 || ay > nvmvbs) continue;
            ancestor = est->mvs[ay] + ax;
            od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
            OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
             "Updated ancestor (%2i, %2i): dd %5i, dr %5i",
             ax, ay, ancestor->dd, ancestor->dr));
          }
          /*Add back in the components that do not apply to the interior
             corner.*/
          ddd = -ddd;
          if (vx & mask) cx = vx;
          else cy = vy;
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Checking interior corner (%2i, %2i): ddd %i", cx, cy, ddd));
          ancestor = est->mvs[cy] + cx;
          od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Updated corner (%2i, %2i): dd %5i, dr %5i",
           cx, cy, ancestor->dd, ancestor->dr));
          /*And update all the interior corner's ancestors, which also, of
             necessity, must contain the affected quadrilateral, and must not
             have been decimated yet.*/
          nanc = OD_NANCESTORS[cy & 3][cx & 3];
          anc = OD_ANCESTORS[cy & 3][cx & 3];
          for (ai = 0; ai < nanc; ai++) {
            ax = cx + anc[ai][0];
            if (ax < 0 || ax > nhmvbs) continue;
            ay = cy + anc[ai][1];
            if (ay < 0 || ay > nvmvbs) continue;
            ancestor = est->mvs[ay] + ax;
            od_mv_dec_update(est, ancestor, ancestor->dd + ddd, ancestor->dr);
            OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
             "Updated ancestor (%2i, %2i): dd %5i, dr %5i",
             ax, ay, ancestor->dd, ancestor->dr));
          }
        }
      }
      /*Otherwise, we eliminated several smaller blocks.
        Update the SAD and block setup for the larger block that took their
         place.*/
      else {
        int oc;
        bx = vx - (1 << log_mvb_sz);
        by = vy - (1 << log_mvb_sz);
        log_mvb_sz++;
        mask = (1 << (log_mvb_sz + 1)) - 1;
        oc = !!(bx & mask);
        if (by & mask) oc = 3 - oc;
        block = est->mvs[by] + bx;
        block->log_mvb_sz = log_mvb_sz;
        block->oc = oc;
        block->s = 3;
        if (log_mvb_sz < 2) {
          block->sad =
           est->sad_cache[log_mvb_sz][by >> log_mvb_sz][bx >> log_mvb_sz][3];
        }
        /*At the top level, we cached the SAD in the dd field.*/
        else block->sad = block->dd;
      }
      /*If we just decimated our target vertex, stop.*/
      if (merge == dec) break;
    }
  }
#if defined(OD_LOGGING_ENABLED)
  if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN)) {
    od_mv_est_check_rd_state(est, ref, 2);
  }
#endif
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Finished merging."));
  /*if (dec != NULL) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Node (%i, %i) dd %i, dr %i, dopt %i: not enough.", dec->vx, dec->vy,
     dec->dd, dec->dr, dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE)));
  }*/
  /*if (state->mv_grid[31][1].valid) {
    dec = est->mvs[31] + 1;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "(%i, %i) remains. dd: %5i, dr: %2i, dopt: %6i.", dec->vx, dec->vy,
     dec->dd, dec->dr, dec->dr*est->lambda + (dec->dd << OD_ERROR_SCALE)));
  }*/
}

/*STAGE 3: Iterated Dynamic Programming.*/

/*The list of MVs that can be predicted by a level 0 MV, excluding those not
   yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED0[17] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  2, -2 }, {  1, -1 }, {  2,  2 }, {  1,  1 }, {  0,  4 }, {  4,  4 },
  /*The remaining ones are not.*/
  { -2, -2 }, {  0, -2 }, { -1, -1 }, {  0, -1 }, { -1,  0 }, { -2,  0 },
  { -1,  1 }, {  0,  1 }, { -2,  2 }, {  0,  2 }, { -4,  4 }
};
/*The list of MVs that can be predicted by a level 1 MV, excluding those
   not yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED1[10] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  1, -1 }, {  1,  1 },
  /*The remaining ones are not.*/
  {  0, -2 }, { -1, -1 }, {  0, -1 }, { -2,  0 }, { -1,  0 }, { -1,  1 },
  {  0,  1 }, {  0,  2 }
};
/*The list of MVs that can be predicted by a level 2 MV, excluding those
   not yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED2[7] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  1, -1 }, {  1,  1 },
  /*The remaining ones are not.*/
  { -1, -1 }, {  0, -1 }, { -1,  0 }, { -1,  1 }, {  0,  1 }
};
/*The list of MVs that can be predicted by a level 3 MV, excluding those
   not yet considered by DP across rows.*/
static const od_offset OD_ROW_PREDICTED3[3] = {
  /*These predicted MVs are NOT changeable by future MVs in the DP path.*/
  {  0, -1 }, { -1,  0 }, {  0,  1 }
};

/*The list of MVs that can be predicted by a level 0 MV, excluding those not
   yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED0[17] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  {  2,  2 }, { -2,  2 }, { -1,  1 }, {  1,  1 }, {  4,  4 },
  /*The remaining ones are not.*/
  { -2, -2 }, {  0, -2 }, {  2, -2 }, { -1, -1 }, {  0, -1 }, {  1, -1 },
  { -2,  0 }, { -1,  0 }, {  1,  0 }, {  2,  0 }, {  4,  0 }, { -4,  4 }
};
/*The list of MVs that can be predicted by a level 1 MV, excluding those
   not yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED1[10] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 },
  /*The remaining ones are not.*/
  {  0, -2 }, { -1, -1 }, {  0, -1 }, {  1, -1 }, { -2,  0 }, { -1,  0 },
  {  1,  0 }, {  2,  0 }
};
/*The list of MVs that can be predicted by a level 2 MV, excluding those
   not yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED2[7] = {
  /*These predicted MVs are changeable by future MVs in the DP path.*/
  { -1,  1 }, {  1,  1 },
  /*The remaining ones are not.*/
  { -1, -1 }, {  0, -1 }, {  1, -1 }, { -1,  0 }, {  1,  0 }
};
/*The list of MVs that can be predicted by a level 3 MV, excluding those
   not yet considered by DP across columns.*/
static const od_offset OD_COL_PREDICTED3[3] = {
  /*These predicted MVs are NOT changeable by future MVs in the DP path.*/
  {  0, -1 }, { -1,  0 }, {  1,  0 }
};

/*The number of predicted MVs in each list.*/
static const int OD_NPREDICTED[5] = { 17, 10, 7, 3, 0 };
/*The number of changeable predicted MVs in each list.*/
static const int OD_NROW_PRED_CHANGEABLE[4] = { 6, 2, 2, 0 };
/*The number of changeable predicted MVs in each list.*/
static const int OD_NCOL_PRED_CHANGEABLE[4] = { 5, 2, 2, 0 };
/*The lists of offsets to predicted MVs for each level.*/
static const od_offset *const OD_ROW_PREDICTED[4] = {
  OD_ROW_PREDICTED0,
  OD_ROW_PREDICTED1,
  OD_ROW_PREDICTED2,
  OD_ROW_PREDICTED3
};
/*The lists of offsets to predicted MVs for each level.*/
static const od_offset *const OD_COL_PREDICTED[4] = {
  OD_COL_PREDICTED0,
  OD_COL_PREDICTED1,
  OD_COL_PREDICTED2,
  OD_COL_PREDICTED3
};

/*This should be the maximum value found in either OD_ROW_PRED_HIST_SIZE or
   OD_COL_PRED_HIST_SIZE.*/
#define OD_PRED_HIST_SIZE_MAX (8)

/*The amount of history to restore in the trellis state to ensure predicted MVs
   are evaluated correctly in row refinement.*/
static const int OD_ROW_PRED_HIST_SIZE[5] = { 8, 4, 2, 2, 1 };
/*The amount of history to restore in the trellis state to ensure predicted MVs
   are evaluated correctly in column refinement.*/
static const int OD_COL_PRED_HIST_SIZE[5] = { 8, 4, 2, 2, 1 };

/*Returns the boundary case indicating which motion vector range edges the
   current motion vector is abutting.
  vx: The horizontal position of the node.
  vy: The vertical position of the node.
  dx: The horizontal component of the motion vector.
  dy: The vertical component of the motion vector.
  dsz: The amount the vector is being adjusted by.
  log_blk_sz: The log base 2 of the maximum size of a block the vector can
               belong to.
  Return: A set of flags indicating the boundary conditions, after the
   documentation at OD_SQUARE_SITES.*/
static int od_mv_est_get_boundary_case(od_state *state,
 int vx, int vy, int dx, int dy, int dsz, int log_blk_sz) {
  int mvxmin;
  int mvxmax;
  int mvymin;
  int mvymax;
  int blk_sz;
  int bx;
  int by;
  blk_sz = 1 << log_blk_sz;
  bx = (vx - 2) << 2;
  by = (vy - 2) << 2;
  mvxmin = (OD_MAXI(bx - blk_sz - 32, -16) - (bx - blk_sz)) << 3;
  mvxmax = ((OD_MINI(bx + blk_sz + 32, state->frame_width + 16)
   - (bx + blk_sz)) << 3) - dsz;
  mvymin = (OD_MAXI(by - blk_sz - 32, -16) - (by - blk_sz)) << 3;
  mvymax = ((OD_MINI(by + blk_sz + 32, state->frame_height + 16)
   - (by + blk_sz)) << 3) - dsz;
  return (dx <= mvxmin) | (dx >= mvxmax) << 1 |
   (dy <= mvymin) << 2 | (dy >= mvymax) << 3;
}

/*Computes the SAD of the specified block.*/
static ogg_int32_t od_mv_est_block_sad8(od_mv_est_ctx *est, int ref,
 od_mv_node *block) {
  ogg_int32_t ret;
  ret = od_mv_est_sad8(est, ref, block->vx, block->vy,
   block->oc, block->s, block->log_mvb_sz);
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Adding SAD (%3i, %3i) [%2ix%2i]: %6i", (block->vx - 2) << 2,
   (block->vy - 2) << 2, 4 << block->log_mvb_sz, 4 << block->log_mvb_sz, ret));
  return ret;
}

/*Gets the change in SAD for the blocks affected by the given DP node, using
   the current state of the grid.*/
static ogg_int32_t od_mv_dp_get_sad_change8(od_mv_est_ctx *est, int ref,
 od_mv_dp_node *dp, ogg_int32_t block_sads[OD_DP_NBLOCKS_MAX]) {
  int bi;
  ogg_int32_t dd;
  dd = 0;
  for (bi = 0; bi < dp->nblocks; bi++) {
    od_mv_node *block;
    block = dp->blocks[bi];
    block_sads[bi] = od_mv_est_block_sad8(est, ref, block);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "SAD change for block (%i, %i) [%ix%i]: %i - %i = %i", block->vx,
     block->vy, 1 << (block->log_mvb_sz + 2), 1 << (block->log_mvb_sz + 2),
     block_sads[bi], block->sad, block_sads[bi] - block->sad));
    dd += block_sads[bi] - block->sad;
  }
  return dd;
}

/*Computes a rate adjustment for the predictors changed by following the given
   trellis path.
  As a side effect, enough of the trellis needed to evaluate that change is
   loaded into the MV grid.
  dp: The current DP node.
  cur_mv_rate: Returns the MV rate for the motion vector set by the current DP
                node.
  pred_mv_rates: Returns the MV rate for the motion vectors predicted by the
                  MV set by the current DP node.
  prevsi: The state index to follow in the previous DP node.
  mv_res: The motion vector resolution (0 = 1/8th pel to 2 = 1/2 pel).
  Return: The change in rate for the preceding MVs.*/
static int od_mv_dp_get_rate_change(od_state *state, od_mv_dp_node *dp,
 int *cur_mv_rate, int pred_mv_rates[OD_DP_NPREDICTED_MAX],
 int prevsi, int mv_res) {
  od_mv_node *mv;
  od_mv_grid_pt *mvg;
  int nhmvbs;
  int nvmvbs;
  int pred[2];
  int pi;
  int dr;
  /*Move the state from the current trellis path into the grid.*/
  if (dp->min_predictor_node != NULL) {
    int pred_sis[OD_PRED_HIST_SIZE_MAX];
    int pred_si;
    int npreds;
    od_mv_dp_node *pred_dp;
    npreds = dp - dp->min_predictor_node;
    OD_ASSERT2(npreds <= OD_PRED_HIST_SIZE_MAX, "Too far back!");
    OD_LOG_PARTIAL((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "Restoring "));
    /*First, follow the trellis path backwards to find the state used in each
       node.*/
    pred_si = pred_sis[npreds - 1] = prevsi;
    for (pi = 2; pi <= npreds; pi++) {
      pred_dp = dp - pi;
      pred_si = pred_sis[npreds - pi] = pred_dp[1].states[pred_si].prevsi;
    }
    /*Then restore that state going FORWARDS.*/
    for (pred_dp = dp->min_predictor_node; pred_dp < dp; pred_dp++) {
      pred_si = pred_sis[pred_dp - dp->min_predictor_node];
      /*Restore the state for this MV itself.*/
      pred_dp->mv->mv_rate = pred_dp->states[pred_si].mv_rate;
      mvg = pred_dp->mvg;
      mvg->mv[0] = pred_dp->states[pred_si].mv[0];
      mvg->mv[1] = pred_dp->states[pred_si].mv[1];
      OD_LOG_PARTIAL((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "(%i, %i: %i)->(%i , %i) ",
       pred_dp->mv->vx, pred_dp->mv->vy, pred_si, mvg->mv[0], mvg->mv[1]));
      /*Restore the state for the MVs this one predicted.*/
      for (pi = 0; pi < pred_dp->npred_changeable; pi++) {
        pred_dp->predicted_mvs[pi]->mv_rate =
         pred_dp->states[pred_si].pred_mv_rates[pi];
      }
    }
    OD_LOG_PARTIAL((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "\n"));
  }
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  /*Compute the new rate for the current MV.*/
  mv = dp->mv;
  if (mv->vx < 2 || mv->vx > nhmvbs - 2 || mv->vy < 2 || mv->vy > nvmvbs - 2) {
    *cur_mv_rate = dr = 0;
  }
  else {
    od_state_get_predictor(state, pred, mv->vx, mv->vy,
     OD_MC_LEVEL[mv->vy & 3][mv->vx & 3], mv_res);
    mvg = dp->mvg;
    *cur_mv_rate = od_mv_est_bits(mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res,
     pred[0], pred[1]);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Current MV rate: %i - %i = %i",
     *cur_mv_rate, mv->mv_rate, *cur_mv_rate - mv->mv_rate));
    dr = *cur_mv_rate - mv->mv_rate;
    /*Compute the new rates for the MVs this one predicts.*/
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Calculating predicted pred_mv_rates for node (%i, %i):",
     dp->mv->vx, dp->mv->vy));
    for (pi = 0; pi < dp->npredicted; pi++) {
      mv = dp->predicted_mvs[pi];
      mvg = dp->predicted_mvgs[pi];
      od_state_get_predictor(state, pred, mv->vx, mv->vy,
       OD_MC_LEVEL[mv->vy & 3][mv->vx & 3], mv_res);
      pred_mv_rates[pi] = od_mv_est_bits(
       mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res, pred[0], pred[1]);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Calculated predicted mv_rate of %i for (%i, %i)",
       pred_mv_rates[pi], mv->vx, mv->vy));
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Predictor was: (%i, %i)   MV was: (%i, %i)",
       pred[0], pred[1], mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res));
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Predicted MV (%i, %i) rate: %i - %i = %i", mv->vx, mv->vy,
       pred_mv_rates[pi], mv->mv_rate, pred_mv_rates[pi] - mv->mv_rate));
      dr += pred_mv_rates[pi] - mv->mv_rate;
    }
  }
  return dr;
}

#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
static const unsigned char OD_YCbCr_EDGE[3] = { 41, 240, 110 };

static void od_mv_dp_animate_state(od_state *state,
 int ref, od_mv_dp_node *dp, int has_gap) {
  od_mv_dp_node *dp0;
  char iter_label[16];
  int active_states[OD_DP_NSTATES_MAX];
  int prev_active_states[OD_DP_NSTATES_MAX];
  int nactive_states;
  int nprev_active_states;
  int dp_state;
  int si;
  int d0vx;
  int d0vy;
  int d1vx;
  int d1vy;
  int x0;
  int y0;
  od_state_mc_predict(state, ref);
  od_state_fill_vis(state);
  /*Now, draw the current state of the DP.*/
  dp0 = dp;
  /*First draw the candidate edge labels for the active trellis paths.*/
  for (si = 0; si < dp0->nstates; si++) prev_active_states[si] = si;
  nprev_active_states = dp0->nstates;
  nactive_states = 0;
  do {
    if (nactive_states > 0) {
      d0vx = dp[0].mv->vx;
      d0vy = dp[0].mv->vy;
      x0 = ((d0vx - 2) << 3) + (OD_UMV_PADDING << 1);
      y0 = ((d0vy - 2) << 3) + (OD_UMV_PADDING << 1);
      if (nprev_active_states > 0) {
        int mvb_sz;
        int x1;
        int y1;
        d1vx = dp[1].mv->vx;
        d1vy = dp[1].mv->vy;
        x1 = ((d1vx - 2) << 3) + (OD_UMV_PADDING << 1);
        y1 = ((d1vy - 2) << 3) + (OD_UMV_PADDING << 1);
        od_img_draw_line(&state->vis_img, x0, y0, x1, y1, OD_YCbCr_EDGE);
        if (d1vx - d0vx > 1) {
          mvb_sz = d1vx - d0vx;
          if (!has_gap || dp + 1 != dp0) mvb_sz >>= 1;
          if (!state->mv_grid[d0vy][d0vx + mvb_sz].valid) {
            if (d0vy >= mvb_sz
             && state->mv_grid[d0vy - mvb_sz][d0vx + mvb_sz].valid) {
              od_img_draw_line(&state->vis_img,
               x0 + (mvb_sz << 3), y0 - (mvb_sz << 3), x0 + (mvb_sz << 3), y1,
               OD_YCbCr_EDGE);
            }
            if (dp[0].mv->vy <= ((state->nvmbs + 1) << 2) - mvb_sz
             && state->mv_grid[d0vy + mvb_sz][d0vx + mvb_sz].valid) {
              od_img_draw_line(&state->vis_img,
               x0 + (mvb_sz << 3), y0 + (mvb_sz << 3),
               x0 + (mvb_sz << 3), y1, OD_YCbCr_EDGE);
            }
          }
        }
        else if (d1vy - d0vy > 1) {
          mvb_sz = d1vy - d0vy;
          if (!has_gap || dp + 1 != dp0) mvb_sz >>= 1;
          if (!state->mv_grid[d0vy + mvb_sz][d0vx].valid) {
            if (d0vx >= mvb_sz
             && state->mv_grid[d0vy + mvb_sz][d0vx - mvb_sz].valid) {
              od_img_draw_line(&state->vis_img,
               x0 - (mvb_sz << 3), y0 + (mvb_sz << 3), x1, y0 + (mvb_sz << 3),
               OD_YCbCr_EDGE);
            }
            if (d0vx <= ((state->nhmbs + 1) << 2) - mvb_sz
             && state->mv_grid[d0vy + mvb_sz][d0vx + mvb_sz].valid) {
              od_img_draw_line(&state->vis_img,
               x0 + (mvb_sz << 3), y0 + (mvb_sz << 3), x1, y0 + (mvb_sz << 3),
               OD_YCbCr_EDGE);
            }
          }
        }
      }
    }
    OD_COPY(active_states, prev_active_states, nprev_active_states);
    nactive_states = nprev_active_states;
    /*Follow the chain backwards to find the new active states.*/
    nprev_active_states = 0;
    for (si = 0; si < nactive_states; si++) {
      int sj;
      dp_state = dp[0].states[active_states[si]].prevsi;
      for (sj = 0;
       sj < nprev_active_states && prev_active_states[sj] != dp_state; sj++);
      if (sj >= nprev_active_states) {
        prev_active_states[nprev_active_states++] = dp_state;
      }
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
  /*Now, draw all the candidate MVs in the active trellis paths.
    These two steps used to be together; now they're apart.
    Sorry for the mess that caused.*/
  /*Redraw the MVs, so they appear over the edge labels above.*/
  od_state_draw_mvs(state);
  for (si = 0; si < dp0->nstates; si++) prev_active_states[si] = si;
  nprev_active_states = dp0->nstates;
  nactive_states = 0;
  dp = dp0;
  do {
    if (nactive_states > 0) {
      d0vx = dp[0].mv->vx;
      d0vy = dp[0].mv->vy;
      x0 = ((d0vx - 2) << 3) + (OD_UMV_PADDING << 1);
      y0 = ((d0vy - 2) << 3) + (OD_UMV_PADDING << 1);
      if (!has_gap || dp + 1 != dp0) {
        x0 = ((d1vx - 2) << 3) + (OD_UMV_PADDING << 1);
        y0 = ((d1vy - 2) << 3) + (OD_UMV_PADDING << 1);
        for (si = 0; si < nactive_states; si++) {
          dp_state = active_states[si];
          od_img_draw_line(&state->vis_img, x0, y0,
           x0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[0], 2, 2),
           y0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[1], 2, 2),
           OD_YCbCr_MVCAND);
        }
      }
    }
    OD_COPY(active_states, prev_active_states, nprev_active_states);
    nactive_states = nprev_active_states;
    /*Follow the chain backwards to find the new active states.*/
    nprev_active_states = 0;
    for (si = 0; si < nactive_states; si++) {
      int sj;
      dp_state = dp[0].states[active_states[si]].prevsi;
      for (sj = 0;
       sj < nprev_active_states && prev_active_states[sj] != dp_state; sj++);
      if (sj >= nprev_active_states) {
        prev_active_states[nprev_active_states++] = dp_state;
      }
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
  /*Draw the first state's MV's.*/
  d1vx = dp[1].mv->vx;
  d1vy = dp[1].mv->vy;
  x0 = ((d1vx - 2) << 3) + (OD_UMV_PADDING << 1);
  y0 = ((d1vy - 2) << 3) + (OD_UMV_PADDING << 1);
  for (si = 0; si < nactive_states; si++) {
    dp_state = active_states[si];
    od_img_draw_line(&state->vis_img, x0, y0,
     x0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[0], 2, 2),
     y0 + OD_DIV_ROUND_POW2(dp[1].states[dp_state].mv[1], 2, 2),
     OD_YCbCr_MVCAND);
  }
  sprintf(iter_label, "ani%08i", state->ani_iter++);
  od_state_dump_img(state, &state->vis_img, iter_label);
}
#endif

/*Row refinement.*/

static void od_mv_dp_row_init(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy, od_mv_dp_node *prev_dp) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  state = &est->enc->state;
  dp->mv = est->mvs[vy] + vx;
  dp->mvg = state->mv_grid[vy] + vx;
  dp->original_mv[0] = dp->mvg->mv[0];
  dp->original_mv[1] = dp->mvg->mv[1];
  dp->original_mv_rate = dp->mv->mv_rate;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  if (vx < 2 || vx > nhmvbs - 2 || vy < 2 || vy > nvmvbs - 2) {
    /*Strictly speaking, we may be used to predict others, but since our MV
       can't possibly change, neither can their rate.*/
    dp->npredicted = dp->npred_changeable = 0;
    /*No one else is used to predict us, or any other MV we predict.
      However, we may still need to load the previous MV into the grid to
       estimate our SADs properly.*/
    dp->min_predictor_node = prev_dp;
  }
  else {
    int level;
    int pred_hist;
    int npred;
    int nchangeable;
    int pi;
    /*Get the list of MVs we help predict.*/
    level = OD_MC_LEVEL[vy & 3][vx & 3];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Initializing node (%i, %i) [%i, %i] at level %i:",
     vx, vy, (vx - 2) << 2, (vy - 2) << 2, level));
    npred = nchangeable = 0;
    for (pi = 0; pi < OD_NPREDICTED[level]; pi++) {
      int px;
      int py;
      px = vx + OD_ROW_PREDICTED[level][pi][0];
      if (px < 2 || px > nhmvbs - 2) continue;
      py = vy + OD_ROW_PREDICTED[level][pi][1];
      if (py < 2 || py > nvmvbs - 2) continue;
      if (state->mv_grid[py][px].valid) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Adding (%i, %i) as a PREDICTED MV.", px, py));
        dp->predicted_mvgs[npred] = state->mv_grid[py] + px;
        dp->predicted_mvs[npred] = est->mvs[py] + px;
        if (pi < OD_NROW_PRED_CHANGEABLE[level]) {
          OD_LOG((OD_LOG_MOTION_ESTIMATION,
           OD_LOG_DEBUG, "It is CHANGEABLE."));
          dp->original_mv_rates[npred] = est->mvs[py][px].mv_rate;
          nchangeable++;
        }
        npred++;
      }
    }
    dp->npredicted = npred;
    dp->npred_changeable = nchangeable;
    /*Now, figure out the earliest DP node that influences our own prediction,
       or that of one of the other MVs we predict.*/
    pred_hist = OD_ROW_PRED_HIST_SIZE[level];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Marking history up to %i back: %i>=%i",
     pred_hist, prev_dp != NULL ? prev_dp->mv->vx : -1, vx-pred_hist));
    if (prev_dp != NULL && prev_dp->mv->vx >= vx - pred_hist) {
      od_mv_dp_node *dp_pred;
      dp_pred = prev_dp;
      while (dp_pred->mv->vx > vx-pred_hist
       && dp_pred->states[0].prevsi >= 0) {
        dp_pred--;
      }
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Stopped at (%i, %i) (%i <= %i? %i) (%i < 0? %i)",
       dp_pred->mv->vx, dp_pred->mv->vy, dp_pred->mv->vx, vx-pred_hist,
       dp_pred->mv->vx <= vx-pred_hist,
       dp_pred->states[0].prevsi, dp_pred->states[0].prevsi < 0));
      if (dp_pred->mv->vx < vx-pred_hist) {
        dp_pred++;
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Too far, incrementing to (%i, %i).",
         dp_pred->mv->vx, dp_pred->mv->vy));
      }
      dp->min_predictor_node = dp_pred;
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "State will be restored back to (%i, %i).",
       dp_pred->mv->vx, dp_pred->mv->vy));
    }
    else dp->min_predictor_node = NULL;
  }
}

static void od_mv_dp_first_row_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nvmvbs = (state->nvmbs + 1) << 2;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  nblocks = 0;
  if (vx > 2) {
    if (level >= 3) {
      if (vy >= mvb_sz) {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      }
      if (vy <= nvmvbs - mvb_sz) {
        dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
      }
    }
    else {
      int half_mvb_sz;
      int mvb_off;
      half_mvb_sz = mvb_sz >> 1;
      if (vy >= mvb_sz) {
        if (state->mv_grid[vy - half_mvb_sz][vx - half_mvb_sz].valid) {
          int ulvx;
          int ulvy;
          ulvx = vx - (half_mvb_sz >> 1);
          ulvy = vy - (half_mvb_sz >> 1);
          if (level > 0 || !state->mv_grid[ulvy][ulvx].valid) {
            mvb_off = half_mvb_sz;
          }
          else mvb_off = half_mvb_sz >> 1;
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - mvb_off;
          if (!state->mv_grid[vy - mvb_off][vx].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx - mvb_off;
          }
          if (!state->mv_grid[vy][vx - mvb_off].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - mvb_off] + vx - (mvb_off << 1);
          }
        }
        else dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      }
      if (vy <= nvmvbs - mvb_sz) {
        if (state->mv_grid[vy + half_mvb_sz][vx - half_mvb_sz].valid) {
          int dlvx;
          int dlvy;
          dlvx = vx - (half_mvb_sz >> 1);
          dlvy = vy + (half_mvb_sz >> 1);
          if (level > 0 || !state->mv_grid[dlvy][dlvx].valid) {
            mvb_off = half_mvb_sz;
          }
          else mvb_off = half_mvb_sz >> 1;
          dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_off;
          if (!state->mv_grid[vy + mvb_off][vx].valid) {
            dp->blocks[nblocks++] = est->mvs[vy + mvb_off] + vx - mvb_off;
          }
          if (!state->mv_grid[vy][vx - mvb_off].valid) {
            dp->blocks[nblocks++] = est->mvs[vy] + vx - (mvb_off << 1);
          }
        }
        else dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
      }
    }
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_prev_row_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nvmvbs;
  int level;
  int prev_level;
  int log_mvb_sz;
  int prev_log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nvmvbs = (state->nvmbs + 1) << 2;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  prev_level = OD_MC_LEVEL[vy & 3][(vx - mvb_sz) & 3];
  prev_log_mvb_sz = (4 - prev_level) >> 1;
  nblocks = 0;
  if (level >= 3) {
    if (vy >= mvb_sz) {
      dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      if (prev_log_mvb_sz > log_mvb_sz
       && !state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid) {
        dp->blocks[nblocks++] = est->mvs[vy - (mvb_sz << 1)] + vx - mvb_sz;
      }
    }
    if (vy <= nvmvbs - mvb_sz) {
      dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
      if (prev_log_mvb_sz > log_mvb_sz
       && !state->mv_grid[vy + mvb_sz][vx - mvb_sz].valid) {
        dp->blocks[nblocks++] = est->mvs[vy + mvb_sz] + vx - mvb_sz;
      }
    }
  }
  else {
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz = mvb_sz >> 1;
    if (vy >= mvb_sz) {
      if (state->mv_grid[vy - half_mvb_sz][vx - half_mvb_sz].valid) {
        int ulvx;
        int ulvy;
        ulvx = vx - (half_mvb_sz >> 1);
        ulvy = vy - (half_mvb_sz >> 1);
        if (level > 0 || !state->mv_grid[ulvy][ulvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - mvb_off;
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx - mvb_off;
        }
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - (mvb_off << 1);
          if (!state->mv_grid[vy - mvb_off][vx - (mvb_off << 1)].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx - (mvb_off << 1);
          }
        }
      }
      else {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_sz << 1)] + vx - mvb_sz;
        }
      }
    }
    if (vy <= nvmvbs - mvb_sz) {
      if (state->mv_grid[vy + half_mvb_sz][vx - half_mvb_sz].valid) {
        int dlvx;
        int dlvy;
        dlvx = vx - (half_mvb_sz >> 1);
        dlvy = vy + (half_mvb_sz >> 1);
        if (level > 0 || !state->mv_grid[dlvy][dlvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_off;
        if (!state->mv_grid[vy + mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy + mvb_off] + vx - mvb_off;
        }
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy] + vx - (mvb_off << 1);
          if (!state->mv_grid[vy + mvb_off][vx - (mvb_off << 1)].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy + mvb_off] + vx - (mvb_off << 1);
          }
        }
      }
      else {
        dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy + mvb_sz][vx - mvb_sz].valid) {
          dp->blocks[nblocks++] = est->mvs[vy + mvb_sz] + vx - mvb_sz;
        }
      }
    }
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_last_row_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nvmvbs = (state->nvmbs + 1) << 2;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  nblocks = 0;
  if (level >= 3) {
    if (vy >= mvb_sz) dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
    if (vy <= nvmvbs - mvb_sz) dp->blocks[nblocks++] = est->mvs[vy] + vx;
  }
  else {
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz = mvb_sz >> 1;
    if (vy >= mvb_sz) {
      if (state->mv_grid[vy - half_mvb_sz][vx + half_mvb_sz].valid) {
        int urvx;
        int urvy;
        urvx = vx + (half_mvb_sz >> 1);
        urvy = vy - (half_mvb_sz >> 1);
        if (level > 0
         || !state->mv_grid[urvy][urvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx;
        if (!state->mv_grid[vy][vx + mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx + mvb_off;
        }
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx;
        }
      }
      else dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
    }
    if (vy <= nvmvbs - mvb_sz) {
      if (state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
        int drvx;
        int drvy;
        drvx = vx + (half_mvb_sz >> 1);
        drvy = vy + (half_mvb_sz >> 1);
        if (level > 0
         || !state->mv_grid[drvy][drvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy] + vx;
        if (!state->mv_grid[vy][vx + mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy] + vx + mvb_off;
        }
        if (!state->mv_grid[vy + mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy + mvb_off] + vx;
        }
      }
      else dp->blocks[nblocks++] = est->mvs[vy] + vx;
    }
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_restore_row_state(od_mv_dp_node *dp) {
  od_mv_grid_pt *mvg;
  int pi;
  do {
    /*Restore the state for this MV itself.*/
    dp->mv->mv_rate = dp->original_mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->original_mv[0];
    mvg->mv[1] = dp->original_mv[1];
    for (pi = 0; pi < dp->npred_changeable; pi++) {
      /*Restore the state for the MVs this one predicted.*/
      dp->predicted_mvs[pi]->mv_rate = dp->original_mv_rates[pi];
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
}

static void od_mv_dp_install_row_state(od_mv_dp_node *dp, int prevsi) {
  od_mv_dp_node *dp0;
  od_mv_grid_pt *mvg;
  int nextsi;
  int si;
  int pi;
  int bi;
  /*We must install the state going FORWARDS, since the pred_mv_rates may have
     changed several times over the course of the trellis.
    Therefore, first we reverse all of the prevsi pointers to make them act
     like nextsi pointers.*/
  nextsi = -1;
  dp0 = dp;
  for (si = prevsi; si >= 0; si = prevsi) {
    dp--;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Node %i, prevsi: %i nextsi: %i", (int)(dp0 - dp), prevsi, nextsi));
    prevsi = dp->states[si].prevsi;
    dp->states[si].prevsi = nextsi;
    nextsi = si;
  }
  /*Now we traverse forward installing the rest of the state.*/
  for (si = nextsi; dp < dp0; dp++) {
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Installing state %i for (%i, %i):", si, dp->mv->vx, dp->mv->vy));
    /*Install the state for this MV itself.*/
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Installing current mv_rate for (%i, %i): %i",
     dp->mv->vx, dp->mv->vy, dp->states[si].mv_rate));
    dp->mv->mv_rate = dp->states[si].mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->states[si].mv[0];
    mvg->mv[1] = dp->states[si].mv[1];
    /*Install the new block SADs.*/
    for (bi = 0; bi < dp->nblocks; bi++) {
      dp->blocks[bi]->sad = dp->states[si].block_sads[bi];
    }
    /*Install the state for the MVs this one predicted.*/
    for (pi = 0; pi < dp->npredicted; pi++) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Installing predicted mv_rate for (%i, %i): %i",
       dp->predicted_mvs[pi]->vx, dp->predicted_mvs[pi]->vy,
       dp->states[si].pred_mv_rates[pi]));
      dp->predicted_mvs[pi]->mv_rate = dp->states[si].pred_mv_rates[pi];
    }
    si = dp->states[si].prevsi;
  }
}

static ogg_int32_t od_mv_est_refine_row(od_mv_est_ctx *est,
 int ref, int vy, int log_dsz, int mv_res, const int *pattern_nsites,
 const od_pattern *pattern) {
  od_state *state;
  od_mv_grid_pt *grid;
  od_mv_grid_pt *pmvg;
  od_mv_grid_pt *mvg;
  od_mv_dp_node *dp_node;
  od_mv_dp_state *cstate;
  od_mv_dp_state *pstate;
  ogg_int32_t dcost;
  int nhmvbs;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int labels_only;
  int nsites;
  int sitei;
  int site;
  int curx;
  int cury;
  int vx;
  int b;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  labels_only = vy < 2 || vy > nvmvbs - 2;
  grid = state->mv_grid[vy];
  dcost = 0;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Refining row %i (%i)...", vy, (vy - 2) << 2));
  for (vx = 0;; vx++) {
    ogg_int32_t block_sads[OD_DP_NSTATES_MAX][OD_DP_NBLOCKS_MAX];
    ogg_int32_t best_cost;
    ogg_int32_t cost;
    ogg_int32_t best_dd;
    ogg_int32_t dd;
    int cur_mv_rates[OD_DP_NSTATES_MAX];
    int pred_mv_rates[OD_DP_NSTATES_MAX][OD_DP_NPREDICTED_MAX];
    int best_dr;
    int dr;
    int best_si;
    int si;
    for (; vx <= nhmvbs && !grid[vx].valid; vx++);
    if (vx > nhmvbs) break;
    level = OD_MC_LEVEL[vy & 3][vx & 3];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Starting DP at vertex %i (%i), level %i", vx, (vx - 2) << 2, level));
    log_mvb_sz = (4 - level) >> 1;
    mvb_sz = 1 << log_mvb_sz;
    mvg = grid + vx;
    curx = mvg->mv[0];
    cury = mvg->mv[1];
    dp_node = est->dp_nodes;
    od_mv_dp_row_init(est, dp_node, vx, vy, NULL);
    od_mv_dp_first_row_block_setup(est, dp_node, vx, vy);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
    od_mv_dp_get_sad_change8(est, ref, dp_node, block_sads[0]);
    /*Compute the set of states for the first node.*/
    if (vx >= 2 && !labels_only) {
      b = od_mv_est_get_boundary_case(state, vx, vy, curx, cury,
       1 << log_dsz, log_mvb_sz + 2);
      nsites = pattern_nsites[b];
    }
    else b = nsites = 0;
    for (sitei = 0, site = 4;; sitei++) {
      cstate = dp_node[0].states + sitei;
      cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
      cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
      cstate->prevsi = -1;
      mvg->mv[0] = cstate->mv[0];
      mvg->mv[1] = cstate->mv[1];
      cstate->dr = od_mv_dp_get_rate_change(state, dp_node,
       &cstate->mv_rate, cstate->pred_mv_rates, -1, mv_res);
      cstate->dd = od_mv_dp_get_sad_change8(est, ref, dp_node,
       cstate->block_sads);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "State: %i (%g, %g)  dr: %i  dd: %i  dopt: %i",
       sitei, 0.125*cstate->mv[0], 0.125*cstate->mv[1], cstate->dr, cstate->dd,
       cstate->dr*est->lambda + (cstate->dd << OD_ERROR_SCALE)));
      if (sitei >= nsites) break;
      site = pattern[b][sitei];
    }
    dp_node[0].nstates = nsites + 1;
    pmvg = mvg;
    while (vx < nhmvbs) {
      /*Find the next available MV to advance to.*/
      if (level & 1) {
        if (!grid[vx + mvb_sz].valid) {
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Gap found at %i (%i), stopping", vx, (vx - 2) << 2));
          break;
        }
        else if (level >= 3) vx++;
        else if (!grid[vx + 1].valid) vx += mvb_sz;
        else vx++;
      }
      else if (level >= 4) vx++;
      else if (!grid[vx + (mvb_sz >> 1)].valid) vx += mvb_sz;
      else if (level >= 2 || !grid[vx + 1].valid) vx += mvb_sz >> 1;
      else vx++;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        od_mv_dp_restore_row_state(dp_node);
        od_mv_dp_animate_state(state, ref, dp_node, 0);
      }
#endif
      level = OD_MC_LEVEL[vy & 3][vx & 3];
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Continuing DP at vertex %i (%i), level %i",
       vx, (vx - 2) << 2, level));
      log_mvb_sz = (4 - level) >> 1;
      mvb_sz = 1 << log_mvb_sz;
      mvg = grid + vx;
      curx = mvg->mv[0];
      cury = mvg->mv[1];
      od_mv_dp_row_init(est, dp_node + 1, vx, vy, dp_node);
      od_mv_dp_prev_row_block_setup(est, dp_node + 1, vx, vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[0] = dp_node[0].original_mv[0];
        od_mv_dp_get_sad_change8(est, ref, dp_node + 1, block_sads[0]);
      }
      /*Compute the set of states for this node.*/
      if (vx <= nhmvbs - 2 && !labels_only) {
        b = od_mv_est_get_boundary_case(state,
         vx, vy, curx, cury, 1 << log_dsz, log_mvb_sz + 2);
        nsites = pattern_nsites[b];
      }
      else nsites = 0;
      for (sitei = 0, site = 4;; sitei++) {
        cstate = dp_node[1].states + sitei;
        cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
        cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
        best_si = 0;
        best_dr = dp_node[0].states[0].dr;
        best_dd = dp_node[0].states[0].dd;
        best_cost = INT_MAX;
        mvg->mv[0] = cstate->mv[0];
        mvg->mv[1] = cstate->mv[1];
        for (si = 0; si < dp_node[0].nstates; si++) {
          pstate = dp_node[0].states + si;
          /*Get the rate change for this state using previous state si.
            This automatically loads the required bits of the trellis path into
             the grid, like the previous MV.*/
          cstate->dr = od_mv_dp_get_rate_change(state, dp_node + 1,
           cur_mv_rates + si, pred_mv_rates[si], si, mv_res);
          /*Test against the previous state.*/
          dr = pstate->dr + cstate->dr;
          dd = pstate->dd + od_mv_dp_get_sad_change8(est,
           ref, dp_node + 1, block_sads[si]);
          cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "State: %2i (%7g, %7g) P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
           sitei, 0.125*cstate->mv[0], 0.125*cstate->mv[1],
           si, dr, dd, cost));
          if (cost < best_cost) {
            best_si = si;
            best_cost = cost;
            best_dd = dd;
            best_dr = dr;
          }
        }
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: %2i  Best P.State: %i", sitei, best_si));
        cstate->prevsi = best_si;
        cstate->dr = best_dr;
        cstate->dd = best_dd;
        OD_COPY(cstate->block_sads, block_sads[best_si], dp_node[1].nblocks);
        cstate->mv_rate = cur_mv_rates[best_si];
        OD_COPY(cstate->pred_mv_rates, pred_mv_rates[best_si],
         dp_node[1].npredicted);
        if (sitei >= nsites) break;
        site = pattern[b][sitei];
      }
      dp_node[1].nstates = nsites + 1;
      dp_node++;
      pmvg = mvg;
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Finished DP at vertex %i (%i)",
     dp_node[0].mv->vx, (dp_node[0].mv->vx - 2) << 2));
    best_si = 0;
    best_cost = INT_MAX;
    /*TODO: Once we stop optimizing at arbitrary places, we'll need to
       compute the rate change of MVs we didn't get to.*/
    dp_node[1].npredicted = dp_node[1].npred_changeable = 0;
    if (dp_node[0].mv->vx < nhmvbs - 2) {
      od_mv_dp_last_row_block_setup(est, dp_node + 1, dp_node[0].mv->vx, vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[1] = dp_node[0].original_mv[1];
        od_mv_dp_get_sad_change8(est, ref, dp_node + 1, block_sads[0]);
      }
      for (si = 0; si < dp_node[0].nstates; si++) {
        pstate = dp_node[0].states + si;
        pmvg->mv[0] = pstate->mv[0];
        pmvg->mv[1] = pstate->mv[1];
        /*Test against the state with a following edge.*/
        dr = pstate->dr;
        dd = pstate->dd + od_mv_dp_get_sad_change8(est,
         ref, dp_node + 1, block_sads[si]);
        cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: --  P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
         si, dr, dd, cost));
        if (cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else {
      dp_node[1].nblocks = 0;
      for (si = 0; si < dp_node[0].nstates; si++) {
        pstate = dp_node[0].states + si;
        cost = pstate->dr*est->lambda + (pstate->dd << OD_ERROR_SCALE);
        if (cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Best P.State: %i dopt: %7i", best_si, best_cost));
    if (best_cost > 0) {
      /*Our optimal path is worse than what we started with!
        Restore the original state and give up.*/
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Best cost (%7i) > 0! Optimization failed.", best_cost));
      od_mv_dp_restore_row_state(dp_node);
    }
    else {
      int bi;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        char iter_label[16];
        od_mv_dp_restore_row_state(dp_node);
        od_mv_dp_animate_state(state, ref, dp_node, 0);
        od_mv_dp_install_row_state(dp_node + 1, best_si);
        od_state_mc_predict(state, ref);
        od_state_fill_vis(state);
        sprintf(iter_label, "ani%08i", state->ani_iter++);
        od_state_dump_img(state, &state->vis_img, iter_label);
      }
#endif
      /*Update the state along the optimal path.*/
      od_mv_dp_install_row_state(dp_node + 1, best_si);
      /*Store the SADs from this last node, too.*/
      for (bi = 0; bi < dp_node[1].nblocks; bi++) {
        dp_node[1].blocks[bi]->sad = block_sads[best_si][bi];
      }
      dcost += best_cost;
    }
  }
#if defined(OD_LOGGING_ENABLED)
  if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN)) {
    od_mv_est_check_rd_state(est, ref, mv_res);
  }
#endif
  return dcost;
}

/*Column refinement.*/

static void od_mv_dp_col_init(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy, od_mv_dp_node *prev_dp) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  state = &est->enc->state;
  dp->mv = est->mvs[vy] + vx;
  dp->mvg = state->mv_grid[vy] + vx;
  dp->original_mv[0] = dp->mvg->mv[0];
  dp->original_mv[1] = dp->mvg->mv[1];
  dp->original_mv_rate = dp->mv->mv_rate;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  if (vx < 2 || vx > nhmvbs - 2 || vy < 2 || vy > nvmvbs - 2) {
    /*Strictly speaking, we may be used to predict others, but since our MV
       can't possibly change, neither can their rate.*/
    dp->npredicted = dp->npred_changeable = 0;
    /*No one else is used to predict us, or any other MV we predict.
      However, we may still need to load the previous MV into the grid to
       estimate our SADs properly.*/
    dp->min_predictor_node = prev_dp;
  }
  else {
    int level;
    int pred_hist;
    int npred;
    int nchangeable;
    int pi;
    /*Get the list of MVs we help predict.*/
    level = OD_MC_LEVEL[vy & 3][vx & 3];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Initializing node (%i, %i) [%i, %i] at level %i:",
     vx, vy, (vx - 2) <<2, (vy - 2) <<2, level));
    npred = nchangeable = 0;
    for (pi = 0; pi < OD_NPREDICTED[level]; pi++) {
      int px;
      int py;
      px = vx + OD_COL_PREDICTED[level][pi][0];
      if (px < 2 || px > nhmvbs - 2) continue;
      py = vy + OD_COL_PREDICTED[level][pi][1];
      if (py < 2 || py > nvmvbs - 2) continue;
      if (state->mv_grid[py][px].valid) {
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "Adding (%i, %i) as a PREDICTED MV.", px, py));
        dp->predicted_mvgs[npred] = state->mv_grid[py] + px;
        dp->predicted_mvs[npred] = est->mvs[py] + px;
        if (pi < OD_NCOL_PRED_CHANGEABLE[level]) {
          OD_LOG((OD_LOG_MOTION_ESTIMATION,
           OD_LOG_DEBUG, "It is CHANGEABLE."));
          dp->original_mv_rates[npred] = est->mvs[py][px].mv_rate;
          nchangeable++;
        }
        npred++;
      }
    }
    dp->npredicted = npred;
    dp->npred_changeable = nchangeable;
    /*Now, figure out the earliest DP node that influences our own prediction,
       or that of one of the other MVs we predict.*/
    pred_hist = OD_COL_PRED_HIST_SIZE[level];
    if (prev_dp != NULL && prev_dp->mv->vy >= vy - pred_hist) {
      od_mv_dp_node *dp_pred;
      for (dp_pred = prev_dp; dp_pred->mv->vy > vy-pred_hist
       && dp_pred->states[0].prevsi >= 0; dp_pred--);
      if (dp_pred->mv->vy < vy - pred_hist) dp_pred++;
      dp->min_predictor_node = dp_pred;
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "State will be restored back to (%i, %i).",
       dp_pred->mv->vx, dp_pred->mv->vy));
    }
    else dp->min_predictor_node = NULL;
  }
}

static void od_mv_dp_first_col_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  nblocks = 0;
  if (vy > 2) {
    if (level >= 3) {
      if (vx >= mvb_sz) {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      }
      if (vx <= nhmvbs - mvb_sz) {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
      }
    }
    else {
      int half_mvb_sz;
      int mvb_off;
      half_mvb_sz = mvb_sz >> 1;
      if (vx >= mvb_sz) {
        if (state->mv_grid[vy - half_mvb_sz][vx - half_mvb_sz].valid) {
          int ulvx;
          int ulvy;
          ulvx = vx - (half_mvb_sz >> 1);
          ulvy = vy - (half_mvb_sz >> 1);
          if (level > 0
           || !state->mv_grid[ulvy][ulvx].valid) {
            mvb_off = half_mvb_sz;
          }
          else mvb_off = half_mvb_sz >> 1;
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - mvb_off;
          if (!state->mv_grid[vy][vx - mvb_off].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - mvb_off] + vx - (mvb_off << 1);
          }
          if (!state->mv_grid[vy - mvb_off][vx].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx - mvb_off;
          }
        }
        else dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      }
      if (vx <= nhmvbs - mvb_sz) {
        if (state->mv_grid[vy - half_mvb_sz][vx + half_mvb_sz].valid) {
          int urvx;
          int urvy;
          urvx = vx + (half_mvb_sz >> 1);
          urvy = vy + (half_mvb_sz >> 1);
          if (level > 0
           || !state->mv_grid[urvy][urvx].valid) {
            mvb_off = half_mvb_sz;
          }
          else mvb_off = half_mvb_sz >> 1;
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx;
          if (!state->mv_grid[vy][vx + mvb_off].valid) {
            dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx + mvb_off;
          }
          if (!state->mv_grid[vy - mvb_off][vx].valid) {
            dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx;
          }
        }
        else dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
      }
    }
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_prev_col_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nhmvbs;
  int level;
  int prev_level;
  int log_mvb_sz;
  int prev_log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  prev_level = OD_MC_LEVEL[(vy - mvb_sz) & 3][vx & 3];
  prev_log_mvb_sz = (4 - prev_level) >> 1;
  nblocks = 0;
  if (level >= 3) {
    if (vx >= mvb_sz) {
      dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
      if (prev_log_mvb_sz > log_mvb_sz
       && !state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid) {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - (mvb_sz << 1);
      }
    }
    if (vx <= nhmvbs - mvb_sz) {
      dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
      if (prev_log_mvb_sz > log_mvb_sz
       && !state->mv_grid[vy - mvb_sz][vx + mvb_sz].valid) {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx + mvb_sz;
      }
    }
  }
  else {
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz = mvb_sz >> 1;
    if (vx >= mvb_sz) {
      if (state->mv_grid[vy - half_mvb_sz][vx - half_mvb_sz].valid) {
        int ulvx;
        int ulvy;
        ulvx = vx - (half_mvb_sz >> 1);
        ulvy = vy - (half_mvb_sz >> 1);
        if (level > 0
         || !state->mv_grid[ulvy][ulvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - mvb_off;
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx - (mvb_off << 1);
        }
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx - mvb_off;
          if (!state->mv_grid[vy - (mvb_off << 1)][vx - mvb_off].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx - (mvb_off << 1);
          }
        }
      }
      else {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - mvb_sz;
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy - mvb_sz][vx - mvb_sz].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx - (mvb_sz << 1);
        }
      }
    }
    if (vx <= nhmvbs - mvb_sz) {
      if (state->mv_grid[vy - half_mvb_sz][vx + half_mvb_sz].valid) {
        int urvx;
        int urvy;
        urvx = vx + (half_mvb_sz >> 1);
        urvy = vy - (half_mvb_sz >> 1);
        if (level > 0
         || !state->mv_grid[urvy][urvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx;
        if (!state->mv_grid[vy][vx + mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_off] + vx + mvb_off;
        }
        if (!state->mv_grid[vy - mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - (mvb_off << 1)] + vx;
          if (!state->mv_grid[vy - (mvb_off << 1)][vx + mvb_off].valid) {
            dp->blocks[nblocks++] =
             est->mvs[vy - (mvb_off << 1)] + vx + mvb_off;
          }
        }
      }
      else {
        dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx;
        if (prev_log_mvb_sz > log_mvb_sz
         && !state->mv_grid[vy - mvb_sz][vx + mvb_sz].valid) {
          dp->blocks[nblocks++] = est->mvs[vy - mvb_sz] + vx + mvb_sz;
        }
      }
    }
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_last_col_block_setup(od_mv_est_ctx *est,
 od_mv_dp_node *dp, int vx, int vy) {
  od_state *state;
  int nhmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int nblocks;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  level = OD_MC_LEVEL[vy & 3][vx & 3];
  log_mvb_sz = (4 - level) >> 1;
  mvb_sz = 1 << log_mvb_sz;
  nblocks = 0;
  if (level >= 3) {
    if (vx >= mvb_sz) dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
    if (vx <= nhmvbs - mvb_sz) dp->blocks[nblocks++] = est->mvs[vy] + vx;
  }
  else {
    int half_mvb_sz;
    int mvb_off;
    half_mvb_sz = mvb_sz >> 1;
    if (vx >= mvb_sz) {
      if (state->mv_grid[vy + half_mvb_sz][vx - half_mvb_sz].valid) {
        int dlvx;
        int dlvy;
        dlvx = vx - (half_mvb_sz >> 1);
        dlvy = vy + (half_mvb_sz >> 1);
        if (level > 0
         || !state->mv_grid[dlvy][dlvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_off;
        if (!state->mv_grid[vy][vx - mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy] + vx - (mvb_off << 1);
        }
        if (!state->mv_grid[vy + mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy + mvb_off] + vx - mvb_off;
        }
      }
      else dp->blocks[nblocks++] = est->mvs[vy] + vx - mvb_sz;
    }
    if (vx <= nhmvbs - mvb_sz) {
      if (state->mv_grid[vy + half_mvb_sz][vx + half_mvb_sz].valid) {
        int drvx;
        int drvy;
        drvx = vx + (half_mvb_sz >> 1);
        drvy = vy + (half_mvb_sz >> 1);
        if (level > 0
         || !state->mv_grid[drvy][drvx].valid) {
          mvb_off = half_mvb_sz;
        }
        else mvb_off = half_mvb_sz >> 1;
        dp->blocks[nblocks++] = est->mvs[vy] + vx;
        if (!state->mv_grid[vy][vx + mvb_off].valid) {
          dp->blocks[nblocks++] = est->mvs[vy] + vx + mvb_off;
        }
        if (!state->mv_grid[vy + mvb_off][vx].valid) {
          dp->blocks[nblocks++] = est->mvs[vy + mvb_off] + vx;
        }
      }
      else dp->blocks[nblocks++] = est->mvs[vy] + vx;
    }
  }
  dp->nblocks = nblocks;
}

static void od_mv_dp_restore_col_state(od_mv_dp_node *dp) {
  od_mv_grid_pt *mvg;
  int pi;
  do {
    /*Restore the state for this MV itself.*/
    dp->mv->mv_rate = dp->original_mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->original_mv[0];
    mvg->mv[1] = dp->original_mv[1];
    for (pi = 0; pi < dp->npred_changeable; pi++) {
      /*Restore the state for the MVs this one predicted.*/
      dp->predicted_mvs[pi]->mv_rate = dp->original_mv_rates[pi];
    }
  }
  while ((dp--)->states[0].prevsi >= 0);
}

static void od_mv_dp_install_col_state(od_mv_dp_node *dp, int prevsi) {
  od_mv_dp_node *dp0;
  od_mv_grid_pt *mvg;
  int            nextsi;
  int            si;
  int            pi;
  int            bi;
  /*We must install the state going FORWARDS, since the pred_mv_rates may have
     changed several times over the course of the trellis.
    Therefore, first we reverse all of the prevsi pointers to make them act
     like nextsi pointers.*/
  nextsi = -1;
  dp0 = dp;
  for (si = prevsi; si >= 0; si = prevsi) {
    dp--;
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Node %i, prevsi: %i nextsi: %i", (int)(dp0 - dp), prevsi, nextsi));
    prevsi = dp->states[si].prevsi;
    dp->states[si].prevsi = nextsi;
    nextsi = si;
  }
  /*Now we traverse forward installing the rest of the state.*/
  for (si = nextsi; dp < dp0; dp++) {
    /*Install the state for this MV itself.*/
    dp->mv->mv_rate = dp->states[si].mv_rate;
    mvg = dp->mvg;
    mvg->mv[0] = dp->states[si].mv[0];
    mvg->mv[1] = dp->states[si].mv[1];
    /*Install the new block SADs.*/
    for (bi = 0; bi < dp->nblocks; bi++) {
      dp->blocks[bi]->sad = dp->states[si].block_sads[bi];
    }
    /*Install the state for the MVs this one predicted.*/
    for (pi = 0; pi < dp->npredicted; pi++) {
      dp->predicted_mvs[pi]->mv_rate = dp->states[si].pred_mv_rates[pi];
    }
    si = dp->states[si].prevsi;
  }
}

static ogg_int32_t od_mv_est_refine_col(od_mv_est_ctx *est,
 int ref, int vx, int log_dsz, int mv_res, const int *pattern_nsites,
 const od_pattern *pattern) {
  od_state *state;
  od_mv_grid_pt **grid;
  od_mv_grid_pt *pmvg;
  od_mv_grid_pt *mvg;
  od_mv_dp_node *dp_node;
  od_mv_dp_state *cstate;
  od_mv_dp_state *pstate;
  ogg_int32_t dcost;
  int nhmvbs;
  int nvmvbs;
  int level;
  int log_mvb_sz;
  int mvb_sz;
  int labels_only;
  int nsites;
  int sitei;
  int site;
  int curx;
  int cury;
  int vy;
  int b;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  labels_only = vx < 2 || vx > nhmvbs - 2;
  grid = state->mv_grid;
  dcost = 0;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Refining column %i (%i)...", vx, (vx - 2) << 2));
  for (vy = 0;; vy++) {
    ogg_int32_t block_sads[OD_DP_NSTATES_MAX][OD_DP_NBLOCKS_MAX];
    ogg_int32_t best_cost;
    ogg_int32_t cost;
    ogg_int32_t best_dd;
    ogg_int32_t dd;
    int cur_mv_rates[OD_DP_NSTATES_MAX];
    int pred_mv_rates[OD_DP_NSTATES_MAX][OD_DP_NPREDICTED_MAX];
    int best_dr;
    int dr;
    int best_si;
    int si;
    for (; vy <= nvmvbs && !grid[vy][vx].valid; vy++);
    if (vy > nvmvbs) break;
    level = OD_MC_LEVEL[vy & 3][vx & 3];
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Starting DP at vertex %i (%i), level %i", vy, (vy - 2) << 2, level));
    log_mvb_sz = (4 - level) >> 1;
    mvb_sz = 1 << log_mvb_sz;
    mvg = grid[vy] + vx;
    curx = mvg->mv[0];
    cury = mvg->mv[1];
    dp_node = est->dp_nodes;
    od_mv_dp_col_init(est, dp_node, vx, vy, NULL);
    od_mv_dp_first_col_block_setup(est, dp_node, vx, vy);
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
    if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
      od_mv_dp_get_sad_change8(est, ref, dp_node, block_sads[0]);
    }
    /*Compute the set of states for the first node.*/
    if (vy >= 2 && !labels_only) {
      b = od_mv_est_get_boundary_case(state,
       vx, vy, curx, cury, 1 << log_dsz, log_mvb_sz + 2);
      nsites = pattern_nsites[b];
    }
    else b = nsites = 0;
    for (sitei = 0, site = 4;; sitei++) {
      cstate = dp_node[0].states + sitei;
      cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
      cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
      cstate->prevsi = -1;
      mvg->mv[0] = cstate->mv[0];
      mvg->mv[1] = cstate->mv[1];
      cstate->dr = od_mv_dp_get_rate_change(state, dp_node,
       &cstate->mv_rate, cstate->pred_mv_rates, -1, mv_res);
      cstate->dd = od_mv_dp_get_sad_change8(est, ref, dp_node,
       cstate->block_sads);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "State: %i  dr: %i  dd: %i  dopt: %i", sitei, cstate->dr, cstate->dd,
       cstate->dr*est->lambda + (cstate->dd << OD_ERROR_SCALE)));
      if (sitei >= nsites) break;
      site = pattern[b][sitei];
    }
    dp_node[0].nstates = nsites + 1;
    pmvg = mvg;
    while (vy < nvmvbs) {
      /*Find the next available MV to advance to.*/
      if (level & 1) {
        if (!grid[vy + mvb_sz][vx].valid) {
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "Gap found at %i (%i), stopping", vy, (vy - 2) << 2));
          break;
        }
        else if (level >= 3) vy++;
        else if (!grid[vy + 1][vx].valid) vy += mvb_sz;
        else vy++;
      }
      else if (level >= 4) vy++;
      else if (!grid[vy + (mvb_sz >> 1)][vx].valid) vy += mvb_sz;
      else if (level >= 2 || !grid[vy + 1][vx].valid) vy += mvb_sz >> 1;
      else vy++;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        od_mv_dp_restore_col_state(dp_node);
        od_mv_dp_animate_state(state, ref, dp_node, 0);
      }
#endif
      level = OD_MC_LEVEL[vy & 3][vx & 3];
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Continuing DP at vertex %i (%i), level %i", vy, (vy - 2) << 2, level));
      log_mvb_sz = (4 - level) >> 1;
      mvb_sz = 1 << log_mvb_sz;
      mvg = grid[vy] + vx;
      curx = mvg->mv[0];
      cury = mvg->mv[1];
      od_mv_dp_col_init(est, dp_node + 1, vx, vy, dp_node);
      od_mv_dp_prev_col_block_setup(est, dp_node + 1, vx, vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[0] = dp_node[0].original_mv[0];
        od_mv_dp_get_sad_change8(est, ref, dp_node + 1, block_sads[0]);
      }
      /*Compute the set of states for this node.*/
      if (vy <= nvmvbs - 2 && !labels_only) {
        b = od_mv_est_get_boundary_case(state,
         vx, vy, curx, cury, 1 << log_dsz, log_mvb_sz + 2);
        nsites = pattern_nsites[b];
      }
      else nsites = 0;
      for (sitei = 0, site = 4;; sitei++) {
        cstate = dp_node[1].states + sitei;
        cstate->mv[0] = curx + (OD_SITE_DX[site] << log_dsz);
        cstate->mv[1] = cury + (OD_SITE_DY[site] << log_dsz);
        best_si = 0;
        best_dr = dp_node[0].states[0].dr;
        best_dd = dp_node[0].states[0].dd;
        best_cost = INT_MAX;
        mvg->mv[0] = cstate->mv[0];
        mvg->mv[1] = cstate->mv[1];
        for (si = 0; si < dp_node[0].nstates; si++) {
          pstate = dp_node[0].states + si;
          /*Get the rate change for this state using previous state si.
            This automatically loads the required bits of the trellis path into
             the grid, like the previous MV.*/
          cstate->dr = od_mv_dp_get_rate_change(state, dp_node + 1,
           cur_mv_rates + si, pred_mv_rates[si], si, mv_res);
          /*Test against the previous state.*/
          dr = pstate->dr + cstate->dr;
          dd = pstate->dd + od_mv_dp_get_sad_change8(est,
           ref, dp_node + 1, block_sads[si]);
          cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
          OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
           "State: %2i  P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
           sitei, si, dr, dd, cost));
          if (cost < best_cost) {
            best_si = si;
            best_cost = cost;
            best_dd = dd;
            best_dr = dr;
          }
        }
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: %2i  Best P.State: %i", sitei, best_si));
        cstate->prevsi = best_si;
        cstate->dr = best_dr;
        cstate->dd = best_dd;
        OD_COPY(cstate->block_sads, block_sads[cstate->prevsi],
         dp_node[1].nblocks);
        cstate->mv_rate = cur_mv_rates[best_si];
        OD_COPY(cstate->pred_mv_rates, pred_mv_rates[best_si],
         dp_node[1].npredicted);
        if (sitei >= nsites) break;
        site = pattern[b][sitei];
      }
      dp_node[1].nstates = nsites + 1;
      dp_node++;
      pmvg = mvg;
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Finished DP at vertex %i (%i)",
     dp_node[0].mv->vy, (dp_node[0].mv->vy - 2) << 2));
    best_si = 0;
    best_cost = INT_MAX;
    /*TODO: Once we stop optimizing at arbitrary places, we'll need to
       compute the rate change of MVs we didn't get to.*/
    dp_node[1].npredicted = dp_node[1].npred_changeable = 0;
    if (dp_node[0].mv->vy < nvmvbs - 2) {
      od_mv_dp_last_col_block_setup(est, dp_node + 1, vx, dp_node[0].mv->vy);
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "TESTING block SADs:"));
      if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG)) {
        pmvg->mv[0] = dp_node[0].original_mv[0];
        pmvg->mv[1] = dp_node[0].original_mv[1];
        od_mv_dp_get_sad_change8(est, ref, dp_node + 1, block_sads[0]);
      }
      for (si = 0; si < dp_node[0].nstates; si++) {
        pstate = dp_node[0].states + si;
        pmvg->mv[0] = pstate->mv[0];
        pmvg->mv[1] = pstate->mv[1];
        /*Test against the state with a following edge.*/
        dr = pstate->dr;
        dd = pstate->dd + od_mv_dp_get_sad_change8(est,
         ref, dp_node + 1, block_sads[si]);
        cost = dr*est->lambda + (dd << OD_ERROR_SCALE);
        OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
         "State: --  P.State: %i  dr: %3i  dd: %6i  dopt: %7i",
         si, dr, dd, cost));
        if (best_si < 0 || cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    /*There are no blocks to accumulate SAD after this one, so pick the best
       state so far.*/
    else {
      dp_node[1].nblocks = 0;
      for (si = 0; si < dp_node[0].nstates; si++) {
        dp_node[1].nblocks = 0;
        pstate = dp_node[0].states + si;
        cost = pstate->dr*est->lambda + (pstate->dd << OD_ERROR_SCALE);
        if (best_si < 0 || cost < best_cost) {
          best_si = si;
          best_cost = cost;
        }
      }
    }
    OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
     "Best P.State: %i dopt: %7i", best_si, best_cost));
    if (best_cost > 0) {
      /*Our optimal path is worse than what we started with!
        Restore the original state and give up.*/
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
       "Best cost (%7i) > 0! Optimization failed.", best_cost));
      od_mv_dp_restore_col_state(dp_node);
    }
    else {
      int bi;
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
      if (daala_granule_basetime(state, state->cur_time) == ANI_FRAME) {
        char iter_label[16];
        od_mv_dp_restore_col_state(dp_node);
        od_mv_dp_animate_state(state, ref, dp_node, 0);
        od_mv_dp_install_col_state(dp_node + 1, best_si);
        od_state_mc_predict(state, ref);
        od_state_fill_vis(state);
        sprintf(iter_label, "ani%08i", state->ani_iter++);
        od_state_dump_img(state, &state->vis_img, iter_label);
      }
#endif
      /*Update the state along the optimal path.*/
      od_mv_dp_install_col_state(dp_node + 1, best_si);
      /*Store the SADs from this last node, too.*/
      for (bi = 0; bi < dp_node[1].nblocks; bi++) {
        dp_node[1].blocks[bi]->sad = block_sads[best_si][bi];
      }
      dcost += best_cost;
    }
  }
#if defined(OD_LOGGING_ENABLED)
  if (od_logging_active(OD_LOG_MOTION_ESTIMATION, OD_LOG_WARN)) {
    od_mv_est_check_rd_state(est, ref, mv_res);
  }
#endif
  return dcost;
}

static ogg_int32_t od_mv_est_refine(od_mv_est_ctx *est, int ref, int log_dsz,
 int mv_res, const int *pattern_nsites, const od_pattern *pattern) {
  od_state *state;
  ogg_int32_t dcost;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG,
   "Refining with displacements of %0g and 1/%i pel MV resolution.",
   (1 << log_dsz)*0.125, 1 << (3 - mv_res)));
  dcost = 0;
  for (vy = 0; vy <= nvmvbs; vy++) {
    if (est->row_counts[vy]) {
      dcost += od_mv_est_refine_row(est, ref, vy, log_dsz, mv_res,
       pattern_nsites, pattern);
    }
  }
  for (vx = 0; vx <= nhmvbs; vx++) {
    if (est->col_counts[vx]) {
      dcost += od_mv_est_refine_col(est, ref, vx, log_dsz, mv_res,
       pattern_nsites, pattern);
    }
  }
  return dcost;
}

/*STAGE 4: Sub-pel Refinement.*/

/*Stores the full-pel MVs for use by EPZS^2 in the next frame before sub-pel
   refinement.*/
void od_mv_est_update_fullpel_mvs(od_mv_est_ctx *est, int ref) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  for (vy = 2; vy <= nvmvbs - 2; vy++) {
    for (vx = 2; vx <= nhmvbs - 2; vx++) {
      od_mv_grid_pt *mvg;
      od_mv_node *mv;
      mvg = state->mv_grid[vy] + vx;
      if (!mvg->valid) continue;
      mv = est->mvs[vy] + vx;
      mv->bma_mvs[0][ref][0] = mvg->mv[0] >> 3;
      mv->bma_mvs[0][ref][1] = mvg->mv[1] >> 3;
    }
  }
}

/*Sets the mv_rate of each node in the mesh, using the given MV resolution.
  Returns the change in rate.*/
int od_mv_est_update_mv_rates(od_mv_est_ctx *est, int mv_res) {
  od_state *state;
  int nhmvbs;
  int nvmvbs;
  int vx;
  int vy;
  int dr;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  dr = 0;
  for (vy = 2; vy <= nvmvbs - 2; vy++) {
    for (vx = 2; vx <= nhmvbs - 2; vx++) {
      od_mv_grid_pt *mvg;
      od_mv_node *mv;
      int pred[2];
      mvg = state->mv_grid[vy] + vx;
      if (!mvg->valid) continue;
      mv = est->mvs[vy] + vx;
      od_state_get_predictor(state, pred,
       vx, vy, OD_MC_LEVEL[vy & 3][vx & 3], mv_res);
      dr -= mv->mv_rate;
      mv->mv_rate = od_mv_est_bits(mvg->mv[0] >> mv_res, mvg->mv[1] >> mv_res,
       pred[0], pred[1]);
      dr += mv->mv_rate;
    }
  }
  return dr;
}

od_mv_est_ctx *od_mv_est_alloc(od_enc_ctx *enc) {
  od_mv_est_ctx *ret;
  ret = (od_mv_est_ctx *)_ogg_malloc(sizeof(*ret));
  od_mv_est_init(ret, enc);
  return ret;
}

void od_mv_est_free(od_mv_est_ctx *est) {
  if (est != NULL) {
    od_mv_est_clear(est);
    _ogg_free(est);
  }
}

void od_mv_subpel_refine(od_mv_est_ctx *est, int ref, int cost_thresh) {
  od_state *state;
  od_mv_grid_pt **grid;
  ogg_int32_t dcost;
  ogg_int32_t subpel_cost;
  int nhmvbs;
  int nvmvbs;
  int mv_res;
  int best_mv_res;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  /*Save the fullpell MVs now for use by EPZS^2 on the next frame.
    We could also try rounding the results after refinement, I guess.
    I'm not sure it makes much difference*/
  od_mv_est_update_fullpel_mvs(est, ref);
  do {
    dcost = od_mv_est_refine(est, ref, 2, 2,
     OD_DIAMOND_NSITES, OD_DIAMOND_SITES);
  }
  while (dcost < cost_thresh);
  for (best_mv_res = mv_res = 2; mv_res-- > est->mv_res_min;) {
    subpel_cost = od_mv_est_update_mv_rates(est, mv_res)*est->lambda;
    /*If the rate penalty for refining is small, bump the termination threshold
       down to make sure we actually get a decent improvement.
      We make sure not to let it get too small, however, so we're not here all
       day (a motion field of all (0, 0)'s would have a rate penalty of 0!).*/
    cost_thresh = OD_MAXI(cost_thresh,
     -OD_MAXI(subpel_cost, 16 << OD_ERROR_SCALE));
    OD_COPY(est->refine_grid[0], state->mv_grid[0],
     (nhmvbs + 1)*(nvmvbs + 1));
    do {
      dcost = od_mv_est_refine(est, ref, mv_res, mv_res,
       OD_DIAMOND_NSITES, OD_DIAMOND_SITES);
      subpel_cost += dcost;
    }
    while (dcost < cost_thresh);
    if (subpel_cost >= 0) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_INFO,
       "1/%i refinement FAILED:    dopt %7i", 1 << (3 - mv_res), subpel_cost));
      grid = est->refine_grid;
      est->refine_grid = state->mv_grid;
      state->mv_grid = grid;
      break;
    }
    else {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_INFO,
       "1/%i refinement SUCCEEDED: dopt %7i", 1 << (3 - mv_res), subpel_cost));
      best_mv_res = mv_res;
    }
  }
  state->mv_res = best_mv_res;
}

void od_mv_est(od_mv_est_ctx *est, int ref, int lambda) {
  od_state *state;
  od_img_plane *iplane;
  ogg_int32_t dcost;
  int cost_thresh;
  int nhmvbs;
  int nvmvbs;
  int pli;
  state = &est->enc->state;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  iplane = state->io_imgs[OD_FRAME_INPUT].planes + 0;
  /*Sanitize user parameters*/
  est->level_min = OD_MINI(est->enc->params.mv_level_min,
   est->enc->params.mv_level_max);
  est->level_max = est->enc->params.mv_level_max;
  /*If the luma plane is decimated for some reason, then our distortions will
     be smaller, so scale lambda appropriately.*/
  est->lambda = lambda >> (iplane->xdec + iplane->ydec);
  /*Compute termination thresholds for EPZS^2.*/
  est->thresh1[0] = 16 >> (iplane->xdec + iplane->ydec);
  est->thresh1[1] = 64 >> (iplane->xdec + iplane->ydec);
  est->thresh1[2] = 256 >> (iplane->xdec + iplane->ydec);
  /*If we're using the chroma planes, then our distortions will be larger.
    Compensate by increasing lambda and the termination thresholds.*/
  if (est->flags & OD_MC_USE_CHROMA) {
    for (pli = 1; pli < state->io_imgs[OD_FRAME_INPUT].nplanes; pli++) {
      iplane = state->io_imgs[OD_FRAME_INPUT].planes + pli;
      est->lambda +=
       lambda >> (iplane->xdec + iplane->ydec + OD_MC_CHROMA_SCALE);
      est->thresh1[0] +=
       16 >> (iplane->xdec + iplane->ydec + OD_MC_CHROMA_SCALE);
      est->thresh1[1] +=
       64 >> (iplane->xdec + iplane->ydec + OD_MC_CHROMA_SCALE);
      est->thresh1[2] +=
       256 >> (iplane->xdec + iplane->ydec + OD_MC_CHROMA_SCALE);
    }
  }
  OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_DEBUG, "lambda: %i", est->lambda));
  est->thresh2_offs[0] = est->thresh1[0] >> 1;
  est->thresh2_offs[1] = est->thresh1[1] >> 1;
  est->thresh2_offs[2] = est->thresh1[2] >> 1;
  /*Accelerated predictor weights.*/
  est->mvapw[ref][0] = 0x20000;
  est->mvapw[ref][1] = 0x10000;
  /*TODO: Constant velocity predictor weight.*/
#if defined(OD_DUMP_IMAGES) && defined(OD_ANIMATE)
  /*Set some initial state.
    This would get reset eventually by the algorithm in a more convenient
     place, but is needed earlier by the visualization.*/
  if (daala_granule_basetime(&est->enc->state, est->enc->state.cur_time) ==
   ANI_FRAME) {
    od_state_mvs_clear(&est->enc->state);
  }
#endif
  od_mv_est_init_mvs(est, ref);
  od_mv_est_decimate(est, ref);
  /*This threshold is somewhat arbitrary.
    Chen and Willson use 6000 (with SSD as an error metric).
    We would like something more dependent on the frame size.
    For CIF, there are a maximum of 6992 vertices in the mesh, which is pretty
     close to 6000.
    With a SAD error metric like we use, the square root of 6000 would be a
     more appropriate value, however that gives a PSNR improvement of less than
     0.01 dB, and requires almost twice as many iterations to achieve.*/
  cost_thresh = -nhmvbs*nvmvbs << OD_ERROR_SCALE;
#if defined(OD_USE_LOGARITHMIC_REFINEMENT)
  /*Logarithmic search.*/
  do {
    dcost = 0;
    dcost += od_mv_est_refine(est, ref, 5, 2,
     OD_SQUARE_NSITES, OD_SQUARE_SITES);
    dcost += od_mv_est_refine(est, ref, 4, 2,
     OD_SQUARE_NSITES, OD_SQUARE_SITES);
    dcost += od_mv_est_refine(est, ref, 3, 2,
     OD_SQUARE_NSITES, OD_SQUARE_SITES);
  }
  while (dcost < cost_thresh);
#else
  /*Diamond search.
    This appears to give the same quality as the logarithmic search, but at
     nearly 10 times the speed.*/
  do {
    dcost = od_mv_est_refine(est, ref, 3, 2,
     OD_DIAMOND_NSITES, OD_DIAMOND_SITES);
  }
  while (dcost < cost_thresh);
#endif
  od_mv_subpel_refine(est, ref, cost_thresh);
  od_restore_fpu(state);
}
