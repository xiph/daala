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

#if !defined(_mc_H)
# define _mc_H (1)

typedef struct od_mv_grid_pt od_mv_grid_pt;

# include "state.h"

/*Indexing of vetices and edges:
  0--0--1
  |     |
  3     1
  |     |
  3--2--2

  Interpolation formulas:
  i = horizontal position in the destination block, 0...blk_sz - 1
  j = vertical position  in the destination block, 0...blk_sz - 1
  blk_sz = length of an edge in the destination block in pixels
  blk_sz2 = blk_sz*blk_sz
  w[k] = bilinear weight of vertex k:
  w[0] = (blk_sz - i)*(blk_sz - j)/blk_sz2
  w[1] = i*(blk_sz - j)/blk_sz2
  w[2] = i*j/blk_sz2
  w[3] = (blk_sz - i)*j/blk_sz2
  sw[k, oc] = BLEND weight of vertex k in a block with two unsplit edges and
   outside corner oc:
  sw[oc + 0, oc] = w[oc] + w[oc + 1]/2 + w[oc + 3]/2
  sw[oc + 1, oc] = w[oc + 1]/2
  sw[oc + 2, oc] = w[oc + 2]
  sw[oc + 3, oc] = w[oc + 3]/2
  m[k] = motion vector k
  src[m] = source image value of pixel (i, j) offset by the given motion
   vector, bilinearly interpolated from the values at half-pel locations.
  Basic OBMC blending is then just
    src[m[0]]*w[0] + src[m[1]]*w[1] + src[m[2]]*w[2] + src[m[3]]*w[3]

  Subpel at resolutions finer than halfpel is done with bilinear interpolation.
  First, horizontal interpolation is done, and the result is truncated to
   an integer, followed by vertical interpolation and truncation to an
   integer.

  Lots of indexing and finite differences where tiny errors may be made here.
  This is also a good candidate for SIMD optimization.*/

/*We very intentionally do NOT use the quadtree data structure proposed by
   Balmelli, as it is encumbered by patents.
  Fortunately, our relatively limited levels of subdivision makes it mostly
   unnecessary.*/
struct od_mv_grid_pt {
  /*The x, y offsets of the motion vector in units of 1/8th pixels.*/
  int mv[2];
  /*Whether or not this MV actually has a valid value.*/
  unsigned valid:1;
};

void od_mc_predict8(od_state *state, unsigned char *dst, int dystride,
 const unsigned char *src, int systride, const ogg_int32_t mvx[4],
 const ogg_int32_t mvy[4], int oc, int s, int log_xblk_sz, int log_yblk_sz);
void od_state_mvs_clear(od_state *state);
int od_state_get_predictor(od_state *state, int pred[2],
 int vx, int vy, int level, int mv_res);
int od_mv_level1_ctx(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level1_probz(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level2_ctx(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level2_probz(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level3_ctx(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level3_probz(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level4_ctx(od_mv_grid_pt **grid, int vx, int vy);
int od_mv_level4_probz(od_mv_grid_pt **grid, int vx, int vy);

#endif
