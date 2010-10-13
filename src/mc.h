/*
    Daala video codec
    Copyright (C) 2010 Timothy B. Terriberry and Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


#if !defined(_mc_H)
# define _mc_H (1)

typedef struct od_mv_grid_pt od_mv_grid_pt;

# include "state.h"

/*Interpolation along an edge by interpolating vectors.*/
#define OD_MC_INTERP_VECTOR (1)
/*Interpolation along an edge by blending samples from different vectors.*/
#define OD_MC_INTERP_BLEND  (0)

#define od_mc_interp_type(_a,_b,_c,_d) ((_a)|((_b)<<1)|((_c)<<2)|((_d)<<3))

#define OD_MC_INTERP_BBBB (0)
#define OD_MC_INTERP_VBBB (1)
#define OD_MC_INTERP_BVBB (2)
#define OD_MC_INTERP_VVBB (3)
#define OD_MC_INTERP_BBVB (4)
#define OD_MC_INTERP_VBVB (5)
#define OD_MC_INTERP_BVVB (6)
#define OD_MC_INTERP_VVVB (7)
#define OD_MC_INTERP_BBBV (8)
#define OD_MC_INTERP_VBBV (9)
#define OD_MC_INTERP_BVBV (10)
#define OD_MC_INTERP_VVBV (11)
#define OD_MC_INTERP_BBVV (12)
#define OD_MC_INTERP_VBVV (13)
#define OD_MC_INTERP_BVVV (14)
#define OD_MC_INTERP_VVVV (15)

/*
  Indexing of vetices and edges:
  0--0--1
  |     |
  3     1
  |     |
  3--2--2

  Interpolation formulas:
  i = horizontal position in the destination block, 0...blk_sz-1
  j = vertical position  in the destination block, 0...blk_sz-1
  blk_sz = length of an edge in the destination block in pixels
  blk_sz2 = blk_sz*blk_sz
  w[k]   = bilinear VECTOR weight of vertex k:
  w[0]   = (blk_sz-i)*(blk_sz-j)/blk_sz2
  w[1]   = i*(blk_sz-j)/blk_sz2
  w[2]   = i*j/blk_sz2
  w[3]   = (blk_sz-i)*j/blk_sz2
  u[k]   = bilinear weight of edge k:
  u[0]   = (blk_sz-j)/blk_sz
  u[1]   = i/blk_sz
  u[2]   = j/blk_sz
  u[3]   = (blk_sz-i)/blk_sz
  s[k]   = bilinear weight of vertex k:
  m[k]   = motion vector k
  s[k,c] = BLEND weight of vertex k in a split block whose outside corner is c:
  s[c+0,c] = w[c]+w[c+1]/2+w[c+3]/2
  s[c+1,c] = w[c+1]/2
  s[c+2,c] = w[c+2]
  s[c+3,c] = w[c+3]/2
  src[m] = source image value of pixel (i,j) offset by the given motion
   vector, bilinearly interpolated from the values at half-pel locations.
   First, horizontal interpolation is done, and the result is truncated to
   an integer, followed by vertical interpolation and truncation to an
   integer.

  Blending and interpolation functions are defined by the type of edges:
  V: vector interpolation edge.
     Along this edge, the motion vector is linearly interpolated from the
      value at one vertex to the other.
  B: blend edge
     Along this edge, the predicted value is linearly interpolated between
      that of the source image offset by the motion vector at one vertex
      to the value of the source image offset by the motion vector at the
      other vertex.
  This allows for C^0 switching between overlapped block motion compensation
   and smooth, per-pixel motion vector fields.
  The various interpolation formulas below are chosen so that:
   1) They have the proper form along each edge,
   2) They still match each other when a block is split to add a new motion
       vector in the middle (subject to the modifications to the blending
       weights listed above).
   2) All the vectors can be computed with simple bilinear interpolation,
   3) All the resulting source image pixels can be blended with simple
       bilinear interpolation.

  There are 6 similar classes of edge configurations:

  VVVV:
     src[m[0]*w[0]+m[1]*w[1]+m[2]*w[2]+m[3]*w[3]]
  BVVV (and rotations thereof):
     src[u[0]*m[0]+w[2]*m[2]+w[3]*m[3]]*w[0]+
      src[u[0]*m[1]+w[2]*m[2]+w[3]*m[3]]*w[1]+
      src[m[0]*w[0]+m[1]*w[1]+m[2]*w[2]+m[3]*w[3]]*u[2]
  BVBV (and rotation thereof):
     src[u[0]*m[0]+u[3]*m[2]]*u[2]+
      src[u[0]*m[1]+u[3]*m[3]]*u[1]
  VVBB (and rotations thereof):
     src[(w[0]+u[2])*m[0]+w[1]*m[1]]*w[0]+
      src[w[0]*m[0]+w[1]*m[1]+w[2]*m[2]+w[3]*m[3]]*w[1]+
      src[w[1]*m[1]+(w[2]+u[3])*m[2]]*w[2]+
      src[m[3]]*w[3]
  VBBB (and rotations thereof):
     src[(u[3]+w[2])*m[0]+w[1]*m[1]]*w[0]+
      src[w[0]*m[0]+(w[1]+u[2])*m[1]]*w[1]+
      src[m[2]]*w[2]+src[m[3]]*w[3]
  BBBB:
     src[m[0]]*w[0]+src[m[1]]*w[1]+src[m[2]]*w[2]+src[m[3]]*w[3]

  The remaining 10 edge configurations can be obtained by rotations of these
   formulas.

  Lots of indexing and finite differences where tiny errors may be made here.
  This is also a good candidate for SIMD optimization.
  Each case can be optimized separately.
  The first and last are expected to be more common than the others.
*/

/*We very intentionally do NOT use the quadtree data structure proposed by
   Balmelli, as it is encumbered by patents.
  Fortunately, our relatively limited levels of subdivision makes it mostly
   unnecessary.*/
struct od_mv_grid_pt{
  int      mv[2];
  unsigned valid:1;
  unsigned right:1;
  unsigned down:1;
};

void od_mc_setup_mvc(ogg_int32_t _dmv[4],const ogg_int32_t _mvs[4],
 const int _m[4],int _r,int _log_xblk_sz,int _log_yblk_sz);
void od_mc_predict8(od_state *_state,unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],int _interp_type,int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz);

#endif
