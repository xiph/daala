/*Daala video codec
Copyright (c) 2003-2013 Daala project contributors.  All rights reserved.
Copyright (c) 2015, Cisco Systems (Thor deblocking filter)

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

#include <stdlib.h>
#include <limits.h>
#include "filter.h"
#include "state.h"
#include "block_size.h"

/*Pre-/post-filter pairs of various sizes.
  For a FIR, PR, LP filter bank, the pre-filter must have the structure:

   0.5*{{I,J},{J,-I}}*{{U,0},{0,V}},{{I,J},{J,-I}}

  To create a pre-/post-filter for an M-sample block transform, I, J, U and V
   are blocks of size floor(M/2).
  I is the identity, J is the reversal matrix, and U and V are arbitrary
   invertible matrices.
  If U and V are orthonormal, then the pre-/post-filters are orthonormal.
  Otherwise they are merely biorthogonal.
  For filters with larger inputs, multiple pre-/post-filter stages can be used.
  However, one level should be sufficient to eliminate blocking artifacts, and
   additional levels would expand the influence of ringing artifacts.

  All the filters here are an even narrower example of the above structure.
  U is taken to be the identity, as it does not help much with energy
   compaction.
  V can be chosen so that the filter pair is (1,2) regular, assuming a block
   filter with 0 DC leakage, such as the DCT.
  This may be particularly important for the way that motion-compensated
   prediction is done.
  A 2-regular synthesis filter can reproduce a linear ramp from just the DC
   coefficient, which matches the linearly interpolated offsets that were
   subtracted out of the motion-compensation phase.

  In order to yield a fast implementation, the V filter is chosen to be of the
   form:
    x0 -s0----------...----+---- y0
            p0 |           | u0
    x1 -s1-----+----...--+------ y1
              p1 |       | u1
    x2 -s2-------+--...+-------- y2
                p2 |   | u2
    x3 -s3---------+...--------- y3
                     .
                     .
                     .
  Here, s(i) is a scaling factor, such that s(i) >= 1, to preserve perfect
   reconstruction given an integer implementation.
  p(i) and u(i) are arbitrary, so long as p(i)u(i) != -1, to keep the transform
   invertible.
  In order to make V (1,2) regular, we have the conditions:
    s0+M*u0==M
    (2*i+1)*s(i)+M*u(i)+M*(1-u(i-1))*p(i-1)==M, i in [1..K-2]
    (M-1)*s(K-1)+M*(1-u(K-2))*p(K-2)==M
   where K=floor(M/2).
  These impose additonal constraints on u(i) and p(i), derived from the
   constraints above, such as u(0) <= (M-1)/M.
  It is desirable to have u(i), p(i) and 1/s(i) be dyadic rationals, the latter
   to provide for a fast inverse transform.
  However, as can be seen by the constraints, it is very easy to make u(i) and
   p(i) be dyadic, or 1/s(i) be dyadic, but solutions with all of them dyadic
   are very sparse, and require at least s(0) to be a power of 2.
  Such solutions do not have a good coding gain, and so we settle for having
   s(i) be dyadic instead of 1/s(i).

  Or we use the form
    x0 -s0---------------+---- y0
                    p0 | | u0
    x1 -s1-----------+-+------ y1
                p1 | | u1
    x2 -s2-------+-+---------- y2
            p2 | | u2
    x3 -s3-----+-------------- y3
                 .
                 .
                 .
   which yields slightly higher coding gains, but the conditions for
    (1,2) regularity
    s0+M*u0==M
    (2*i+1)*s(i)+M*u(i)+(2*i-1)*s(i-1)*p(i-1)==M, i in [1..K-2]
    (M-1)*s(K-1)+(M-3)*s(K-2)*p(K-2)==M
   make dyadic rational approximations more sparse, such that the coding gains
   of the approximations are actually lower. This is selected when the TYPE3
   defines are set.

  The maximum denominator for all coefficients was allowed to be 64.*/

const od_filter_func OD_PRE_FILTER[OD_NBSIZES] = {
  od_pre_filter4,
  od_pre_filter8,
  od_pre_filter16,
  od_pre_filter32
};

const od_filter_func OD_POST_FILTER[OD_NBSIZES] = {
  od_post_filter4,
  od_post_filter8,
  od_post_filter16,
  od_post_filter32
};

const od_filter_dering_direction_func
 OD_DERING_DIRECTION_C[OD_DERINGSIZES] = {
  od_filter_dering_direction_4x4_c,
  od_filter_dering_direction_8x8_c
};

const od_filter_dering_orthogonal_func
 OD_DERING_ORTHOGONAL_C[OD_DERINGSIZES] = {
  od_filter_dering_orthogonal_4x4_c,
  od_filter_dering_orthogonal_8x8_c
};

/* Generated from gen_filter_tables.c. */
const int direction_offsets_table[8][3] = {
  {-1*OD_FILT_BSTRIDE + 1, -2*OD_FILT_BSTRIDE + 2, -3*OD_FILT_BSTRIDE + 3  },
  { 0*OD_FILT_BSTRIDE + 1, -1*OD_FILT_BSTRIDE + 2, -1*OD_FILT_BSTRIDE + 3  },
  { 0*OD_FILT_BSTRIDE + 1,  0*OD_FILT_BSTRIDE + 2,  0*OD_FILT_BSTRIDE + 3  },
  { 0*OD_FILT_BSTRIDE + 1,  1*OD_FILT_BSTRIDE + 2,  1*OD_FILT_BSTRIDE + 3  },
  { 1*OD_FILT_BSTRIDE + 1,  2*OD_FILT_BSTRIDE + 2,  3*OD_FILT_BSTRIDE + 3  },
  { 1*OD_FILT_BSTRIDE + 0,  2*OD_FILT_BSTRIDE + 1,  3*OD_FILT_BSTRIDE + 1  },
  { 1*OD_FILT_BSTRIDE + 0,  2*OD_FILT_BSTRIDE + 0,  3*OD_FILT_BSTRIDE + 0  },
  { 1*OD_FILT_BSTRIDE + 0,  2*OD_FILT_BSTRIDE - 1,  3*OD_FILT_BSTRIDE - 1  },
};

/** Strength of the bilinear smoothing for each plane. */
static const int OD_BILINEAR_STRENGTH[OD_NPLANES_MAX] = {5, 20, 20, 5};

/*Filter parameters for the pre/post filters.
  When changing these the intra-predictors in
  initdata.c must be updated.*/

/*R=f
  6-bit s0=1.328125, s1=1.171875, p0=-0.234375, u0=0.515625
  Ar95_Cg = 8.55232 dB, SBA = 23.3660, Width = 4.6896
  BiOrth = 0.004968, Subset1_Cg = 9.62133 */
#define OD_FILTER_PARAMS4_0 (85)
#define OD_FILTER_PARAMS4_1 (75)
#define OD_FILTER_PARAMS4_2 (-15)
#define OD_FILTER_PARAMS4_3 (33)

const int OD_FILTER_PARAMS4[4] = {
  OD_FILTER_PARAMS4_0, OD_FILTER_PARAMS4_1, OD_FILTER_PARAMS4_2,
  OD_FILTER_PARAMS4_3
};

void od_pre_filter4(od_coeff _y[4], const od_coeff _x[4]) {
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 4; i++) {
    _y[i] = _x[i];
  }
#else
  int t[4];
  /*+1/-1 butterflies (required for FIR, PR, LP).*/
  t[3] = _x[0]-_x[3];
  t[2] = _x[1]-_x[2];
  t[1] = _x[1]-(t[2]>>1);
  t[0] = _x[0]-(t[3]>>1);
  /*U filter (arbitrary invertible, omitted).*/
  /*V filter (arbitrary invertible).*/
  /*Scaling factors: the biorthogonal part.*/
  /*Note: t[i]+=t[i]>>(OD_COEFF_BITS-1)&1 is equivalent to: if(t[i]>0)t[i]++
    This step ensures that the scaling is trivially invertible on the decoder's
     side, with perfect reconstruction.*/
  OD_DCT_OVERFLOW_CHECK(t[2], OD_FILTER_PARAMS4_0, 32, 51);
  /*s0*/
# if OD_FILTER_PARAMS4_0 != 64
  t[2] = t[2]*OD_FILTER_PARAMS4_0>>6;
  t[2] += -t[2]>>(OD_COEFF_BITS-1)&1;
# endif
  OD_DCT_OVERFLOW_CHECK(t[3], OD_FILTER_PARAMS4_1, 32, 52);
  /*s1*/
# if OD_FILTER_PARAMS4_1 != 64
  t[3] = t[3]*OD_FILTER_PARAMS4_1>>6;
  t[3] += -t[3]>>(OD_COEFF_BITS-1)&1;
# endif
  /*Rotation:*/
  OD_DCT_OVERFLOW_CHECK(t[2], OD_FILTER_PARAMS4_2, 32, 53);
  /*p0*/
  t[3] += (t[2]*OD_FILTER_PARAMS4_2+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[3], OD_FILTER_PARAMS4_3, 32, 86);
  /*u0*/
  t[2] += (t[3]*OD_FILTER_PARAMS4_3+32)>>6;
  /*More +1/-1 butterflies (required for FIR, PR, LP).*/
  t[0] += t[3]>>1;
  _y[0] = (od_coeff)t[0];
  t[1] += t[2]>>1;
  _y[1] = (od_coeff)t[1];
  _y[2] = (od_coeff)(t[1]-t[2]);
  _y[3] = (od_coeff)(t[0]-t[3]);
#endif
}

void od_post_filter4(od_coeff _x[4], const od_coeff _y[4]) {
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 4; i++) {
    _x[i] = _y[i];
  }
#else
  int t[4];
  t[3] = _y[0]-_y[3];
  t[2] = _y[1]-_y[2];
  t[1] = _y[1]-(t[2]>>1);
  t[0] = _y[0]-(t[3]>>1);
  t[2] -= (t[3]*OD_FILTER_PARAMS4_3+32)>>6;
  t[3] -= (t[2]*OD_FILTER_PARAMS4_2+32)>>6;
#if OD_FILTER_PARAMS4_1 != 64
  t[3] = (t[3]*(1 << 6))/OD_FILTER_PARAMS4_1;
#endif
#if OD_FILTER_PARAMS4_0 != 64
  t[2] = (t[2]*(1 << 6))/OD_FILTER_PARAMS4_0;
#endif
  t[0] += t[3]>>1;
  _x[0] = (od_coeff)t[0];
  t[1] += t[2]>>1;
  _x[1] = (od_coeff)t[1];
  _x[2] = (od_coeff)(t[1]-t[2]);
  _x[3] = (od_coeff)(t[0]-t[3]);
#endif
}

/*R=f
  6-bit
  Subset3_2d_Cg = 16.7948051638528391
  This filter has some of the scale factors
   forced to one to reduce complexity. Without
   this constraint we get a filter with Subset3_2d_Cg
   of 16.8035257369844686. The small cg loss is likely
   worth the reduction in multiplies and adds.*/
#if 0
# define OD_FILTER8_TYPE3 (1)
# define OD_FILTER_PARAMS8_0 (86)
# define OD_FILTER_PARAMS8_1 (64)
# define OD_FILTER_PARAMS8_2 (64)
# define OD_FILTER_PARAMS8_3 (70)
# define OD_FILTER_PARAMS8_4 (-29)
# define OD_FILTER_PARAMS8_5 (-25)
# define OD_FILTER_PARAMS8_6 (-13)
# define OD_FILTER_PARAMS8_7 (45)
# define OD_FILTER_PARAMS8_8 (29)
# define OD_FILTER_PARAMS8_9 (17)
#elif 1
/*Optimal 1d subset3 Cg*/
# define OD_FILTER8_TYPE3 (1)
# define OD_FILTER_PARAMS8_0 (93)
# define OD_FILTER_PARAMS8_1 (72)
# define OD_FILTER_PARAMS8_2 (73)
# define OD_FILTER_PARAMS8_3 (78)
# define OD_FILTER_PARAMS8_4 (-28)
# define OD_FILTER_PARAMS8_5 (-23)
# define OD_FILTER_PARAMS8_6 (-10)
# define OD_FILTER_PARAMS8_7 (50)
# define OD_FILTER_PARAMS8_8 (37)
# define OD_FILTER_PARAMS8_9 (23)
#else
/*1,2 regular*/
# define OD_FILTER8_TYPE3 (0)
# define OD_FILTER_PARAMS8_0 (88)
# define OD_FILTER_PARAMS8_1 (75)
# define OD_FILTER_PARAMS8_2 (76)
# define OD_FILTER_PARAMS8_3 (76)
# define OD_FILTER_PARAMS8_4 (-24)
# define OD_FILTER_PARAMS8_5 (-20)
# define OD_FILTER_PARAMS8_6 (-4)
# define OD_FILTER_PARAMS8_7 (53)
# define OD_FILTER_PARAMS8_8 (40)
# define OD_FILTER_PARAMS8_9 (24)
#endif

const int OD_FILTER_PARAMS8[10] = {
  OD_FILTER_PARAMS8_0, OD_FILTER_PARAMS8_1, OD_FILTER_PARAMS8_2,
  OD_FILTER_PARAMS8_3, OD_FILTER_PARAMS8_4, OD_FILTER_PARAMS8_5,
  OD_FILTER_PARAMS8_6, OD_FILTER_PARAMS8_7, OD_FILTER_PARAMS8_8,
  OD_FILTER_PARAMS8_9
};

void od_pre_filter8(od_coeff _y[8], const od_coeff _x[8]) {
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 8; i++) {
    _y[i] = _x[i];
  }
#else
  int t[8];
  /*+1/-1 butterflies (required for FIR, PR, LP).*/
  t[7] = _x[0]-_x[7];
  t[6] = _x[1]-_x[6];
  t[5] = _x[2]-_x[5];
  t[4] = _x[3]-_x[4];
  t[3] = _x[3]-(t[4]>>1);
  t[2] = _x[2]-(t[5]>>1);
  t[1] = _x[1]-(t[6]>>1);
  t[0] = _x[0]-(t[7]>>1);
  /*U filter (arbitrary invertible, omitted).*/
  /*V filter (arbitrary invertible, can be optimized for speed or accuracy).*/
  /*Scaling factors: the biorthogonal part.*/
  /*Note: t[i]+=t[i]>>(OD_COEFF_BITS-1)&1; is equivalent to: if(t[i]>0)t[i]++;
    This step ensures that the scaling is trivially invertible on the
     decoder's side, with perfect reconstruction.*/
#if OD_FILTER_PARAMS8_0 != 64
  OD_DCT_OVERFLOW_CHECK(t[4], OD_FILTER_PARAMS8_0, 0, 54);
  t[4] = (t[4]*OD_FILTER_PARAMS8_0)>>6;
  t[4] += -t[4]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS8_1 != 64
  OD_DCT_OVERFLOW_CHECK(t[5], OD_FILTER_PARAMS8_1, 0, 55);
  t[5] = (t[5]*OD_FILTER_PARAMS8_1)>>6;
  t[5] += -t[5]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS8_2 != 64
  OD_DCT_OVERFLOW_CHECK(t[6], OD_FILTER_PARAMS8_2, 0, 56);
  t[6] = (t[6]*OD_FILTER_PARAMS8_2)>>6;
  t[6] += -t[6]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS8_3 != 64
  OD_DCT_OVERFLOW_CHECK(t[7], OD_FILTER_PARAMS8_3, 0, 57);
  t[7] = (t[7]*OD_FILTER_PARAMS8_3)>>6;
  t[7] += -t[7]>>(OD_COEFF_BITS-1)&1;
#endif
  /*Rotations:*/
#if OD_FILTER8_TYPE3
  OD_DCT_OVERFLOW_CHECK(t[6], OD_FILTER_PARAMS8_6, 32, 58);
  t[7] += (t[6]*OD_FILTER_PARAMS8_6+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[7], OD_FILTER_PARAMS8_9, 32, 59);
  t[6] += (t[7]*OD_FILTER_PARAMS8_9+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[5], OD_FILTER_PARAMS8_5, 32, 60);
  t[6] += (t[5]*OD_FILTER_PARAMS8_5+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[6], OD_FILTER_PARAMS8_8, 32, 61);
  t[5] += (t[6]*OD_FILTER_PARAMS8_8+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[4], OD_FILTER_PARAMS8_4, 32, 62);
  t[5] += (t[4]*OD_FILTER_PARAMS8_4+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[5], OD_FILTER_PARAMS8_7, 32, 63);
  t[4] += (t[5]*OD_FILTER_PARAMS8_7+32)>>6;
#else
  OD_DCT_OVERFLOW_CHECK(t[4], OD_FILTER_PARAMS8_4, 32, 58);
  t[5] += (t[4]*OD_FILTER_PARAMS8_4+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[5], OD_FILTER_PARAMS8_5, 32, 59);
  t[6] += (t[5]*OD_FILTER_PARAMS8_5+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[6], OD_FILTER_PARAMS8_6, 32, 60);
  t[7] += (t[6]*OD_FILTER_PARAMS8_6+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[7], OD_FILTER_PARAMS8_7, 32, 61);
  t[6] += (t[7]*OD_FILTER_PARAMS8_9+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[6], OD_FILTER_PARAMS8_8, 32, 62);
  t[5] += (t[6]*OD_FILTER_PARAMS8_8+32)>>6;
  OD_DCT_OVERFLOW_CHECK(t[5], OD_FILTER_PARAMS8_9, 32, 63);
  t[4] += (t[5]*OD_FILTER_PARAMS8_7+32)>>6;
#endif
  /*More +1/-1 butterflies (required for FIR, PR, LP).*/
  t[0] += t[7]>>1;
  _y[0] = (od_coeff)t[0];
  t[1] += t[6]>>1;
  _y[1] = (od_coeff)t[1];
  t[2] += t[5]>>1;
  _y[2] = (od_coeff)t[2];
  t[3] += t[4]>>1;
  _y[3] = (od_coeff)t[3];
  _y[4] = (od_coeff)(t[3]-t[4]);
  _y[5] = (od_coeff)(t[2]-t[5]);
  _y[6] = (od_coeff)(t[1]-t[6]);
  _y[7] = (od_coeff)(t[0]-t[7]);
#endif
}

void od_post_filter8(od_coeff _x[8], const od_coeff _y[8]) {
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 8; i++) {
    _x[i] = _y[i];
  }
#else
  int t[8];
  t[7] = _y[0]-_y[7];
  t[6] = _y[1]-_y[6];
  t[5] = _y[2]-_y[5];
  t[4] = _y[3]-_y[4];
  t[3] = _y[3]-(t[4]>>1);
  t[2] = _y[2]-(t[5]>>1);
  t[1] = _y[1]-(t[6]>>1);
  t[0] = _y[0]-(t[7]>>1);
#if OD_FILTER8_TYPE3
  t[4] -= (t[5]*OD_FILTER_PARAMS8_7+32)>>6;
  t[5] -= (t[4]*OD_FILTER_PARAMS8_4+32)>>6;
  t[5] -= (t[6]*OD_FILTER_PARAMS8_8+32)>>6;
  t[6] -= (t[5]*OD_FILTER_PARAMS8_5+32)>>6;
  t[6] -= (t[7]*OD_FILTER_PARAMS8_9+32)>>6;
  t[7] -= (t[6]*OD_FILTER_PARAMS8_6+32)>>6;
#else
  t[4] -= (t[5]*OD_FILTER_PARAMS8_7+32)>>6;
  t[5] -= (t[6]*OD_FILTER_PARAMS8_8+32)>>6;
  t[6] -= (t[7]*OD_FILTER_PARAMS8_9+32)>>6;
  t[7] -= (t[6]*OD_FILTER_PARAMS8_6+32)>>6;
  t[6] -= (t[5]*OD_FILTER_PARAMS8_5+32)>>6;
  t[5] -= (t[4]*OD_FILTER_PARAMS8_4+32)>>6;
#endif
#if OD_FILTER_PARAMS8_3 != 64
  t[7] = (t[7]<<6)/OD_FILTER_PARAMS8_3;
#endif
#if OD_FILTER_PARAMS8_2 != 64
  t[6] = (t[6]<<6)/OD_FILTER_PARAMS8_2;
#endif
#if OD_FILTER_PARAMS8_1 != 64
  t[5] = (t[5]<<6)/OD_FILTER_PARAMS8_1;
#endif
#if OD_FILTER_PARAMS8_0 != 64
  t[4] = (t[4]<<6)/OD_FILTER_PARAMS8_0;
#endif
  t[0] += t[7]>>1;
  _x[0] = (od_coeff)t[0];
  t[1] += t[6]>>1;
  _x[1] = (od_coeff)t[1];
  t[2] += t[5]>>1;
  _x[2] = (od_coeff)t[2];
  t[3] += t[4]>>1;
  _x[3] = (od_coeff)t[3];
  _x[4] = (od_coeff)(t[3]-t[4]);
  _x[5] = (od_coeff)(t[2]-t[5]);
  _x[6] = (od_coeff)(t[1]-t[6]);
  _x[7] = (od_coeff)(t[0]-t[7]);
#endif
}

/*R=f Type-3
  6-bit
  Subset3_2d_Cg = 17.0007008994465814
  This filter has some of the scale factors
   forced to 1, without this we get a filter
   with a Cg of 17.0010659559160686. The tiny
   CG loss is likely worth the reduction in
   multiplies and adds.*/
#if 0
# define OD_FILTER16_TYPE3 (1)
# define OD_FILTER_PARAMS16_0 (90)
# define OD_FILTER_PARAMS16_1 (67)
# define OD_FILTER_PARAMS16_2 (64)
# define OD_FILTER_PARAMS16_3 (64)
# define OD_FILTER_PARAMS16_4 (64)
# define OD_FILTER_PARAMS16_5 (64)
# define OD_FILTER_PARAMS16_6 (66)
# define OD_FILTER_PARAMS16_7 (67)
# define OD_FILTER_PARAMS16_8 (32)
# define OD_FILTER_PARAMS16_9 (37)
# define OD_FILTER_PARAMS16_10 (35)
# define OD_FILTER_PARAMS16_11 (30)
# define OD_FILTER_PARAMS16_12 (23)
# define OD_FILTER_PARAMS16_13 (16)
# define OD_FILTER_PARAMS16_14 (7)
# define OD_FILTER_PARAMS16_15 (52)
# define OD_FILTER_PARAMS16_16 (42)
# define OD_FILTER_PARAMS16_17 (36)
# define OD_FILTER_PARAMS16_18 (31)
# define OD_FILTER_PARAMS16_19 (25)
# define OD_FILTER_PARAMS16_20 (18)
# define OD_FILTER_PARAMS16_21 (9)
#elif 1
/*Optimal 1D subset3 Cg*/
# define OD_FILTER16_TYPE3 (1)
# define OD_FILTER_PARAMS16_0 (94)
# define OD_FILTER_PARAMS16_1 (71)
# define OD_FILTER_PARAMS16_2 (68)
# define OD_FILTER_PARAMS16_3 (68)
# define OD_FILTER_PARAMS16_4 (68)
# define OD_FILTER_PARAMS16_5 (69)
# define OD_FILTER_PARAMS16_6 (70)
# define OD_FILTER_PARAMS16_7 (73)
# define OD_FILTER_PARAMS16_8 (-32)
# define OD_FILTER_PARAMS16_9 (-37)
# define OD_FILTER_PARAMS16_10 (-36)
# define OD_FILTER_PARAMS16_11 (-32)
# define OD_FILTER_PARAMS16_12 (-26)
# define OD_FILTER_PARAMS16_13 (-17)
# define OD_FILTER_PARAMS16_14 (-7)
# define OD_FILTER_PARAMS16_15 (56)
# define OD_FILTER_PARAMS16_16 (49)
# define OD_FILTER_PARAMS16_17 (45)
# define OD_FILTER_PARAMS16_18 (40)
# define OD_FILTER_PARAMS16_19 (34)
# define OD_FILTER_PARAMS16_20 (26)
# define OD_FILTER_PARAMS16_21 (15)
#else
/*1,2-regular*/
# define OD_FILTER16_TYPE3 (0)
# define OD_FILTER_PARAMS16_0 (80)
# define OD_FILTER_PARAMS16_1 (72)
# define OD_FILTER_PARAMS16_2 (73)
# define OD_FILTER_PARAMS16_3 (68)
# define OD_FILTER_PARAMS16_4 (72)
# define OD_FILTER_PARAMS16_5 (74)
# define OD_FILTER_PARAMS16_6 (70)
# define OD_FILTER_PARAMS16_7 (73)
# define OD_FILTER_PARAMS16_8 (-32)
# define OD_FILTER_PARAMS16_9 (-28)
# define OD_FILTER_PARAMS16_10 (-24)
# define OD_FILTER_PARAMS16_11 (-32)
# define OD_FILTER_PARAMS16_12 (-24)
# define OD_FILTER_PARAMS16_13 (-13)
# define OD_FILTER_PARAMS16_14 (-2)
# define OD_FILTER_PARAMS16_15 (59)
# define OD_FILTER_PARAMS16_16 (53)
# define OD_FILTER_PARAMS16_17 (46)
# define OD_FILTER_PARAMS16_18 (41)
# define OD_FILTER_PARAMS16_19 (35)
# define OD_FILTER_PARAMS16_20 (24)
# define OD_FILTER_PARAMS16_21 (12)
#endif

const int OD_FILTER_PARAMS16[22]={
  OD_FILTER_PARAMS16_0,OD_FILTER_PARAMS16_1,OD_FILTER_PARAMS16_2,
  OD_FILTER_PARAMS16_3,OD_FILTER_PARAMS16_4,OD_FILTER_PARAMS16_5,
  OD_FILTER_PARAMS16_6,OD_FILTER_PARAMS16_7,OD_FILTER_PARAMS16_8,
  OD_FILTER_PARAMS16_9,OD_FILTER_PARAMS16_10,OD_FILTER_PARAMS16_11,
  OD_FILTER_PARAMS16_12,OD_FILTER_PARAMS16_13,OD_FILTER_PARAMS16_14,
  OD_FILTER_PARAMS16_15,OD_FILTER_PARAMS16_16,OD_FILTER_PARAMS16_17,
  OD_FILTER_PARAMS16_18,OD_FILTER_PARAMS16_19,OD_FILTER_PARAMS16_20,
  OD_FILTER_PARAMS16_21
};

void od_pre_filter16(od_coeff _y[16],const od_coeff _x[16]){
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 16; i++) {
    _y[i] = _x[i];
  }
#else
   int t[16];
   /*+1/-1 butterflies (required for FIR, PR, LP).*/
   t[15]=_x[0]-_x[15];
   t[14]=_x[1]-_x[14];
   t[13]=_x[2]-_x[13];
   t[12]=_x[3]-_x[12];
   t[11]=_x[4]-_x[11];
   t[10]=_x[5]-_x[10];
   t[9]=_x[6]-_x[9];
   t[8]=_x[7]-_x[8];
   t[7]=_x[7]-(t[8]>>1);
   t[6]=_x[6]-(t[9]>>1);
   t[5]=_x[5]-(t[10]>>1);
   t[4]=_x[4]-(t[11]>>1);
   t[3]=_x[3]-(t[12]>>1);
   t[2]=_x[2]-(t[13]>>1);
   t[1]=_x[1]-(t[14]>>1);
   t[0]=_x[0]-(t[15]>>1);
   /*U filter (arbitrary invertible, omitted).*/
   /*V filter (arbitrary invertible, can be optimized for speed or accuracy).*/
   /*Scaling factors: the biorthogonal part.*/
   /*Note: t[i]+=t[i]>>(OD_COEFF_BITS-1)&1; is equivalent to: if(t[i]>0)t[i]++;
     This step ensures that the scaling is trivially invertible on the
      decoder's side, with perfect reconstruction.*/
#if OD_FILTER_PARAMS16_0!=64
   OD_DCT_OVERFLOW_CHECK(t[8], OD_FILTER_PARAMS16_0, 0, 64);
   t[8]=(t[8]*OD_FILTER_PARAMS16_0)>>6;
   t[8]+=-t[8]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_1!=64
   OD_DCT_OVERFLOW_CHECK(t[9], OD_FILTER_PARAMS16_1, 0, 65);
   t[9]=(t[9]*OD_FILTER_PARAMS16_1)>>6;
   t[9]+=-t[9]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_2!=64
   OD_DCT_OVERFLOW_CHECK(t[10], OD_FILTER_PARAMS16_2, 0, 66);
   t[10]=(t[10]*OD_FILTER_PARAMS16_2)>>6;
   t[10]+=-t[10]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_3!=64
   OD_DCT_OVERFLOW_CHECK(t[11], OD_FILTER_PARAMS16_3, 0, 67);
   t[11]=(t[11]*OD_FILTER_PARAMS16_3)>>6;
   t[11]+=-t[11]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_4!=64
   OD_DCT_OVERFLOW_CHECK(t[12], OD_FILTER_PARAMS16_4, 0, 68);
   t[12]=(t[12]*OD_FILTER_PARAMS16_4)>>6;
   t[12]+=-t[12]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_5!=64
   OD_DCT_OVERFLOW_CHECK(t[13], OD_FILTER_PARAMS16_5, 0, 69);
   t[13]=(t[13]*OD_FILTER_PARAMS16_5)>>6;
   t[13]+=-t[13]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_6!=64
   OD_DCT_OVERFLOW_CHECK(t[14], OD_FILTER_PARAMS16_6, 0, 70);
   t[14]=t[14]*OD_FILTER_PARAMS16_6>>6;
   t[14]+=-t[14]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_7!=64
   OD_DCT_OVERFLOW_CHECK(t[15], OD_FILTER_PARAMS16_7, 0, 71);
   t[15]=t[15]*OD_FILTER_PARAMS16_7>>6;
   t[15]+=-t[15]>>(OD_COEFF_BITS-1)&1;
#endif
   /*Rotations:*/
#if OD_FILTER16_TYPE3
   OD_DCT_OVERFLOW_CHECK(t[14], OD_FILTER_PARAMS16_14, 32, 72);
   t[15]+=(t[14]*OD_FILTER_PARAMS16_14+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[15], OD_FILTER_PARAMS16_21, 32, 73);
   t[14]+=(t[15]*OD_FILTER_PARAMS16_21+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[13], OD_FILTER_PARAMS16_13, 32, 74);
   t[14]+=(t[13]*OD_FILTER_PARAMS16_13+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[14], OD_FILTER_PARAMS16_20, 32, 75);
   t[13]+=(t[14]*OD_FILTER_PARAMS16_20+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[12], OD_FILTER_PARAMS16_12, 32, 76);
   t[13]+=(t[12]*OD_FILTER_PARAMS16_12+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[13], OD_FILTER_PARAMS16_19, 32, 77);
   t[12]+=(t[13]*OD_FILTER_PARAMS16_19+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[11], OD_FILTER_PARAMS16_11, 32, 78);
   t[12]+=(t[11]*OD_FILTER_PARAMS16_11+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[12], OD_FILTER_PARAMS16_18, 32, 79);
   t[11]+=(t[12]*OD_FILTER_PARAMS16_18+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[10], OD_FILTER_PARAMS16_10, 32, 80);
   t[11]+=(t[10]*OD_FILTER_PARAMS16_10+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[11], OD_FILTER_PARAMS16_17, 32, 81);
   t[10]+=(t[11]*OD_FILTER_PARAMS16_17+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[9], OD_FILTER_PARAMS16_9, 32, 82);
   t[10]+=(t[9]*OD_FILTER_PARAMS16_9+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[10], OD_FILTER_PARAMS16_16, 32, 83);
   t[9]+=(t[10]*OD_FILTER_PARAMS16_16+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[8], OD_FILTER_PARAMS16_8, 32, 84);
   t[9]+=(t[8]*OD_FILTER_PARAMS16_8+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[9], OD_FILTER_PARAMS16_15, 32, 85);
   t[8]+=(t[9]*OD_FILTER_PARAMS16_15+32)>>6;
#else
   OD_DCT_OVERFLOW_CHECK(t[8], OD_FILTER_PARAMS16_8, 32, 72);
   t[9]+=(t[8]*OD_FILTER_PARAMS16_8+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[9], OD_FILTER_PARAMS16_9, 32, 73);
   t[10]+=(t[9]*OD_FILTER_PARAMS16_9+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[10], OD_FILTER_PARAMS16_10, 32, 74);
   t[11]+=(t[10]*OD_FILTER_PARAMS16_10+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[11], OD_FILTER_PARAMS16_11, 32, 75);
   t[12]+=(t[11]*OD_FILTER_PARAMS16_11+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[12], OD_FILTER_PARAMS16_12, 32, 76);
   t[13]+=(t[12]*OD_FILTER_PARAMS16_12+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[13], OD_FILTER_PARAMS16_13, 32, 77);
   t[14]+=(t[13]*OD_FILTER_PARAMS16_13+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[14], OD_FILTER_PARAMS16_14, 32, 78);
   t[15]+=(t[14]*OD_FILTER_PARAMS16_14+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[15], OD_FILTER_PARAMS16_21, 32, 79);
   t[14]+=(t[15]*OD_FILTER_PARAMS16_21+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[14], OD_FILTER_PARAMS16_20, 32, 80);
   t[13]+=(t[14]*OD_FILTER_PARAMS16_20+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[13], OD_FILTER_PARAMS16_19, 32, 81);
   t[12]+=(t[13]*OD_FILTER_PARAMS16_19+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[12], OD_FILTER_PARAMS16_18, 32, 82);
   t[11]+=(t[12]*OD_FILTER_PARAMS16_18+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[11], OD_FILTER_PARAMS16_17, 32, 83);
   t[10]+=(t[11]*OD_FILTER_PARAMS16_17+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[10], OD_FILTER_PARAMS16_16, 32, 84);
   t[9]+=(t[10]*OD_FILTER_PARAMS16_16+32)>>6;
   OD_DCT_OVERFLOW_CHECK(t[9], OD_FILTER_PARAMS16_15, 32, 85);
   t[8]+=(t[9]*OD_FILTER_PARAMS16_15+32)>>6;
#endif
   /*More +1/-1 butterflies (required for FIR, PR, LP).*/
   t[0]+=t[15]>>1;
   _y[0]=(od_coeff)t[0];
   t[1]+=t[14]>>1;
   _y[1]=(od_coeff)t[1];
   t[2]+=t[13]>>1;
   _y[2]=(od_coeff)t[2];
   t[3]+=t[12]>>1;
   _y[3]=(od_coeff)t[3];
   t[4]+=t[11]>>1;
   _y[4]=(od_coeff)t[4];
   t[5]+=t[10]>>1;
   _y[5]=(od_coeff)t[5];
   t[6]+=t[9]>>1;
   _y[6]=(od_coeff)t[6];
   t[7]+=t[8]>>1;
   _y[7]=(od_coeff)t[7];
   _y[8]=(od_coeff)(t[7]-t[8]);
   _y[9]=(od_coeff)(t[6]-t[9]);
   _y[10]=(od_coeff)(t[5]-t[10]);
   _y[11]=(od_coeff)(t[4]-t[11]);
   _y[12]=(od_coeff)(t[3]-t[12]);
   _y[13]=(od_coeff)(t[2]-t[13]);
   _y[14]=(od_coeff)(t[1]-t[14]);
   _y[15]=(od_coeff)(t[0]-t[15]);
#endif
}

void od_post_filter16(od_coeff _x[16],const od_coeff _y[16]){
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 16; i++) {
    _x[i] = _y[i];
  }
#else
   int t[16];
   t[15]=_y[0]-_y[15];
   t[14]=_y[1]-_y[14];
   t[13]=_y[2]-_y[13];
   t[12]=_y[3]-_y[12];
   t[11]=_y[4]-_y[11];
   t[10]=_y[5]-_y[10];
   t[9]=_y[6]-_y[9];
   t[8]=_y[7]-_y[8];
   t[7]=_y[7]-(t[8]>>1);
   t[6]=_y[6]-(t[9]>>1);
   t[5]=_y[5]-(t[10]>>1);
   t[4]=_y[4]-(t[11]>>1);
   t[3]=_y[3]-(t[12]>>1);
   t[2]=_y[2]-(t[13]>>1);
   t[1]=_y[1]-(t[14]>>1);
   t[0]=_y[0]-(t[15]>>1);
#if OD_FILTER16_TYPE3
   t[8]-=(t[9]*OD_FILTER_PARAMS16_15+32)>>6;
   t[9]-=(t[8]*OD_FILTER_PARAMS16_8+32)>>6;
   t[9]-=(t[10]*OD_FILTER_PARAMS16_16+32)>>6;
   t[10]-=(t[9]*OD_FILTER_PARAMS16_9+32)>>6;
   t[10]-=(t[11]*OD_FILTER_PARAMS16_17+32)>>6;
   t[11]-=(t[10]*OD_FILTER_PARAMS16_10+32)>>6;
   t[11]-=(t[12]*OD_FILTER_PARAMS16_18+32)>>6;
   t[12]-=(t[11]*OD_FILTER_PARAMS16_11+32)>>6;
   t[12]-=(t[13]*OD_FILTER_PARAMS16_19+32)>>6;
   t[13]-=(t[12]*OD_FILTER_PARAMS16_12+32)>>6;
   t[13]-=(t[14]*OD_FILTER_PARAMS16_20+32)>>6;
   t[14]-=(t[13]*OD_FILTER_PARAMS16_13+32)>>6;
   t[14]-=(t[15]*OD_FILTER_PARAMS16_21+32)>>6;
   t[15]-=(t[14]*OD_FILTER_PARAMS16_14+32)>>6;
#else
   t[8]-=(t[9]*OD_FILTER_PARAMS16_15+32)>>6;
   t[9]-=(t[10]*OD_FILTER_PARAMS16_16+32)>>6;
   t[10]-=(t[11]*OD_FILTER_PARAMS16_17+32)>>6;
   t[11]-=(t[12]*OD_FILTER_PARAMS16_18+32)>>6;
   t[12]-=(t[13]*OD_FILTER_PARAMS16_19+32)>>6;
   t[13]-=(t[14]*OD_FILTER_PARAMS16_20+32)>>6;
   t[14]-=(t[15]*OD_FILTER_PARAMS16_21+32)>>6;
   t[15]-=(t[14]*OD_FILTER_PARAMS16_14+32)>>6;
   t[14]-=(t[13]*OD_FILTER_PARAMS16_13+32)>>6;
   t[13]-=(t[12]*OD_FILTER_PARAMS16_12+32)>>6;
   t[12]-=(t[11]*OD_FILTER_PARAMS16_11+32)>>6;
   t[11]-=(t[10]*OD_FILTER_PARAMS16_10+32)>>6;
   t[10]-=(t[9]*OD_FILTER_PARAMS16_9+32)>>6;
   t[9]-=(t[8]*OD_FILTER_PARAMS16_8+32)>>6;
#endif
#if OD_FILTER_PARAMS16_7!=64
   t[15]=(t[15]<<6)/OD_FILTER_PARAMS16_7;
#endif
#if OD_FILTER_PARAMS16_6!=64
   t[14]=(t[14]<<6)/OD_FILTER_PARAMS16_6;
#endif
#if OD_FILTER_PARAMS16_5!=64
   t[13]=(t[13]<<6)/OD_FILTER_PARAMS16_5;
#endif
#if OD_FILTER_PARAMS16_4!=64
   t[12]=(t[12]<<6)/OD_FILTER_PARAMS16_4;
#endif
#if OD_FILTER_PARAMS16_3!=64
   t[11]=(t[11]<<6)/OD_FILTER_PARAMS16_3;
#endif
#if OD_FILTER_PARAMS16_2!=64
   t[10]=(t[10]<<6)/OD_FILTER_PARAMS16_2;
#endif
#if OD_FILTER_PARAMS16_1!=64
   t[9]=(t[9]<<6)/OD_FILTER_PARAMS16_1;
#endif
#if OD_FILTER_PARAMS16_0!=64
   t[8]=(t[8]<<6)/OD_FILTER_PARAMS16_0;
#endif
   t[0]+=t[15]>>1;
   _x[0]=(od_coeff)t[0];
   t[1]+=t[14]>>1;
   _x[1]=(od_coeff)t[1];
   t[2]+=t[13]>>1;
   _x[2]=(od_coeff)t[2];
   t[3]+=t[12]>>1;
   _x[3]=(od_coeff)t[3];
   t[4]+=t[11]>>1;
   _x[4]=(od_coeff)t[4];
   t[5]+=t[10]>>1;
   _x[5]=(od_coeff)t[5];
   t[6]+=t[9]>>1;
   _x[6]=(od_coeff)t[6];
   t[7]+=t[8]>>1;
   _x[7]=(od_coeff)t[7];
   _x[8]=(od_coeff)(t[7]-t[8]);
   _x[9]=(od_coeff)(t[6]-t[9]);
   _x[10]=(od_coeff)(t[5]-t[10]);
   _x[11]=(od_coeff)(t[4]-t[11]);
   _x[12]=(od_coeff)(t[3]-t[12]);
   _x[13]=(od_coeff)(t[2]-t[13]);
   _x[14]=(od_coeff)(t[1]-t[14]);
   _x[15]=(od_coeff)(t[0]-t[15]);
#endif
}

/* "Optimal" 1D AR95 cg = 10.0537088081722175 */
# define OD_FILTER32_TYPE3 (1)
# define OD_FILTER_PARAMS32_0 (91)
# define OD_FILTER_PARAMS32_1 (70)
# define OD_FILTER_PARAMS32_2 (68)
# define OD_FILTER_PARAMS32_3 (67)
# define OD_FILTER_PARAMS32_4 (67)
# define OD_FILTER_PARAMS32_5 (67)
# define OD_FILTER_PARAMS32_6 (67)
# define OD_FILTER_PARAMS32_7 (66)
# define OD_FILTER_PARAMS32_8 (66)
# define OD_FILTER_PARAMS32_9 (67)
# define OD_FILTER_PARAMS32_10 (67)
# define OD_FILTER_PARAMS32_11 (66)
# define OD_FILTER_PARAMS32_12 (67)
# define OD_FILTER_PARAMS32_13 (67)
# define OD_FILTER_PARAMS32_14 (67)
# define OD_FILTER_PARAMS32_15 (70)
# define OD_FILTER_PARAMS32_16 (-32)
# define OD_FILTER_PARAMS32_17 (-41)
# define OD_FILTER_PARAMS32_18 (-42)
# define OD_FILTER_PARAMS32_19 (-41)
# define OD_FILTER_PARAMS32_20 (-40)
# define OD_FILTER_PARAMS32_21 (-38)
# define OD_FILTER_PARAMS32_22 (-36)
# define OD_FILTER_PARAMS32_23 (-34)
# define OD_FILTER_PARAMS32_24 (-32)
# define OD_FILTER_PARAMS32_25 (-29)
# define OD_FILTER_PARAMS32_26 (-24)
# define OD_FILTER_PARAMS32_27 (-19)
# define OD_FILTER_PARAMS32_28 (-14)
# define OD_FILTER_PARAMS32_29 (-9)
# define OD_FILTER_PARAMS32_30 (-5)
# define OD_FILTER_PARAMS32_31 (58)
# define OD_FILTER_PARAMS32_32 (52)
# define OD_FILTER_PARAMS32_33 (50)
# define OD_FILTER_PARAMS32_34 (48)
# define OD_FILTER_PARAMS32_35 (45)
# define OD_FILTER_PARAMS32_36 (43)
# define OD_FILTER_PARAMS32_37 (40)
# define OD_FILTER_PARAMS32_38 (38)
# define OD_FILTER_PARAMS32_39 (35)
# define OD_FILTER_PARAMS32_40 (32)
# define OD_FILTER_PARAMS32_41 (29)
# define OD_FILTER_PARAMS32_42 (24)
# define OD_FILTER_PARAMS32_43 (18)
# define OD_FILTER_PARAMS32_44 (13)
# define OD_FILTER_PARAMS32_45 (8)

const int OD_FILTER_PARAMS32[46] = {
  OD_FILTER_PARAMS32_0, OD_FILTER_PARAMS32_1, OD_FILTER_PARAMS32_2,
  OD_FILTER_PARAMS32_3, OD_FILTER_PARAMS32_4, OD_FILTER_PARAMS32_5,
  OD_FILTER_PARAMS32_6, OD_FILTER_PARAMS32_7, OD_FILTER_PARAMS32_8,
  OD_FILTER_PARAMS32_9, OD_FILTER_PARAMS32_10, OD_FILTER_PARAMS32_11,
  OD_FILTER_PARAMS32_12, OD_FILTER_PARAMS32_13, OD_FILTER_PARAMS32_14,
  OD_FILTER_PARAMS32_15, OD_FILTER_PARAMS32_16, OD_FILTER_PARAMS32_17,
  OD_FILTER_PARAMS32_18, OD_FILTER_PARAMS32_19, OD_FILTER_PARAMS32_20,
  OD_FILTER_PARAMS32_21, OD_FILTER_PARAMS32_22, OD_FILTER_PARAMS32_23,
  OD_FILTER_PARAMS32_24, OD_FILTER_PARAMS32_25, OD_FILTER_PARAMS32_26,
  OD_FILTER_PARAMS32_27, OD_FILTER_PARAMS32_28, OD_FILTER_PARAMS32_29,
  OD_FILTER_PARAMS32_30, OD_FILTER_PARAMS32_31, OD_FILTER_PARAMS32_32,
  OD_FILTER_PARAMS32_33, OD_FILTER_PARAMS32_34, OD_FILTER_PARAMS32_35,
  OD_FILTER_PARAMS32_36, OD_FILTER_PARAMS32_37, OD_FILTER_PARAMS32_38,
  OD_FILTER_PARAMS32_39, OD_FILTER_PARAMS32_40, OD_FILTER_PARAMS32_41,
  OD_FILTER_PARAMS32_42, OD_FILTER_PARAMS32_43, OD_FILTER_PARAMS32_44,
  OD_FILTER_PARAMS32_45
};

void od_pre_filter32(od_coeff _y[32], const od_coeff _x[32]) {
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 32; i++) {
    _y[i] = _x[i];
  }
#else
  int t[32];
  t[31] = _x[0] - _x[31];
  t[30] = _x[1] - _x[30];
  t[29] = _x[2] - _x[29];
  t[28] = _x[3] - _x[28];
  t[27] = _x[4] - _x[27];
  t[26] = _x[5] - _x[26];
  t[25] = _x[6] - _x[25];
  t[24] = _x[7] - _x[24];
  t[23] = _x[8] - _x[23];
  t[22] = _x[9] - _x[22];
  t[21] = _x[10] - _x[21];
  t[20] = _x[11] - _x[20];
  t[19] = _x[12] - _x[19];
  t[18] = _x[13] - _x[18];
  t[17] = _x[14] - _x[17];
  t[16] = _x[15] - _x[16];
  t[15] = _x[15] - (t[16] >> 1);
  t[14] = _x[14] - (t[17] >> 1);
  t[13] = _x[13] - (t[18] >> 1);
  t[12] = _x[12] - (t[19] >> 1);
  t[11] = _x[11] - (t[20] >> 1);
  t[10] = _x[10] - (t[21] >> 1);
  t[9] = _x[9] - (t[22] >> 1);
  t[8] = _x[8] - (t[23] >> 1);
  t[7] = _x[7] - (t[24] >> 1);
  t[6] = _x[6] - (t[25] >> 1);
  t[5] = _x[5] - (t[26] >> 1);
  t[4] = _x[4] - (t[27] >> 1);
  t[3] = _x[3] - (t[28] >> 1);
  t[2] = _x[2] - (t[29] >> 1);
  t[1] = _x[1] - (t[30] >> 1);
  t[0] = _x[0] - (t[31] >> 1);
#if OD_FILTER_PARAMS32_0 != 64
  OD_DCT_OVERFLOW_CHECK(t[16], OD_FILTER_PARAMS32_0, 0, 87);
  t[16] = (t[16]*OD_FILTER_PARAMS32_0) >> 6;
  t[16] += -t[16] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_1 != 64
  OD_DCT_OVERFLOW_CHECK(t[17], OD_FILTER_PARAMS32_1, 0, 88);
  t[17] = (t[17]*OD_FILTER_PARAMS32_1) >> 6;
  t[17] += -t[17] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_2 != 64
  OD_DCT_OVERFLOW_CHECK(t[18], OD_FILTER_PARAMS32_2, 0, 89);
  t[18] = (t[18]*OD_FILTER_PARAMS32_2) >> 6;
  t[18] += -t[18] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_3 != 64
  OD_DCT_OVERFLOW_CHECK(t[19], OD_FILTER_PARAMS32_3, 0, 90);
  t[19] = (t[19]*OD_FILTER_PARAMS32_3) >> 6;
  t[19] += -t[19] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_4 != 64
  OD_DCT_OVERFLOW_CHECK(t[20], OD_FILTER_PARAMS32_4, 0, 91);
  t[20] = (t[20]*OD_FILTER_PARAMS32_4) >> 6;
  t[20] += -t[20] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_5 != 64
  OD_DCT_OVERFLOW_CHECK(t[21], OD_FILTER_PARAMS32_5, 0, 92);
  t[21] = (t[21]*OD_FILTER_PARAMS32_5) >> 6;
  t[21] += -t[21] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_6 != 64
  OD_DCT_OVERFLOW_CHECK(t[22], OD_FILTER_PARAMS32_6, 0, 93);
  t[22] = (t[22]*OD_FILTER_PARAMS32_6) >> 6;
  t[22] += -t[22] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_7 != 64
  OD_DCT_OVERFLOW_CHECK(t[23], OD_FILTER_PARAMS32_7, 0, 94);
  t[23] = (t[23]*OD_FILTER_PARAMS32_7) >> 6;
  t[23] += -t[23] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_8 != 64
  OD_DCT_OVERFLOW_CHECK(t[24], OD_FILTER_PARAMS32_8, 0, 95);
  t[24] = (t[24]*OD_FILTER_PARAMS32_8) >> 6;
  t[24] += -t[24] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_9 != 64
  OD_DCT_OVERFLOW_CHECK(t[25], OD_FILTER_PARAMS32_9, 0, 96);
  t[25] = (t[25]*OD_FILTER_PARAMS32_9) >> 6;
  t[25] += -t[25] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_10 != 64
  OD_DCT_OVERFLOW_CHECK(t[26], OD_FILTER_PARAMS32_10, 0, 97);
  t[26] = (t[26]*OD_FILTER_PARAMS32_10) >> 6;
  t[26] += -t[26] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_11 != 64
  OD_DCT_OVERFLOW_CHECK(t[27], OD_FILTER_PARAMS32_11, 0, 98);
  t[27] = (t[27]*OD_FILTER_PARAMS32_11) >> 6;
  t[27] += -t[27] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_12 != 64
  OD_DCT_OVERFLOW_CHECK(t[28], OD_FILTER_PARAMS32_12, 0, 99);
  t[28] = (t[28]*OD_FILTER_PARAMS32_12) >> 6;
  t[28] += -t[28] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_13 != 64
  OD_DCT_OVERFLOW_CHECK(t[29], OD_FILTER_PARAMS32_13, 0, 100);
  t[29] = (t[29]*OD_FILTER_PARAMS32_13) >> 6;
  t[29] += -t[29] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_14 != 64
  OD_DCT_OVERFLOW_CHECK(t[30], OD_FILTER_PARAMS32_14, 0, 101);
  t[30] = (t[30]*OD_FILTER_PARAMS32_14) >> 6;
  t[30] += -t[30] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER_PARAMS32_15 != 64
  OD_DCT_OVERFLOW_CHECK(t[31], OD_FILTER_PARAMS32_15, 0, 102);
  t[31] = (t[31]*OD_FILTER_PARAMS32_15) >> 6;
  t[31] += -t[31] >> (OD_COEFF_BITS - 1)&1;
#endif
#if OD_FILTER32_TYPE3
  OD_DCT_OVERFLOW_CHECK(t[30], OD_FILTER_PARAMS32_30, 32, 103);
  t[31] += (t[30]*OD_FILTER_PARAMS32_30 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[31], OD_FILTER_PARAMS32_45, 32, 104);
  t[30] += (t[31]*OD_FILTER_PARAMS32_45 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[29], OD_FILTER_PARAMS32_29, 32, 105);
  t[30] += (t[29]*OD_FILTER_PARAMS32_29 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[30], OD_FILTER_PARAMS32_44, 32, 106);
  t[29] += (t[30]*OD_FILTER_PARAMS32_44 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[28], OD_FILTER_PARAMS32_28, 32, 107);
  t[29] += (t[28]*OD_FILTER_PARAMS32_28 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[29], OD_FILTER_PARAMS32_43, 32, 108);
  t[28] += (t[29]*OD_FILTER_PARAMS32_43 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[27], OD_FILTER_PARAMS32_27, 32, 109);
  t[28] += (t[27]*OD_FILTER_PARAMS32_27 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[28], OD_FILTER_PARAMS32_42, 32, 110);
  t[27] += (t[28]*OD_FILTER_PARAMS32_42 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[26], OD_FILTER_PARAMS32_26, 32, 111);
  t[27] += (t[26]*OD_FILTER_PARAMS32_26 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[27], OD_FILTER_PARAMS32_41, 32, 112);
  t[26] += (t[27]*OD_FILTER_PARAMS32_41 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[25], OD_FILTER_PARAMS32_25, 32, 113);
  t[26] += (t[25]*OD_FILTER_PARAMS32_25 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[26], OD_FILTER_PARAMS32_40, 32, 114);
  t[25] += (t[26]*OD_FILTER_PARAMS32_40 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[24], OD_FILTER_PARAMS32_24, 32, 115);
  t[25] += (t[24]*OD_FILTER_PARAMS32_24 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[25], OD_FILTER_PARAMS32_39, 32, 116);
  t[24] += (t[25]*OD_FILTER_PARAMS32_39 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[23], OD_FILTER_PARAMS32_23, 32, 117);
  t[24] += (t[23]*OD_FILTER_PARAMS32_23 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[24], OD_FILTER_PARAMS32_38, 32, 118);
  t[23] += (t[24]*OD_FILTER_PARAMS32_38 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[22], OD_FILTER_PARAMS32_22, 32, 119);
  t[23] += (t[22]*OD_FILTER_PARAMS32_22 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[23], OD_FILTER_PARAMS32_37, 32, 120);
  t[22] += (t[23]*OD_FILTER_PARAMS32_37 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[21], OD_FILTER_PARAMS32_21, 32, 121);
  t[22] += (t[21]*OD_FILTER_PARAMS32_21 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[22], OD_FILTER_PARAMS32_36, 32, 122);
  t[21] += (t[22]*OD_FILTER_PARAMS32_36 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[20], OD_FILTER_PARAMS32_20, 32, 123);
  t[21] += (t[20]*OD_FILTER_PARAMS32_20 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[21], OD_FILTER_PARAMS32_35, 32, 124);
  t[20] += (t[21]*OD_FILTER_PARAMS32_35 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[19], OD_FILTER_PARAMS32_19, 32, 125);
  t[20] += (t[19]*OD_FILTER_PARAMS32_19 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[20], OD_FILTER_PARAMS32_34, 32, 126);
  t[19] += (t[20]*OD_FILTER_PARAMS32_34 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[18], OD_FILTER_PARAMS32_18, 32, 127);
  t[19] += (t[18]*OD_FILTER_PARAMS32_18 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[19], OD_FILTER_PARAMS32_33, 32, 128);
  t[18] += (t[19]*OD_FILTER_PARAMS32_33 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[17], OD_FILTER_PARAMS32_17, 32, 129);
  t[18] += (t[17]*OD_FILTER_PARAMS32_17 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[18], OD_FILTER_PARAMS32_32, 32, 130);
  t[17] += (t[18]*OD_FILTER_PARAMS32_32 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[16], OD_FILTER_PARAMS32_16, 32, 131);
  t[17] += (t[16]*OD_FILTER_PARAMS32_16 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[17], OD_FILTER_PARAMS32_31, 32, 132);
  t[16] += (t[17]*OD_FILTER_PARAMS32_31 + 32) >> 6;
#else
  OD_DCT_OVERFLOW_CHECK(t[16], OD_FILTER_PARAMS32_16, 32, 103);
  t[17] += (t[16]*OD_FILTER_PARAMS32_16 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[17], OD_FILTER_PARAMS32_17, 32, 104);
  t[18] += (t[17]*OD_FILTER_PARAMS32_17 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[18], OD_FILTER_PARAMS32_18, 32, 105);
  t[19] += (t[18]*OD_FILTER_PARAMS32_18 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[19], OD_FILTER_PARAMS32_19, 32, 106);
  t[20] += (t[19]*OD_FILTER_PARAMS32_19 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[20], OD_FILTER_PARAMS32_20, 32, 107);
  t[21] += (t[20]*OD_FILTER_PARAMS32_20 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[21], OD_FILTER_PARAMS32_21, 32, 108);
  t[22] += (t[21]*OD_FILTER_PARAMS32_21 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[22], OD_FILTER_PARAMS32_22, 32, 109);
  t[23] += (t[22]*OD_FILTER_PARAMS32_22 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[23], OD_FILTER_PARAMS32_23, 32, 110);
  t[24] += (t[23]*OD_FILTER_PARAMS32_23 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[24], OD_FILTER_PARAMS32_24, 32, 111);
  t[25] += (t[24]*OD_FILTER_PARAMS32_24 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[25], OD_FILTER_PARAMS32_25, 32, 112);
  t[26] += (t[25]*OD_FILTER_PARAMS32_25 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[26], OD_FILTER_PARAMS32_26, 32, 113);
  t[27] += (t[26]*OD_FILTER_PARAMS32_26 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[27], OD_FILTER_PARAMS32_27, 32, 114);
  t[28] += (t[27]*OD_FILTER_PARAMS32_27 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[28], OD_FILTER_PARAMS32_28, 32, 115);
  t[29] += (t[28]*OD_FILTER_PARAMS32_28 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[29], OD_FILTER_PARAMS32_29, 32, 116);
  t[30] += (t[29]*OD_FILTER_PARAMS32_29 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[30], OD_FILTER_PARAMS32_30, 32, 117);
  t[31] += (t[30]*OD_FILTER_PARAMS32_30 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[31], OD_FILTER_PARAMS32_45, 32, 118);
  t[30] += (t[31]*OD_FILTER_PARAMS32_45 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[30], OD_FILTER_PARAMS32_44, 32, 119);
  t[29] += (t[30]*OD_FILTER_PARAMS32_44 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[29], OD_FILTER_PARAMS32_43, 32, 120);
  t[28] += (t[29]*OD_FILTER_PARAMS32_43 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[28], OD_FILTER_PARAMS32_42, 32, 121);
  t[27] += (t[28]*OD_FILTER_PARAMS32_42 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[27], OD_FILTER_PARAMS32_41, 32, 122);
  t[26] += (t[27]*OD_FILTER_PARAMS32_41 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[26], OD_FILTER_PARAMS32_40, 32, 123);
  t[25] += (t[26]*OD_FILTER_PARAMS32_40 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[25], OD_FILTER_PARAMS32_39, 32, 124);
  t[24] += (t[25]*OD_FILTER_PARAMS32_39 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[24], OD_FILTER_PARAMS32_38, 32, 125);
  t[23] += (t[24]*OD_FILTER_PARAMS32_38 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[23], OD_FILTER_PARAMS32_37, 32, 126);
  t[22] += (t[23]*OD_FILTER_PARAMS32_37 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[22], OD_FILTER_PARAMS32_36, 32, 127);
  t[21] += (t[22]*OD_FILTER_PARAMS32_36 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[21], OD_FILTER_PARAMS32_35, 32, 128);
  t[20] += (t[21]*OD_FILTER_PARAMS32_35 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[20], OD_FILTER_PARAMS32_34, 32, 129);
  t[19] += (t[20]*OD_FILTER_PARAMS32_34 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[19], OD_FILTER_PARAMS32_33, 32, 130);
  t[18] += (t[19]*OD_FILTER_PARAMS32_33 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[18], OD_FILTER_PARAMS32_32, 32, 131);
  t[17] += (t[18]*OD_FILTER_PARAMS32_32 + 32) >> 6;
  OD_DCT_OVERFLOW_CHECK(t[17], OD_FILTER_PARAMS32_31, 32, 132);
  t[16] += (t[17]*OD_FILTER_PARAMS32_31 + 32) >> 6;
#endif
  t[0] += t[31] >> 1;
  _y[0] = (od_coeff)t[0];
  t[1] += t[30] >> 1;
  _y[1] = (od_coeff)t[1];
  t[2] += t[29] >> 1;
  _y[2] = (od_coeff)t[2];
  t[3] += t[28] >> 1;
  _y[3] = (od_coeff)t[3];
  t[4] += t[27] >> 1;
  _y[4] = (od_coeff)t[4];
  t[5] += t[26] >> 1;
  _y[5] = (od_coeff)t[5];
  t[6] += t[25] >> 1;
  _y[6] = (od_coeff)t[6];
  t[7] += t[24] >> 1;
  _y[7] = (od_coeff)t[7];
  t[8] += t[23] >> 1;
  _y[8] = (od_coeff)t[8];
  t[9] += t[22] >> 1;
  _y[9] = (od_coeff)t[9];
  t[10] += t[21] >> 1;
  _y[10] = (od_coeff)t[10];
  t[11] += t[20] >> 1;
  _y[11] = (od_coeff)t[11];
  t[12] += t[19] >> 1;
  _y[12] = (od_coeff)t[12];
  t[13] += t[18] >> 1;
  _y[13] = (od_coeff)t[13];
  t[14] += t[17] >> 1;
  _y[14] = (od_coeff)t[14];
  t[15] += t[16] >> 1;
  _y[15] = (od_coeff)t[15];
  _y[16] = (od_coeff)(t[15] - t[16]);
  _y[17] = (od_coeff)(t[14] - t[17]);
  _y[18] = (od_coeff)(t[13] - t[18]);
  _y[19] = (od_coeff)(t[12] - t[19]);
  _y[20] = (od_coeff)(t[11] - t[20]);
  _y[21] = (od_coeff)(t[10] - t[21]);
  _y[22] = (od_coeff)(t[9] - t[22]);
  _y[23] = (od_coeff)(t[8] - t[23]);
  _y[24] = (od_coeff)(t[7] - t[24]);
  _y[25] = (od_coeff)(t[6] - t[25]);
  _y[26] = (od_coeff)(t[5] - t[26]);
  _y[27] = (od_coeff)(t[4] - t[27]);
  _y[28] = (od_coeff)(t[3] - t[28]);
  _y[29] = (od_coeff)(t[2] - t[29]);
  _y[30] = (od_coeff)(t[1] - t[30]);
  _y[31] = (od_coeff)(t[0] - t[31]);
#endif
}

void od_post_filter32(od_coeff _x[32], const od_coeff _y[32]) {
#if OD_DISABLE_FILTER
  int i;
  for (i = 0; i < 32; i++) {
    _x[i] = _y[i];
  }
#else
  int t[32];
  t[31] = _y[0] - _y[31];
  t[30] = _y[1] - _y[30];
  t[29] = _y[2] - _y[29];
  t[28] = _y[3] - _y[28];
  t[27] = _y[4] - _y[27];
  t[26] = _y[5] - _y[26];
  t[25] = _y[6] - _y[25];
  t[24] = _y[7] - _y[24];
  t[23] = _y[8] - _y[23];
  t[22] = _y[9] - _y[22];
  t[21] = _y[10] - _y[21];
  t[20] = _y[11] - _y[20];
  t[19] = _y[12] - _y[19];
  t[18] = _y[13] - _y[18];
  t[17] = _y[14] - _y[17];
  t[16] = _y[15] - _y[16];
  t[15] = _y[15] - (t[16] >> 1);
  t[14] = _y[14] - (t[17] >> 1);
  t[13] = _y[13] - (t[18] >> 1);
  t[12] = _y[12] - (t[19] >> 1);
  t[11] = _y[11] - (t[20] >> 1);
  t[10] = _y[10] - (t[21] >> 1);
  t[9] = _y[9] - (t[22] >> 1);
  t[8] = _y[8] - (t[23] >> 1);
  t[7] = _y[7] - (t[24] >> 1);
  t[6] = _y[6] - (t[25] >> 1);
  t[5] = _y[5] - (t[26] >> 1);
  t[4] = _y[4] - (t[27] >> 1);
  t[3] = _y[3] - (t[28] >> 1);
  t[2] = _y[2] - (t[29] >> 1);
  t[1] = _y[1] - (t[30] >> 1);
  t[0] = _y[0] - (t[31] >> 1);
#if OD_FILTER32_TYPE3
  t[16] -= (t[17]*OD_FILTER_PARAMS32_31 + 32) >> 6;
  t[17] -= (t[16]*OD_FILTER_PARAMS32_16 + 32) >> 6;
  t[17] -= (t[18]*OD_FILTER_PARAMS32_32 + 32) >> 6;
  t[18] -= (t[17]*OD_FILTER_PARAMS32_17 + 32) >> 6;
  t[18] -= (t[19]*OD_FILTER_PARAMS32_33 + 32) >> 6;
  t[19] -= (t[18]*OD_FILTER_PARAMS32_18 + 32) >> 6;
  t[19] -= (t[20]*OD_FILTER_PARAMS32_34 + 32) >> 6;
  t[20] -= (t[19]*OD_FILTER_PARAMS32_19 + 32) >> 6;
  t[20] -= (t[21]*OD_FILTER_PARAMS32_35 + 32) >> 6;
  t[21] -= (t[20]*OD_FILTER_PARAMS32_20 + 32) >> 6;
  t[21] -= (t[22]*OD_FILTER_PARAMS32_36 + 32) >> 6;
  t[22] -= (t[21]*OD_FILTER_PARAMS32_21 + 32) >> 6;
  t[22] -= (t[23]*OD_FILTER_PARAMS32_37 + 32) >> 6;
  t[23] -= (t[22]*OD_FILTER_PARAMS32_22 + 32) >> 6;
  t[23] -= (t[24]*OD_FILTER_PARAMS32_38 + 32) >> 6;
  t[24] -= (t[23]*OD_FILTER_PARAMS32_23 + 32) >> 6;
  t[24] -= (t[25]*OD_FILTER_PARAMS32_39 + 32) >> 6;
  t[25] -= (t[24]*OD_FILTER_PARAMS32_24 + 32) >> 6;
  t[25] -= (t[26]*OD_FILTER_PARAMS32_40 + 32) >> 6;
  t[26] -= (t[25]*OD_FILTER_PARAMS32_25 + 32) >> 6;
  t[26] -= (t[27]*OD_FILTER_PARAMS32_41 + 32) >> 6;
  t[27] -= (t[26]*OD_FILTER_PARAMS32_26 + 32) >> 6;
  t[27] -= (t[28]*OD_FILTER_PARAMS32_42 + 32) >> 6;
  t[28] -= (t[27]*OD_FILTER_PARAMS32_27 + 32) >> 6;
  t[28] -= (t[29]*OD_FILTER_PARAMS32_43 + 32) >> 6;
  t[29] -= (t[28]*OD_FILTER_PARAMS32_28 + 32) >> 6;
  t[29] -= (t[30]*OD_FILTER_PARAMS32_44 + 32) >> 6;
  t[30] -= (t[29]*OD_FILTER_PARAMS32_29 + 32) >> 6;
  t[30] -= (t[31]*OD_FILTER_PARAMS32_45 + 32) >> 6;
  t[31] -= (t[30]*OD_FILTER_PARAMS32_30 + 32) >> 6;
#else
  t[16] -= (t[17]*OD_FILTER_PARAMS32_31 + 32) >> 6;
  t[17] -= (t[18]*OD_FILTER_PARAMS32_32 + 32) >> 6;
  t[18] -= (t[19]*OD_FILTER_PARAMS32_33 + 32) >> 6;
  t[19] -= (t[20]*OD_FILTER_PARAMS32_34 + 32) >> 6;
  t[20] -= (t[21]*OD_FILTER_PARAMS32_35 + 32) >> 6;
  t[21] -= (t[22]*OD_FILTER_PARAMS32_36 + 32) >> 6;
  t[22] -= (t[23]*OD_FILTER_PARAMS32_37 + 32) >> 6;
  t[23] -= (t[24]*OD_FILTER_PARAMS32_38 + 32) >> 6;
  t[24] -= (t[25]*OD_FILTER_PARAMS32_39 + 32) >> 6;
  t[25] -= (t[26]*OD_FILTER_PARAMS32_40 + 32) >> 6;
  t[26] -= (t[27]*OD_FILTER_PARAMS32_41 + 32) >> 6;
  t[27] -= (t[28]*OD_FILTER_PARAMS32_42 + 32) >> 6;
  t[28] -= (t[29]*OD_FILTER_PARAMS32_43 + 32) >> 6;
  t[29] -= (t[30]*OD_FILTER_PARAMS32_44 + 32) >> 6;
  t[30] -= (t[31]*OD_FILTER_PARAMS32_45 + 32) >> 6;
  t[31] -= (t[30]*OD_FILTER_PARAMS32_30 + 32) >> 6;
  t[30] -= (t[29]*OD_FILTER_PARAMS32_29 + 32) >> 6;
  t[29] -= (t[28]*OD_FILTER_PARAMS32_28 + 32) >> 6;
  t[28] -= (t[27]*OD_FILTER_PARAMS32_27 + 32) >> 6;
  t[27] -= (t[26]*OD_FILTER_PARAMS32_26 + 32) >> 6;
  t[26] -= (t[25]*OD_FILTER_PARAMS32_25 + 32) >> 6;
  t[25] -= (t[24]*OD_FILTER_PARAMS32_24 + 32) >> 6;
  t[24] -= (t[23]*OD_FILTER_PARAMS32_23 + 32) >> 6;
  t[23] -= (t[22]*OD_FILTER_PARAMS32_22 + 32) >> 6;
  t[22] -= (t[21]*OD_FILTER_PARAMS32_21 + 32) >> 6;
  t[21] -= (t[20]*OD_FILTER_PARAMS32_20 + 32) >> 6;
  t[20] -= (t[19]*OD_FILTER_PARAMS32_19 + 32) >> 6;
  t[19] -= (t[18]*OD_FILTER_PARAMS32_18 + 32) >> 6;
  t[18] -= (t[17]*OD_FILTER_PARAMS32_17 + 32) >> 6;
  t[17] -= (t[16]*OD_FILTER_PARAMS32_16 + 32) >> 6;
#endif
#if OD_FILTER_PARAMS32_15 != 64
  t[31] = (t[31] << 6)/OD_FILTER_PARAMS32_15;
#endif
#if OD_FILTER_PARAMS32_14 != 64
  t[30] = (t[30] << 6)/OD_FILTER_PARAMS32_14;
#endif
#if OD_FILTER_PARAMS32_13 != 64
  t[29] = (t[29] << 6)/OD_FILTER_PARAMS32_13;
#endif
#if OD_FILTER_PARAMS32_12 != 64
  t[28] = (t[28] << 6)/OD_FILTER_PARAMS32_12;
#endif
#if OD_FILTER_PARAMS32_11 != 64
  t[27] = (t[27] << 6)/OD_FILTER_PARAMS32_11;
#endif
#if OD_FILTER_PARAMS32_10 != 64
  t[26] = (t[26] << 6)/OD_FILTER_PARAMS32_10;
#endif
#if OD_FILTER_PARAMS32_9 != 64
  t[25] = (t[25] << 6)/OD_FILTER_PARAMS32_9;
#endif
#if OD_FILTER_PARAMS32_8 != 64
  t[24] = (t[24] << 6)/OD_FILTER_PARAMS32_8;
#endif
#if OD_FILTER_PARAMS32_7 != 64
  t[23] = (t[23] << 6)/OD_FILTER_PARAMS32_7;
#endif
#if OD_FILTER_PARAMS32_6 != 64
  t[22] = (t[22] << 6)/OD_FILTER_PARAMS32_6;
#endif
#if OD_FILTER_PARAMS32_5 != 64
  t[21] = (t[21] << 6)/OD_FILTER_PARAMS32_5;
#endif
#if OD_FILTER_PARAMS32_4 != 64
  t[20] = (t[20] << 6)/OD_FILTER_PARAMS32_4;
#endif
#if OD_FILTER_PARAMS32_3 != 64
  t[19] = (t[19] << 6)/OD_FILTER_PARAMS32_3;
#endif
#if OD_FILTER_PARAMS32_2 != 64
  t[18] = (t[18] << 6)/OD_FILTER_PARAMS32_2;
#endif
#if OD_FILTER_PARAMS32_1 != 64
  t[17] = (t[17] << 6)/OD_FILTER_PARAMS32_1;
#endif
#if OD_FILTER_PARAMS32_0 != 64
  t[16] = (t[16] << 6)/OD_FILTER_PARAMS32_0;
#endif
  t[0] += t[31] >> 1;
  _x[0] = (od_coeff)t[0];
  t[1] += t[30] >> 1;
  _x[1] = (od_coeff)t[1];
  t[2] += t[29] >> 1;
  _x[2] = (od_coeff)t[2];
  t[3] += t[28] >> 1;
  _x[3] = (od_coeff)t[3];
  t[4] += t[27] >> 1;
  _x[4] = (od_coeff)t[4];
  t[5] += t[26] >> 1;
  _x[5] = (od_coeff)t[5];
  t[6] += t[25] >> 1;
  _x[6] = (od_coeff)t[6];
  t[7] += t[24] >> 1;
  _x[7] = (od_coeff)t[7];
  t[8] += t[23] >> 1;
  _x[8] = (od_coeff)t[8];
  t[9] += t[22] >> 1;
  _x[9] = (od_coeff)t[9];
  t[10] += t[21] >> 1;
  _x[10] = (od_coeff)t[10];
  t[11] += t[20] >> 1;
  _x[11] = (od_coeff)t[11];
  t[12] += t[19] >> 1;
  _x[12] = (od_coeff)t[12];
  t[13] += t[18] >> 1;
  _x[13] = (od_coeff)t[13];
  t[14] += t[17] >> 1;
  _x[14] = (od_coeff)t[14];
  t[15] += t[16] >> 1;
  _x[15] = (od_coeff)t[15];
  _x[16] = (od_coeff)(t[15] - t[16]);
  _x[17] = (od_coeff)(t[14] - t[17]);
  _x[18] = (od_coeff)(t[13] - t[18]);
  _x[19] = (od_coeff)(t[12] - t[19]);
  _x[20] = (od_coeff)(t[11] - t[20]);
  _x[21] = (od_coeff)(t[10] - t[21]);
  _x[22] = (od_coeff)(t[9] - t[22]);
  _x[23] = (od_coeff)(t[8] - t[23]);
  _x[24] = (od_coeff)(t[7] - t[24]);
  _x[25] = (od_coeff)(t[6] - t[25]);
  _x[26] = (od_coeff)(t[5] - t[26]);
  _x[27] = (od_coeff)(t[4] - t[27]);
  _x[28] = (od_coeff)(t[3] - t[28]);
  _x[29] = (od_coeff)(t[2] - t[29]);
  _x[30] = (od_coeff)(t[1] - t[30]);
  _x[31] = (od_coeff)(t[0] - t[31]);
#endif
}

#if OD_DEBLOCKING

/* This is a slightly modified implementation of the Thor deblocking filter
   adapted for use in Daala. See draft-fuldseth-netvc-thor-00 or the Thor
   codebase for more details. */

int beta_table[52] = {
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
     16, 17, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64
};

static const int tc_table[56] =
{
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,5,5,6,6,7,8,9,9,10,10,11,11,12,12,13,13,14,14
};

#define OD_DEBLOCK_BETA(q) (beta_table[OD_CLAMPI(0, (q) - 8, 51)] << OD_COEFF_SHIFT)
#define OD_DEBLOCK_TC(q) (tc_table[OD_CLAMPI(0, (q) - 8, 51)] << OD_COEFF_SHIFT)

void od_thor_deblock_col8(od_coeff *c0, int stride, int q) {
  od_coeff p12;
  od_coeff p02;
  od_coeff q12;
  od_coeff q02;
  od_coeff p15;
  od_coeff p05;
  od_coeff q15;
  od_coeff q05;
  od_coeff d;

  p12 = c0[2*stride - 2];
  p02 = c0[2*stride - 1];
  q02 = c0[2*stride + 0];
  q12 = c0[2*stride + 1];
  p15 = c0[5*stride - 2];
  p05 = c0[5*stride - 1];
  q05 = c0[5*stride + 0];
  q15 = c0[5*stride + 1];
  d = abs(p12-p02) + abs(q12-q02) + abs(p15-p05) + abs(q15-q05);
  if (d < OD_DEBLOCK_BETA(q)) {
    int k;
    od_coeff tc;
    tc = OD_DEBLOCK_TC(q);
    for (k = 0; k < 8; k++) {
      od_coeff p0;
      od_coeff p1;
      od_coeff q0;
      od_coeff q1;
      od_coeff delta;
      p1 = c0[k*stride - 2];
      p0 = c0[k*stride - 1];
      q0 = c0[k*stride + 0];
      q1 = c0[k*stride + 1];
      delta = (18*(q0-p0) - 6*(q1-p1) + 16)>>5;
      delta = OD_CLAMPI(-tc, delta,tc);
      c0[k*stride - 2] = p1 + delta/2;
      c0[k*stride - 1] = p0 + delta;
      c0[k*stride + 0] = q0 - delta;
      c0[k*stride + 1] = q1 - delta/2;
    }
  }
}

void od_thor_deblock_row8(od_coeff *c0, int stride, int q) {
  od_coeff p12;
  od_coeff p02;
  od_coeff q12;
  od_coeff q02;
  od_coeff p15;
  od_coeff p05;
  od_coeff q15;
  od_coeff q05;
  od_coeff d;

  p12 = c0[2 - 2*stride];
  p02 = c0[2 - 1*stride];
  q02 = c0[2 + 0*stride];
  q12 = c0[2 + 1*stride];
  p15 = c0[5 - 2*stride];
  p05 = c0[5 - 1*stride];
  q05 = c0[5 + 0*stride];
  q15 = c0[5 + 1*stride];
  d = abs(p12-p02) + abs(q12-q02) + abs(p15-p05) + abs(q15-q05);
  if (d < OD_DEBLOCK_BETA(q)) {
    int k;
    od_coeff tc;
    tc = OD_DEBLOCK_TC(q);
    for (k = 0; k < 8; k++) {
      od_coeff p0;
      od_coeff p1;
      od_coeff q0;
      od_coeff q1;
      od_coeff delta;
      p1 = c0[k - 2*stride];
      p0 = c0[k - 1*stride];
      q0 = c0[k + 0*stride];
      q1 = c0[k + 1*stride];
      delta = (18*(q0-p0) - 6*(q1-p1) + 16)>>5;
      delta = OD_CLAMPI(-tc, delta,tc);
      c0[k - 2*stride] = p1 + delta/2;
      c0[k - 1*stride] = p0 + delta;
      c0[k + 0*stride] = q0 - delta;
      c0[k + 1*stride] = q1 - delta/2;
    }
  }
}
#endif

#define OD_BLOCK_SIZE4x4_DEC(bsize, bstride, bx, by, dec) \
 OD_MAXI(OD_BLOCK_SIZE4x4(bsize, bstride, bx, by), dec)

void od_prefilter_split(od_coeff *c0, int stride, int bs, int f) {
#if OD_DEBLOCKING
#else
  int i;
  int j;
  od_coeff *c;
  c = c0 + ((2 << bs) - (2 << f))*stride;
  for (j = 0; j < 4 << bs; j++) {
    int k;
    od_coeff t[4 << OD_NBSIZES];
    for (k = 0; k < 4 << f; k++) t[k] = c[stride*k + j];
    (*OD_PRE_FILTER[f])(t, t);
    for (k = 0; k < 4 << f; k++) c[stride*k + j] = t[k];
  }
  c = c0 + (2 << bs) - (2 << f);
  for (i = 0; i < 4 << bs; i++) {
    (*OD_PRE_FILTER[f])(c + i*stride, c + i*stride);
  }
#endif
}

void od_postfilter_split(od_coeff *c0, int stride, int bs, int f, int q,
 unsigned char *skip, int skip_stride) {
  int i;
  od_coeff *c;
#if OD_DEBLOCKING
  if (bs==0) return;
  c = c0 + (2 << bs);
  for (i = 0; i < 4 << bs; i += 8) {
    if (!skip[(i >> 2)*skip_stride + (1 << bs >> 1) - 1]
     || !skip[(i >> 2)*skip_stride + (1 << bs >> 1)]) {
      od_thor_deblock_col8(c + i*stride, stride, q);
    }
  }
  c = c0 + (2 << bs)*stride;
  for (i = 0; i < 4 << bs; i += 8) {
    if (!skip[((1 << bs >> 1) - 1)*skip_stride + (i >> 2)]
     || !skip[(1 << bs >> 1)*skip_stride + (i >> 2)]) {
      od_thor_deblock_row8(c + i, stride, q);
    }
  }
#else
  int j;
  (void)q;
  (void)skip;
  (void)skip_stride;
  c = c0 + (2 << bs) - (2 << f);
  for (i = 0; i < 4 << bs; i++) {
    (*OD_POST_FILTER[f])(c + i*stride, c + i*stride);
  }
  c = c0 + ((2 << bs) - (2 << f))*stride;
  for (j = 0; j < 4 << bs; j++) {
    int k;
    od_coeff t[4 << OD_NBSIZES];
    for (k = 0; k < 4 << f; k++) t[k] = c[stride*k + j];
    (*OD_POST_FILTER[f])(t, t);
    for (k = 0; k < 4 << f; k++) c[stride*k + j] = t[k];
  }
#endif
}

void od_apply_prefilter_frame_sbs(od_coeff *c0, int stride, int nhsb, int nvsb,
 int xdec, int ydec) {
#if OD_DEBLOCKING
#else
  int sbx;
  int sby;
  int i;
  int j;
  int f;
  od_coeff *c;
  f = OD_FILT_SIZE(OD_NBSIZES - 1, xdec);
  c = c0 + ((OD_BSIZE_MAX >> ydec) - (2 << f))*stride;
  for (sby = 1; sby < nvsb; sby++) {
    for (j = 0; j < nhsb << OD_LOG_BSIZE_MAX >> xdec; j++) {
      int k;
      od_coeff t[4 << OD_NBSIZES];
      for (k = 0; k < 4 << f; k++) t[k] = c[stride*k + j];
      (*OD_PRE_FILTER[f])(t, t);
      for (k = 0; k < 4 << f; k++) c[stride*k + j] = t[k];
    }
    c += OD_BSIZE_MAX*stride >> ydec;
  }
  c = c0 + (OD_BSIZE_MAX >> ydec) - (2 << f);
  for (sbx = 1; sbx < nhsb; sbx++) {
    for (i = 0; i < nvsb << OD_LOG_BSIZE_MAX >> ydec; i++) {
      (*OD_PRE_FILTER[f])(c + i*stride, c + i*stride);
    }
    c += OD_BSIZE_MAX >> xdec;
  }
#endif
}

void od_apply_postfilter_frame_sbs(od_coeff *c0, int stride, int nhsb,
 int nvsb, int xdec, int ydec, int q, unsigned char *skip, int skip_stride) {
#if OD_DEBLOCKING
  od_coeff *c;
  int sbx;
  int sby;
  c = c0 + (OD_BSIZE_MAX >> ydec);
  for (sbx = 1; sbx < nhsb; sbx++) {
    int i;
    for (i = 0; i < nvsb << OD_LOG_BSIZE_MAX >> ydec; i+=8) {
      if (!skip[(i >> 2)*skip_stride + (sbx << 3 >> xdec) - 1]
       || !skip[(i >> 2)*skip_stride + (sbx << 3 >> xdec)]) {
        od_thor_deblock_col8(c + i*stride, stride, q);
      }
    }
    c += OD_BSIZE_MAX >> xdec;
  }
  c = c0 + (OD_BSIZE_MAX >> ydec)*stride;
  for (sby = 1; sby < nvsb; sby++) {
    int j;
    for (j = 0; j < nhsb << OD_LOG_BSIZE_MAX >> xdec; j+=8) {
      if (!skip[((sby << 3 >> xdec) - 1)*skip_stride + (j >> 2)]
       || !skip[(sby << 3 >> xdec)*skip_stride + (j >> 2)]) {
        od_thor_deblock_row8(c + j, stride, q);
      }
    }
    c += OD_BSIZE_MAX*stride >> ydec;
  }
#else
  int sbx;
  int sby;
  int i;
  int j;
  int f;
  od_coeff *c;
  (void)q;
  (void)skip;
  (void)skip_stride;
  f = OD_FILT_SIZE(OD_NBSIZES - 1, xdec);
  c = c0 + (OD_BSIZE_MAX >> ydec) - (2 << f);
  for (sbx = 1; sbx < nhsb; sbx++) {
    for (i = 0; i < nvsb << OD_LOG_BSIZE_MAX >> ydec; i++) {
      (*OD_POST_FILTER[f])(c + i*stride, c + i*stride);
    }
    c += OD_BSIZE_MAX >> xdec;
  }
  c = c0 + ((OD_BSIZE_MAX >> ydec) - (2 << f))*stride;
  for (sby = 1; sby < nvsb; sby++) {
    for (j = 0; j < nhsb << OD_LOG_BSIZE_MAX >> xdec; j++) {
      int k;
      od_coeff t[4 << OD_NBSIZES];
      for (k = 0; k < 4 << f; k++) t[k] = c[stride*k + j];
      (*OD_POST_FILTER[f])(t, t);
      for (k = 0; k < 4 << f; k++) c[stride*k + j] = t[k];
    }
    c += OD_BSIZE_MAX*stride >> ydec;
  }
#endif
}

/* Detect direction. 0 means 45-degree up-right, 2 is horizontal, and so on.
   The search minimizes the weighted variance along all the lines in a
   particular direction, i.e. the squared error between the input and a
   "predicted" block where each pixel is replaced by the average along a line
   in a particular direction. Since each direction have the same sum(x^2) term,
   that term is never computed. See Section 2, step 2, of:
   http://jmvalin.ca/notes/intra_paint.pdf */
static int od_dir_find8(const int16_t *img, int stride, int32_t *var) {
  int i;
  int cost[8] = {0};
  int partial[8][15] = {{0}};
  int best_cost = 0;
  int best_dir = 0;
  for (i = 0; i < 8; i++) {
    int j;
    for (j = 0; j < 8; j++) {
      int x;
      x = img[i*stride + j] >> OD_COEFF_SHIFT;
      partial[0][i + j] += x;
      partial[1][i + j/2] += x;
      partial[2][i] += x;
      partial[3][3 + i - j/2] += x;
      partial[4][7 + i - j] += x;
      partial[5][3 - i/2 + j] += x;
      partial[6][j] += x;
      partial[7][i/2 + j] += x;
    }
  }
  for (i = 0; i < 8; i++) {
    cost[2] += partial[2][i]*partial[2][i] >> 3;
    cost[6] += partial[6][i]*partial[6][i] >> 3;
  }
  for (i = 0; i < 7; i++) {
    cost[0] += OD_DIVU_SMALL(partial[0][i]*partial[0][i], i + 1)
     + OD_DIVU_SMALL(partial[0][14 - i]*partial[0][14 - i], i + 1);
    cost[4] += OD_DIVU_SMALL(partial[4][i]*partial[4][i], i + 1)
     + OD_DIVU_SMALL(partial[4][14 - i]*partial[4][14 - i], i + 1);
  }
  cost[0] += partial[0][7]*partial[0][8 - 1] >> 3;
  cost[4] += partial[4][7]*partial[4][8 - 1] >> 3;
  for (i = 1; i < 8; i += 2) {
    int j;
    for (j = 0; j < 4 + 1; j++) {
      cost[i] += partial[i][3 + j]*partial[i][3 + j] >> 3;
    }
    for (j = 0; j < 4 - 1; j++) {
      cost[i] += OD_DIVU_SMALL(partial[i][j]*partial[i][j], 2*j + 2)
       + OD_DIVU_SMALL(partial[i][10 - j]*partial[i][10 - j], 2*j + 2);
    }
  }
  for (i = 0; i < 8; i++) {
    if (cost[i] > best_cost) {
      best_cost = cost[i];
      best_dir = i;
    }
  }
  /* Difference between the optimal variance and the variance along the
     orthogonal direction. Again, the sum(x^2) terms cancel out. */
  *var = best_cost - cost[(best_dir + 4) & 7];
  return best_dir;
}

#define OD_DERING_VERY_LARGE (30000)
#define OD_DERING_INBUF_SIZE ((OD_BSIZE_MAX + 2*OD_FILT_BORDER)*\
 (OD_BSIZE_MAX + 2*OD_FILT_BORDER))

/* Smooth in the direction detected. */
void od_filter_dering_direction_c(int16_t *y, int ystride, int16_t *in,
 int ln, int threshold, int dir) {
  int i;
  int j;
  int k;
  static const int taps[3] = {3, 2, 2};
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      od_coeff sum;
      od_coeff xx;
      od_coeff yy;
      xx = in[i*OD_FILT_BSTRIDE + j];
      sum= 0;
      for (k = 0; k < 3; k++) {
        od_coeff p0;
        od_coeff p1;
        p0 = in[i*OD_FILT_BSTRIDE + j + direction_offsets_table[dir][k]] - xx;
        p1 = in[i*OD_FILT_BSTRIDE + j - direction_offsets_table[dir][k]] - xx;
        if (abs(p0) < threshold) sum += taps[k]*p0;
        if (abs(p1) < threshold) sum += taps[k]*p1;
      }
      yy = xx + ((sum + 8) >> 4);
      y[i*ystride + j] = yy;
    }
  }
}

void od_filter_dering_direction_4x4_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  od_filter_dering_direction_c(y, ystride, in, 2, threshold, dir);
}

void od_filter_dering_direction_8x8_c(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  od_filter_dering_direction_c(y, ystride, in, 3, threshold, dir);
}

/* Smooth in the direction orthogonal to what was detected. */
void od_filter_dering_orthogonal_c(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int ln, int threshold, int dir) {
  int i;
  int j;
  int offset;
  if (dir <= 4) offset = OD_FILT_BSTRIDE;
  else offset = 1;
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      od_coeff athresh;
      od_coeff yy;
      od_coeff sum;
      od_coeff p;
      /* Deringing orthogonal to the direction uses a tighter threshold
         because we want to be conservative. We've presumably already
         achieved some deringing, so the amount of change is expected
         to be low. Also, since we might be filtering across an edge, we
         want to make sure not to blur it. That being said, we might want
         to be a little bit more aggressive on pure horizontal/vertical
         since the ringing there tends to be directional, so it doesn't
         get removed by the directional filtering. */
      athresh = OD_MINI(threshold, threshold/3
       + abs(in[i*OD_FILT_BSTRIDE + j] - x[i*xstride + j]));
      yy = in[i*OD_FILT_BSTRIDE + j];
      sum = 0;
      p = in[i*OD_FILT_BSTRIDE + j + offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j - offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j + 2*offset] - yy;
      if (abs(p) < athresh) sum += p;
      p = in[i*OD_FILT_BSTRIDE + j - 2*offset] - yy;
      if (abs(p) < athresh) sum += p;
      y[i*ystride + j] = yy + ((3*sum + 8) >> 4);
    }
  }
}

void od_filter_dering_orthogonal_4x4_c(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int threshold, int dir) {
  od_filter_dering_orthogonal_c(y, ystride, in, x, xstride, 2, threshold, dir);
}

void od_filter_dering_orthogonal_8x8_c(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int threshold, int dir) {
  od_filter_dering_orthogonal_c(y, ystride, in, x, xstride, 3, threshold, dir);
}

/* This table approximates x^0.16 with the index being log2(x). It is clamped
   to [-.5, 3]. The table is computed as:
   round(256*min(3, max(.5, 1.08*(sqrt(2)*2.^([0:17]+8)/256/256).^.16))) */
static int16_t od_thresh_table_q8[18] = {
  128, 134, 150, 168, 188, 210, 234, 262,
  292, 327, 365, 408, 455, 509, 569, 635,
  710, 768,
};

/* Compute deringing filter threshold for each 8x8 block based on the
   directional variance difference. A high variance difference means that we
   have a highly directional pattern (e.g. a high contrast edge), so we can
   apply more deringing. A low variance means that we either have a low
   contrast edge, or a non-directional texture, so we want to be careful not
   to blur. */
static void od_compute_thresh(int thresh[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int threshold, int32_t var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int32_t sb_var, int nhb, int nvb) {
  int bx;
  int by;
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      int v1;
      int v2;
      /* We use both the variance of 8x8 blocks and the variance of the
         entire superblock to determine the threshold. */
      v1 = OD_MINI(32767, var[by][bx] >> 6);
      v2 = OD_MINI(32767, sb_var/(OD_BSIZE_MAX*OD_BSIZE_MAX));
      thresh[by][bx] = threshold*od_thresh_table_q8[OD_CLAMPI(0,
       OD_ILOG(v1*v2) - 9, 17)] >> 8;
    }
  }
}

void od_dering(od_state *state, int16_t *y, int ystride, int16_t *x, int
 xstride, int ln, int sbx, int sby, int nhsb, int nvsb, int q, int xdec,
 int dir[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS],
 int pli, unsigned char *bskip, int skip_stride) {
  int i;
  int j;
  int n;
  int threshold;
  int bx;
  int by;
  int16_t inbuf[OD_DERING_INBUF_SIZE];
  int16_t *in;
  int nhb;
  int nvb;
  int bsize;
  int varsum = 0;
  int32_t var[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
  int thresh[OD_DERING_NBLOCKS][OD_DERING_NBLOCKS];
  n = 1 << ln;
  bsize = 3 - xdec;
  nhb = nvb = n >> bsize;
  in = inbuf + OD_FILT_BORDER*OD_FILT_BSTRIDE + OD_FILT_BORDER;
  /* We avoid filtering the pixels for which some of the pixels to average
     are outside the frame. We could change the filter instead, but it would
     add special cases for any future vectorization. */
  for (i = 0; i < OD_DERING_INBUF_SIZE; i++) inbuf[i] = OD_DERING_VERY_LARGE;
  for (i = -OD_FILT_BORDER*(sby != 0); i < n
   + OD_FILT_BORDER*(sby != nvsb - 1); i++) {
    for (j = -OD_FILT_BORDER*(sbx != 0); j < n
     + OD_FILT_BORDER*(sbx != nhsb - 1); j++) {
      in[i*OD_FILT_BSTRIDE + j] = x[i*xstride + j];
    }
  }
  /* The threshold is meant to be the estimated amount of ringing for a given
     quantizer. Ringing is mostly proportional to the quantizer, but we
     use an exponent slightly smaller than unity because as quantization
     becomes coarser, the relative effect of quantization becomes slightly
     smaller as many unquantized coefficients are already close to zero. The
     value here comes from observing that on ntt-short, the best threshold for
     -v 5 appeared to be around 0.5*q, while the best threshold for -v 400
     was 0.25*q, i.e. 1-log(.5/.25)/log(400/5) = 0.84182 */
  threshold = 1.0*pow(q, 0.84182);
  if (pli == 0) {
    for (by = 0; by < nvb; by++) {
      for (bx = 0; bx < nhb; bx++) {
        dir[by][bx] = od_dir_find8(&x[8*by*xstride + 8*bx], xstride,
         &var[by][bx]);
        varsum += var[by][bx];
      }
    }
    od_compute_thresh(thresh, threshold, var, varsum, nhb, nvb);
  }
  else {
    for (by = 0; by < nvb; by++) {
      for (bx = 0; bx < nhb; bx++) {
        thresh[by][bx] = threshold;
      }
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      int xstart;
      int ystart;
      int xend;
      int yend;
      int skip;
      xstart = (sbx == 0) ? 0 : -1;
      ystart = (sby == 0) ? 0 : -1;
      xend = (2 >> xdec) + (sbx != nhsb - 1);
      yend = (2 >> xdec) + (sby != nvsb - 1);
      skip = 1;
      /* We look at whether the current block and its 4x4 surrounding (due to
         lapping) are skipped to avoid filtering the same content multiple
         times. */
      for (i = ystart; i < yend; i++) {
        for (j = xstart; j < xend; j++) {
          skip = skip && bskip[((by << 1 >> xdec) + i)*skip_stride
           + (bx << 1 >> xdec) + j];
        }
      }
      if (skip) thresh[by][bx] = 0;
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      (*state->opt_vtbl.filter_dering_direction[bsize - OD_LOG_BSIZE0])(
       &y[(by*ystride << bsize) + (bx << bsize)], ystride,
       &in[(by*OD_FILT_BSTRIDE << bsize) + (bx << bsize)],
       thresh[by][bx], dir[by][bx]);
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      in[i*OD_FILT_BSTRIDE + j] = y[i*ystride + j];
    }
  }
  for (by = 0; by < nvb; by++) {
    for (bx = 0; bx < nhb; bx++) {
      (*state->opt_vtbl.filter_dering_orthogonal[bsize - OD_LOG_BSIZE0])(
       &y[(by*ystride << bsize) + (bx << bsize)], ystride,
       &in[(by*OD_FILT_BSTRIDE << bsize) + (bx << bsize)],
       &x[(by*xstride << bsize) + (bx << bsize)], xstride,
       thresh[by][bx], dir[by][bx]);
    }
  }
}

/** Smoothes a block using bilinear interpolation from its four corners.
 *  The interpolation is applied using a weight that depends on the amount
 *  amount of distortion it causes to the signal compared to the quantization
 *  noise.
 * @param [in,out] x      block pixels
 * @param [in]     ln     log2 of block size
 * @param [in]     stride stride of x
 * @param [in]     q      quantizer
 * @param [in]     pli    plane index
 */
void od_bilinear_smooth(od_coeff *x, int ln, int stride, int q, int pli) {
  od_coeff x00;
  od_coeff x01;
  od_coeff x10;
  od_coeff x11;
  od_coeff a00;
  od_coeff a01;
  od_coeff a10;
  od_coeff a11;
  od_coeff y[OD_BSIZE_MAX][OD_BSIZE_MAX];
  int32_t dist;
  int shift;
  int w;
  int i;
  int j;
  int n;
  n = 1 << ln;
  x00 = x[0];
  x01 = x[n - 1];
  x10 = x[(n - 1)*stride];
  x11 = x[(n - 1)*stride + (n - 1)];
  a00 = x00;
  a01 = x01 - x00;
  a10 = x10 - x00;
  a11 = x11 + x00 - x10 - x01;
  /* Multiply by 1+1/n (approximation of n/(n-1)) here so that we can divide
     by n in the loop instead of dividing by n-1. */
  a01 += (a01 + n/2) >> ln;
  a10 += (a10 + n/2) >> ln;
  a11 += (2*a10 + n/2) >> ln;
  /* Minimal amount of shifting to avoid an overflow. */
  shift = OD_MAXI(0, 2*OD_COEFF_SHIFT + 2*ln - 16);
  dist = 0;
  /* Bilinear interpolation with non-linear x*y term. */
  for (i = 0; i < n; i++) {
    int32_t partial;
    partial = 0;
    for (j = 0; j < n; j++) {
      y[i][j] = a00 + ((j*a01 + i*a10 + (j*i*a11 >> ln) + n/2) >> ln);
      partial += (y[i][j] - x[i*stride + j])*(y[i][j] - x[i*stride + j]);
    }
    dist += partial >> shift;
  }
  /* Compensates for truncating above. */
  dist += n/2;
  dist >>= 2*ln - shift;
  /* Compute 1 - Wiener filter gain = strength * (q^2/12) / dist. */
  w = OD_MINI(1024, OD_BILINEAR_STRENGTH[pli]*q*q/(1 + 12*dist));
  /* Square the theoretical gain to attenuate the effect when we're unsure
     whether it's useful. */
  w = w*w >> 12;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      x[i*stride + j] -= (w*(x[i*stride + j]-y[i][j]) + 128) >> 8;
    }
  }
}

void od_smooth_recursive(od_coeff *c, unsigned char *bsize, int bstride,
 int bx, int by, int bsi, int w, int xdec, int ydec, int min_bs,
 int quantizer, int pli) {
  int obs;
  int bs;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  obs = OD_BLOCK_SIZE4x4(bsize, bstride, bx << bsi, by << bsi);
  bs = OD_MAXI(obs, xdec);
  OD_ASSERT(bs <= bsi);
  if (bs == bsi) {
    if (bs >= min_bs) {
      bs += OD_LOG_BSIZE0 - xdec;
      od_bilinear_smooth(&c[(by << bs)*w + (bx << bs)], bs, w, quantizer, pli);
    }
  }
  else {
    bx <<= 1;
    by <<= 1;
    bsi--;
    od_smooth_recursive(c, bsize, bstride, bx + 0, by + 0, bsi, w, xdec, ydec,
     min_bs, quantizer, pli);
    od_smooth_recursive(c, bsize, bstride, bx + 1, by + 0, bsi, w, xdec, ydec,
     min_bs, quantizer, pli);
    od_smooth_recursive(c, bsize, bstride, bx + 0, by + 1, bsi, w, xdec, ydec,
     min_bs, quantizer, pli);
    od_smooth_recursive(c, bsize, bstride, bx + 1, by + 1, bsi, w, xdec, ydec,
     min_bs, quantizer, pli);
  }
}

#if defined(TEST)
# include <stdio.h>

int minv[32];
int maxv[32];

int main(void) {
  int32_t min[16];
  int32_t max[16];
  int mini[16];
  int maxi[16];
  int dims;
  int i;
  int j;
  int err;
  err = 0;
  for (dims = 4; dims <= 16; dims <<= 1) {
    printf("filter%d:\n", dims);
    for (j = 0; j < dims; j++) {
      min[j] = mini[j] = INT_MAX;
      max[j] = maxi[j] = INT_MIN;
    }
    for (i = 0; i < (1<<dims); i++) {
      od_coeff x[16];
      od_coeff y[16];
      od_coeff x2[16];
      for (j = 0; j < dims; j++) x[j] = i>>j&1 ? 676 : -676;
      switch (dims) {
        case 16: od_pre_filter16(y, x);
          break;
        case  8: od_pre_filter8(y, x);
          break;
        case  4: od_pre_filter4(y, x);
          break;
      }
      for (j = 0; j < dims; j++) {
        if (y[j] < min[j]) {
          min[j] = y[j];
          mini[j] = i;
        }
        else if (y[j] > max[j]) {
          max[j] = y[j];
          maxi[j] = i;
        }
      }
      switch (dims) {
        case 16: od_post_filter16(x2, y);
          break;
        case  8: od_post_filter8(x2, y);
          break;
        case  4: od_post_filter4(x2, y);
          break;
      }
      for (j = 0; j < dims; j++)
        if (x[j] != x2[j]) {
          printf("Mismatch:\n");
          printf("in:    ");
          for (j = 0; j < dims; j++) printf(" %i", x[j]);
          printf("\nxform: ");
          for (j = 0; j < dims; j++) printf(" %i", y[j]);
          printf("\nout:    ");
          for (j = 0; j < dims; j++) printf(" %i", x2[j]);
          printf("\n\n");
          err |= 1;
        }
    }
    printf(" Min: ");
    for (j = 0; j < dims; j++)
      printf("%5i%c", min[j], j == dims-1 ? '\n' : ' ');
    printf(" Max: ");
    for (j = 0; j < dims; j++)
      printf("%5i%c", max[j], j == dims-1 ? '\n' : ' ');
  }
  return err;
}
#endif
