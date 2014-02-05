/*Daala video codec
Copyright (c) 2003-2013 Daala project contributors.  All rights reserved.

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

#include "filter.h"
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

  Or we use the form
    x0 ---------------------------------+------ y0
                                p0 | u0 | r0 |
    x1 ---------------------+------+---------+- y1
                    p1 | u1 | r1 |
    x2 ---------+------+---------+------------- y2
        p2 | u2 | r2 |
    x3 ----+---------+------------------------- y3
                           .
                           .
                           .

   which is a construction using three lifting steps to build a filter with
    an overall determinant of one.
   This is a necessary condition for having perfect
    reconstruction in the reverse directions:

   The scaling steps in the other factorizations give perfect reconstruction
    in only the forward direction: post(pre(image)) = image, but are
    not sufficient to make prediction for static scenes instantly convergent.
   For that we need pre(post(image)) = image as well.
  */

const od_filter_func OD_PRE_FILTER[OD_NBSIZES] = {
  od_pre_filter4,
  od_pre_filter8,
  od_pre_filter16
};

const od_filter_func OD_POST_FILTER[OD_NBSIZES] = {
  od_post_filter4,
  od_post_filter8,
  od_post_filter16
};


#if defined(TEST)
extern int minv[32];
extern int maxv[32];
# define _Check(_val, _idx) \
  if ((_val) < minv[(_idx)]) minv[(_idx)] = (_val); \
  if ((_val) > maxv[(_idx)]) maxv[(_idx)] = (_val);
#else
# define _Check(_val, _idx)
#endif

/*Filter parameters for the pre/post filters.
  When changing these the intra-predictors in
  initdata.c must be updated.*/
#define OD_FILTER_PARAMS4_0 (64)
#define OD_FILTER_PARAMS4_1 (64)
#define OD_FILTER_PARAMS4_2 (-42)

const int OD_FILTER_PARAMS4[3] = {
  OD_FILTER_PARAMS4_0, OD_FILTER_PARAMS4_1, OD_FILTER_PARAMS4_2
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
  /*Rotation:*/
#if OD_FILTER_PARAMS4_0 != 64
  t[3] += (t[2]*OD_FILTER_PARAMS4_0 + 32) >> 6;
#else
  t[3] += t[2];
#endif
#if OD_FILTER_PARAMS4_1 != 64
  t[2] = ((t[3]*OD_FILTER_PARAMS4_1 + 32) >> 6) - t[2];
#else
  t[2] = t[3] - t[2];
#endif
  t[3] += (t[2]*OD_FILTER_PARAMS4_2 + 32) >> 6;
  /*More +1/-1 butterflies (required for FIR, PR, LP).*/
  t[0] += t[2]>>1;
  _y[0] = (od_coeff)t[0];
  t[1] += t[3]>>1;
  _y[1] = (od_coeff)t[1];
  _y[2] = (od_coeff)(t[1]-t[3]);
  _y[3] = (od_coeff)(t[0]-t[2]);
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
  t[2] = _y[0]-_y[3];
  t[3] = _y[1]-_y[2];
  t[1] = _y[1]-(t[3]>>1);
  t[0] = _y[0]-(t[2]>>1);
  t[3] -= (t[2]*OD_FILTER_PARAMS4_2 + 32) >> 6;
#if OD_FILTER_PARAMS4_1 != 64
  t[2] = ((t[3]*OD_FILTER_PARAMS4_1 + 32) >> 6) - t[2];
#else
  t[2] = t[3] - t[2];
#endif
#if OD_FILTER_PARAMS4_0 != 64
  t[3] -= (t[2]*OD_FILTER_PARAMS4_0 + 32) >> 6;
#else
  t[3] -= t[2];
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
  t[4] = (t[4]*OD_FILTER_PARAMS8_0)>>6;
  t[4] += -t[4]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS8_1 != 64
  t[5] = (t[5]*OD_FILTER_PARAMS8_1)>>6;
  t[5] += -t[5]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS8_2 != 64
  t[6] = (t[6]*OD_FILTER_PARAMS8_2)>>6;
  t[6] += -t[6]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS8_3 != 64
  t[7] = (t[7]*OD_FILTER_PARAMS8_3)>>6;
  t[7] += -t[7]>>(OD_COEFF_BITS-1)&1;
#endif
  /*Rotations:*/
#if OD_FILTER8_TYPE3
  t[7] += (t[6]*OD_FILTER_PARAMS8_6+32)>>6;
  t[6] += (t[7]*OD_FILTER_PARAMS8_9+32)>>6;
  t[6] += (t[5]*OD_FILTER_PARAMS8_5+32)>>6;
  t[5] += (t[6]*OD_FILTER_PARAMS8_8+32)>>6;
  t[5] += (t[4]*OD_FILTER_PARAMS8_4+32)>>6;
  t[4] += (t[5]*OD_FILTER_PARAMS8_7+32)>>6;
#else
  t[5] += (t[4]*OD_FILTER_PARAMS8_4+32)>>6;
  t[6] += (t[5]*OD_FILTER_PARAMS8_5+32)>>6;
  t[7] += (t[6]*OD_FILTER_PARAMS8_6+32)>>6;
  t[6] += (t[7]*OD_FILTER_PARAMS8_9+32)>>6;
  t[5] += (t[6]*OD_FILTER_PARAMS8_8+32)>>6;
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
   t[8]=(t[8]*OD_FILTER_PARAMS16_0)>>6;
   t[8]+=-t[8]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_1!=64
   t[9]=(t[9]*OD_FILTER_PARAMS16_1)>>6;
   t[9]+=-t[9]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_2!=64
   t[10]=(t[10]*OD_FILTER_PARAMS16_2)>>6;
   t[10]+=-t[10]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_3!=64
   t[11]=(t[11]*OD_FILTER_PARAMS16_3)>>6;
   t[11]+=-t[11]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_4!=64
   t[12]=(t[12]*OD_FILTER_PARAMS16_4)>>6;
   t[12]+=-t[12]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_5!=64
   t[13]=(t[13]*OD_FILTER_PARAMS16_5)>>6;
   t[13]+=-t[13]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_6!=64
   t[14]=t[14]*OD_FILTER_PARAMS16_6>>6;
   t[14]+=-t[14]>>(OD_COEFF_BITS-1)&1;
#endif
#if OD_FILTER_PARAMS16_7!=64
   t[15]=t[15]*OD_FILTER_PARAMS16_7>>6;
   t[15]+=-t[15]>>(OD_COEFF_BITS-1)&1;
#endif
   /*Rotations:*/
#if OD_FILTER16_TYPE3
   t[15]+=(t[14]*OD_FILTER_PARAMS16_14+32)>>6;
   t[14]+=(t[15]*OD_FILTER_PARAMS16_21+32)>>6;
   t[14]+=(t[13]*OD_FILTER_PARAMS16_13+32)>>6;
   t[13]+=(t[14]*OD_FILTER_PARAMS16_20+32)>>6;
   t[13]+=(t[12]*OD_FILTER_PARAMS16_12+32)>>6;
   t[12]+=(t[13]*OD_FILTER_PARAMS16_19+32)>>6;
   t[12]+=(t[11]*OD_FILTER_PARAMS16_11+32)>>6;
   t[11]+=(t[12]*OD_FILTER_PARAMS16_18+32)>>6;
   t[11]+=(t[10]*OD_FILTER_PARAMS16_10+32)>>6;
   t[10]+=(t[11]*OD_FILTER_PARAMS16_17+32)>>6;
   t[10]+=(t[9]*OD_FILTER_PARAMS16_9+32)>>6;
   t[9]+=(t[10]*OD_FILTER_PARAMS16_16+32)>>6;
   t[9]+=(t[8]*OD_FILTER_PARAMS16_8+32)>>6;
   t[8]+=(t[9]*OD_FILTER_PARAMS16_15+32)>>6;
#else
   t[9]+=(t[8]*OD_FILTER_PARAMS16_8+32)>>6;
   t[10]+=(t[9]*OD_FILTER_PARAMS16_9+32)>>6;
   t[11]+=(t[10]*OD_FILTER_PARAMS16_10+32)>>6;
   t[12]+=(t[11]*OD_FILTER_PARAMS16_11+32)>>6;
   t[13]+=(t[12]*OD_FILTER_PARAMS16_12+32)>>6;
   t[14]+=(t[13]*OD_FILTER_PARAMS16_13+32)>>6;
   t[15]+=(t[14]*OD_FILTER_PARAMS16_14+32)>>6;
   t[14]+=(t[15]*OD_FILTER_PARAMS16_21+32)>>6;
   t[13]+=(t[14]*OD_FILTER_PARAMS16_20+32)>>6;
   t[12]+=(t[13]*OD_FILTER_PARAMS16_19+32)>>6;
   t[11]+=(t[12]*OD_FILTER_PARAMS16_18+32)>>6;
   t[10]+=(t[11]*OD_FILTER_PARAMS16_17+32)>>6;
   t[9]+=(t[10]*OD_FILTER_PARAMS16_16+32)>>6;
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

#define ZERO_FILTERS (0)

/*Remove OD_MINI if we ever support 32-point filters.*/
#define OD_BLOCK_SIZE4x4_DEC(bsize, bstride, bx, by, dec) \
 OD_MAXI(OD_MINI(OD_BLOCK_SIZE4x4(bsize, bstride, bx, by), 2 + (dec)), dec)

static void od_apply_filter_cols(od_coeff *c, int stride, int bx, int by,
 int out, unsigned char l, const unsigned char *bsize, int bstride, int xdec,
 int ydec, int inv) {
  int nbx;
  int nby;
  unsigned char n;
  int f;
  int i;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  /*Compute the neighbor's block size.*/
  nbx = bx << l;
  nby = ((by + out) << l) + (out - 1);
  n = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, nbx, nby, xdec);
  /*Assume we use the filter for the current block's size.*/
  f = l;
  /*If this block is larger, search *my* neighbors for smallest block size.*/
  if (l > n) {
    for (i = 0; i < 1 << l; i++) {
      f = OD_MINI(f, OD_BLOCK_SIZE4x4_DEC(bsize, bstride, nbx + i, nby, ydec));
    }
  }
  /*If neighbor is larger, search *its* neighbors for smallest block size.*/
  if (l < n) {
    /*Compute the start of the edge in the neighboring block.*/
    nbx &= ~((1 << n) - 1);
    nby += out ? -1 : 1;
    for (i = 0; i < 1 << n; i++) {
      f = OD_MINI(f, OD_BLOCK_SIZE4x4_DEC(bsize, bstride, nbx + i, nby, ydec));
    }
  }
  /*Block size in c is l - xdec.*/
  l -= xdec;
  /*Filter size in c is f - xdec.*/
  f -= xdec;
  OD_ASSERT(0 <= f && f <= OD_NBSIZES);
  c += (((by + out) << (l + 2)) - (2 << f))*stride + (bx << (l + 2));
  /*Apply the column filter across the edge.*/
  for (i = 0; i < 4 << l; i++) {
    int j;
    od_coeff t[4 << OD_NBSIZES];
    for (j = 0; j < 4 << f; j++) {
      t[j] = c[stride*j + i];
    }
    (*(inv ? OD_POST_FILTER : OD_PRE_FILTER)[f])(t, t);
    for (j = 0; j < 4 << f; j++) {
      c[stride*j + i] = t[j];
#if ZERO_FILTERS
      c[stride*j + i] = 0;
#endif
    }
  }
}

static void od_apply_filter_rows(od_coeff *c, int stride, int bx, int by,
 int out, unsigned char l, const unsigned char *bsize, int bstride, int xdec,
 int ydec, int inv) {
  int nbx;
  int nby;
  unsigned char n;
  int f;
  int i;
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  /*Compute the neighbor's block size.*/
  nbx = ((bx + out) << l) + (out - 1);
  nby = by << l;
  n = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, nbx, nby, xdec);
  /*Assume we use the filter for the current block's size.*/
  f = l;
  /*If this block is larger, search *my* neighbors for smallest block size.*/
  if (l > n) {
    for (i = 0; i < 1 << l; i++) {
      f = OD_MINI(f, OD_BLOCK_SIZE4x4_DEC(bsize, bstride, nbx, nby + i, ydec));
    }
  }
  /*If neighbor is larger, search *its* neighbors for smallest block size.*/
  if (l < n) {
    /*Compute the start of the edge in the neighboring block.*/
    nbx += out ? -1 : 1;
    nby &= ~((1 << n) - 1);
    for (i = 0; i < 1 << n; i++) {
      f = OD_MINI(f, OD_BLOCK_SIZE4x4_DEC(bsize, bstride, nbx, nby + i, ydec));
    }
  }
  /*Block size in c is l - xdec.*/
  l -= xdec;
  /*Filter size in c is f - xdec.*/
  f -= xdec;
  OD_ASSERT(0 <= f && f <= OD_NBSIZES);
  c += (by << (l + 2))*stride + ((bx + out) << (l + 2)) - (2 << f);
  /* Apply the row filter down the edge. */
  for (i = 0; i < 4 << l; i++) {
    (*(inv ? OD_POST_FILTER : OD_PRE_FILTER)[f])(c, c);
#if ZERO_FILTERS
    {
      int j;
      for (j = 0; j < 4 << f; j++) {
        c[j] = 0;
      }
    }
#endif
    c += stride;
  }
}

void od_apply_prefilter(od_coeff *c, int w, int bx, int by, unsigned char l,
 const unsigned char *bsize, int bstride, int xdec, int ydec, int edge) {
  unsigned d;
  OD_ASSERT((edge & ~(OD_BOTTOM_EDGE | OD_LEFT_EDGE)) == 0);
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  d = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, bx << l, by << l, xdec);
  OD_ASSERT(d <= l);
  if (d == l) {
    if (edge & OD_BOTTOM_EDGE) {
      od_apply_filter_cols(c, w, bx, by, 1, l, bsize, bstride, xdec, ydec, 0);
    }
    if (edge & OD_LEFT_EDGE) {
      od_apply_filter_rows(c, w, bx, by, 0, l, bsize, bstride, xdec, ydec, 0);
    }
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_apply_prefilter(c, w, bx + 0, by + 0, l, bsize, bstride, xdec, ydec,
     edge | OD_BOTTOM_EDGE);
    od_apply_prefilter(c, w, bx + 1, by + 0, l, bsize, bstride, xdec, ydec,
     OD_BOTTOM_EDGE | OD_LEFT_EDGE);
    od_apply_prefilter(c, w, bx + 0, by + 1, l, bsize, bstride, xdec, ydec,
     edge);
    od_apply_prefilter(c, w, bx + 1, by + 1, l, bsize, bstride, xdec, ydec,
     edge | OD_LEFT_EDGE);
  }
}

void od_apply_postfilter(od_coeff *c, int w, int bx, int by, unsigned char l,
 const unsigned char *bsize, int bstride, int xdec, int ydec, int edge) {
  unsigned char d;
  OD_ASSERT((edge & ~(OD_RIGHT_EDGE | OD_TOP_EDGE)) == 0);
  /*This code assumes 4:4:4 or 4:2:0 input.*/
  OD_ASSERT(xdec == ydec);
  d = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, bx << l, by << l, xdec);
  OD_ASSERT(d <= l);
  if (d == l) {
    if (edge & OD_RIGHT_EDGE) {
      od_apply_filter_rows(c, w, bx, by, 1, l, bsize, bstride, xdec, ydec, 1);
    }
    if (edge & OD_TOP_EDGE) {
      od_apply_filter_cols(c, w, bx, by, 0, l, bsize, bstride, xdec, ydec, 1);
    }
  }
  else {
    l--;
    bx <<= 1;
    by <<= 1;
    od_apply_postfilter(c, w, bx + 0, by + 0, l, bsize, bstride, xdec, ydec,
     edge | OD_RIGHT_EDGE);
    od_apply_postfilter(c, w, bx + 1, by + 0, l, bsize, bstride, xdec, ydec,
     edge);
    od_apply_postfilter(c, w, bx + 0, by + 1, l, bsize, bstride, xdec, ydec,
     OD_TOP_EDGE | OD_RIGHT_EDGE);
    od_apply_postfilter(c, w, bx + 1, by + 1, l, bsize, bstride, xdec, ydec,
     edge | OD_TOP_EDGE);
  }
}

#if defined(TEST)
# include <stdio.h>

int minv[32];
int maxv[32];

int main(void) {
  ogg_int32_t min[16];
  ogg_int32_t max[16];
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
      min[j] = mini[j] = 2147483647;
      max[j] = maxi[j] = -2147483647-1;
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
