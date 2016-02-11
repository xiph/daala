/*Daala video codec
Copyright (c) 2002-2013 Daala project contributors.  All rights reserved.

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

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "block_size.h"
#include "dct.h"
#include "tf.h"

/*This flag controls the use of embedded DCTs.
  Each n-point orthonormal Type-II DCT can be decomposed via a stage of
   "butterfly" operations into an n/2-point asymmetric Type-II DCT and an
   n/2-point asymmetric Type-IV DST.
  Similarly, each n-point asymmetric Type-II DCT can be decomposed via a stage
   of "butterfly" operations into an n/2-point orthonormal Type-II DCT and an
   n/2-point orthonormal Type-IV DST.
  When this flag is enabled, the DCTs used are exactly the embedded 4-point and
   16-point orthonormal Type-II DCTs from the 64-point transform, and the
   embedded 8-point orthonormal Type-II DCT from the 32-point transform.*/
#define OD_USE_EMBEDDED_DCTS (0)

/*Making function pointer tables at least one entry
   longer than needed makes it highly likely that an
   off-by-one will result in a null-pointer rather than
   another otherwise compatible function pointer.
  This can help avoid difficult to diagnose misbehavior.*/

/*Function tables suffixed with _C are for generic implementations.
  Code should use the tables in od_state.opt_vtbl to get optimized
   implementations when they are available.*/
const od_dct_func_2d OD_FDCT_2D_C[OD_NBSIZES + 1] = {
  od_bin_fdct4x4,
  od_bin_fdct8x8,
  od_bin_fdct16x16,
  od_bin_fdct32x32,
  od_bin_fdct64x64
};

const od_dct_func_2d OD_IDCT_2D_C[OD_NBSIZES + 1] = {
  od_bin_idct4x4,
  od_bin_idct8x8,
  od_bin_idct16x16,
  od_bin_idct32x32,
  od_bin_idct64x64
};

const od_fdct_func_1d OD_FDCT_1D[OD_NBSIZES + 1] = {
  od_bin_fdct4,
  od_bin_fdct8,
  od_bin_fdct16,
  od_bin_fdct32,
  od_bin_fdct64
};

const od_idct_func_1d OD_IDCT_1D[OD_NBSIZES + 1] = {
  od_bin_idct4,
  od_bin_idct8,
  od_bin_idct16,
  od_bin_idct32,
  od_bin_idct64
};

#if !OD_USE_EMBEDDED_DCTS
void od_bin_fdct4(od_coeff y[4], const od_coeff *x, int xstride) {
  /*9 adds, 2 shifts, 3 "muls".*/
  int t0;
  int t1;
  int t2;
  int t2h;
  int t3;
  /*Initial permutation:*/
  t0 = *(x + 0*xstride);
  t2 = *(x + 1*xstride);
  t1 = *(x + 2*xstride);
  t3 = *(x + 3*xstride);
  /*+1/-1 butterflies:*/
  t3 = t0 - t3;
  t2 += t1;
  t2h = OD_DCT_RSHIFT(t2, 1);
  t1 = t2h - t1;
  t0 -= OD_DCT_RSHIFT(t3, 1);
  /*+ Embedded 2-point type-II DCT.*/
  t0 += t2h;
  t2 = t0 - t2;
  /*+ Embedded 2-point type-IV DST.*/
  /*23013/32768 ~= 4*sin(\frac{\pi}{8}) - 2*tan(\frac{\pi}{8}) ~=
     0.70230660471416898931046248770220*/
  OD_DCT_OVERFLOW_CHECK(t1, 23013, 16384, 0);
  t3 -= (t1*23013 + 16384) >> 15;
  /*21407/32768~=\sqrt{1/2}*cos(\frac{\pi}{8}))
     ~=0.65328148243818826392832158671359*/
  OD_DCT_OVERFLOW_CHECK(t3, 21407, 16384, 1);
  t1 += (t3*21407 + 16384) >> 15;
  /*18293/16384 ~= 4*sin(\frac{\pi}{8}) - tan(\frac{\pi}{8}) ~=
     1.1165201670872640381121512119119*/
  OD_DCT_OVERFLOW_CHECK(t1, 18293, 8192, 2);
  t3 -= (t1*18293 + 8192) >> 14;
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
}

void od_bin_idct4(od_coeff *x, int xstride, const od_coeff y[4]) {
  int t0;
  int t1;
  int t2;
  int t2h;
  int t3;
  t0 = y[0];
  t1 = y[1];
  t2 = y[2];
  t3 = y[3];
  t3 += (t1*18293 + 8192) >> 14;
  t1 -= (t3*21407 + 16384) >> 15;
  t3 += (t1*23013 + 16384) >> 15;
  t2 = t0 - t2;
  t2h = OD_DCT_RSHIFT(t2, 1);
  t0 -= t2h - OD_DCT_RSHIFT(t3, 1);
  t1 = t2h - t1;
  *(x + 0*xstride) = (od_coeff)t0;
  *(x + 1*xstride) = (od_coeff)(t2 - t1);
  *(x + 2*xstride) = (od_coeff)t1;
  *(x + 3*xstride) = (od_coeff)(t0 - t3);
}
#endif

void od_bin_fdct4x4(od_coeff *y, int ystride, const od_coeff *x, int xstride) {
  od_coeff z[4*4];
  int i;
  for (i = 0; i < 4; i++) od_bin_fdct4(z + 4*i, x + i, xstride);
  for (i = 0; i < 4; i++) od_bin_fdct4(y + ystride*i, z + i, 4);
}

void od_bin_idct4x4(od_coeff *x, int xstride, const od_coeff *y, int ystride) {
  od_coeff z[4*4];
  int i;
  for (i = 0; i < 4; i++) od_bin_idct4(z + i, 4, y + ystride*i);
  for (i = 0; i < 4; i++) od_bin_idct4(x + i, xstride, z + 4*i);
}

#if !OD_USE_EMBEDDED_DCTS
void od_bin_fdct8(od_coeff y[8], const od_coeff *x, int xstride) {
  /*31 adds, 5 shifts, 15 "muls".*/
  /*The minimum theoretical number of multiplies for a uniformly-scaled 8-point
     transform is 11, but the best I've been able to come up with for a
     reversible version with orthonormal scaling is 15.
    We pick up 3 multiplies when computing the DC, since we have an odd number
     of summation stages, leaving no chance to cancel the asymmetry in the last
     one.
    Instead, we have to implement it as a rotation by \frac{\pi}{4} using
     lifting steps.
    We pick up one more multiply when computing the Type IV DCT in the odd
     half.
    This comes from using 3 lifting steps to implement another rotation by
     \frac{\pi}{4} (with asymmetrically scaled inputs and outputs) instead of
     simply scaling two values by \sqrt{2}.*/
  int t0;
  int t1;
  int t1h;
  int t2;
  int t3;
  int t4;
  int t4h;
  int t5;
  int t6;
  int t6h;
  int t7;
  /*Initial permutation:*/
  t0 = *(x + 0*xstride);
  t4 = *(x + 1*xstride);
  t2 = *(x + 2*xstride);
  t6 = *(x + 3*xstride);
  t7 = *(x + 4*xstride);
  t3 = *(x + 5*xstride);
  t5 = *(x + 6*xstride);
  t1 = *(x + 7*xstride);
  /*+1/-1 butterflies:*/
  t1 = t0 - t1;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t0 -= t1h;
  t4 += t5;
  t4h = OD_DCT_RSHIFT(t4, 1);
  t5 -= t4h;
  t3 = t2 - t3;
  t2 -= OD_DCT_RSHIFT(t3, 1);
  t6 += t7;
  t6h = OD_DCT_RSHIFT(t6, 1);
  t7 = t6h - t7;
  /*+ Embedded 4-point type-II DCT.*/
  t0 += t6h;
  t6 = t0 - t6;
  t2 = t4h - t2;
  t4 = t2 - t4;
  /*|-+ Embedded 2-point type-II DCT.*/
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(t4, 13573, 16384, 3);
  t0 -= (t4*13573 + 16384) >> 15;
  /*11585/16384 ~= \sqrt{\frac{1}{2}} ~= 0.70710678118654752440084436210485*/
  OD_DCT_OVERFLOW_CHECK(t0, 11585, 8192, 4);
  t4 += (t0*11585 + 8192) >> 14;
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(t4, 13573, 16384, 5);
  t0 -= (t4*13573 + 16384) >> 15;
  /*|-+ Embedded 2-point type-IV DST.*/
  /*21895/32768 ~= \frac{1 - cos(\frac{3\pi}{8})}{\sin(\frac{3\pi}{8})} ~=
     0.66817863791929891999775768652308*/
  OD_DCT_OVERFLOW_CHECK(t2, 21895, 16384, 6);
  t6 -= (t2*21895 + 16384) >> 15;
  /*15137/16384~=sin(\frac{3\pi}{8})~=0.92387953251128675612818318939679*/
  OD_DCT_OVERFLOW_CHECK(t6, 15137, 8192, 7);
  t2 += (t6*15137 + 8192) >> 14;
  /*21895/32768 ~= \frac{1 - cos(\frac{3\pi}{8})}{\sin(\frac{3\pi}{8})}~=
     0.66817863791929891999775768652308*/
  OD_DCT_OVERFLOW_CHECK(t2, 21895, 16384, 8);
  t6 -= (t2*21895 + 16384) >> 15;
  /*+ Embedded 4-point type-IV DST.*/
  /*19195/32768 ~= 2 - \sqrt{2} ~= 0.58578643762690495119831127579030*/
  OD_DCT_OVERFLOW_CHECK(t5, 19195, 16384, 9);
  t3 += (t5*19195 + 16384) >> 15;
  /*11585/16384 ~= \sqrt{\frac{1}{2}} ~= 0.70710678118654752440084436210485*/
  OD_DCT_OVERFLOW_CHECK(t3, 11585, 8192, 10);
  t5 += (t3*11585 + 8192) >> 14;
  /*7489/8192 ~= \sqrt{2}-\frac{1}{2} ~= 0.91421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(t5, 7489, 4096, 11);
  t3 -= (t5*7489 + 4096) >> 13;
  t7 = OD_DCT_RSHIFT(t5, 1) - t7;
  t5 -= t7;
  t3 = t1h - t3;
  t1 -= t3;
  /*3227/32768 ~= \frac{1 - cos(\frac{\pi}{16})}{sin(\frac{\pi}{16})} ~=
     0.098491403357164253077197521291327*/
  OD_DCT_OVERFLOW_CHECK(t1, 3227, 16384, 12);
  t7 += (t1*3227 + 16384) >> 15;
  /*6393/32768 ~= sin(\frac{\pi}{16}) ~= 0.19509032201612826784828486847702*/
  OD_DCT_OVERFLOW_CHECK(t7, 6393, 16384, 13);
  t1 -= (t7*6393 + 16384) >> 15;
  /*3227/32768 ~= \frac{1 - cos(\frac{\pi}{16})}{sin(\frac{\pi}{16})} ~=
     0.098491403357164253077197521291327*/
  OD_DCT_OVERFLOW_CHECK(t1, 3227, 16384, 14);
  t7 += (t1*3227 + 16384) >> 15;
  /*2485/8192 ~= \frac{1 - cos(\frac{3\pi}{16})}{sin(\frac{3\pi}{16})} ~=
     0.30334668360734239167588394694130*/
  OD_DCT_OVERFLOW_CHECK(t3, 2485, 4096, 15);
  t5 += (t3*2485 + 4096) >> 13;
  /*18205/32768 ~= sin(\frac{3\pi}{16}) ~= 0.55557023301960222474283081394853*/
  OD_DCT_OVERFLOW_CHECK(t5, 18205, 16384, 16);
  t3 -= (t5*18205 + 16384) >> 15;
  /*2485/8192 ~= \frac{1 - cos(\frac{3\pi}{16})}{sin(\frac{3\pi}{16})} ~=
     0.30334668360734239167588394694130*/
  OD_DCT_OVERFLOW_CHECK(t3, 2485, 4096, 17);
  t5 += (t3*2485 + 4096) >> 13;
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
  y[4] = (od_coeff)t4;
  y[5] = (od_coeff)t5;
  y[6] = (od_coeff)t6;
  y[7] = (od_coeff)t7;
}

void od_bin_idct8(od_coeff *x, int xstride, const od_coeff y[16]) {
  int t0;
  int t1;
  int t1h;
  int t2;
  int t3;
  int t4;
  int t4h;
  int t5;
  int t6;
  int t6h;
  int t7;
  t0 = y[0];
  t1 = y[1];
  t2 = y[2];
  t3 = y[3];
  t4 = y[4];
  t5 = y[5];
  t6 = y[6];
  t7 = y[7];
  t5 -= (t3*2485 + 4096) >> 13;
  t3 += (t5*18205 + 16384) >> 15;
  t5 -= (t3*2485 + 4096) >> 13;
  t7 -= (t1*3227 + 16384) >> 15;
  t1 += (t7*6393 + 16384) >> 15;
  t7 -= (t1*3227 + 16384) >> 15;
  t1 += t3;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t3 = t1h - t3;
  t5 += t7;
  t7 = OD_DCT_RSHIFT(t5, 1) - t7;
  t3 += (t5*7489 + 4096) >> 13;
  t5 -= (t3*11585 + 8192) >> 14;
  t3 -= (t5*19195 + 16384) >> 15;
  t6 += (t2*21895 + 16384) >> 15;
  t2 -= (t6*15137 + 8192) >> 14;
  t6 += (t2*21895 + 16384) >> 15;
  t0 += (t4*13573 + 16384) >> 15;
  t4 -= (t0*11585 + 8192) >> 14;
  t0 += (t4*13573 + 16384) >> 15;
  t4 = t2 - t4;
  t4h = OD_DCT_RSHIFT(t4, 1);
  t2 = t4h - t2;
  t6 = t0 - t6;
  t6h = OD_DCT_RSHIFT(t6, 1);
  t0 -= t6h;
  t7 = t6h - t7;
  t6 -= t7;
  t2 += OD_DCT_RSHIFT(t3, 1);
  t3 = t2 - t3;
  t5 += t4h;
  t4 -= t5;
  t0 += t1h;
  t1 = t0 - t1;
  *(x + 0*xstride) = (od_coeff)t0;
  *(x + 1*xstride) = (od_coeff)t4;
  *(x + 2*xstride) = (od_coeff)t2;
  *(x + 3*xstride) = (od_coeff)t6;
  *(x + 4*xstride) = (od_coeff)t7;
  *(x + 5*xstride) = (od_coeff)t3;
  *(x + 6*xstride) = (od_coeff)t5;
  *(x + 7*xstride) = (od_coeff)t1;
}
#endif

void od_bin_fdct8x8(od_coeff *y, int ystride, const od_coeff *x, int xstride) {
  od_coeff z[8*8];
  int i;
  for (i = 0; i < 8; i++) od_bin_fdct8(z + 8*i, x + i, xstride);
  for (i = 0; i < 8; i++) od_bin_fdct8(y + ystride*i, z + i, 8);
}

void od_bin_idct8x8(od_coeff *x, int xstride, const od_coeff *y, int ystride) {
  od_coeff z[8*8];
  int i;
  for (i = 0; i < 8; i++) od_bin_idct8(z + i, 8, y + ystride*i);
  for (i = 0; i < 8; i++) od_bin_idct8(x + i, xstride, z + 8*i);
}

#if !OD_USE_EMBEDDED_DCTS
void od_bin_fdct16(od_coeff y[16], const od_coeff *x, int xstride) {
  /*83 adds, 16 shifts, 33 "muls".*/
  /*The minimum theoretical number of multiplies is 26~\cite{DH87}, but the
     best practical algorithm I know is 31~\cite{LLM89}.
    This is a modification of the Loeffler et al. factorization that allows us
     to have a reversible integer transform with true orthonormal scaling.
    This required some major reworking of the odd quarter of the even half
     (the 4-point Type IV DST), and added two multiplies in order to implement
     two rotations by \frac{\pi}{4} with lifting steps (requiring 3 multiplies
     instead of the normal 2).
    @INPROCEEDINGS{DH87,
      author={Pierre Duhamel and Hedi H'Mida},
      title="New 2^n {DCT} Algorithms Suitable for {VLSI} Implementation",
      booktitle="Proc. $12^\textrm{th}$ International Conference on Acoustics,
       Speech, and Signal Processing (ICASSP'87)",
      volume=12,
      pages="1805--1808",
      address="Issy-les-Moulineaux, France",
      month=Apr,
      year=1987
    }
    @INPROCEEDINGS{LLM89,
      author="Christoph Loeffler and Adriaan Lightenberg and
       George S. Moschytz",
      title="Practical Fast {1-D} {DCT} Algorithms with 11 Multiplications",
      booktitle="Proc. $14^\textrm{th}$ International Conference on Acoustics,
       Speech, and Signal Processing (ICASSP'89)",
      volume=2,
      pages="988--991",
      address="Zhurich, Switzerland",
      month=May,
      year=1989
    }*/
  int t0;
  int t1;
  int t1h;
  int t2;
  int t2h;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  int t8;
  int t8h;
  int t9;
  int ta;
  int tah;
  int tb;
  int tbh;
  int tc;
  int tch;
  int td;
  int tdh;
  int te;
  int tf;
  int tfh;
  /*Initial permutation:*/
  t0 = *(x + 0*xstride);
  t8 = *(x + 1*xstride);
  t4 = *(x + 2*xstride);
  tc = *(x + 3*xstride);
  te = *(x + 4*xstride);
  ta = *(x + 5*xstride);
  t6 = *(x + 6*xstride);
  t2 = *(x + 7*xstride);
  t3 = *(x + 8*xstride);
  td = *(x + 9*xstride);
  t9 = *(x + 10*xstride);
  tf = *(x + 11*xstride);
  t1 = *(x + 12*xstride);
  t7 = *(x + 13*xstride);
  tb = *(x + 14*xstride);
  t5 = *(x + 15*xstride);
  /*+1/-1 butterflies:*/
  t5 = t0 - t5;
  t8 += tb;
  t7 = t4 - t7;
  tc += t1;
  tf = te - tf;
  ta += t9;
  td = t6 - td;
  t2 += t3;
  t0 -= OD_DCT_RSHIFT(t5, 1);
  t8h = OD_DCT_RSHIFT(t8, 1);
  tb = t8h - tb;
  t4 -= OD_DCT_RSHIFT(t7, 1);
  tch = OD_DCT_RSHIFT(tc, 1);
  t1 = tch - t1;
  te -= OD_DCT_RSHIFT(tf, 1);
  tah = OD_DCT_RSHIFT(ta, 1);
  t9 = tah - t9;
  t6 -= OD_DCT_RSHIFT(td, 1);
  t2h = OD_DCT_RSHIFT(t2, 1);
  t3 = t2h - t3;
  /*+ Embedded 8-point type-II DCT.*/
  t0 += t2h;
  t6 = t8h - t6;
  t4 += tah;
  te = tch - te;
  t2 = t0 - t2;
  t8 -= t6;
  ta = t4 - ta;
  tc -= te;
  /*|-+ Embedded 4-point type-II DCT.*/
  tc = t0 - tc;
  t8 += t4;
  t8h = OD_DCT_RSHIFT(t8, 1);
  t4 = t8h - t4;
  t0 -= OD_DCT_RSHIFT(tc, 1);
  /*|-|-+ Embedded 2-point type-II DCT.*/
  t0 += t8h;
  t8 = t0 - t8;
  /*|-|-+ Embedded 2-point type-IV DST.*/
  /*32013/32768 ~= 4*sin(\frac{\pi}{8}) - 2*tan(\frac{\pi}{8}) ~=
     0.70230660471416898931046248770220*/
  OD_DCT_OVERFLOW_CHECK(t4, 23013, 16384, 18);
  tc -= (t4*23013 + 16384) >> 15;
  /*10703/16384 ~= \sqrt{1/2}*cos(\frac{\pi}{8})) ~=
     0.65328148243818826392832158671359*/
  OD_DCT_OVERFLOW_CHECK(tc, 10703, 8192, 19);
  t4 += (tc*10703 + 8192) >> 14;
  /*9147/8192 ~= 4*sin(\frac{\pi}{8}) - tan(\frac{\pi}{8}) ~=
     1.1165201670872640381121512119119*/
  OD_DCT_OVERFLOW_CHECK(t4, 9147, 4096, 20);
  tc -= (t4*9147 + 4096) >> 13;
  /*|-+ Embedded 4-point type-IV DST.*/
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(ta, 13573, 16384, 21);
  t6 += (ta*13573 + 16384) >> 15;
  /*11585/16384 ~= \sqrt{\frac{1}{2}} ~= 0.70710678118654752440084436210485*/
  OD_DCT_OVERFLOW_CHECK(t6, 11585, 8192, 22);
  ta -= (t6*11585 + 8192) >> 14;
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(ta, 13573, 16384, 23);
  t6 += (ta*13573 + 16384) >> 15;
  ta += te;
  t2 += t6;
  te = OD_DCT_RSHIFT(ta, 1) - te;
  t6 = OD_DCT_RSHIFT(t2, 1) - t6;
  /*2775/2048 ~= \frac{\sqrt{2} - cos(\frac{\pi}{16})}{2sin(\frac{\pi}{16})}
     ~= 1.1108400393486273201524536919723*/
  OD_DCT_OVERFLOW_CHECK(t2, 2275, 1024, 24);
  te += (t2*2275 + 1024) >> 11;
  /*9041/32768 ~= \sqrt{2}sin(\frac{\pi}{16}) ~=
     0.27589937928294301233595756366937*/
  OD_DCT_OVERFLOW_CHECK(te, 9041, 16384, 25);
  t2 -= (te*9041 + 16384) >> 15;
  /*2873/2048 ~=
     \frac{cos(\frac{\pi}{16}) - \sqrt{\frac{1}{2}}}{sin(\frac{\pi}{16})} ~=
     1.4028297067142967321050338435598*/
  OD_DCT_OVERFLOW_CHECK(t2, 2873, 1024, 26);
  te -= (t2*2873 + 1024) >> 11;
  /*8593/16384 ~=
    \frac{\sqrt{2} - cos(\frac{3\pi}{16})}{2sin(\frac{3\pi}{16})} ~=
    0.52445569924008942966043945081053*/
  OD_DCT_OVERFLOW_CHECK(ta, 8593, 8192, 27);
  t6 -= (ta*8593 + 8192) >> 14;
  /*12873/16384 ~= \sqrt{2}sin(\frac{3\pi}{16}) ~=
     0.78569495838710218127789736765722*/
  OD_DCT_OVERFLOW_CHECK(t6, 12873, 8192, 28);
  ta += (t6*12873 + 8192) >> 14;
  /*7335/32768
     ~=\frac{cos(\frac{3\pi}{16})-\sqrt{\frac{1}{2}}}{sin(\frac{3\pi}{16})}
     ~=0.22384718209265507914012811666071*/
  OD_DCT_OVERFLOW_CHECK(ta, 7335, 16384, 29);
  t6 += (ta*7335 + 16384) >> 15;
  /*+ Embedded 8-point type-IV DST.*/
  /*1035/2048 ~=
    \frac{\sqrt{2} - cos(\frac{7\pi}{32})}{2sin(\frac{7\pi}{32})} ~=
    0.50536719493782972897642806316664*/
  OD_DCT_OVERFLOW_CHECK(t5, 1035, 1024, 30);
  t3 += (t5*1035 + 1024) >> 11;
  /*14699/16384 ~= \sqrt{2}sin(\frac{7\pi}{32}) ~=
     0.89716758634263628425064138922084*/
  OD_DCT_OVERFLOW_CHECK(t3, 14699, 8192, 31);
  t5 -= (t3*14699 + 8192) >> 14;
  /*851/8192 ~=
    \frac{cos(\frac{7\pi}{32}) - \sqrt{\frac{1}{2}}}{sin(\frac{7\pi}{32}} ~=
    0.10388456785615844342131055214354*/
  OD_DCT_OVERFLOW_CHECK(t5, 851, 4096, 32);
  t3 -= (t5*851 + 4096) >> 13;
  /*17515/32768 ~=
     \frac{\sqrt{2} - cos(\frac{11\pi}{32})}{2sin(\frac{11\pi}{32})} ~=
     0.53452437516842143578098634302964*/
  OD_DCT_OVERFLOW_CHECK(td, 17515, 16384, 33);
  tb += (td*17515 + 16384) >> 15;
  /*20435/16384 ~= \sqrt{2}sin(\frac{11\pi}{32}) ~=
     1.2472250129866712325719902627977*/
  OD_DCT_OVERFLOW_CHECK(tb, 20435, 8192, 34);
  td -= (tb*20435 + 8192) >> 14;
  /*4379/16384 ~= \frac{\sqrt{\frac{1}{2}}
     - cos(\frac{11\pi}{32})}{sin(\frac{11\pi}{32})} ~=
     0.26726880719302561523614336238196*/
  OD_DCT_OVERFLOW_CHECK(td, 4379, 8192, 35);
  tb += (td*4379 + 8192) >> 14;
  /*12905/16384 ~=
     \frac{\sqrt{2} - cos(\frac{3\pi}{32})}{2sin(\frac{3\pi}{32})} ~=
     0.78762894232967441973847776796517*/
  OD_DCT_OVERFLOW_CHECK(t7, 12905, 8192, 36);
  t9 += (t7*12905 + 8192) >> 14;
  /*3363/8192 ~= \sqrt{2}sin(\frac{3\pi}{32}) ~=
     0.41052452752235738115636923775513*/
  OD_DCT_OVERFLOW_CHECK(t9, 3363, 4096, 37);
  t7 -= (t9*3363 + 4096) >> 13;
  /*14101/16384 ~=
     \frac{cos(\frac{3\pi}{32}) - \sqrt{\frac{1}{2}}}{sin(\frac{3\pi}{32})} ~=
     0.86065016213948579370059934044795*/
  OD_DCT_OVERFLOW_CHECK(t7, 14101, 8192, 38);
  t9 -= (t7*14101 + 8192) >> 14;
  /*5417/8192 ~=
     \frac{\sqrt{2} - cos(\frac{15\pi}{32})}{2sin(\frac{15\pi}{32})} ~=
     0.66128246684651710406296283785232*/
  OD_DCT_OVERFLOW_CHECK(tf, 5417, 4096, 39);
  t1 += (tf*5417 + 4096) >> 13;
  /*23059/16384 ~= \sqrt{2}sin(\frac{15\pi}{32}) ~=
     1.4074037375263824590260782229840*/
  OD_DCT_OVERFLOW_CHECK(t1, 23059, 8192, 40);
  tf -= (t1*23059 + 8192) >> 14;
  /*20055/32768 ~=
    \frac{\sqrt{\frac{1}{2}} - cos(\frac{15\pi}{32})}{sin(\frac{15\pi}{32})} ~=
    0.61203676516793497752436407720666*/
  OD_DCT_OVERFLOW_CHECK(tf, 20055, 16384, 41);
  t1 += (tf*20055 + 16384) >> 15;
  tf = t3 - tf;
  td += t9;
  tfh = OD_DCT_RSHIFT(tf, 1);
  t3 -= tfh;
  tdh = OD_DCT_RSHIFT(td, 1);
  t9 = tdh - t9;
  t1 += t5;
  tb = t7 - tb;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t5 = t1h - t5;
  tbh = OD_DCT_RSHIFT(tb, 1);
  t7 -= tbh;
  t3 += tbh;
  t5 = tdh - t5;
  t9 += tfh;
  t7 = t1h - t7;
  tb -= t3;
  td -= t5;
  tf = t9 - tf;
  t1 -= t7;
  /*10947/16384 ~= \frac{1 - cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})} ~=
     0.66817863791929891999775768652308*/
  OD_DCT_OVERFLOW_CHECK(tb, 10947, 8192, 42);
  t5 -= (tb*10947 + 8192) >> 14;
  /*15137/16384 ~= sin(\frac{3\pi}{8}) ~= 0.92387953251128675612818318939679*/
  OD_DCT_OVERFLOW_CHECK(t5, 15137, 8192, 43);
  tb += (t5*15137 + 8192) >> 14;
  /*10947/16384 ~= \frac{1 - cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})} ~=
     0.66817863791929891999775768652308*/
  OD_DCT_OVERFLOW_CHECK(tb, 10947, 8192, 44);
  t5 -= (tb*10947 + 8192) >> 14;
  /*21895/32768 ~= \frac{1 - cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})} ~=
     0.66817863791929891999775768652308*/
  OD_DCT_OVERFLOW_CHECK(t3, 21895, 16384, 45);
  td += (t3*21895 + 16384) >> 15;
  /*15137/16384 ~= sin(\frac{3\pi}{8}) ~= 0.92387953251128675612818318939679*/
  OD_DCT_OVERFLOW_CHECK(td, 15137, 8192, 46);
  t3 -= (td*15137 + 8192) >> 14;
  /*10947/16384 ~= \frac{1 - cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})} ~=
     0.66817863791929891999775768652308*/
  OD_DCT_OVERFLOW_CHECK(t3, 10947, 8192, 47);
  td += (t3*10947 + 8192) >> 14;
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(tf, 13573, 16384, 48);
  t1 -= (tf*13573 + 16384) >> 15;
  /*11585/16384 ~= \sqrt{\frac{1}{2}} ~= 0.70710678118654752440084436210485*/
  OD_DCT_OVERFLOW_CHECK(t1, 11585, 8192, 49);
  tf += (t1*11585 + 8192) >> 14;
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  OD_DCT_OVERFLOW_CHECK(tf, 13573, 16384, 50);
  t1 -= (tf*13573 + 16384) >> 15;
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
  y[4] = (od_coeff)t4;
  y[5] = (od_coeff)t5;
  y[6] = (od_coeff)t6;
  y[7] = (od_coeff)t7;
  y[8] = (od_coeff)t8;
  y[9] = (od_coeff)t9;
  y[10] = (od_coeff)ta;
  y[11] = (od_coeff)tb;
  y[12] = (od_coeff)tc;
  y[13] = (od_coeff)td;
  y[14] = (od_coeff)te;
  y[15] = (od_coeff)tf;
}

void od_bin_idct16(od_coeff *x, int xstride, const od_coeff y[16]) {
  int t0;
  int t1;
  int t1h;
  int t2;
  int t2h;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  int t8;
  int t8h;
  int t9;
  int ta;
  int tah;
  int tb;
  int tbh;
  int tc;
  int tch;
  int td;
  int tdh;
  int te;
  int tf;
  int tfh;
  t0 = y[0];
  t1 = y[1];
  t2 = y[2];
  t3 = y[3];
  t4 = y[4];
  t5 = y[5];
  t6 = y[6];
  t7 = y[7];
  t8 = y[8];
  t9 = y[9];
  ta = y[10];
  tb = y[11];
  tc = y[12];
  td = y[13];
  te = y[14];
  tf = y[15];
  t1 += (tf*13573 + 16384) >> 15;
  tf -= (t1*11585 + 8192) >> 14;
  t1 += ((tf*13573 + 16384) >> 15)+t7;
  td -= (t3*10947 + 8192) >> 14;
  t3 += (td*15137 + 8192) >> 14;
  t5 += (tb*10947 + 8192) >> 14;
  tb -= (t5*15137 + 8192) >> 14;
  t5 += (tb*10947 + 8192) >> 14;
  td += t5 - ((t3*21895 + 16384) >> 15);
  tf = t9 - tf;
  tb += t3;
  tfh = OD_DCT_RSHIFT(tf, 1);
  t9 -= tfh;
  tbh = OD_DCT_RSHIFT(tb, 1);
  t3 += tfh - tbh;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t7 = t1h - t7 + tbh;
  tdh = OD_DCT_RSHIFT(td, 1);
  t5 += t1h - tdh;
  t9 = tdh - t9;
  td -= t9;
  tf = t3 - tf;
  t1 -= t5 + ((tf*20055 + 16384) >> 15);
  tf += (t1*23059 + 8192) >> 14;
  t1 -= (tf*5417 + 4096) >> 13;
  tb = t7 - tb;
  t9 += (t7*14101 + 8192) >> 14;
  t7 += (t9*3363 + 4096) >> 13;
  t9 -= (t7*12905 + 8192) >> 14;
  tb -= (td*4379 + 8192) >> 14;
  td += (tb*20435 + 8192) >> 14;
  tb -= (td*17515 + 16384) >> 15;
  t3 += (t5*851 + 4096) >> 13;
  t5 += (t3*14699 + 8192) >> 14;
  t3 -= (t5*1035 + 1024) >> 11;
  t6 -= (ta*7335 + 16384) >> 15;
  ta -= (t6*12873 + 8192) >> 14;
  te += (t2*2873 + 1024) >> 11;
  t2 += (te*9041 + 16384) >> 15;
  t6 = OD_DCT_RSHIFT(t2, 1) - t6 - ((ta*8593 + 8192) >> 14);
  te = OD_DCT_RSHIFT(ta, 1) - te + ((t2*2275 + 1024) >> 11);
  t2 -= t6;
  ta -= te;
  t6 -= (ta*13573 + 16384) >> 15;
  ta += (t6*11585 + 8192) >> 14;
  t6 -= (ta*13573 + 16384) >> 15;
  tc += (t4*9147 + 4096) >> 13;
  t4 -= (tc*10703 + 8192) >> 14;
  tc += (t4*23013 + 16384) >> 15;
  t8 = t0 - t8;
  t8h = OD_DCT_RSHIFT(t8, 1);
  t0 -= t8h - OD_DCT_RSHIFT(tc, 1);
  t4 = t8h - t4;
  t8 += t6 - t4;
  tc = t0 - tc + te;
  ta = t4 - ta;
  t2 = t0 - t2;
  tch = OD_DCT_RSHIFT(tc, 1);
  te = tch - te;
  tah = OD_DCT_RSHIFT(ta, 1);
  t4 -= tah;
  t8h = OD_DCT_RSHIFT(t8, 1);
  t6 = t8h - t6;
  t2h = OD_DCT_RSHIFT(t2, 1);
  t0 -= t2h;
  t3 = t2h - t3;
  t6 += OD_DCT_RSHIFT(td, 1);
  t9 = tah - t9;
  te += OD_DCT_RSHIFT(tf, 1);
  t1 = tch - t1;
  t4 += OD_DCT_RSHIFT(t7, 1);
  tb = t8h - tb;
  t0 += OD_DCT_RSHIFT(t5, 1);
  *(x + 0*xstride) = (od_coeff)t0;
  *(x + 1*xstride) = (od_coeff)(t8 - tb);
  *(x + 2*xstride) = (od_coeff)t4;
  *(x + 3*xstride) = (od_coeff)(tc - t1);
  *(x + 4*xstride) = (od_coeff)te;
  *(x + 5*xstride) = (od_coeff)(ta - t9);
  *(x + 6*xstride) = (od_coeff)t6;
  *(x + 7*xstride) = (od_coeff)(t2 - t3);
  *(x + 8*xstride) = (od_coeff)t3;
  *(x + 9*xstride) = (od_coeff)(t6 - td);
  *(x + 10*xstride) = (od_coeff)t9;
  *(x + 11*xstride) = (od_coeff)(te - tf);
  *(x + 12*xstride) = (od_coeff)t1;
  *(x + 13*xstride) = (od_coeff)(t4 - t7);
  *(x + 14*xstride) = (od_coeff)tb;
  *(x + 15*xstride) = (od_coeff)(t0 - t5);
}
#endif

void od_bin_fdct16x16(od_coeff *y, int ystride,
 const od_coeff *x, int xstride) {
  od_coeff z[16*16];
  int i;
  for (i = 0; i < 16; i++) od_bin_fdct16(z + 16*i, x + i, xstride);
  for (i = 0; i < 16; i++) od_bin_fdct16(y + ystride*i, z + i, 16);
}

void od_bin_idct16x16(od_coeff *x, int xstride,
 const od_coeff *y, int ystride) {
  od_coeff z[16*16];
  int i;
  for (i = 0; i < 16; i++) od_bin_idct16(z + i, 16, y + ystride*i);
  for (i = 0; i < 16; i++) od_bin_idct16(x + i, xstride, z + 16*i);
}

#define OD_FDCT_2(t0, t1) \
  /* Embedded 2-point orthonormal Type-II fDCT. */ \
  do { \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 13573, 16384, 100); \
    t0 -= (t1*13573 + 16384) >> 15; \
    /* 5793/8192 ~= Sin[pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 5793, 4096, 101); \
    t1 += (t0*5793 + 4096) >> 13; \
    /* 3393/8192 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 3393, 4096, 102); \
    t0 -= (t1*3393 + 4096) >> 13; \
  } \
  while (0)

#define OD_IDCT_2(t0, t1) \
  /* Embedded 2-point orthonormal Type-II iDCT. */ \
  do { \
    /* 3393/8192 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    t0 += (t1*3393 + 4096) >> 13; \
    /* 5793/8192 ~= Sin[pi/4] ~= 0.707106781186547 */ \
    t1 -= (t0*5793 + 4096) >> 13; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    t0 += (t1*13573 + 16384) >> 15; \
  } \
  while (0)

#define OD_FDCT_2_ASYM(p0, p1, p1h) \
  /* Embedded 2-point asymmetric Type-II fDCT. */ \
  do { \
    p0 += p1h; \
    p1 = p0 - p1; \
  } \
  while (0)

#define OD_IDCT_2_ASYM(p0, p1, p1h) \
  /* Embedded 2-point asymmetric Type-II iDCT. */ \
  do { \
    p1 = p0 - p1; \
    p1h = OD_DCT_RSHIFT(p1, 1); \
    p0 -= p1h; \
  } \
  while (0)

#define OD_FDST_2(t0, t1) \
  /* Embedded 2-point orthonormal Type-IV fDST. */ \
  do { \
    /* 10947/16384 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 10947, 8192, 103); \
    t0 -= (t1*10947 + 8192) >> 14; \
    /* 473/512 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 473, 256, 104); \
    t1 += (t0*473 + 256) >> 9; \
    /* 10947/16384 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 10947, 8192, 105); \
    t0 -= (t1*10947 + 8192) >> 14; \
  } \
  while (0)

#define OD_IDST_2(t0, t1) \
  /* Embedded 2-point orthonormal Type-IV iDST. */ \
  do { \
    /* 10947/16384 ~= Tan[3*Pi/16]) ~= 0.668178637919299 */ \
    t0 += (t1*10947 + 8192) >> 14; \
    /* 473/512 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t1 -= (t0*473 + 256) >> 9; \
    /* 10947/16384 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    t0 += (t1*10947 + 8192) >> 14; \
  } \
  while (0)

#define OD_FDST_2_ASYM(p0, p1) \
  /* Embedded 2-point asymmetric Type-IV fDST. */ \
  do { \
    /* 11507/16384 ~= 4*Sin[Pi/8] - 2*Tan[Pi/8] ~= 0.702306604714169 */ \
    OD_DCT_OVERFLOW_CHECK(p1, 11507, 8192, 187); \
    p0 -= (p1*11507 + 8192) >> 14; \
    /* 669/1024 ~= Cos[Pi/8]/Sqrt[2] ~= 0.653281482438188 */ \
    OD_DCT_OVERFLOW_CHECK(p0, 669, 512, 188); \
    p1 += (p0*669 + 512) >> 10; \
    /* 4573/4096 ~= 4*Sin[Pi/8] - Tan[Pi/8] ~= 1.11652016708726 */ \
    OD_DCT_OVERFLOW_CHECK(p1, 4573, 2048, 189); \
    p0 -= (p1*4573 + 2048) >> 12; \
  } \
  while (0)

#define OD_IDST_2_ASYM(p0, p1) \
  /* Embedded 2-point asymmetric Type-IV iDST. */ \
  do { \
    /* 4573/4096 ~= 4*Sin[Pi/8] - Tan[Pi/8] ~= 1.11652016708726 */ \
    p0 += (p1*4573 + 2048) >> 12; \
    /* 669/1024 ~= Cos[Pi/8]/Sqrt[2] ~= 0.653281482438188 */ \
    p1 -= (p0*669 + 512) >> 10; \
    /* 11507/16384 ~= 4*Sin[Pi/8] - 2*Tan[Pi/8] ~= 0.702306604714169 */ \
    p0 += (p1*11507 + 8192) >> 14; \
  } \
  while (0)

#define OD_FDCT_4(q0, q2, q1, q3) \
  /* Embedded 4-point orthonormal Type-II fDCT. */ \
  do { \
    int q2h; \
    int q3h; \
    q3 = q0 - q3; \
    q3h = OD_DCT_RSHIFT(q3, 1); \
    q0 -= q3h; \
    q2 += q1; \
    q2h = OD_DCT_RSHIFT(q2, 1); \
    q1 = q2h - q1; \
    OD_FDCT_2_ASYM(q0, q2, q2h); \
    OD_FDST_2_ASYM(q3, q1); \
  } \
  while (0)

#define OD_IDCT_4(q0, q2, q1, q3) \
  /* Embedded 4-point orthonormal Type-II iDCT. */ \
  do { \
    int q1h; \
    int q3h; \
    OD_IDST_2_ASYM(q3, q2); \
    OD_IDCT_2_ASYM(q0, q1, q1h); \
    q3h = OD_DCT_RSHIFT(q3, 1); \
    q0 += q3h; \
    q3 = q0 - q3; \
    q2 = q1h - q2; \
    q1 -= q2; \
  } \
  while (0)

#define OD_FDCT_4_ASYM(t0, t2, t2h, t1, t3, t3h) \
  /* Embedded 4-point asymmetric Type-II fDCT. */ \
  do { \
    t0 += t3h; \
    t3 = t0 - t3; \
    t1 = t2h - t1; \
    t2 = t1 - t2; \
    OD_FDCT_2(t0, t2); \
    OD_FDST_2(t3, t1); \
  } \
  while (0)

#define OD_IDCT_4_ASYM(t0, t2, t1, t1h, t3, t3h) \
  /* Embedded 4-point asymmetric Type-II iDCT. */ \
  do { \
    OD_IDST_2(t3, t2); \
    OD_IDCT_2(t0, t1); \
    t1 = t2 - t1; \
    t1h = OD_DCT_RSHIFT(t1, 1); \
    t2 = t1h - t2; \
    t3 = t0 - t3; \
    t3h = OD_DCT_RSHIFT(t3, 1); \
    t0 -= t3h; \
  } \
  while (0)

#define OD_FDST_4(q0, q2, q1, q3) \
  /* Embedded 4-point orthonormal Type-IV fDST. */ \
  do { \
    int q0h; \
    int q1h; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(q1, 13573, 16384, 190); \
    q2 += (q1*13573 + 16384) >> 15; \
    /* 5793/8192 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(q2, 5793, 4096, 191); \
    q1 -= (q2*5793 + 4096) >> 13; \
    /* 3393/8192 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(q1, 3393, 4096, 192); \
    q2 += (q1*3393 + 4096) >> 13; \
    q0 += q2; \
    q0h = OD_DCT_RSHIFT(q0, 1); \
    q2 = q0h - q2; \
    q1 += q3; \
    q1h = OD_DCT_RSHIFT(q1, 1); \
    q3 -= q1h; \
    /* 537/1024 ~= (1/Sqrt[2] - Cos[3*Pi/16]/2)/Sin[3*Pi/16] ~=
        0.524455699240090 */ \
    OD_DCT_OVERFLOW_CHECK(q1, 537, 512, 193); \
    q2 -= (q1*537 + 512) >> 10; \
    /* 1609/2048 ~= Sqrt[2]*Sin[3*Pi/16] ~= 0.785694958387102 */ \
    OD_DCT_OVERFLOW_CHECK(q2, 1609, 1024, 194); \
    q1 += (q2*1609 + 1024) >> 11; \
    /* 7335/32768 ~= (1/Sqrt[2] - Cos[3*Pi/16])/Sin[3*Pi/16] ~=
        0.223847182092655 */ \
    OD_DCT_OVERFLOW_CHECK(q1, 7335, 16384, 195); \
    q2 += (q1*7335 + 16384) >> 15; \
    /* 5091/8192 ~= (1/Sqrt[2] - Cos[7*Pi/16]/2)/Sin[7*Pi/16] ~=
        0.6215036383171189 */ \
    OD_DCT_OVERFLOW_CHECK(q0, 5091, 4096, 196); \
    q3 += (q0*5091 + 4096) >> 13; \
    /* 5681/4096 ~= Sqrt[2]*Sin[7*Pi/16] ~= 1.38703984532215 */ \
    OD_DCT_OVERFLOW_CHECK(q3, 5681, 2048, 197); \
    q0 -= (q3*5681 + 2048) >> 12; \
    /* 4277/8192 ~= (1/Sqrt[2] - Cos[7*Pi/16])/Sin[7*Pi/16] ~=
        0.52204745462729 */ \
    OD_DCT_OVERFLOW_CHECK(q0, 4277, 4096, 198); \
    q3 += (q0*4277 + 4096) >> 13; \
  } \
  while (0)

#define OD_IDST_4(q0, q2, q1, q3) \
  /* Embedded 4-point orthonormal Type-IV iDST. */ \
  do { \
    int q0h; \
    int q2h; \
    /* 4277/8192 ~= (1/Sqrt[2] - Cos[7*Pi/16])/Sin[7*Pi/16] ~=
        0.52204745462729 */ \
    q3 -= (q0*4277 + 4096) >> 13; \
    /* 5681/4096 ~= Sqrt[2]*Sin[7*Pi/16] ~= 1.38703984532215 */ \
    q0 += (q3*5681 + 2048) >> 12; \
    /* 5091/8192 ~= (1/Sqrt[2] - Cos[7*Pi/16]/2)/Sin[7*Pi/16] ~=
        0.6215036383171189 */ \
    q3 -= (q0*5091 + 4096) >> 13; \
    /* 7335/32768 ~= (1/Sqrt[2] - Cos[3*Pi/16])/Sin[3*Pi/16] ~=
        0.223847182092655 */ \
    q1 -= (q2*7335 + 16384) >> 15; \
    /* 1609/2048 ~= Sqrt[2]*Sin[3*Pi/16] ~= 0.785694958387102 */ \
    q2 -= (q1*1609 + 1024) >> 11; \
    /* 537/1024 ~= (1/Sqrt[2] - Cos[3*Pi/16]/2)/Sin[3*Pi/16] ~=
        0.524455699240090 */ \
    q1 += (q2*537 + 512) >> 10; \
    q2h = OD_DCT_RSHIFT(q2, 1); \
    q3 += q2h; \
    q2 -= q3; \
    q0h = OD_DCT_RSHIFT(q0, 1); \
    q1 = q0h - q1; \
    q0 -= q1; \
    /* 3393/8192 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    q1 -= (q2*3393 + 4096) >> 13; \
    /* 5793/8192 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    q2 += (q1*5793 + 4096) >> 13; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    q1 -= (q2*13573 + 16384) >> 15; \
  } \
  while (0)

#define OD_FDST_4_ASYM(t0, t0h, t2, t1, t3) \
  /* Embedded 4-point asymmetric Type-IV fDST. */ \
  do { \
    /* 7489/8192 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 7489, 4096, 106); \
    t2 -= (t1*7489 + 4096) >> 13; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 11585, 8192, 107); \
    t1 += (t2*11585 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 19195, 16384, 108); \
    t2 += (t1*19195 + 16384) >> 15; \
    t3 += OD_DCT_RSHIFT(t2, 1); \
    t2 -= t3; \
    t1 = t0h - t1; \
    t0 -= t1; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 6723, 4096, 109); \
    t3 += (t0*6723 + 4096) >> 13; \
    /* 8035/8192 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 8035, 4096, 110); \
    t0 -= (t3*8035 + 4096) >> 13; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 6723, 4096, 111); \
    t3 += (t0*6723 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 8757, 8192, 112); \
    t2 += (t1*8757 + 8192) >> 14; \
    /* 6811/8192 ~= Sin[5*Pi/16] ~= 0.831469612302545 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 6811, 4096, 113); \
    t1 -= (t2*6811 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 8757, 8192, 114); \
    t2 += (t1*8757 + 8192) >> 14; \
  } \
  while (0)

#define OD_IDST_4_ASYM(t0, t0h, t2, t1, t3) \
  /* Embedded 4-point asymmetric Type-IV iDST. */ \
  do { \
    /* 8757/16384 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    t1 -= (t2*8757 + 8192) >> 14; \
    /* 6811/8192 ~= Sin[5*Pi/16] ~= 0.831469612302545 */ \
    t2 += (t1*6811 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    t1 -= (t2*8757 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    t3 -= (t0*6723 + 4096) >> 13; \
    /* 8035/8192 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    t0 += (t3*8035 + 4096) >> 13; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    t3 -= (t0*6723 + 4096) >> 13; \
    t0 += t2; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t2 = t0h - t2; \
    t1 += t3; \
    t3 -= OD_DCT_RSHIFT(t1, 1); \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    t1 -= (t2*19195 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    t2 -= (t1*11585 + 8192) >> 14; \
    /* 7489/8192 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    t1 += (t2*7489 + 4096) >> 13; \
  } \
  while (0)

#define OD_FDCT_8(r0, r4, r2, r6, r1, r5, r3, r7) \
  /* Embedded 8-point orthonormal Type-II fDCT. */ \
  do { \
    int r4h; \
    int r5h; \
    int r6h; \
    int r7h; \
    r7 = r0 - r7; \
    r7h = OD_DCT_RSHIFT(r7, 1); \
    r0 -= r7h; \
    r6 += r1; \
    r6h = OD_DCT_RSHIFT(r6, 1); \
    r1 = r6h - r1; \
    r5 = r2 - r5; \
    r5h = OD_DCT_RSHIFT(r5, 1); \
    r2 -= r5h; \
    r4 += r3; \
    r4h = OD_DCT_RSHIFT(r4, 1); \
    r3 = r4h - r3; \
    OD_FDCT_4_ASYM(r0, r4, r4h, r2, r6, r6h); \
    OD_FDST_4_ASYM(r7, r7h, r3, r5, r1); \
  } \
  while (0)

#define OD_IDCT_8(r0, r4, r2, r6, r1, r5, r3, r7) \
  /* Embedded 8-point orthonormal Type-II iDCT. */ \
  do { \
    int r1h; \
    int r3h; \
    int r5h; \
    int r7h; \
    OD_IDST_4_ASYM(r7, r7h, r5, r6, r4); \
    OD_IDCT_4_ASYM(r0, r2, r1, r1h, r3, r3h); \
    r0 += r7h; \
    r7 = r0 - r7; \
    r6 = r1h - r6; \
    r1 -= r6; \
    r5h = OD_DCT_RSHIFT(r5, 1); \
    r2 += r5h; \
    r5 = r2 - r5; \
    r4 = r3h - r4; \
    r3 -= r4; \
  } \
  while (0)

#define OD_FDCT_8_ASYM(r0, r4, r4h, r2, r6, r6h, r1, r5, r5h, r3, r7, r7h) \
  /* Embedded 8-point asymmetric Type-II fDCT. */ \
  do { \
    r0 += r7h; \
    r7 = r0 - r7; \
    r1 = r6h - r1; \
    r6 -= r1; \
    r2 += r5h; \
    r5 = r2 - r5; \
    r3 = r4h - r3; \
    r4 -= r3; \
    OD_FDCT_4(r0, r4, r2, r6); \
    OD_FDST_4(r7, r3, r5, r1); \
  } \
  while (0)

#define OD_IDCT_8_ASYM(r0, r4, r2, r6, r1, r1h, r5, r5h, r3, r3h, r7, r7h) \
  /* Embedded 8-point asymmetric Type-II iDCT. */ \
  do { \
    OD_IDST_4(r7, r5, r6, r4); \
    OD_IDCT_4(r0, r2, r1, r3); \
    r7 = r0 - r7; \
    r7h = OD_DCT_RSHIFT(r7, 1); \
    r0 -= r7h; \
    r1 += r6; \
    r1h = OD_DCT_RSHIFT(r1, 1); \
    r6 = r1h - r6; \
    r5 = r2 - r5; \
    r5h = OD_DCT_RSHIFT(r5, 1); \
    r2 -= r5h; \
    r3 += r4; \
    r3h = OD_DCT_RSHIFT(r3, 1); \
    r4 = r3h - r4; \
  } \
  while (0)

#define OD_FDST_8(t0, t4, t2, t6, t1, t5, t3, t7) \
  /* Embedded 8-point orthonormal Type-IV fDST. */ \
  do { \
    int t0h; \
    int t2h; \
    int t5h; \
    int t7h; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 13573, 16384, 115); \
    t6 -= (t1*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 11585, 8192, 116); \
    t1 += (t6*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 13573, 16384, 117); \
    t6 -= (t1*13573 + 16384) >> 15; \
    /* 21895/32768 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 21895, 16384, 118); \
    t5 -= (t2*21895 + 16384) >> 15; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 15137, 8192, 119); \
    t2 += (t5*15137 + 8192) >> 14; \
    /* 10947/16384 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 10947, 8192, 120); \
    t5 -= (t2*10947 + 8192) >> 14; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 3259, 8192, 121); \
    t4 -= (t3*3259 + 8192) >> 14; \
    /* 3135/8192 ~= Sin[Pi/8] ~= 0.382683432365090 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 3135, 4096, 122); \
    t3 += (t4*3135 + 4096) >> 13; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 3259, 8192, 123); \
    t4 -= (t3*3259 + 8192) >> 14; \
    t7 += t1; \
    t7h = OD_DCT_RSHIFT(t7, 1); \
    t1 -= t7h; \
    t2 = t3 - t2; \
    t2h = OD_DCT_RSHIFT(t2, 1); \
    t3 -= t2h; \
    t0 -= t6; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t6 += t0h; \
    t5 = t4 - t5; \
    t5h = OD_DCT_RSHIFT(t5, 1); \
    t4 -= t5h; \
    t1 += t5h; \
    t5 = t1 - t5; \
    t4 += t0h; \
    t0 -= t4; \
    t6 -= t2h; \
    t2 += t6; \
    t3 -= t7h; \
    t7 += t3; \
    /* TODO: Can we move this into another operation */ \
    t7 = -t7; \
    /* 7425/8192 ~= Tan[15*Pi/64] ~= 0.906347169019147 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 7425, 4096, 124); \
    t0 -= (t7*7425 + 4096) >> 13; \
    /* 8153/8192 ~= Sin[15*Pi/32] ~= 0.995184726672197 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 8153, 4096, 125); \
    t7 += (t0*8153 + 4096) >> 13; \
    /* 7425/8192 ~= Tan[15*Pi/64] ~= 0.906347169019147 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 7425, 4096, 126); \
    t0 -= (t7*7425 + 4096) >> 13; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.148335987538347 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 4861, 16384, 127); \
    t6 -= (t1*4861 + 16384) >> 15; \
    /* 1189/4096 ~= Sin[3*Pi/32] ~= 0.290284677254462 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 1189, 2048, 128); \
    t1 += (t6*1189 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.148335987538347 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 4861, 16384, 129); \
    t6 -= (t1*4861 + 16384) >> 15; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.599376933681924 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 2455, 2048, 130); \
    t2 -= (t5*2455 + 2048) >> 12; \
    /* 7225/8192 ~= Sin[11*Pi/32] ~= 0.881921264348355 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 7225, 4096, 131); \
    t5 += (t2*7225 + 4096) >> 13; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.599376933681924 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 2455, 2048, 132); \
    t2 -= (t5*2455 + 2048) >> 12; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.357805721314524 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 11725, 16384, 133); \
    t4 -= (t3*11725 + 16384) >> 15; \
    /* 5197/8192 ~= Sin[7*Pi/32] ~= 0.634393284163645 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 5197, 4096, 134); \
    t3 += (t4*5197 + 4096) >> 13; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.357805721314524 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 11725, 16384, 135); \
    t4 -= (t3*11725 + 16384) >> 15; \
  } \
  while (0)

#define OD_IDST_8(t0, t4, t2, t6, t1, t5, t3, t7) \
  /* Embedded 8-point orthonormal Type-IV iDST. */ \
  do { \
    int t0h; \
    int t2h; \
    int t5h_; \
    int t7h_; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.357805721314524 */ \
    t1 += (t6*11725 + 16384) >> 15; \
    /* 5197/8192 ~= Sin[7*Pi/32] ~= 0.634393284163645 */ \
    t6 -= (t1*5197 + 4096) >> 13; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.357805721314524 */ \
    t1 += (t6*11725 + 16384) >> 15; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.599376933681924 */ \
    t2 += (t5*2455 + 2048) >> 12; \
    /* 7225/8192 ~= Sin[11*Pi/32] ~= 0.881921264348355 */ \
    t5 -= (t2*7225 + 4096) >> 13; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.599376933681924 */ \
    t2 += (t5*2455 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.148335987538347 */ \
    t3 += (t4*4861 + 16384) >> 15; \
    /* 1189/4096 ~= Sin[3*Pi/32] ~= 0.290284677254462 */ \
    t4 -= (t3*1189 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.148335987538347 */ \
    t3 += (t4*4861 + 16384) >> 15; \
    /* 7425/8192 ~= Tan[15*Pi/64] ~= 0.906347169019147 */ \
    t0 += (t7*7425 + 4096) >> 13; \
    /* 8153/8192 ~= Sin[15*Pi/32] ~= 0.995184726672197 */ \
    t7 -= (t0*8153 + 4096) >> 13; \
    /* 7425/8192 ~= Tan[15*Pi/64] ~= 0.906347169019147 */ \
    t0 += (t7*7425 + 4096) >> 13; \
    /* TODO: Can we move this into another operation */ \
    t7 = -t7; \
    t7 -= t6; \
    t7h_ = OD_DCT_RSHIFT(t7, 1); \
    t6 += t7h_; \
    t2 -= t3; \
    t2h = OD_DCT_RSHIFT(t2, 1); \
    t3 += t2h; \
    t0 += t1; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t1 -= t0h; \
    t5 = t4 - t5; \
    t5h_ = OD_DCT_RSHIFT(t5, 1); \
    t4 -= t5h_; \
    t1 += t5h_; \
    t5 = t1 - t5; \
    t3 -= t0h; \
    t0 += t3; \
    t6 += t2h; \
    t2 = t6 - t2; \
    t4 += t7h_; \
    t7 -= t4; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    t1 += (t6*3259 + 8192) >> 14; \
    /* 3135/8192 ~= Sin[Pi/8] ~= 0.382683432365090 */ \
    t6 -= (t1*3135 + 4096) >> 13; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    t1 += (t6*3259 + 8192) >> 14; \
    /* 10947/16384 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    t5 += (t2*10947 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t2 -= (t5*15137 + 8192) >> 14; \
    /* 21895/32768 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    t5 += (t2*21895 + 16384) >> 15; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    t3 += (t4*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    t4 -= (t3*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    t3 += (t4*13573 + 16384) >> 15; \
  } \
  while (0)

/* Rewrite this so that t0h can be passed in. */
#define OD_FDST_8_ASYM(t0, t4, t2, t6, t1, t5, t3, t7) \
  /* Embedded 8-point asymmetric Type-IV fDST. */ \
  do { \
    int t0h; \
    int t2h; \
    int t5h; \
    int t7h; \
    /* 1035/2048 ~= (Sqrt[2] - Cos[7*Pi/32])/(2*Sin[7*Pi/32]) ~=
        0.505367194937830 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 1035, 1024, 199); \
    t6 += (t1*1035 + 1024) >> 11; \
    /* 3675/4096 ~= Sqrt[2]*Sin[7*Pi/32] ~= 0.897167586342636 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 3675, 2048, 200); \
    t1 -= (t6*3675 + 2048) >> 12; \
    /* 851/8192 ~= (Cos[7*Pi/32] - 1/Sqrt[2])/Sin[7*Pi/32] ~=
        0.103884567856159 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 851, 4096, 201); \
    t6 -= (t1*851 + 4096) >> 13; \
    /* 4379/8192 ~= (Sqrt[2] - Sin[5*Pi/32])/(2*Cos[5*Pi/32]) ~=
        0.534524375168421 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 4379, 4096, 202); \
    t5 += (t2*4379 + 4096) >> 13; \
    /* 10217/8192 ~= Sqrt[2]*Cos[5*Pi/32] ~= 1.24722501298667 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 10217, 4096, 203); \
    t2 -= (t5*10217 + 4096) >> 13; \
    /* 4379/16384 ~= (1/Sqrt[2] - Sin[5*Pi/32])/Cos[5*Pi/32] ~=
        0.267268807193026 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 4379, 8192, 204); \
    t5 += (t2*4379 + 8192) >> 14; \
    /* 12905/16384 ~= (Sqrt[2] - Cos[3*Pi/32])/(2*Sin[3*Pi/32]) ~=
        0.787628942329675 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 12905, 8192, 205); \
    t4 += (t3*12905 + 8192) >> 14; \
    /* 3363/8192 ~= Sqrt[2]*Sin[3*Pi/32] ~= 0.410524527522357 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 3363, 4096, 206); \
    t3 -= (t4*3363 + 4096) >> 13; \
    /* 3525/4096 ~= (Cos[3*Pi/32] - 1/Sqrt[2])/Sin[3*Pi/32] ~=
        0.860650162139486 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 3525, 2048, 207); \
    t4 -= (t3*3525 + 2048) >> 12; \
    /* 5417/8192 ~= (Sqrt[2] - Sin[Pi/32])/(2*Cos[Pi/32]) ~=
        0.661282466846517 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 5417, 4096, 208); \
    t7 += (t0*5417 + 4096) >> 13; \
    /* 5765/4096 ~= Sqrt[2]*Cos[Pi/32] ~= 1.40740373752638 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 5765, 2048, 209); \
    t0 -= (t7*5765 + 2048) >> 12; \
    /* 2507/4096 ~= (1/Sqrt[2] - Sin[Pi/32])/Cos[Pi/32] ~=
        0.612036765167935 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 2507, 2048, 210); \
    t7 += (t0*2507 + 2048) >> 12; \
    t0 += t1; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t1 -= t0h; \
    t2 -= t3; \
    t2h = OD_DCT_RSHIFT(t2, 1); \
    t3 += t2h; \
    t5 -= t4; \
    t5h = OD_DCT_RSHIFT(t5, 1); \
    t4 += t5h; \
    t7 += t6; \
    t7h = OD_DCT_RSHIFT(t7, 1); \
    t6 = t7h - t6; \
    t4 = t7h - t4; \
    t7 -= t4; \
    t1 += t5h; \
    t5 = t1 - t5; \
    t6 += t2h; \
    t2 = t6 - t2; \
    t3 -= t0h; \
    t0 += t3; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 3259, 8192, 211); \
    t1 += (t6*3259 + 8192) >> 14; \
    /* 3135/8192 ~= Sin[Pi/8] ~= 0.382683432365090 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 3135, 4096, 212); \
    t6 -= (t1*3135 + 4096) >> 13; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 3259, 8192, 213); \
    t1 += (t6*3259 + 8192) >> 14; \
    /* 2737/4096 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 2737, 2048, 214); \
    t5 += (t2*2737 + 2048) >> 12; \
    /* 473/512 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 473, 256, 215); \
    t2 -= (t5*473 + 256) >> 9; \
    /* 2737/4096 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 2737, 2048, 216); \
    t5 += (t2*2737 + 2048) >> 12; \
    /* 3393/8192 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 3393, 4096, 217); \
    t3 += (t4*3393 + 4096) >> 13; \
    /* 5793/8192 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 5793, 4096, 218); \
    t4 -= (t3*5793 + 4096) >> 13; \
    /* 3393/8192 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 3393, 4096, 219); \
    t3 += (t4*3393 + 4096) >> 13; \
  } \
  while (0)

#define OD_IDST_8_ASYM(t0, t4, t2, t6, t1, t5, t3, t7) \
  /* Embedded 8-point asymmetric Type-IV iDST. */ \
  do { \
    int t0h; \
    int t2h; \
    int t5h__; \
    int t7h__; \
    /* 3393/8192 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    t6 -= (t1*3393 + 4096) >> 13; \
    /* 5793/8192 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    t1 += (t6*5793 + 4096) >> 13; \
    /* 3393/8192 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    t6 -= (t1*3393 + 4096) >> 13; \
    /* 2737/4096 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    t5 -= (t2*2737 + 2048) >> 12; \
    /* 473/512 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t2 += (t5*473 + 256) >> 9; \
    /* 2737/4096 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    t5 -= (t2*2737 + 2048) >> 12; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    t4 -= (t3*3259 + 8192) >> 14; \
    /* 3135/8192 ~= Sin[Pi/8] ~= 0.382683432365090 */ \
    t3 += (t4*3135 + 4096) >> 13; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    t4 -= (t3*3259 + 8192) >> 14; \
    t0 -= t6; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t6 += t0h; \
    t2 = t3 - t2; \
    t2h = OD_DCT_RSHIFT(t2, 1); \
    t3 -= t2h; \
    t5 = t4 - t5; \
    t5h__ = OD_DCT_RSHIFT(t5, 1); \
    t4 -= t5h__; \
    t7 += t1; \
    t7h__ = OD_DCT_RSHIFT(t7, 1); \
    t1 = t7h__ - t1; \
    t3 = t7h__ - t3; \
    t7 -= t3; \
    t1 -= t5h__; \
    t5 += t1; \
    t6 -= t2h; \
    t2 += t6; \
    t4 += t0h; \
    t0 -= t4; \
    /* 2507/4096 ~= (1/Sqrt[2] - Sin[Pi/32])/Cos[Pi/32] ~=
        0.612036765167935 */ \
    t7 -= (t0*2507 + 2048) >> 12; \
    /* 5765/4096 ~= Sqrt[2]*Cos[Pi/32] ~= 1.40740373752638 */ \
    t0 += (t7*5765 + 2048) >> 12; \
    /* 5417/8192 ~= (Sqrt[2] - Sin[Pi/32])/(2*Cos[Pi/32]) ~=
        0.661282466846517 */ \
    t7 -= (t0*5417 + 4096) >> 13; \
    /* 3525/4096 ~= (Cos[3*Pi/32] - 1/Sqrt[2])/Sin[3*Pi/32] ~=
        0.860650162139486 */ \
    t1 += (t6*3525 + 2048) >> 12; \
    /* 3363/8192 ~= Sqrt[2]*Sin[3*Pi/32] ~= 0.410524527522357 */ \
    t6 += (t1*3363 + 4096) >> 13; \
    /* 12905/16384 ~= (1/Sqrt[2] - Cos[3*Pi/32]/1)/Sin[3*Pi/32] ~=
        0.787628942329675 */ \
    t1 -= (t6*12905 + 8192) >> 14; \
    /* 4379/16384 ~= (1/Sqrt[2] - Sin[5*Pi/32])/Cos[5*Pi/32] ~=
        0.267268807193026 */ \
    t5 -= (t2*4379 + 8192) >> 14; \
    /* 10217/8192 ~= Sqrt[2]*Cos[5*Pi/32] ~= 1.24722501298667 */ \
    t2 += (t5*10217 + 4096) >> 13; \
    /* 4379/8192 ~= (Sqrt[2] - Sin[5*Pi/32])/(2*Cos[5*Pi/32]) ~=
        0.534524375168421 */ \
    t5 -= (t2*4379 + 4096) >> 13; \
    /* 851/8192 ~= (Cos[7*Pi/32] - 1/Sqrt[2])/Sin[7*Pi/32] ~=
        0.103884567856159 */ \
    t3 += (t4*851 + 4096) >> 13; \
    /* 3675/4096 ~= Sqrt[2]*Sin[7*Pi/32] ~= 0.897167586342636 */ \
    t4 += (t3*3675 + 2048) >> 12; \
    /* 1035/2048 ~= (Sqrt[2] - Cos[7*Pi/32])/(2*Sin[7*Pi/32]) ~=
        0.505367194937830 */ \
    t3 -= (t4*1035 + 1024) >> 11; \
  } \
  while (0)

#define OD_FDCT_16(s0, s8, s4, sc, s2, sa, s6, se, \
 s1, s9, s5, sd, s3, sb, s7, sf) \
  /* Embedded 16-point orthonormal Type-II fDCT. */ \
  do { \
    int s8h; \
    int sah; \
    int sch; \
    int seh; \
    int sfh; \
    sf = s0 - sf; \
    sfh = OD_DCT_RSHIFT(sf, 1); \
    s0 -= sfh; \
    se += s1; \
    seh = OD_DCT_RSHIFT(se, 1); \
    s1 = seh - s1; \
    sd = s2 - sd; \
    s2 -= OD_DCT_RSHIFT(sd, 1); \
    sc += s3; \
    sch = OD_DCT_RSHIFT(sc, 1); \
    s3 = sch - s3; \
    sb = s4 - sb; \
    s4 -= OD_DCT_RSHIFT(sb, 1); \
    sa += s5; \
    sah = OD_DCT_RSHIFT(sa, 1); \
    s5 = sah - s5; \
    s9 = s6 - s9; \
    s6 -= OD_DCT_RSHIFT(s9, 1); \
    s8 += s7; \
    s8h = OD_DCT_RSHIFT(s8, 1); \
    s7 = s8h - s7; \
    OD_FDCT_8_ASYM(s0, s8, s8h, s4, sc, sch, s2, sa, sah, s6, se, seh); \
    OD_FDST_8_ASYM(sf, s7, sb, s3, sd, s5, s9, s1); \
  } \
  while (0)

#define OD_IDCT_16(s0, s8, s4, sc, s2, sa, s6, se, \
 s1, s9, s5, sd, s3, sb, s7, sf) \
  /* Embedded 16-point orthonormal Type-II iDCT. */ \
  do { \
    int s1h; \
    int s3h; \
    int s5h; \
    int s7h; \
    int sfh; \
    OD_IDST_8_ASYM(sf, sb, sd, s9, se, sa, sc, s8); \
    OD_IDCT_8_ASYM(s0, s4, s2, s6, s1, s1h, s5, s5h, s3, s3h, s7, s7h); \
    sfh = OD_DCT_RSHIFT(sf, 1); \
    s0 += sfh; \
    sf = s0 - sf; \
    se = s1h - se; \
    s1 -= se; \
    s2 += OD_DCT_RSHIFT(sd, 1); \
    sd = s2 - sd; \
    sc = s3h - sc; \
    s3 -= sc; \
    s4 += OD_DCT_RSHIFT(sb, 1); \
    sb = s4 - sb; \
    sa = s5h - sa; \
    s5 -= sa; \
    s6 += OD_DCT_RSHIFT(s9, 1); \
    s9 = s6 - s9; \
    s8 = s7h - s8; \
    s7 -= s8; \
  } \
  while (0)

#define OD_FDCT_16_ASYM(t0, t8, t8h, t4, tc, tch, t2, ta, tah, t6, te, teh, \
 t1, t9, t9h, t5, td, tdh, t3, tb, tbh, t7, tf, tfh) \
  /* Embedded 16-point asymmetric Type-II fDCT. */ \
  do { \
    t0 += tfh; \
    tf = t0 - tf; \
    t1 -= teh; \
    te += t1; \
    t2 += tdh; \
    td = t2 - td; \
    t3 -= tch; \
    tc += t3; \
    t4 += tbh; \
    tb = t4 - tb; \
    t5 -= tah; \
    ta += t5; \
    t6 += t9h; \
    t9 = t6 - t9; \
    t7 -= t8h; \
    t8 += t7; \
    OD_FDCT_8(t0, t8, t4, tc, t2, ta, t6, te); \
    OD_FDST_8(tf, t7, tb, t3, td, t5, t9, t1); \
  } \
  while (0)

#define OD_IDCT_16_ASYM(t0, t8, t4, tc, t2, ta, t6, te, \
 t1, t1h, t9, t9h, t5, t5h, td, tdh, t3, t3h, tb, tbh, t7, t7h, tf, tfh) \
  /* Embedded 16-point asymmetric Type-II iDCT. */ \
  do { \
    OD_IDST_8(tf, tb, td, t9, te, ta, tc, t8); \
    OD_IDCT_8(t0, t4, t2, t6, t1, t5, t3, t7); \
    t1 -= te; \
    t1h = OD_DCT_RSHIFT(t1, 1); \
    te += t1h; \
    t9 = t6 - t9; \
    t9h = OD_DCT_RSHIFT(t9, 1); \
    t6 -= t9h; \
    t5 -= ta; \
    t5h = OD_DCT_RSHIFT(t5, 1); \
    ta += t5h; \
    td = t2 - td; \
    tdh = OD_DCT_RSHIFT(td, 1); \
    t2 -= tdh; \
    t3 -= tc; \
    t3h = OD_DCT_RSHIFT(t3, 1); \
    tc += t3h; \
    tb = t4 - tb; \
    tbh = OD_DCT_RSHIFT(tb, 1); \
    t4 -= tbh; \
    t7 -= t8; \
    t7h = OD_DCT_RSHIFT(t7, 1); \
    t8 += t7h; \
    tf = t0 - tf; \
    tfh = OD_DCT_RSHIFT(tf, 1); \
    t0 -= tfh; \
  } \
  while (0)

#define OD_FDST_16(s0, s8, s4, sc, s2, sa, s6, se, \
 s1, s9, s5, sd, s3, sb, s7, sf) \
  /* Embedded 16-point orthonormal Type-IV fDST. */ \
  do { \
    int s0h; \
    int s2h; \
    int sdh; \
    int sfh; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(s3, 13573, 16384, 220); \
    s1 += (se*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(s1, 11585, 8192, 221); \
    se -= (s1*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(s3, 13573, 16384, 222); \
    s1 += (se*13573 + 16384) >> 15; \
    /* 21895/32768 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(s2, 21895, 16384, 223); \
    sd += (s2*21895 + 16384) >> 15; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(sd, 15137, 16384, 224); \
    s2 -= (sd*15137 + 8192) >> 14; \
    /* 21895/32768 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    OD_DCT_OVERFLOW_CHECK(s2, 21895, 16384, 225); \
    sd += (s2*21895 + 16384) >> 15; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    OD_DCT_OVERFLOW_CHECK(s3, 3259, 8192, 226); \
    sc += (s3*3259 + 8192) >> 14; \
    /* 3135/8192 ~= Sin[Pi/8] ~= 0.382683432365090 */ \
    OD_DCT_OVERFLOW_CHECK(sc, 3135, 4096, 227); \
    s3 -= (sc*3135 + 4096) >> 13; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    OD_DCT_OVERFLOW_CHECK(s3, 3259, 8192, 228); \
    sc += (s3*3259 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(s5, 13573, 16384, 229); \
    sa += (s5*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(sa, 11585, 8192, 230); \
    s5 -= (sa*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[Pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(s5, 13573, 16384, 231); \
    sa += (s5*13573 + 16384) >> 15; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(s9, 13573, 16384, 232); \
    s6 += (s9*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[pi/4] ~= 0.707106781186547 */ \
    OD_DCT_OVERFLOW_CHECK(s6, 11585, 8192, 233); \
    s9 -= (s6*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(s9, 13573, 16384, 234); \
    s6 += (s9*13573 + 16384) >> 15; \
    sf += se; \
    sfh = OD_DCT_RSHIFT(sf, 1); \
    se = sfh - se; \
    s0 += s1; \
    s0h = OD_DCT_RSHIFT(s0, 1); \
    s1 = s0h - s1; \
    s2 = s3 - s2; \
    s2h = OD_DCT_RSHIFT(s2, 1); \
    s3 -= s2h; \
    sd -= sc; \
    sdh = OD_DCT_RSHIFT(sd, 1); \
    sc += sdh; \
    sa = s4 - sa; \
    s4 -= OD_DCT_RSHIFT(sa, 1); \
    s5 += sb; \
    sb = OD_DCT_RSHIFT(s5, 1) - sb; \
    s8 += s6; \
    s6 -= OD_DCT_RSHIFT(s8, 1); \
    s7 = s9 - s7; \
    s9 -= OD_DCT_RSHIFT(s7, 1); \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(sb, 6723, 4096, 235); \
    s4 += (sb*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    OD_DCT_OVERFLOW_CHECK(s4, 16069, 8192, 236); \
    sb -= (s4*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32]) ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(sb, 6723, 4096, 237); \
    s4 += (sb*6723 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32]) ~= 0.534511135950792 */ \
    OD_DCT_OVERFLOW_CHECK(s5, 8757, 8192, 238); \
    sa += (s5*8757 + 8192) >> 14; \
    /* 6811/8192 ~= Sin[5*Pi/16] ~= 0.831469612302545 */ \
    OD_DCT_OVERFLOW_CHECK(sa, 6811, 4096, 239); \
    s5 -= (sa*6811 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    OD_DCT_OVERFLOW_CHECK(s5, 8757, 8192, 240); \
    sa += (s5*8757 + 8192) >> 14; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    OD_DCT_OVERFLOW_CHECK(s9, 2485, 4096, 241); \
    s6 += (s9*2485 + 4096) >> 13; \
    /* 4551/8192 ~= Sin[3*Pi/16] ~= 0.555570233019602 */ \
    OD_DCT_OVERFLOW_CHECK(s6, 4551, 4096, 242); \
    s9 -= (s6*4551 + 4096) >> 13; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    OD_DCT_OVERFLOW_CHECK(s9, 2485, 4096, 243); \
    s6 += (s9*2485 + 4096) >> 13; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    OD_DCT_OVERFLOW_CHECK(s8, 3227, 16384, 244); \
    s7 += (s8*3227 + 16384) >> 15; \
    /* 6393/32768 ~= Sin[Pi/16] ~= 0.19509032201612825 */ \
    OD_DCT_OVERFLOW_CHECK(s7, 6393, 16384, 245); \
    s8 -= (s7*6393 + 16384) >> 15; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    OD_DCT_OVERFLOW_CHECK(s8, 3227, 16384, 246); \
    s7 += (s8*3227 + 16384) >> 15; \
    s1 -= s2h; \
    s2 += s1; \
    se += sdh; \
    sd = se - sd; \
    s3 += sfh; \
    sf -= s3; \
    sc = s0h - sc; \
    s0 -= sc; \
    sb += OD_DCT_RSHIFT(s8, 1); \
    s8 = sb - s8; \
    s4 += OD_DCT_RSHIFT(s7, 1); \
    s7 -= s4; \
    s6 += OD_DCT_RSHIFT(s5, 1); \
    s5 = s6 - s5; \
    s9 -= OD_DCT_RSHIFT(sa, 1); \
    sa += s9; \
    s8 += s0; \
    s0 -= OD_DCT_RSHIFT(s8, 1); \
    sf += s7; \
    s7 = OD_DCT_RSHIFT(sf, 1) - s7; \
    s1 -= s6; \
    s6 += OD_DCT_RSHIFT(s1, 1); \
    s9 += se; \
    se = OD_DCT_RSHIFT(s9, 1) - se; \
    s2 += sa; \
    sa = OD_DCT_RSHIFT(s2, 1) - sa; \
    s5 += sd; \
    sd -= OD_DCT_RSHIFT(s5, 1); \
    s4 = sc - s4; \
    sc -= OD_DCT_RSHIFT(s4, 1); \
    s3 -= sb; \
    sb += OD_DCT_RSHIFT(s3, 1); \
    /* 2799/4096 ~= (1/Sqrt[2] - Cos[31*Pi/64]/2)/Sin[31*Pi/64] ~=
        0.6833961245841154 */ \
    OD_DCT_OVERFLOW_CHECK(sf, 2799, 2048, 247); \
    s0 -= (sf*2799 + 2048) >> 12; \
    /* 2893/2048 ~= Sqrt[2]*Sin[31*Pi/64] ~= 1.4125100802019777 */ \
    OD_DCT_OVERFLOW_CHECK(s0, 2893, 1024, 248); \
    sf += (s0*2893 + 1024) >> 11; \
    /* 5397/8192 ~= (Cos[Pi/4] - Cos[31*Pi/64])/Sin[31*Pi/64] ~=
        0.6588326996993819 */ \
    OD_DCT_OVERFLOW_CHECK(sf, 5397, 4096, 249); \
    s0 -= (sf*5397 + 4096) >> 13; \
    /* 41/64 ~= (1/Sqrt[2] - Cos[29*Pi/64]/2)/Sin[29*Pi/64] ~=
        0.6406758931036793 */ \
    OD_DCT_OVERFLOW_CHECK(s1, 41, 32, 250); \
    se += (s1*41 + 32) >> 6; \
    /* 2865/2048 ~= Sqrt[2]*Sin[29*Pi/64] ~= 1.3989068359730783 */ \
    OD_DCT_OVERFLOW_CHECK(se, 2865, 1024, 251); \
    s1 -= (se*2865 + 1024) >> 11; \
    /* 4641/8192 ~= (1/Sqrt[2] - Cos[29*Pi/64])/Sin[29*Pi/64] ~=
        0.5665078993345056 */ \
    OD_DCT_OVERFLOW_CHECK(s1, 4641, 4096, 252); \
    se += (s1*4641 + 4096) >> 13; \
    /* 2473/4096 ~= (1/Sqrt[2] - Cos[27*Pi/64]/2)/Sin[27*Pi/64] ~=
        0.603709096285651 */ \
    OD_DCT_OVERFLOW_CHECK(s2, 2473, 2048, 253); \
    sd += (s2*2473 + 2048) >> 12; \
    /* 5619/4096 ~= Sqrt[2]*Sin[27*Pi/64] ~= 1.371831354193494 */ \
    OD_DCT_OVERFLOW_CHECK(sd, 5619, 2048, 254); \
    s2 -= (sd*5619 + 2048) >> 12; \
    /* 7839/16384 ~= (1/Sqrt[2] - Cos[27*Pi/64])/Sin[27*Pi/64] ~=
        0.47846561618999817 */ \
    OD_DCT_OVERFLOW_CHECK(s2, 7839, 8192, 255); \
    sd += (s2*7839 + 8192) >> 14; \
    /* 5747/8192 ~= (1/Sqrt[2] - Cos[7*Pi/64]/2)/Sin[7*Pi/64] ~=
        0.7015193429405162 */ \
    OD_DCT_OVERFLOW_CHECK(s3, 5747, 4096, 256); \
    sc -= (s3*5747 + 4096) >> 13; \
    /* 3903/8192 ~= Sqrt[2]*Sin[7*Pi/64] ~= 0.47643419969316125 */ \
    OD_DCT_OVERFLOW_CHECK(sc, 3903, 4096, 257); \
    s3 += (sc*3903 + 4096) >> 13; \
    /* 5701/8192 ~= (1/Sqrt[2] - Cos[7*Pi/64])/Sin[7*Pi/64] ~=
        0.6958870433047222 */ \
    OD_DCT_OVERFLOW_CHECK(s3, 5701, 4096, 258); \
    sc += (s3*5701 + 4096) >> 13; \
    /* 4471/8192 ~= (1/Sqrt[2] - Cos[23*Pi/64]/2)/Sin[23*Pi/64] ~=
        0.5457246432276498 */ \
    OD_DCT_OVERFLOW_CHECK(s4, 4471, 4096, 259); \
    sb += (s4*4471 + 4096) >> 13; \
    /* 1309/1024 ~= Sqrt[2]*Sin[23*Pi/64] ~= 1.278433918575241 */ \
    OD_DCT_OVERFLOW_CHECK(sb, 1309, 512, 260); \
    s4 -= (sb*1309 + 512) >> 10; \
    /* 5067/16384 ~= (1/Sqrt[2] - Cos[23*Pi/64])/Sin[23*Pi/64] ~=
        0.30924225528198984 */ \
    OD_DCT_OVERFLOW_CHECK(s4, 5067, 8192, 261); \
    sb += (s4*5067 + 8192) >> 14; \
    /* 2217/4096 ~= (1/Sqrt[2] - Cos[11*Pi/64]/2)/Sin[11*Pi/64] ~=
        0.5412195895259334 */ \
    OD_DCT_OVERFLOW_CHECK(s5, 2217, 2048, 262); \
    sa -= (s5*2217 + 2048) >> 12; \
    /* 1489/2048 ~= Sqrt[2]*Sin[11*Pi/64] ~= 0.72705107329128 */ \
    OD_DCT_OVERFLOW_CHECK(sa, 1489, 1024, 263); \
    s5 += (sa*1489 + 1024) >> 11; \
    /* 75/256 ~= (1/Sqrt[2] - Cos[11*Pi/64])/Sin[11*Pi/64] ~=
        0.2929800132658202 */ \
    OD_DCT_OVERFLOW_CHECK(s5, 75, 128, 264); \
    sa += (s5*75 + 128) >> 8; \
    /* 2087/4096 ~= (1/Sqrt[2] - Cos[19*Pi/64]/2)/Sin[19*Pi/64] ~=
        0.5095285002941893 */ \
    OD_DCT_OVERFLOW_CHECK(s9, 2087, 2048, 265); \
    s6 -= (s9*2087 + 2048) >> 12; \
    /* 4653/4096 ~= Sqrt[2]*Sin[19*Pi/64] ~= 1.1359069844201428 */ \
    OD_DCT_OVERFLOW_CHECK(s6, 4653, 2048, 266); \
    s9 += (s6*4653 + 2048) >> 12; \
    /* 4545/32768 ~= (1/Sqrt[2] - Cos[19*Pi/64])/Sin[19*Pi/64] ~=
        0.13870322715817154 */ \
    OD_DCT_OVERFLOW_CHECK(s9, 4545, 16384, 267); \
    s6 -= (s9*4545 + 16384) >> 15; \
    /* 2053/4096 ~= (1/Sqrt[2] - Cos[15*Pi/64]/2)/Sin[15*Pi/64] ~=
        0.5012683042634027 */ \
    OD_DCT_OVERFLOW_CHECK(s8, 2053, 2048, 268); \
    s7 += (s8*2053 + 2048) >> 12; \
    /* 1945/2048 ~= Sqrt[2]*Sin[15*Pi/64] ~= 0.9497277818777543 */ \
    OD_DCT_OVERFLOW_CHECK(s7, 1945, 1024, 269); \
    s8 -= (s7*1945 + 1024) >> 11; \
    /* 1651/32768 ~= (1/Sqrt[2] - Cos[15*Pi/64])/Sin[15*Pi/64] ~=
        0.05039668360333519 */ \
    OD_DCT_OVERFLOW_CHECK(s8, 1651, 16384, 270); \
    s7 -= (s8*1651 + 16384) >> 15; \
  } \
  while (0)

#define OD_IDST_16(s0, s8, s4, sc, s2, sa, s6, se, \
 s1, s9, s5, sd, s3, sb, s7, sf) \
  /* Embedded 16-point orthonormal Type-IV iDST. */ \
  do { \
    int s0h; \
    int s4h; \
    int sbh; \
    int sfh; \
    /* 1651/32768 ~= (1/Sqrt[2] - Cos[15*Pi/64])/Sin[15*Pi/64] ~=
        0.05039668360333519 */ \
    se += (s1*1651 + 16384) >> 15; \
    /* 1945/2048 ~= Sqrt[2]*Sin[15*Pi/64] ~= 0.9497277818777543 */ \
    s1 += (se*1945 + 1024) >> 11; \
    /* 2053/4096 ~= (1/Sqrt[2] - Cos[15*Pi/64]/2)/Sin[15*Pi/64] ~=
        0.5012683042634027 */ \
    se -= (s1*2053 + 2048) >> 12; \
    /* 4545/32768 ~= (1/Sqrt[2] - Cos[19*Pi/64])/Sin[19*Pi/64] ~=
        0.13870322715817154 */ \
    s6 += (s9*4545 + 16384) >> 15; \
    /* 4653/32768 ~= Sqrt[2]*Sin[19*Pi/64] ~= 1.1359069844201428 */ \
    s9 -= (s6*4653 + 2048) >> 12; \
    /* 2087/4096 ~= (1/Sqrt[2] - Cos[19*Pi/64]/2)/Sin[19*Pi/64] ~=
        0.5095285002941893 */ \
    s6 += (s9*2087 + 2048) >> 12; \
    /* 75/256 ~= (1/Sqrt[2] - Cos[11*Pi/64])/Sin[11*Pi/64] ~=
        0.2929800132658202 */ \
    s5 -= (sa*75 + 128) >> 8; \
    /* 1489/2048 ~= Sqrt[2]*Sin[11*Pi/64] ~= 0.72705107329128 */ \
    sa -= (s5*1489 + 1024) >> 11; \
    /* 2217/4096 ~= (1/Sqrt[2] - Cos[11*Pi/64]/2)/Sin[11*Pi/64] ~=
        0.5412195895259334 */ \
    s5 += (sa*2217 + 2048) >> 12; \
    /* 5067/16384 ~= (1/Sqrt[2] - Cos[23*Pi/64])/Sin[23*Pi/64] ~=
        0.30924225528198984 */ \
    sd -= (s2*5067 + 8192) >> 14; \
    /* 1309/1024 ~= Sqrt[2]*Sin[23*Pi/64] ~= 1.278433918575241 */ \
    s2 += (sd*1309 + 512) >> 10; \
    /* 4471/8192 ~= (1/Sqrt[2] - Cos[23*Pi/64]/2)/Sin[23*Pi/64] ~=
        0.5457246432276498 */ \
    sd -= (s2*4471 + 4096) >> 13; \
    /* 5701/8192 ~= (1/Sqrt[2] - Cos[7*Pi/64])/Sin[7*Pi/64] ~=
        0.6958870433047222 */ \
    s3 -= (sc*5701 + 4096) >> 13; \
    /* 3903/8192 ~= Sqrt[2]*Sin[7*Pi/64] ~= 0.47643419969316125 */ \
    sc -= (s3*3903 + 4096) >> 13; \
    /* 5747/8192 ~= (1/Sqrt[2] - Cos[7*Pi/64]/2)/Sin[7*Pi/64] ~=
        0.7015193429405162 */ \
    s3 += (sc*5747 + 4096) >> 13; \
    /* 7839/16384 ~= (1/Sqrt[2] - Cos[27*Pi/64])/Sin[27*Pi/64] ~=
        0.47846561618999817 */ \
    sb -= (s4*7839 + 8192) >> 14; \
    /* 5619/4096 ~= Sqrt[2]*Sin[27*Pi/64] ~= 1.371831354193494 */ \
    s4 += (sb*5619 + 2048) >> 12; \
    /* 2473/4096 ~= (1/Sqrt[2] - Cos[27*Pi/64]/2)/Sin[27*Pi/64] ~=
        0.603709096285651 */ \
    sb -= (s4*2473 + 2048) >> 12; \
    /* 4641/8192 ~= (1/Sqrt[2] - Cos[29*Pi/64])/Sin[29*Pi/64] ~=
        0.5665078993345056 */ \
    s7 -= (s8*4641 + 4096) >> 13; \
    /* 2865/2048 ~= Sqrt[2]*Sin[29*Pi/64] ~= 1.3989068359730783 */ \
    s8 += (s7*2865 + 1024) >> 11; \
    /* 41/64 ~= (1/Sqrt[2] - Cos[29*Pi/64]/2)/Sin[29*Pi/64] ~=
        0.6406758931036793 */ \
    s7 -= (s8*41 + 32) >> 6; \
    /* 5397/8192 ~= (Cos[Pi/4] - Cos[31*Pi/64])/Sin[31*Pi/64] ~=
        0.6588326996993819 */ \
    s0 += (sf*5397 + 4096) >> 13; \
    /* 2893/2048 ~= Sqrt[2]*Sin[31*Pi/64] ~= 1.4125100802019777 */ \
    sf -= (s0*2893 + 1024) >> 11; \
    /* 2799/4096 ~= (1/Sqrt[2] - Cos[31*Pi/64]/2)/Sin[31*Pi/64] ~=
        0.6833961245841154 */ \
    s0 += (sf*2799 + 2048) >> 12; \
    sd -= OD_DCT_RSHIFT(sc, 1); \
    sc += sd; \
    s3 += OD_DCT_RSHIFT(s2, 1); \
    s2 = s3 - s2; \
    sb += OD_DCT_RSHIFT(sa, 1); \
    sa -= sb; \
    s5 = OD_DCT_RSHIFT(s4, 1) - s5; \
    s4 -= s5; \
    s7 = OD_DCT_RSHIFT(s9, 1) - s7; \
    s9 -= s7; \
    s6 -= OD_DCT_RSHIFT(s8, 1); \
    s8 += s6; \
    se = OD_DCT_RSHIFT(sf, 1) - se; \
    sf -= se; \
    s0 += OD_DCT_RSHIFT(s1, 1); \
    s1 -= s0; \
    s5 -= s9; \
    s9 += OD_DCT_RSHIFT(s5, 1); \
    sa = s6 - sa; \
    s6 -= OD_DCT_RSHIFT(sa, 1); \
    se += s2; \
    s2 -= OD_DCT_RSHIFT(se, 1); \
    s1 = sd - s1; \
    sd -= OD_DCT_RSHIFT(s1, 1); \
    s0 += s3; \
    s0h = OD_DCT_RSHIFT(s0, 1); \
    s3 = s0h - s3; \
    sf += sc; \
    sfh = OD_DCT_RSHIFT(sf, 1); \
    sc -= sfh; \
    sb = s7 - sb; \
    sbh = OD_DCT_RSHIFT(sb, 1); \
    s7 -= sbh; \
    s4 -= s8; \
    s4h = OD_DCT_RSHIFT(s4, 1); \
    s8 += s4h; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    se -= (s1*3227 + 16384) >> 15; \
    /* 6393/32768 ~= Sin[Pi/16] ~= 0.19509032201612825 */ \
    s1 += (se*6393 + 16384) >> 15; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    se -= (s1*3227 + 16384) >> 15; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    s6 -= (s9*2485 + 4096) >> 13; \
    /* 4551/8192 ~= Sin[3*Pi/16] ~= 0.555570233019602 */ \
    s9 += (s6*4551 + 4096) >> 13; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    s6 -= (s9*2485 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    s5 -= (sa*8757 + 8192) >> 14; \
    /* 6811/8192 ~= Sin[5*Pi/16] ~= 0.831469612302545 */ \
    sa += (s5*6811 + 4096) >> 13; \
    /* 8757/16384 ~= Tan[5*Pi/32]) ~= 0.534511135950792 */ \
    s5 -= (sa*8757 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32]) ~= 0.820678790828660 */ \
    s2 -= (sd*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    sd += (s2*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    s2 -= (sd*6723 + 4096) >> 13; \
    s9 += OD_DCT_RSHIFT(se, 1); \
    se = s9 - se; \
    s6 += OD_DCT_RSHIFT(s1, 1); \
    s1 -= s6; \
    sd = OD_DCT_RSHIFT(sa, 1) - sd; \
    sa -= sd; \
    s2 += OD_DCT_RSHIFT(s5, 1); \
    s5 = s2 - s5; \
    s3 -= sbh; \
    sb += s3; \
    sc += s4h; \
    s4 = sc - s4; \
    s8 = s0h - s8; \
    s0 -= s8; \
    s7 = sfh - s7; \
    sf -= s7; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    s6 -= (s9*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[pi/4] ~= 0.707106781186547 */ \
    s9 += (s6*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    s6 -= (s9*13573 + 16384) >> 15; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    s5 -= (sa*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[pi/4] ~= 0.707106781186547 */ \
    sa += (s5*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    s5 -= (sa*13573 + 16384) >> 15; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    s3 -= (sc*3259 + 8192) >> 14; \
    /* 3135/8192 ~= Sin[Pi/8] ~= 0.382683432365090 */ \
    sc += (s3*3135 + 4096) >> 13; \
    /* 3259/16384 ~= Tan[Pi/16] ~= 0.198912367379658 */ \
    s3 -= (sc*3259 + 8192) >> 14; \
    /* 21895/32768 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    sb -= (s4*21895 + 16384) >> 15; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    s4 += (sb*15137 + 8192) >> 14; \
    /* 21895/32768 ~= Tan[3*Pi/16] ~= 0.668178637919299 */ \
    sb -= (s4*21895 + 16384) >> 15; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    s8 -= (s7*13573 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[pi/4] ~= 0.707106781186547 */ \
    s7 += (s8*11585 + 8192) >> 14; \
    /* 13573/32768 ~= Tan[pi/8] ~= 0.414213562373095 */ \
    s8 -= (s7*13573 + 16384) >> 15; \
  } \
  while (0)

/* TODO: rewrite this to match OD_FDST_16. */
#define OD_FDST_16_ASYM(t0, t0h, t8, t4, t4h, tc, t2, ta, t6, te, \
 t1, t9, t5, td, t3, tb, t7, t7h, tf) \
  /* Embedded 16-point asymmetric Type-IV fDST. */ \
  do { \
    int t2h; \
    int t3h; \
    int t6h; \
    int t8h; \
    int t9h; \
    int tch; \
    int tdh; \
    /* TODO: Can we move these into another operation */ \
    t8 = -t8; \
    t9 = -t9; \
    ta = -ta; \
    tb = -tb; \
    td = -td; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    OD_DCT_OVERFLOW_CHECK(te, 13573, 8192, 136); \
    t1 -= (te*13573 + 8192) >> 14; \
    /* 11585/32768 ~= Sin[Pi/4]/2 ~= 0.353553390593274 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 11585, 16384, 137); \
    te += (t1*11585 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    OD_DCT_OVERFLOW_CHECK(te, 13573, 8192, 138); \
    t1 -= (te*13573 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(td, 4161, 8192, 139); \
    t2 += (td*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 15137, 8192, 140); \
    td -= (t2*15137 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(td, 14341, 8192, 141); \
    t2 += (td*14341 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 14341, 8192, 142); \
    tc -= (t3*14341 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 15137, 8192, 143); \
    t3 += (tc*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 4161, 8192, 144); \
    tc -= (t3*4161 + 8192) >> 14; \
    te = t0h - te; \
    t0 -= te; \
    tf = OD_DCT_RSHIFT(t1, 1) - tf; \
    t1 -= tf; \
    /* TODO: Can we move this into another operation */ \
    tc = -tc; \
    t2 = OD_DCT_RSHIFT(tc, 1) - t2; \
    tc -= t2; \
    t3 = OD_DCT_RSHIFT(td, 1) - t3; \
    td = t3 - td; \
    /* 7489/8192 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 7489, 4096, 145); \
    t9 -= (t6*7489 + 4096) >> 13; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 11585, 8192, 146); \
    t6 += (t9*11585 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 19195, 16384, 147); \
    t9 += (t6*19195 + 16384) >> 15; \
    t8 += OD_DCT_RSHIFT(t9, 1); \
    t9 -= t8; \
    t6 = t7h - t6; \
    t7 -= t6; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 6723, 4096, 148); \
    t8 += (t7*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    OD_DCT_OVERFLOW_CHECK(t8, 16069, 8192, 149); \
    t7 -= (t8*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32]) ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 6723, 4096, 150); \
    t8 += (t7*6723 + 4096) >> 13; \
    /* 17515/32768 ~= Tan[5*Pi/32]) ~= 0.534511135950792 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 17515, 16384, 151); \
    t9 += (t6*17515 + 16384) >> 15; \
    /* 13623/16384 ~= Sin[5*Pi/16] ~= 0.831469612302545 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 13623, 8192, 152); \
    t6 -= (t9*13623 + 8192) >> 14; \
    /* 17515/32768 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 17515, 16384, 153); \
    t9 += (t6*17515 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 13573, 8192, 154); \
    t5 += (ta*13573 + 8192) >> 14; \
    /* 11585/32768 ~= Sin[Pi/4]/2 ~= 0.353553390593274 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 11585, 16384, 155); \
    ta -= (t5*11585 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 13573, 8192, 156); \
    t5 += (ta*13573 + 8192) >> 14; \
    tb += OD_DCT_RSHIFT(t5, 1); \
    t5 = tb - t5; \
    ta += t4h; \
    t4 -= ta; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 2485, 4096, 157); \
    ta += (t5*2485 + 4096) >> 13; \
    /* 18205/32768 ~= Sin[3*Pi/16] ~= 0.555570233019602 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 18205, 16384, 158); \
    t5 -= (ta*18205 + 16384) >> 15; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 2485, 4096, 159); \
    ta += (t5*2485 + 4096) >> 13; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 6723, 4096, 160); \
    tb -= (t4*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 16069, 8192, 161); \
    t4 += (tb*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 6723, 4096, 162); \
    tb -= (t4*6723 + 4096) >> 13; \
    /* TODO: Can we move this into another operation */ \
    t5 = -t5; \
    tc -= tf; \
    tch = OD_DCT_RSHIFT(tc, 1); \
    tf += tch; \
    t3 += t0; \
    t3h = OD_DCT_RSHIFT(t3, 1); \
    t0 -= t3h; \
    td -= t1; \
    tdh = OD_DCT_RSHIFT(td, 1); \
    t1 += tdh; \
    t2 += te; \
    t2h = OD_DCT_RSHIFT(t2, 1); \
    te -= t2h; \
    t8 += t4; \
    t8h = OD_DCT_RSHIFT(t8, 1); \
    t4 = t8h - t4; \
    t7 = tb - t7; \
    t7h = OD_DCT_RSHIFT(t7, 1); \
    tb = t7h - tb; \
    t6 -= ta; \
    t6h = OD_DCT_RSHIFT(t6, 1); \
    ta += t6h; \
    t9 = t5 - t9; \
    t9h = OD_DCT_RSHIFT(t9, 1); \
    t5 -= t9h; \
    t0 -= t7h; \
    t7 += t0; \
    tf += t8h; \
    t8 -= tf; \
    te -= t6h; \
    t6 += te; \
    t1 += t9h; \
    t9 -= t1; \
    tb -= tch; \
    tc += tb; \
    t4 += t3h; \
    t3 -= t4; \
    ta -= tdh; \
    td += ta; \
    t5 = t2h - t5; \
    t2 -= t5; \
    /* TODO: Can we move these into another operation */ \
    t8 = -t8; \
    t9 = -t9; \
    ta = -ta; \
    tb = -tb; \
    tc = -tc; \
    td = -td; \
    tf = -tf; \
    /* 7799/8192 ~= Tan[31*Pi/128] ~= 0.952079146700925 */ \
    OD_DCT_OVERFLOW_CHECK(tf, 7799, 4096, 163); \
    t0 -= (tf*7799 + 4096) >> 13; \
    /* 4091/4096 ~= Sin[31*Pi/64] ~= 0.998795456205172 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 4091, 2048, 164); \
    tf += (t0*4091 + 2048) >> 12; \
    /* 7799/8192 ~= Tan[31*Pi/128] ~= 0.952079146700925 */ \
    OD_DCT_OVERFLOW_CHECK(tf, 7799, 4096, 165); \
    t0 -= (tf*7799 + 4096) >> 13; \
    /* 2417/32768 ~= Tan[3*Pi/128] ~= 0.0737644315224493 */ \
    OD_DCT_OVERFLOW_CHECK(te, 2417, 16384, 166); \
    t1 += (te*2417 + 16384) >> 15; \
    /* 601/4096 ~= Sin[3*Pi/64] ~= 0.146730474455362 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 601, 2048, 167); \
    te -= (t1*601 + 2048) >> 12; \
    /* 2417/32768 ~= Tan[3*Pi/128] ~= 0.0737644315224493 */ \
    OD_DCT_OVERFLOW_CHECK(te, 2417, 16384, 168); \
    t1 += (te*2417 + 16384) >> 15; \
    /* 14525/32768 ~= Tan[17*Pi/128] ~= 0.443269513890864 */ \
    OD_DCT_OVERFLOW_CHECK(t8, 14525, 16384, 169); \
    t7 -= (t8*14525 + 16384) >> 15; \
    /* 3035/4096 ~= Sin[17*Pi/64] ~= 0.740951125354959 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 3035, 2048, 170); \
    t8 += (t7*3035 + 2048) >> 12; \
    /* 7263/16384 ~= Tan[17*Pi/128] ~= 0.443269513890864 */ \
    OD_DCT_OVERFLOW_CHECK(t8, 7263, 8192, 171); \
    t7 -= (t8*7263 + 8192) >> 14; \
    /* 6393/8192 ~= Tan[27*Pi/128] ~= 0.780407659653944 */ \
    OD_DCT_OVERFLOW_CHECK(td, 6393, 4096, 172); \
    t2 -= (td*6393 + 4096) >> 13; \
    /* 3973/4096 ~= Sin[27*Pi/64] ~= 0.970031253194544 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 3973, 2048, 173); \
    td += (t2*3973 + 2048) >> 12; \
    /* 6393/8192 ~= Tan[27*Pi/128] ~= 0.780407659653944 */ \
    OD_DCT_OVERFLOW_CHECK(td, 6393, 4096, 174); \
    t2 -= (td*6393 + 4096) >> 13; \
    /* 9281/16384 ~= Tan[21*Pi/128] ~= 0.566493002730344 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 9281, 8192, 175); \
    t5 -= (ta*9281 + 8192) >> 14; \
    /* 7027/8192 ~= Sin[21*Pi/64] ~= 0.857728610000272 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 7027, 4096, 176); \
    ta += (t5*7027 + 4096) >> 13; \
    /* 9281/16384 ~= Tan[21*Pi/128] ~= 0.566493002730344 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 9281, 8192, 177); \
    t5 -= (ta*9281 + 8192) >> 14; \
    /* 11539/16384 ~= Tan[25*Pi/128] ~= 0.704279460865044 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 11539, 8192, 178); \
    t3 -= (tc*11539 + 8192) >> 14; \
    /* 7713/8192 ~= Sin[25*Pi/64] ~= 0.941544065183021 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 7713, 4096, 179); \
    tc += (t3*7713 + 4096) >> 13; \
    /* 11539/16384 ~= Tan[25*Pi/128] ~= 0.704279460865044 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 11539, 8192, 180); \
    t3 -= (tc*11539 + 8192) >> 14; \
    /* 10375/16384 ~= Tan[23*Pi/128] ~= 0.633243016177569 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 10375, 8192, 181); \
    t4 -= (tb*10375 + 8192) >> 14; \
    /* 7405/8192 ~= Sin[23*Pi/64] ~= 0.903989293123443 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 7405, 4096, 182); \
    tb += (t4*7405 + 4096) >> 13; \
    /* 10375/16384 ~= Tan[23*Pi/128] ~= 0.633243016177569 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 10375, 8192, 183); \
    t4 -= (tb*10375 + 8192) >> 14; \
    /* 8247/16384 ~= Tan[19*Pi/128] ~= 0.503357699799294 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 8247, 8192, 184); \
    t6 -= (t9*8247 + 8192) >> 14; \
    /* 1645/2048 ~= Sin[19*Pi/64] ~= 0.803207531480645 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 1645, 1024, 185); \
    t9 += (t6*1645 + 1024) >> 11; \
    /* 8247/16384 ~= Tan[19*Pi/128] ~= 0.503357699799294 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 8247, 8192, 186); \
    t6 -= (t9*8247 + 8192) >> 14; \
  } \
  while (0)

#define OD_IDST_16_ASYM(t0, t0h, t8, t4, tc, t2, t2h, ta, t6, te, teh, \
 t1, t9, t5, td, t3, tb, t7, tf) \
  /* Embedded 16-point asymmetric Type-IV iDST. */ \
  do { \
    int t1h_; \
    int t3h_; \
    int t4h; \
    int t6h; \
    int t9h_; \
    int tbh_; \
    int tch; \
    /* 8247/16384 ~= Tan[19*Pi/128] ~= 0.503357699799294 */ \
    t6 += (t9*8247 + 8192) >> 14; \
    /* 1645/2048 ~= Sin[19*Pi/64] ~= 0.803207531480645 */ \
    t9 -= (t6*1645 + 1024) >> 11; \
    /* 8247/16384 ~= Tan[19*Pi/128] ~= 0.503357699799294 */ \
    t6 += (t9*8247 + 8192) >> 14; \
    /* 10375/16384 ~= Tan[23*Pi/128] ~= 0.633243016177569 */ \
    t2 += (td*10375 + 8192) >> 14; \
    /* 7405/8192 ~= Sin[23*Pi/64] ~= 0.903989293123443 */ \
    td -= (t2*7405 + 4096) >> 13; \
    /* 10375/16384 ~= Tan[23*Pi/128] ~= 0.633243016177569 */ \
    t2 += (td*10375 + 8192) >> 14; \
    /* 11539/16384 ~= Tan[25*Pi/128] ~= 0.704279460865044 */ \
    tc += (t3*11539 + 8192) >> 14; \
    /* 7713/8192 ~= Sin[25*Pi/64] ~= 0.941544065183021 */ \
    t3 -= (tc*7713 + 4096) >> 13; \
    /* 11539/16384 ~= Tan[25*Pi/128] ~= 0.704279460865044 */ \
    tc += (t3*11539 + 8192) >> 14; \
    /* 9281/16384 ~= Tan[21*Pi/128] ~= 0.566493002730344 */ \
    ta += (t5*9281 + 8192) >> 14; \
    /* 7027/8192 ~= Sin[21*Pi/64] ~= 0.857728610000272 */ \
    t5 -= (ta*7027 + 4096) >> 13; \
    /* 9281/16384 ~= Tan[21*Pi/128] ~= 0.566493002730344 */ \
    ta += (t5*9281 + 8192) >> 14; \
    /* 6393/8192 ~= Tan[27*Pi/128] ~= 0.780407659653944 */ \
    t4 += (tb*6393 + 4096) >> 13; \
    /* 3973/4096 ~= Sin[27*Pi/64] ~= 0.970031253194544 */ \
    tb -= (t4*3973 + 2048) >> 12; \
    /* 6393/8192 ~= Tan[27*Pi/128] ~= 0.780407659653944 */ \
    t4 += (tb*6393 + 4096) >> 13; \
    /* 7263/16384 ~= Tan[17*Pi/128] ~= 0.443269513890864 */ \
    te += (t1*7263 + 8192) >> 14; \
    /* 3035/4096 ~= Sin[17*Pi/64] ~= 0.740951125354959 */ \
    t1 -= (te*3035 + 2048) >> 12; \
    /* 14525/32768 ~= Tan[17*Pi/128] ~= 0.443269513890864 */ \
    te += (t1*14525 + 16384) >> 15; \
    /* 2417/32768 ~= Tan[3*Pi/128] ~= 0.0737644315224493 */ \
    t8 -= (t7*2417 + 16384) >> 15; \
    /* 601/4096 ~= Sin[3*Pi/64] ~= 0.146730474455362 */ \
    t7 += (t8*601 + 2048) >> 12; \
    /* 2417/32768 ~= Tan[3*Pi/128] ~= 0.0737644315224493 */ \
    t8 -= (t7*2417 + 16384) >> 15; \
    /* 7799/8192 ~= Tan[31*Pi/128] ~= 0.952079146700925 */ \
    t0 += (tf*7799 + 4096) >> 13; \
    /* 4091/4096 ~= Sin[31*Pi/64] ~= 0.998795456205172 */ \
    tf -= (t0*4091 + 2048) >> 12; \
    /* 7799/8192 ~= Tan[31*Pi/128] ~= 0.952079146700925 */ \
    t0 += (tf*7799 + 4096) >> 13; \
    /* TODO: Can we move these into another operation */ \
    t1 = -t1; \
    t3 = -t3; \
    t5 = -t5; \
    t9 = -t9; \
    tb = -tb; \
    td = -td; \
    tf = -tf; \
    t4 += ta; \
    t4h = OD_DCT_RSHIFT(t4, 1); \
    ta = t4h - ta; \
    tb -= t5; \
    tbh_ = OD_DCT_RSHIFT(tb, 1); \
    t5 += tbh_; \
    tc += t2; \
    tch = OD_DCT_RSHIFT(tc, 1); \
    t2 -= tch; \
    t3 -= td; \
    t3h_ = OD_DCT_RSHIFT(t3, 1); \
    td += t3h_; \
    t9 += t8; \
    t9h_ = OD_DCT_RSHIFT(t9, 1); \
    t8 -= t9h_; \
    t6 -= t7; \
    t6h = OD_DCT_RSHIFT(t6, 1); \
    t7 += t6h; \
    t1 += tf; \
    t1h_ = OD_DCT_RSHIFT(t1, 1); \
    tf -= t1h_; \
    te -= t0; \
    teh = OD_DCT_RSHIFT(te, 1); \
    t0 += teh; \
    ta += t9h_; \
    t9 = ta - t9; \
    t5 -= t6h; \
    t6 += t5; \
    td = teh - td; \
    te = td - te; \
    t2 = t1h_ - t2; \
    t1 -= t2; \
    t7 += t4h; \
    t4 -= t7; \
    t8 -= tbh_; \
    tb += t8; \
    t0 += tch; \
    tc -= t0; \
    tf -= t3h_; \
    t3 += tf; \
    /* TODO: Can we move this into another operation */ \
    ta = -ta; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    td += (t2*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    t2 -= (td*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.820678790828660 */ \
    td += (t2*6723 + 4096) >> 13; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    t5 -= (ta*2485 + 4096) >> 13; \
    /* 18205/32768 ~= Sin[3*Pi/16] ~= 0.555570233019602 */ \
    ta += (t5*18205 + 16384) >> 15; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    t5 -= (ta*2485 + 4096) >> 13; \
    t2 += t5; \
    t2h = OD_DCT_RSHIFT(t2, 1); \
    t5 -= t2h; \
    ta = td - ta; \
    td -= OD_DCT_RSHIFT(ta, 1); \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    ta -= (t5*13573 + 8192) >> 14; \
    /* 11585/32768 ~= Sin[Pi/4]/2 ~= 0.353553390593274 */ \
    t5 += (ta*11585 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    ta -= (t5*13573 + 8192) >> 14; \
    /* 17515/32768 ~= Tan[5*Pi/32] ~= 0.534511135950792 */ \
    t9 -= (t6*17515 + 16384) >> 15; \
    /* 13623/16384 ~= Sin[5*Pi/16] ~= 0.831469612302545 */ \
    t6 += (t9*13623 + 8192) >> 14; \
    /* 17515/32768 ~= Tan[5*Pi/32]) ~= 0.534511135950792 */ \
    t9 -= (t6*17515 + 16384) >> 15; \
    /* 6723/8192 ~= Tan[7*Pi/32]) ~= 0.820678790828660 */ \
    t1 -= (te*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.980785280403230 */ \
    te += (t1*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32]) ~= 0.820678790828660 */ \
    t1 -= (te*6723 + 4096) >> 13; \
    te += t6; \
    teh = OD_DCT_RSHIFT(te, 1); \
    t6 = teh - t6; \
    t9 += t1; \
    t1 -= OD_DCT_RSHIFT(t9, 1); \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    t9 -= (t6*19195 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    t6 -= (t9*11585 + 8192) >> 14; \
    /* 7489/8192 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    t9 += (t6*7489 + 4096) >> 13; \
    tb = tc - tb; \
    tc = OD_DCT_RSHIFT(tb, 1) - tc; \
    t3 += t4; \
    t4 = OD_DCT_RSHIFT(t3, 1) - t4; \
    /* TODO: Can we move this into another operation */ \
    t3 = -t3; \
    t8 += tf; \
    tf = OD_DCT_RSHIFT(t8, 1) - tf; \
    t0 += t7; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t7 = t0h - t7; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    t3 += (tc*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    tc -= (t3*15137 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    t3 += (tc*14341 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    t4 -= (tb*14341 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    tb += (t4*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    t4 -= (tb*4161 + 8192) >> 14; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    t8 += (t7*13573 + 8192) >> 14; \
    /* 11585/32768 ~= Sin[Pi/4]/2 ~= 0.353553390593274 */ \
    t7 -= (t8*11585 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    t8 += (t7*13573 + 8192) >> 14; \
    /* TODO: Can we move these into another operation */ \
    t1 = -t1; \
    t5 = -t5; \
    t9 = -t9; \
    tb = -tb; \
    td = -td; \
  } \
  while (0)

#define OD_FDCT_32(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, t6, tm, \
 te, tu, t1, th, t9, tp, t5, tl, td, tt, t3, tj, tb, tr, t7, tn, tf, tv) \
  /* Embedded 32-point orthonormal Type-II fDCT. */ \
  do { \
    int tgh; \
    int thh; \
    int tih; \
    int tkh; \
    int tmh; \
    int tnh; \
    int toh; \
    int tqh; \
    int tsh; \
    int tuh; \
    int tvh; \
    tv = t0 - tv; \
    tvh = OD_DCT_RSHIFT(tv, 1); \
    t0 -= tvh; \
    tu += t1; \
    tuh = OD_DCT_RSHIFT(tu, 1); \
    t1 = tuh - t1; \
    tt = t2 - tt; \
    t2 -= OD_DCT_RSHIFT(tt, 1); \
    ts += t3; \
    tsh = OD_DCT_RSHIFT(ts, 1); \
    t3 = tsh - t3; \
    tr = t4 - tr; \
    t4 -= OD_DCT_RSHIFT(tr, 1); \
    tq += t5; \
    tqh = OD_DCT_RSHIFT(tq, 1); \
    t5 = tqh - t5; \
    tp = t6 - tp; \
    t6 -= OD_DCT_RSHIFT(tp, 1); \
    to += t7; \
    toh = OD_DCT_RSHIFT(to, 1); \
    t7 = toh - t7; \
    tn = t8 - tn; \
    tnh = OD_DCT_RSHIFT(tn, 1); \
    t8 -= tnh; \
    tm += t9; \
    tmh = OD_DCT_RSHIFT(tm, 1); \
    t9 = tmh - t9; \
    tl = ta - tl; \
    ta -= OD_DCT_RSHIFT(tl, 1); \
    tk += tb; \
    tkh = OD_DCT_RSHIFT(tk, 1); \
    tb = tkh - tb; \
    tj = tc - tj; \
    tc -= OD_DCT_RSHIFT(tj, 1); \
    ti += td; \
    tih = OD_DCT_RSHIFT(ti, 1); \
    td = tih - td; \
    th = te - th; \
    thh = OD_DCT_RSHIFT(th, 1); \
    te -= thh; \
    tg += tf; \
    tgh = OD_DCT_RSHIFT(tg, 1); \
    tf = tgh - tf; \
    OD_FDCT_16_ASYM(t0, tg, tgh, t8, to, toh, t4, tk, tkh, tc, ts, tsh, \
     t2, ti, tih, ta, tq, tqh, t6, tm, tmh, te, tu, tuh); \
    OD_FDST_16_ASYM(tv, tvh, tf, tn, tnh, t7, tr, tb, tj, t3, \
     tt, td, tl, t5, tp, t9, th, thh, t1); \
  } \
  while (0)

#define OD_IDCT_32(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, t6, tm, \
 te, tu, t1, th, t9, tp, t5, tl, td, tt, t3, tj, tb, tr, t7, tn, tf, tv) \
  /* Embedded 32-point orthonormal Type-II iDCT. */ \
  do { \
    int t1h; \
    int t3h; \
    int t5h; \
    int t7h; \
    int t9h; \
    int tbh; \
    int tdh; \
    int tfh; \
    int thh; \
    int tth; \
    int tvh; \
    OD_IDST_16_ASYM(tv, tvh, tn, tr, tj, tt, tth, tl, tp, th, thh, \
     tu, tm, tq, ti, ts, tk, to, tg); \
    OD_IDCT_16_ASYM(t0, t8, t4, tc, t2, ta, t6, te, \
     t1, t1h, t9, t9h, t5, t5h, td, tdh, t3, t3h, tb, tbh, t7, t7h, tf, tfh); \
    tu = t1h - tu; \
    t1 -= tu; \
    te += thh; \
    th = te - th; \
    tm = t9h - tm; \
    t9 -= tm; \
    t6 += OD_DCT_RSHIFT(tp, 1); \
    tp = t6 - tp; \
    tq = t5h - tq; \
    t5 -= tq; \
    ta += OD_DCT_RSHIFT(tl, 1); \
    tl = ta - tl; \
    ti = tdh - ti; \
    td -= ti; \
    t2 += tth; \
    tt = t2 - tt; \
    ts = t3h - ts; \
    t3 -= ts; \
    tc += OD_DCT_RSHIFT(tj, 1); \
    tj = tc - tj; \
    tk = tbh - tk; \
    tb -= tk; \
    t4 += OD_DCT_RSHIFT(tr, 1); \
    tr = t4 - tr; \
    to = t7h - to; \
    t7 -= to; \
    t8 += OD_DCT_RSHIFT(tn, 1); \
    tn = t8 - tn; \
    tg = tfh - tg; \
    tf -= tg; \
    t0 += tvh; \
    tv = t0 - tv; \
  } \
  while (0)

#define OD_FDCT_32_ASYM(t0, tg, tgh, t8, to, toh, t4, tk, tkh, tc, ts, tsh, \
 t2, ti, tih, ta, tq, tqh, t6, tm, tmh, te, tu, tuh, t1, th, thh, \
 t9, tp, tph, t5, tl, tlh, td, tt, tth, t3, tj, tjh, tb, tr, trh, \
 t7, tn, tnh, tf, tv, tvh) \
  /* Embedded 32-point asymmetric Type-II fDCT. */ \
  do { \
    t0 += tvh; \
    tv = t0 - tv; \
    t1 = tuh - t1; \
    tu -= t1; \
    t2 += tth; \
    tt = t2 - tt; \
    t3 = tsh - t3; \
    ts -= t3; \
    t4 += trh; \
    tr = t4 - tr; \
    t5 = tqh - t5; \
    tq -= t5; \
    t6 += tph; \
    tp = t6 - tp; \
    t7 = toh - t7; \
    to -= t7; \
    t8 += tnh; \
    tn = t8 - tn; \
    t9 = tmh - t9; \
    tm -= t9; \
    ta += tlh; \
    tl = ta - tl; \
    tb = tkh - tb; \
    tk -= tb; \
    tc += tjh; \
    tj = tc - tj; \
    td = tih - td; \
    ti -= td; \
    te += thh; \
    th = te - th; \
    tf = tgh - tf; \
    tg -= tf; \
    OD_FDCT_16(t0, tg, t8, to, t4, tk, tc, ts, \
     t2, ti, ta, tq, t6, tm, te, tu); \
    OD_FDST_16(tv, tf, tn, t7, tr, tb, tj, t3, \
     tt, td, tl, t5, tp, t9, th, t1); \
  } \
  while (0)

#define OD_IDCT_32_ASYM(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, \
 t6, tm, te, tu, t1, t1h, th, thh, t9, t9h, tp, tph, t5, t5h, tl, tlh, \
 td, tdh, tt, tth, t3, t3h, tj, tjh, tb, tbh, tr, trh, t7, t7h, tn, tnh, \
 tf, tfh, tv, tvh) \
  /* Embedded 32-point asymmetric Type-II iDCT. */ \
  do { \
    OD_IDST_16(tv, tn, tr, tj, tt, tl, tp, th, \
     tu, tm, tq, ti, ts, tk, to, tg); \
    OD_IDCT_16(t0, t8, t4, tc, t2, ta, t6, te, \
     t1, t9, t5, td, t3, tb, t7, tf); \
    tv = t0 - tv; \
    tvh = OD_DCT_RSHIFT(tv, 1); \
    t0 -= tvh; \
    t1 += tu; \
    t1h = OD_DCT_RSHIFT(t1, 1); \
    tu = t1h - tu; \
    tt = t2 - tt; \
    tth = OD_DCT_RSHIFT(tt, 1); \
    t2 -= tth; \
    t3 += ts; \
    t3h = OD_DCT_RSHIFT(t3, 1); \
    ts = t3h - ts; \
    tr = t4 - tr; \
    trh = OD_DCT_RSHIFT(tr, 1); \
    t4 -= trh; \
    t5 += tq; \
    t5h = OD_DCT_RSHIFT(t5, 1); \
    tq = t5h - tq; \
    tp = t6 - tp; \
    tph = OD_DCT_RSHIFT(tp, 1); \
    t6 -= tph; \
    t7 += to; \
    t7h = OD_DCT_RSHIFT(t7, 1); \
    to = t7h - to; \
    tn = t8 - tn; \
    tnh = OD_DCT_RSHIFT(tn, 1); \
    t8 -= tnh; \
    t9 += tm; \
    t9h = OD_DCT_RSHIFT(t9, 1); \
    tm = t9h - tm; \
    tl = ta - tl; \
    tlh = OD_DCT_RSHIFT(tl, 1); \
    ta -= tlh; \
    tb += tk; \
    tbh = OD_DCT_RSHIFT(tb, 1); \
    tk = tbh - tk; \
    tj = tc - tj; \
    tjh = OD_DCT_RSHIFT(tj, 1); \
    tc -= tjh; \
    td += ti; \
    tdh = OD_DCT_RSHIFT(td, 1); \
    ti = tdh - ti; \
    th = te - th; \
    thh = OD_DCT_RSHIFT(th, 1); \
    te -= thh; \
    tf += tg; \
    tfh = OD_DCT_RSHIFT(tf, 1); \
    tg = tfh - tg; \
  } \
  while (0)

#define OD_FDST_32_ASYM(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, t6, \
 tm, te, tu, t1, th, t9, tp, t5, tl, td, tt, t3, tj, tb, tr, t7, tn, tf, tv) \
  /* Embedded 32-point asymmetric Type-IV fDST. */ \
  do { \
    int t0h; \
    int t1h; \
    int t4h; \
    int t5h; \
    int tqh; \
    int trh; \
    int tuh; \
    int tvh; \
    \
    tu = -tu; \
    \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    OD_DCT_OVERFLOW_CHECK(tq, 13573, 8192, 271); \
    t5 -= (tq*13573 + 8192) >> 14; \
    /* 11585/32768 ~= Sin[Pi/4]/2 ~= 0.353553390593274 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 11585, 16384, 272); \
    tq += (t5*11585 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    OD_DCT_OVERFLOW_CHECK(tq, 13573, 8192, 273); \
    t5 -= (tq*13573 + 8192) >> 14; \
    /* 29957/32768 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 29957, 16384, 274); \
    tp += (t6*29957 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    OD_DCT_OVERFLOW_CHECK(tp, 11585, 8192, 275); \
    t6 -= (tp*11585 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 19195, 16384, 276); \
    tp -= (t6*19195 + 16384) >> 15; \
    /* 29957/32768 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 29957, 16384, 277); \
    tu += (t1*29957 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    OD_DCT_OVERFLOW_CHECK(tu, 11585, 8192, 278); \
    t1 -= (tu*11585 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 19195, 16384, 279); \
    tu -= (t1*19195 + 16384) >> 15; \
    /* 28681/32768 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 28681, 16384, 280); \
    tt += (t2*28681 + 16384) >> 15; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(tt, 15137, 8192, 281); \
    t2 -= (tt*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 4161, 8192, 282); \
    tt += (t2*4161 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(ts, 4161, 8192, 283); \
    t3 += (ts*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 15137, 8192, 284); \
    ts -= (t3*15137 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(ts, 14341, 8192, 285); \
    t3 += (ts*14341 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    OD_DCT_OVERFLOW_CHECK(tm, 19195, 16384, 286); \
    t9 -= (tm*19195 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 11585, 8192, 287); \
    tm -= (t9*11585 + 8192) >> 14; \
    /* 7489/8192 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(tm, 7489, 4096, 288); \
    t9 += (tm*7489 + 4096) >> 13; \
    /* 3259/8192 ~= 2*Tan[Pi/16] ~= 0.397824734759316 */ \
    OD_DCT_OVERFLOW_CHECK(tl, 3259, 4096, 289); \
    ta += (tl*3259 + 4096) >> 13; \
    /* 3135/16384 ~= Sin[Pi/8]/2 ~= 0.1913417161825449 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 3135, 8192, 290); \
    tl -= (ta*3135 + 8192) >> 14; \
    /* 3259/8192 ~= 2*Tan[Pi/16] ~= 0.397824734759316 */ \
    OD_DCT_OVERFLOW_CHECK(tl, 3259, 4096, 291); \
    ta += (tl*3259 + 4096) >> 13; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(tk, 4161, 8192, 292); \
    tb += (tk*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 15137, 8192, 293); \
    tk -= (tb*15137 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(tk, 14341, 8192, 294); \
    tb += (tk*14341 + 8192) >> 14; \
    /* 29957/32768 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    OD_DCT_OVERFLOW_CHECK(te, 29957, 16384, 295); \
    th += (te*29957 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    OD_DCT_OVERFLOW_CHECK(th, 11585, 8192, 296); \
    te -= (th*11585 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    OD_DCT_OVERFLOW_CHECK(te, 19195, 16384, 297); \
    th -= (te*19195 + 16384) >> 15; \
    /* 28681/32768 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 28681, 16384, 298); \
    tj += (tc*28681 + 16384) >> 15; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(tj, 15137, 8192, 299); \
    tc -= (tj*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 4161, 8192, 300); \
    tj += (tc*4161 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    OD_DCT_OVERFLOW_CHECK(ti, 4161, 8192, 301); \
    td += (ti*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    OD_DCT_OVERFLOW_CHECK(td, 15137, 8192, 302); \
    ti -= (td*15137 + 8192) >> 14; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    OD_DCT_OVERFLOW_CHECK(ti, 14341, 8192, 303); \
    td += (ti*14341 + 8192) >> 14; \
    \
    t1 = -t1; \
    t2 = -t2; \
    t3 = -t3; \
    td = -td; \
    tg = -tg; \
    to = -to; \
    ts = -ts; \
    \
    tr -= OD_DCT_RSHIFT(t5, 1); \
    t5 += tr; \
    tq -= OD_DCT_RSHIFT(t4, 1); /* pass */ \
    t4 += tq; \
    t6 -= OD_DCT_RSHIFT(t7, 1); \
    t7 += t6; \
    to -= OD_DCT_RSHIFT(tp, 1); /* pass */ \
    tp += to; \
    t1 += OD_DCT_RSHIFT(t0, 1); /* pass */ \
    t0 -= t1; \
    tv -= OD_DCT_RSHIFT(tu, 1); \
    tu += tv; \
    t3 -= OD_DCT_RSHIFT(tt, 1); \
    tt += t3; \
    t2 += OD_DCT_RSHIFT(ts, 1); \
    ts -= t2; \
    t9 -= OD_DCT_RSHIFT(t8, 1); /* pass */ \
    t8 += t9; \
    tn += OD_DCT_RSHIFT(tm, 1); \
    tm -= tn; \
    tb += OD_DCT_RSHIFT(ta, 1); \
    ta -= tb; \
    tl -= OD_DCT_RSHIFT(tk, 1); \
    tk += tl; \
    te -= OD_DCT_RSHIFT(tf, 1); /* pass */ \
    tf += te; \
    tg -= OD_DCT_RSHIFT(th, 1); \
    th += tg; \
    tc -= OD_DCT_RSHIFT(ti, 1); \
    ti += tc; \
    td += OD_DCT_RSHIFT(tj, 1); \
    tj -= td; \
    \
    t4 = -t4; \
    \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.8206787908286602 */ \
    OD_DCT_OVERFLOW_CHECK(tr, 6723, 4096, 304); \
    t4 += (tr*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.9807852804032304 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 16069, 8192, 305); \
    tr -= (t4*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.8206787908286602 */ \
    OD_DCT_OVERFLOW_CHECK(tr, 6723, 4096, 306); \
    t4 += (tr*6723 + 4096) >> 13; \
    /* 17515/32768 ~= Tan[5*Pi/32] ~= 0.5345111359507916 */ \
    OD_DCT_OVERFLOW_CHECK(tq, 17515, 16384, 307); \
    t5 += (tq*17515 + 16384) >> 15; \
    /* 13623/16384 ~= Sin[5*Pi/16] ~= 0.8314696123025452 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 13623, 8192, 308); \
    tq -= (t5*13623 + 8192) >> 14; \
    /* 17515/32768 ~= Tan[5*Pi/32] ~= 0.5345111359507916 */ \
    OD_DCT_OVERFLOW_CHECK(tq, 17515, 16384, 309); \
    t5 += (tq*17515 + 16384) >> 15; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    OD_DCT_OVERFLOW_CHECK(to, 3227, 16384, 310); \
    t7 += (to*3227 + 16384) >> 15; \
    /* 6393/32768 ~= Sin[Pi/16] ~= 0.19509032201612825 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 6393, 16384, 311); \
    to -= (t7*6393 + 16384) >> 15; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    OD_DCT_OVERFLOW_CHECK(to, 3227, 16384, 312); \
    t7 += (to*3227 + 16384) >> 15; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    OD_DCT_OVERFLOW_CHECK(tp, 2485, 4096, 313); \
    t6 += (tp*2485 + 4096) >> 13; \
    /* 18205/32768 ~= Sin[3*Pi/16] ~= 0.555570233019602 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 18205, 16384, 314); \
    tp -= (t6*18205 + 16384) >> 15; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    OD_DCT_OVERFLOW_CHECK(tp, 2485, 4096, 315); \
    t6 += (tp*2485 + 4096) >> 13; \
    \
    t5 = -t5; \
    \
    tr += to; \
    trh = OD_DCT_RSHIFT(tr, 1); \
    to -= trh; \
    t4 += t7; \
    t4h = OD_DCT_RSHIFT(t4, 1); \
    t7 -= t4h; \
    t5 += tp; \
    t5h = OD_DCT_RSHIFT(t5, 1); \
    tp -= t5h; \
    tq += t6; \
    tqh = OD_DCT_RSHIFT(tq, 1); \
    t6 -= tqh; \
    t0 -= t3; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t3 += t0h; \
    tv -= ts; \
    tvh = OD_DCT_RSHIFT(tv, 1); \
    ts += tvh; \
    tu += tt; \
    tuh = OD_DCT_RSHIFT(tu, 1); \
    tt -= tuh; \
    t1 -= t2; \
    t1h = OD_DCT_RSHIFT(t1, 1); \
    t2 += t1h; \
    t8 += tb; \
    tb -= OD_DCT_RSHIFT(t8, 1); \
    tn += tk; \
    tk -= OD_DCT_RSHIFT(tn, 1); \
    t9 += tl; \
    tl -= OD_DCT_RSHIFT(t9, 1); \
    tm -= ta; \
    ta += OD_DCT_RSHIFT(tm, 1); \
    tc -= tf; \
    tf += OD_DCT_RSHIFT(tc, 1); \
    tj += tg; \
    tg -= OD_DCT_RSHIFT(tj, 1); \
    td -= te; \
    te += OD_DCT_RSHIFT(td, 1); \
    ti += th; \
    th -= OD_DCT_RSHIFT(ti, 1); \
    \
    t9 = -t9; \
    tl = -tl; \
    \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    OD_DCT_OVERFLOW_CHECK(tn, 805, 8192, 316); \
    t8 += (tn*805 + 8192) >> 14; \
    /* 803/8192 ~= Sin[Pi/32] ~= 0.0980171403295606 */ \
    OD_DCT_OVERFLOW_CHECK(t8, 803, 4096, 317); \
    tn -= (t8*803 + 4096) >> 13; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    OD_DCT_OVERFLOW_CHECK(tn, 805, 8192, 318); \
    t8 += (tn*805 + 8192) >> 14; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 11725, 16384, 319); \
    tk += (tb*11725 + 16384) >> 15; \
    /* 5197/8192 ~= Sin[7*Pi/32] ~= 0.6343932841636455 */ \
    OD_DCT_OVERFLOW_CHECK(tk, 5197, 4096, 320); \
    tb -= (tk*5197 + 4096) >> 13; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 11725, 16384, 321); \
    tk += (tb*11725 + 16384) >> 15; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.5993769336819237 */ \
    OD_DCT_OVERFLOW_CHECK(tl, 2455, 2048, 322); \
    ta += (tl*2455 + 2048) >> 12; \
    /* 14449/16384 ~= Sin[11*Pi/32] ~= 0.881921264348355 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 14449, 8192, 323); \
    tl -= (ta*14449 + 8192) >> 14; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.5993769336819237 */ \
    OD_DCT_OVERFLOW_CHECK(tl, 2455, 2048, 324); \
    ta += (tl*2455 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    OD_DCT_OVERFLOW_CHECK(tm, 4861, 16384, 325); \
    t9 += (tm*4861 + 16384) >> 15; \
    /* 1189/4096 ~= Sin[3*Pi/32] ~= 0.29028467725446233 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 1189, 2048, 326); \
    tm -= (t9*1189 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    OD_DCT_OVERFLOW_CHECK(tm, 4861, 16384, 327); \
    t9 += (tm*4861 + 16384) >> 15; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    OD_DCT_OVERFLOW_CHECK(tg, 805, 8192, 328); \
    tf += (tg*805 + 8192) >> 14; \
    /* 803/8192 ~= Sin[Pi/32] ~= 0.0980171403295606 */ \
    OD_DCT_OVERFLOW_CHECK(tf, 803, 4096, 329); \
    tg -= (tf*803 + 4096) >> 13; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    OD_DCT_OVERFLOW_CHECK(tg, 805, 8192, 330); \
    tf += (tg*805 + 8192) >> 14; \
    /* 2931/8192 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    OD_DCT_OVERFLOW_CHECK(tj, 2931, 4096, 331); \
    tc += (tj*2931 + 4096) >> 13; \
    /* 5197/8192 ~= Sin[7*Pi/32] ~= 0.6343932841636455 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 5197, 4096, 332); \
    tj -= (tc*5197 + 4096) >> 13; \
    /* 2931/8192 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    OD_DCT_OVERFLOW_CHECK(tj, 2931, 4096, 333); \
    tc += (tj*2931 + 4096) >> 13; \
    /* 513/2048 ~= Tan[5*Pi/64] ~= 0.25048696019130545 */ \
    OD_DCT_OVERFLOW_CHECK(ti, 513, 1024, 334); \
    td += (ti*513 + 1024) >> 11; \
    /* 7723/16384 ~= Sin[5*Pi/32] ~= 0.47139673682599764 */ \
    OD_DCT_OVERFLOW_CHECK(td, 7723, 8192, 335); \
    ti -= (td*7723 + 8192) >> 14; \
    /* 513/2048 ~= Tan[5*Pi/64] ~= 0.25048696019130545 */ \
    OD_DCT_OVERFLOW_CHECK(ti, 513, 1024, 336); \
    td += (ti*513 + 1024) >> 11; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    OD_DCT_OVERFLOW_CHECK(th, 4861, 16384, 337); \
    te += (th*4861 + 16384) >> 15; \
    /* 1189/4096 ~= Sin[3*Pi/32] ~= 0.29028467725446233 */ \
    OD_DCT_OVERFLOW_CHECK(te, 1189, 2048, 338); \
    th -= (te*1189 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    OD_DCT_OVERFLOW_CHECK(th, 4861, 16384, 339); \
    te += (th*4861 + 16384) >> 15; \
    \
    ta = -ta; \
    tb = -tb; \
    \
    tt += t5h; \
    t5 -= tt; \
    t2 -= tqh; \
    tq += t2; \
    tp += t1h; \
    t1 -= tp; \
    t6 -= tuh; \
    tu += t6; \
    t7 += tvh; \
    tv -= t7; \
    to += t0h; \
    t0 -= to; \
    t3 -= t4h; \
    t4 += t3; \
    ts += trh; \
    tr -= ts; \
    tf -= OD_DCT_RSHIFT(tn, 1); \
    tn += tf; \
    tg -= OD_DCT_RSHIFT(t8, 1); \
    t8 += tg; \
    tk += OD_DCT_RSHIFT(tc, 1); \
    tc -= tk; \
    tb += OD_DCT_RSHIFT(tj, 1); \
    tj -= tb; \
    ta += OD_DCT_RSHIFT(ti, 1); \
    ti -= ta; \
    tl += OD_DCT_RSHIFT(td, 1); \
    td -= tl; \
    te -= OD_DCT_RSHIFT(tm, 1); \
    tm += te; \
    th -= OD_DCT_RSHIFT(t9, 1); \
    t9 += th; \
    ta -= t5; \
    t5 += OD_DCT_RSHIFT(ta, 1); \
    tq -= tl; \
    tl += OD_DCT_RSHIFT(tq, 1); \
    t2 -= ti; \
    ti += OD_DCT_RSHIFT(t2, 1); \
    td -= tt; \
    tt += OD_DCT_RSHIFT(td, 1); \
    tm += tp; \
    tp -= OD_DCT_RSHIFT(tm, 1); \
    t6 += t9; \
    t9 -= OD_DCT_RSHIFT(t6, 1); \
    te -= tu; \
    tu += OD_DCT_RSHIFT(te, 1); \
    t1 -= th; \
    th += OD_DCT_RSHIFT(t1, 1); \
    t0 -= tg; \
    tg += OD_DCT_RSHIFT(t0, 1); \
    tf += tv; \
    tv -= OD_DCT_RSHIFT(tf, 1); \
    t8 -= t7; \
    t7 += OD_DCT_RSHIFT(t8, 1); \
    to -= tn; \
    tn += OD_DCT_RSHIFT(to, 1); \
    t4 -= tk; \
    tk += OD_DCT_RSHIFT(t4, 1); \
    tb -= tr; \
    tr += OD_DCT_RSHIFT(tb, 1); \
    t3 -= tj; \
    tj += OD_DCT_RSHIFT(t3, 1); \
    tc -= ts; \
    ts += OD_DCT_RSHIFT(tc, 1); \
    \
    tr = -tr; \
    ts = -ts; \
    tt = -tt; \
    tu = -tu; \
    \
    /* 2847/4096 ~= (1/Sqrt[2] - Cos[63*Pi/128]/2)/Sin[63*Pi/128] ~=
        0.6950455016354713 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 2847, 2048, 340); \
    tv += (t0*2847 + 2048) >> 12; \
    /* 5791/4096 ~= Sqrt[2]*Sin[63*Pi/128] ~= 1.413787627688534 */ \
    OD_DCT_OVERFLOW_CHECK(tv, 5791, 2048, 341); \
    t0 -= (tv*5791 + 2048) >> 12; \
    /* 5593/8192 ~= (1/Sqrt[2] - Cos[63*Pi/128])/Sin[63*Pi/128] ~=
        0.6827711905810085 */ \
    OD_DCT_OVERFLOW_CHECK(t0, 5593, 4096, 342); \
    tv += (t0*5593 + 4096) >> 13; \
    /* 4099/8192 ~= (1/Sqrt[2] - Cos[31*Pi/128]/2)/Sin[31*Pi/128] ~=
        0.5003088539809675 */ \
    OD_DCT_OVERFLOW_CHECK(tf, 4099, 4096, 343); \
    tg -= (tf*4099 + 4096) >> 13; \
    /* 1997/2048 ~= Sqrt[2]*Sin[31*Pi/128] ~= 0.9751575901732918 */ \
    OD_DCT_OVERFLOW_CHECK(tg, 1997, 1024, 344); \
    tf += (tg*1997 + 1024) >> 11; \
    /* -815/32768 ~= (1/Sqrt[2] - Cos[31*Pi/128])/Sin[31*Pi/128] ~=
        -0.02485756913896231 */ \
    OD_DCT_OVERFLOW_CHECK(tf, 815, 16384, 345); \
    tg += (tf*815 + 16384) >> 15; \
    /* 2527/4096 ~= (1/Sqrt[2] - Cos[17*Pi/128]/2)/Sin[17*Pi/128] ~=
        0.6169210657818165 */ \
    OD_DCT_OVERFLOW_CHECK(t8, 2527, 2048, 346); \
    tn -= (t8*2527 + 2048) >> 12; \
    /* 4695/8192 ~= Sqrt[2]*Sin[17*Pi/128] ~= 0.5730977622997507 */ \
    OD_DCT_OVERFLOW_CHECK(tn, 4695, 4096, 347); \
    t8 += (tn*4695 + 4096) >> 13; \
    /* -4187/8192 ~= (1/Sqrt[2] - Cos[17*Pi/128])/Sin[17*Pi/128] ~=
        -0.5110608601827629 */ \
    OD_DCT_OVERFLOW_CHECK(t8, 4187, 4096, 348); \
    tn += (t8*4187 + 4096) >> 13; \
    /* 5477/8192 ~= (1/Sqrt[2] - Cos[15*Pi/128]/2)/Sin[15*Pi/128] ~=
        0.6685570995525147 */ \
    OD_DCT_OVERFLOW_CHECK(to, 5477, 4096, 349); \
    t7 += (to*5477 + 4096) >> 13; \
    /* 4169/8192 ~= Sqrt[2]*Sin[15*Pi/128] ~= 0.5089684416985407 */ \
    OD_DCT_OVERFLOW_CHECK(t7, 4169, 4096, 350); \
    to -= (t7*4169 + 4096) >> 13; \
    /* -2571/4096 ~= (1/Sqrt[2] - Cos[15*Pi/128])/Sin[15*Pi/128] ~=
        -0.6276441593165217 */ \
    OD_DCT_OVERFLOW_CHECK(to, 2571, 2048, 351); \
    t7 -= (to*2571 + 2048) >> 12; \
    /* 5331/8192 ~= (1/Sqrt[2] - Cos[59*Pi/128]/2)/Sin[59*Pi/128] ~=
        0.6507957303604222 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 5331, 4096, 352); \
    tt += (t2*5331 + 4096) >> 13; \
    /* 5749/4096 ~= Sqrt[2]*Sin[59*Pi/128] ~= 1.4035780182072333 */ \
    OD_DCT_OVERFLOW_CHECK(tt, 5749, 2048, 353); \
    t2 -= (tt*5749 + 2048) >> 12; \
    /* 2413/4096 ~= (1/Sqrt[2] - Cos[59*Pi/128])/Sin[59*Pi/128] ~=
        0.5891266122920528 */ \
    OD_DCT_OVERFLOW_CHECK(t2, 2413, 2048, 354); \
    tt += (t2*2413 + 2048) >> 12; \
    /* 4167/8192 ~= (1/Sqrt[2] - Cos[27*Pi/128]/2)/Sin[27*Pi/128] ~=
        0.5086435289805458 */ \
    OD_DCT_OVERFLOW_CHECK(td, 4167, 4096, 355); \
    ti -= (td*4167 + 4096) >> 13; \
    /* 891/1024 ~= Sqrt[2]*Sin[27*Pi/128] ~= 0.8700688593994939 */ \
    OD_DCT_OVERFLOW_CHECK(ti, 891, 512, 356); \
    td += (ti*891 + 512) >> 10; \
    /* -4327/32768 ~= (1/Sqrt[2] - Cos[27*Pi/128])/Sin[27*Pi/128] ~=
        -0.13204726103773165 */ \
    OD_DCT_OVERFLOW_CHECK(td, 4327, 16384, 357); \
    ti += (td*4327 + 16384) >> 15; \
    /* 2261/4096 ~= (1/Sqrt[2] - Cos[21*Pi/128]/2)/Sin[21*Pi/128] ~=
        0.5519664910950994 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 2261, 2048, 358); \
    tl -= (ta*2261 + 2048) >> 12; \
    /* 2855/4096 ~= Sqrt[2]*Sin[21*Pi/128] ~= 0.6970633083205415 */ \
    OD_DCT_OVERFLOW_CHECK(tl, 2855, 2048, 359); \
    ta += (tl*2855 + 2048) >> 12; \
    /* -5417/16384 ~= (1/Sqrt[2] - Cos[21*Pi/128])/Sin[21*Pi/128] ~=
        -0.3306569439519963 */ \
    OD_DCT_OVERFLOW_CHECK(ta, 5417, 8192, 360); \
    tl += (ta*5417 + 8192) >> 14; \
    /* 3459/4096 ~= (1/Sqrt[2] - Cos[11*Pi/128]/2)/Sin[11*Pi/128] ~=
        0.8444243553292501 */ \
    OD_DCT_OVERFLOW_CHECK(tq, 3459, 2048, 361); \
    t5 += (tq*3459 + 2048) >> 12; \
    /* 1545/4096 ~= Sqrt[2]*Sin[11*Pi/128] ~= 0.37718879887892737 */ \
    OD_DCT_OVERFLOW_CHECK(t5, 1545, 2048, 362); \
    tq -= (t5*1545 + 2048) >> 12; \
    /* -1971/2048 ~= (1/Sqrt[2] - Cos[11*Pi/128])/Sin[11*Pi/128] ~=
        -0.9623434853244648 */ \
    OD_DCT_OVERFLOW_CHECK(tq, 1971, 1024, 363); \
    t5 -= (tq*1971 + 1024) >> 11; \
    /* 323/512 ~= (1/Sqrt[2] - Cos[57*Pi/128]/2)/Sin[57*Pi/128] ~=
        0.6309143839894504 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 323, 256, 364); \
    ts += (t3*323 + 256) >> 9; \
    /* 5707/4096 ~= Sqrt[2]*Sin[57*Pi/128] ~= 1.3933930045694292 */ \
    OD_DCT_OVERFLOW_CHECK(ts, 5707, 2048, 365); \
    t3 -= (ts*5707 + 2048) >> 12; \
    /* 2229/4096 ~= (1/Sqrt[2] - Cos[57*Pi/128])/Sin[57*Pi/128] ~=
        0.5441561539205226 */ \
    OD_DCT_OVERFLOW_CHECK(t3, 2229, 2048, 366); \
    ts += (t3*2229 + 2048) >> 12; \
    /* 1061/2048 ~= (1/Sqrt[2] - Cos[25*Pi/128]/2)/Sin[25*Pi/128] ~=
        0.5180794213368158 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 1061, 1024, 367); \
    tj -= (tc*1061 + 1024) >> 11; \
    /* 6671/8192 ~= Sqrt[2]*Sin[25*Pi/128] ~= 0.8143157536286402 */ \
    OD_DCT_OVERFLOW_CHECK(tj, 6671, 4096, 368); \
    tc += (tj*6671 + 4096) >> 13; \
    /* -6287/32768 ~= (1/Sqrt[2] - Cos[25*Pi/128])/Sin[25*Pi/128] ~=
        -0.19186603041023065 */ \
    OD_DCT_OVERFLOW_CHECK(tc, 6287, 16384, 369); \
    tj += (tc*6287 + 16384) >> 15; \
    /* 4359/8192 ~= (1/Sqrt[2] - Cos[23*Pi/128]/2)/Sin[23*Pi/128] ~=
        0.5321145141202145 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 4359, 4096, 370); \
    tk -= (tb*4359 + 4096) >> 13; \
    /* 3099/4096 ~= Sqrt[2]*Sin[23*Pi/128] ~= 0.7566008898816587 */ \
    OD_DCT_OVERFLOW_CHECK(tk, 3099, 2048, 371); \
    tb += (tk*3099 + 2048) >> 12; \
    /* -2109/8192 ~= (1/Sqrt[2] - Cos[23*Pi/128])/Sin[23*Pi/128] ~=
        -0.2574717698598901 */ \
    OD_DCT_OVERFLOW_CHECK(tb, 2109, 4096, 372); \
    tk += (tb*2109 + 4096) >> 13; \
    /* 5017/8192 ~= (1/Sqrt[2] - Cos[55*Pi/128]/2)/Sin[55*Pi/128] ~=
        0.6124370775787037 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 5017, 4096, 373); \
    tr += (t4*5017 + 4096) >> 13; \
    /* 1413/1024 ~= Sqrt[2]*Sin[55*Pi/128] ~= 1.3798511851368045 */ \
    OD_DCT_OVERFLOW_CHECK(tr, 1413, 512, 374); \
    t4 -= (tr*1413 + 512) >> 10; \
    /* 8195/16384 ~= (1/Sqrt[2] - Cos[55*Pi/128])/Sin[55*Pi/128] ~=
        0.5001583229201391 */ \
    OD_DCT_OVERFLOW_CHECK(t4, 8195, 8192, 375); \
    tr += (t4*8195 + 8192) >> 14; \
    /* 2373/4096 ~= (1/Sqrt[2] - Cos[19*Pi/128]/2)/Sin[19*Pi/128] ~=
        0.5793773719823809 */ \
    OD_DCT_OVERFLOW_CHECK(tm, 2373, 2048, 376); \
    t9 += (tm*2373 + 2048) >> 12; \
    /* 5209/8192 ~= Sqrt[2]*Sin[19*Pi/128] ~= 0.6358464401941452 */ \
    OD_DCT_OVERFLOW_CHECK(t9, 5209, 4096, 377); \
    tm -= (t9*5209 + 4096) >> 13; \
    /* -3391/8192 ~= (1/Sqrt[2] - Cos[19*Pi/128])/Sin[19*Pi/128] ~=
        -0.41395202418930155 */ \
    OD_DCT_OVERFLOW_CHECK(tm, 3391, 4096, 378); \
    t9 -= (tm*3391 + 4096) >> 13; \
    /* 1517/2048 ~= (1/Sqrt[2] - Cos[13*Pi/128]/2)/Sin[13*Pi/128] ~=
        0.7406956190518837 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 1517, 1024, 379); \
    tp -= (t6*1517 + 1024) >> 11; \
    /* 1817/4096 ~= Sqrt[2]*Sin[13*Pi/128] ~= 0.4436129715409088 */ \
    OD_DCT_OVERFLOW_CHECK(tp, 1817, 2048, 380); \
    t6 += (tp*1817 + 2048) >> 12; \
    /* -6331/8192 ~= (1/Sqrt[2] - Cos[13*Pi/128])/Sin[13*Pi/128] ~=
        -0.772825983107003 */ \
    OD_DCT_OVERFLOW_CHECK(t6, 6331, 4096, 381); \
    tp += (t6*6331 + 4096) >> 13; \
    /* 515/1024 ~= (1/Sqrt[2] - Cos[29*Pi/128]/2)/Sin[29*Pi/128] ~=
        0.5029332763556925 */ \
    OD_DCT_OVERFLOW_CHECK(te, 515, 512, 382); \
    th -= (te*515 + 512) >> 10; \
    /* 7567/8192 ~= Sqrt[2]*Sin[29*Pi/128] ~= 0.9237258930790229 */ \
    OD_DCT_OVERFLOW_CHECK(th, 7567, 4096, 383); \
    te += (th*7567 + 4096) >> 13; \
    /* -2513/32768 ~= (1/Sqrt[2] - Cos[29*Pi/128])/Sin[29*Pi/128] ~=
        -0.07670567731102484 */ \
    OD_DCT_OVERFLOW_CHECK(te, 2513, 16384, 384); \
    th += (te*2513 + 16384) >> 15; \
    /* 2753/4096 ~= (1/Sqrt[2] - Cos[61*Pi/128]/2)/Sin[61*Pi/128] ~=
        0.6721457072988726 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 2753, 2048, 385); \
    tu += (t1*2753 + 2048) >> 12; \
    /* 5777/4096 ~= Sqrt[2]*Sin[61*Pi/128] ~= 1.4103816894602614 */ \
    OD_DCT_OVERFLOW_CHECK(tu, 5777, 2048, 386); \
    t1 -= (tu*5777 + 2048) >> 12; \
    /* 1301/2048 ~= (1/Sqrt[2] - Cos[61*Pi/128])/Sin[61*Pi/128] ~=
        0.6352634915376478 */ \
    OD_DCT_OVERFLOW_CHECK(t1, 1301, 1024, 387); \
    tu += (t1*1301 + 1024) >> 11; \
  } \
  while (0)

#define OD_IDST_32_ASYM(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, t6, \
 tm, te, tu, t1, th, t9, tp, t5, tl, td, tt, t3, tj, tb, tr, t7, tn, tf, tv) \
  /* Embedded 32-point asymmetric Type-IV iDST. */ \
  do { \
    int t0h; \
    int t4h; \
    int tbh; \
    int tfh; \
    int tgh; \
    int tkh; \
    int trh; \
    int tvh; \
    /* 1301/2048 ~= (1/Sqrt[2] - Cos[61*Pi/128])/Sin[61*Pi/128] ~=
        0.6352634915376478 */ \
    tf -= (tg*1301 + 1024) >> 11; \
    /* 5777/4096 ~= Sqrt[2]*Sin[61*Pi/128] ~= 1.4103816894602614 */ \
    tg += (tf*5777 + 2048) >> 12; \
    /* 2753/4096 ~= (1/Sqrt[2] - Cos[61*Pi/128]/2)/Sin[61*Pi/128] ~=
        0.6721457072988726 */ \
    tf -= (tg*2753 + 2048) >> 12; \
    /* -2513/32768 ~= (1/Sqrt[2] - Cos[29*Pi/128])/Sin[29*Pi/128] ~=
        -0.07670567731102484 */ \
    th -= (te*2513 + 16384) >> 15; \
    /* 7567/8192 ~= Sqrt[2]*Sin[29*Pi/128] ~= 0.9237258930790229 */ \
    te -= (th*7567 + 4096) >> 13; \
    /* 515/1024 ~= (1/Sqrt[2] - Cos[29*Pi/128]/2)/Sin[29*Pi/128] ~=
        0.5029332763556925 */ \
    th += (te*515 + 512) >> 10; \
    /* -6331/8192 ~= (1/Sqrt[2] - Cos[13*Pi/128])/Sin[13*Pi/128] ~=
        -0.772825983107003 */ \
    tj -= (tc*6331 + 4096) >> 13; \
    /* 1817/4096 ~= Sqrt[2]*Sin[13*Pi/128] ~= 0.4436129715409088 */ \
    tc -= (tj*1817 + 2048) >> 12; \
    /* 1517/2048 ~= (1/Sqrt[2] - Cos[13*Pi/128]/2)/Sin[13*Pi/128] ~=
        0.7406956190518837 */ \
    tj += (tc*1517 + 1024) >> 11; \
    /* -3391/8192 ~= (1/Sqrt[2] - Cos[19*Pi/128])/Sin[19*Pi/128] ~=
        -0.41395202418930155 */ \
    ti += (td*3391 + 4096) >> 13; \
    /* 5209/8192 ~= Sqrt[2]*Sin[19*Pi/128] ~= 0.6358464401941452 */ \
    td += (ti*5209 + 4096) >> 13; \
    /* 2373/4096 ~= (1/Sqrt[2] - Cos[19*Pi/128]/2)/Sin[19*Pi/128] ~=
        0.5793773719823809 */ \
    ti -= (td*2373 + 2048) >> 12; \
    /* 8195/16384 ~= (1/Sqrt[2] - Cos[55*Pi/128])/Sin[55*Pi/128] ~=
        0.5001583229201391 */ \
    tr -= (t4*8195 + 8192) >> 14; \
    /* 1413/1024 ~= Sqrt[2]*Sin[55*Pi/128] ~= 1.3798511851368045 */ \
    t4 += (tr*1413 + 512) >> 10; \
    /* 5017/8192 ~= (1/Sqrt[2] - Cos[55*Pi/128]/2)/Sin[55*Pi/128] ~=
        0.6124370775787037 */ \
    tr -= (t4*5017 + 4096) >> 13; \
    /* -2109/8192 ~= (1/Sqrt[2] - Cos[23*Pi/128])/Sin[23*Pi/128] ~=
        -0.2574717698598901 */ \
    t5 -= (tq*2109 + 4096) >> 13; \
    /* 3099/4096 ~= Sqrt[2]*Sin[23*Pi/128] ~= 0.7566008898816587 */ \
    tq -= (t5*3099 + 2048) >> 12; \
    /* 4359/8192 ~= (1/Sqrt[2] - Cos[23*Pi/128]/2)/Sin[23*Pi/128] ~=
        0.5321145141202145 */ \
    t5 += (tq*4359 + 4096) >> 13; \
    /* -6287/32768 ~= (1/Sqrt[2] - Cos[25*Pi/128])/Sin[25*Pi/128] ~=
        -0.19186603041023065 */ \
    tp -= (t6*6287 + 16384) >> 15; \
    /* 6671/8192 ~= Sqrt[2]*Sin[25*Pi/128] ~= 0.8143157536286402 */ \
    t6 -= (tp*6671 + 4096) >> 13; \
    /* 1061/2048 ~= (1/Sqrt[2] - Cos[25*Pi/128]/2)/Sin[25*Pi/128] ~=
        0.5180794213368158 */ \
    tp += (t6*1061 + 1024) >> 11; \
    /* 2229/4096 ~= (1/Sqrt[2] - Cos[57*Pi/128])/Sin[57*Pi/128] ~=
        0.5441561539205226 */ \
    t7 -= (to*2229 + 2048) >> 12; \
    /* 5707/4096 ~= Sqrt[2]*Sin[57*Pi/128] ~= 1.3933930045694292 */ \
    to += (t7*5707 + 2048) >> 12; \
    /* 323/512 ~= (1/Sqrt[2] - Cos[57*Pi/128]/2)/Sin[57*Pi/128] ~=
        0.6309143839894504 */ \
    t7 -= (to*323 + 256) >> 9; \
    /* -1971/2048 ~= (1/Sqrt[2] - Cos[11*Pi/128])/Sin[11*Pi/128] ~=
        -0.9623434853244648 */ \
    tk += (tb*1971 + 1024) >> 11; \
    /* 1545/4096 ~= Sqrt[2]*Sin[11*Pi/128] ~= 0.37718879887892737 */ \
    tb += (tk*1545 + 2048) >> 12; \
    /* 3459/4096 ~= (1/Sqrt[2] - Cos[11*Pi/128]/2)/Sin[11*Pi/128] ~=
        0.8444243553292501 */ \
    tk -= (tb*3459 + 2048) >> 12; \
    /* -5417/16384 ~= (1/Sqrt[2] - Cos[21*Pi/128])/Sin[21*Pi/128] ~=
        -0.3306569439519963 */ \
    tl -= (ta*5417 + 8192) >> 14; \
    /* 2855/4096 ~= Sqrt[2]*Sin[21*Pi/128] ~= 0.6970633083205415 */ \
    ta -= (tl*2855 + 2048) >> 12; \
    /* 2261/4096 ~= (1/Sqrt[2] - Cos[21*Pi/128]/2)/Sin[21*Pi/128] ~=
        0.5519664910950994 */ \
    tl += (ta*2261 + 2048) >> 12; \
    /* -4327/32768 ~= (1/Sqrt[2] - Cos[27*Pi/128])/Sin[27*Pi/128] ~=
        -0.13204726103773165 */ \
    t9 -= (tm*4327 + 16384) >> 15; \
    /* 891/1024 ~= Sqrt[2]*Sin[27*Pi/128] ~= 0.8700688593994939 */ \
    tm -= (t9*891 + 512) >> 10; \
    /* 4167/8192 ~= (1/Sqrt[2] - Cos[27*Pi/128]/2)/Sin[27*Pi/128] ~=
        0.5086435289805458 */ \
    t9 += (tm*4167 + 4096) >> 13; \
    /* 2413/4096 ~= (1/Sqrt[2] - Cos[59*Pi/128])/Sin[59*Pi/128] ~=
        0.5891266122920528 */ \
    tn -= (t8*2413 + 2048) >> 12; \
    /* 5749/4096 ~= Sqrt[2]*Sin[59*Pi/128] ~= 1.4035780182072333 */ \
    t8 += (tn*5749 + 2048) >> 12; \
    /* 5331/8192 ~= (1/Sqrt[2] - Cos[59*Pi/128]/2)/Sin[59*Pi/128] ~=
        0.6507957303604222 */ \
    tn -= (t8*5331 + 4096) >> 13; \
    /* -2571/4096 ~= (1/Sqrt[2] - Cos[15*Pi/128])/Sin[15*Pi/128] ~=
        -0.6276441593165217 */ \
    ts += (t3*2571 + 2048) >> 12; \
    /* 4169/8192 ~= Sqrt[2]*Sin[15*Pi/128] ~= 0.5089684416985407 */ \
    t3 += (ts*4169 + 4096) >> 13; \
    /* 5477/8192 ~= (1/Sqrt[2] - Cos[15*Pi/128]/2)/Sin[15*Pi/128] ~=
        0.6685570995525147 */ \
    ts -= (t3*5477 + 4096) >> 13; \
    /* -4187/8192 ~= (1/Sqrt[2] - Cos[17*Pi/128])/Sin[17*Pi/128] ~=
        -0.5110608601827629 */ \
    tt -= (t2*4187 + 4096) >> 13; \
    /* 4695/8192 ~= Sqrt[2]*Sin[17*Pi/128] ~= 0.5730977622997507 */ \
    t2 -= (tt*4695 + 4096) >> 13; \
    /* 2527/4096 ~= (1/Sqrt[2] - Cos[17*Pi/128]/2)/Sin[17*Pi/128] ~=
        0.6169210657818165 */ \
    tt += (t2*2527 + 2048) >> 12; \
    /* -815/32768 ~= (1/Sqrt[2] - Cos[31*Pi/128])/Sin[31*Pi/128] ~=
        -0.02485756913896231 */ \
    t1 -= (tu*815 + 16384) >> 15; \
    /* 1997/2048 ~= Sqrt[2]*Sin[31*Pi/128] ~= 0.9751575901732918 */ \
    tu -= (t1*1997 + 1024) >> 11; \
    /* 4099/8192 ~= (1/Sqrt[2] - Cos[31*Pi/128]/2)/Sin[31*Pi/128] ~=
        0.5003088539809675 */ \
    t1 += (tu*4099 + 4096) >> 13; \
    /* 5593/8192 ~= (1/Sqrt[2] - Cos[63*Pi/128])/Sin[63*Pi/128] ~=
        0.6827711905810085 */ \
    tv -= (t0*5593 + 4096) >> 13; \
    /* 5791/4096 ~= Sqrt[2]*Sin[63*Pi/128] ~= 1.413787627688534 */ \
    t0 += (tv*5791 + 2048) >> 12; \
    /* 2847/4096 ~= (1/Sqrt[2] - Cos[63*Pi/128]/2)/Sin[63*Pi/128] ~=
        0.6950455016354713 */ \
    tv -= (t0*2847 + 2048) >> 12; \
    \
    t7 = -t7; \
    tf = -tf; \
    tn = -tn; \
    tr = -tr; \
    \
    t7 -= OD_DCT_RSHIFT(t6, 1); \
    t6 += t7; \
    tp -= OD_DCT_RSHIFT(to, 1); \
    to += tp; \
    tr -= OD_DCT_RSHIFT(tq, 1); \
    tq += tr; \
    t5 -= OD_DCT_RSHIFT(t4, 1); \
    t4 += t5; \
    tt -= OD_DCT_RSHIFT(t3, 1); \
    t3 += tt; \
    ts -= OD_DCT_RSHIFT(t2, 1); \
    t2 += ts; \
    tv += OD_DCT_RSHIFT(tu, 1); \
    tu -= tv; \
    t1 -= OD_DCT_RSHIFT(t0, 1); \
    t0 += t1; \
    th -= OD_DCT_RSHIFT(tg, 1); \
    tg += th; \
    tf -= OD_DCT_RSHIFT(te, 1); \
    te += tf; \
    ti += OD_DCT_RSHIFT(tc, 1); \
    tc -= ti; \
    tj += OD_DCT_RSHIFT(td, 1); \
    td -= tj; \
    tn -= OD_DCT_RSHIFT(tm, 1); \
    tm += tn; \
    t9 -= OD_DCT_RSHIFT(t8, 1); \
    t8 += t9; \
    tl -= OD_DCT_RSHIFT(tb, 1); \
    tb += tl; \
    tk -= OD_DCT_RSHIFT(ta, 1); \
    ta += tk; \
    \
    ti -= th; \
    th += OD_DCT_RSHIFT(ti, 1); \
    td -= te; \
    te += OD_DCT_RSHIFT(td, 1); \
    tm += tl; \
    tl -= OD_DCT_RSHIFT(tm, 1); \
    t9 += ta; \
    ta -= OD_DCT_RSHIFT(t9, 1); \
    tp += tq; \
    tq -= OD_DCT_RSHIFT(tp, 1); \
    t6 += t5; \
    t5 -= OD_DCT_RSHIFT(t6, 1); \
    t2 -= t1; \
    t1 += OD_DCT_RSHIFT(t2, 1); \
    tt -= tu; \
    tu += OD_DCT_RSHIFT(tt, 1); \
    tr += t7; \
    trh = OD_DCT_RSHIFT(tr, 1); \
    t7 -= trh; \
    t4 -= to; \
    t4h = OD_DCT_RSHIFT(t4, 1); \
    to += t4h; \
    t0 += t3; \
    t0h = OD_DCT_RSHIFT(t0, 1); \
    t3 -= t0h; \
    tv += ts; \
    tvh = OD_DCT_RSHIFT(tv, 1); \
    ts -= tvh; \
    tf -= tc; \
    tfh = OD_DCT_RSHIFT(tf, 1); \
    tc += tfh; \
    tg += tj; \
    tgh = OD_DCT_RSHIFT(tg, 1); \
    tj -= tgh; \
    tb -= t8; \
    tbh = OD_DCT_RSHIFT(tb, 1); \
    t8 += tbh; \
    tk += tn; \
    tkh = OD_DCT_RSHIFT(tk, 1); \
    tn -= tkh; \
    \
    ta = -ta; \
    tq = -tq; \
    \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    te -= (th*4861 + 16384) >> 15; \
    /* 1189/4096 ~= Sin[3*Pi/32] ~= 0.29028467725446233 */ \
    th += (te*1189 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    te -= (th*4861 + 16384) >> 15; \
    /* 513/2048 ~= Tan[5*Pi/64] ~= 0.25048696019130545 */ \
    tm -= (t9*513 + 1024) >> 11; \
    /* 7723/16384 ~= Sin[5*Pi/32] ~= 0.47139673682599764 */ \
    t9 += (tm*7723 + 8192) >> 14; \
    /* 513/2048 ~= Tan[5*Pi/64] ~= 0.25048696019130545 */ \
    tm -= (t9*513 + 1024) >> 11; \
    /* 2931/8192 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    t6 -= (tp*2931 + 4096) >> 13; \
    /* 5197/8192 ~= Sin[7*Pi/32] ~= 0.6343932841636455 */ \
    tp += (t6*5197 + 4096) >> 13; \
    /* 2931/8192 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    t6 -= (tp*2931 + 4096) >> 13; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    tu -= (t1*805 + 8192) >> 14; \
    /* 803/8192 ~= Sin[Pi/32] ~= 0.0980171403295606 */ \
    t1 += (tu*803 + 4096) >> 13; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    tu -= (t1*805 + 8192) >> 14; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    ti -= (td*4861 + 16384) >> 15; \
    /* 1189/4096 ~= Sin[3*Pi/32] ~= 0.29028467725446233 */ \
    td += (ti*1189 + 2048) >> 12; \
    /* 4861/32768 ~= Tan[3*Pi/64] ~= 0.14833598753834742 */ \
    ti -= (td*4861 + 16384) >> 15; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.5993769336819237 */ \
    ta -= (tl*2455 + 2048) >> 12; \
    /* 14449/16384 ~= Sin[11*Pi/32] ~= 0.881921264348355 */ \
    tl += (ta*14449 + 8192) >> 14; \
    /* 2455/4096 ~= Tan[11*Pi/64] ~= 0.5993769336819237 */ \
    ta -= (tl*2455 + 2048) >> 12; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    t5 -= (tq*11725 + 16384) >> 15; \
    /* 5197/8192 ~= Sin[7*Pi/32] ~= 0.6343932841636455 */ \
    tq += (t5*5197 + 4096) >> 13; \
    /* 11725/32768 ~= Tan[7*Pi/64] ~= 0.3578057213145241 */ \
    t5 -= (tq*11725 + 16384) >> 15; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    t2 -= (tt*805 + 8192) >> 14; \
    /* 803/8192 ~= Sin[Pi/32] ~= 0.0980171403295606 */ \
    tt += (t2*803 + 4096) >> 13; \
    /* 805/16384 ~= Tan[Pi/64] ~= 0.04912684976946793 */ \
    t2 -= (tt*805 + 8192) >> 14; \
    \
    tl = -tl; \
    ti = -ti; \
    \
    th += OD_DCT_RSHIFT(t9, 1); \
    t9 -= th; \
    te -= OD_DCT_RSHIFT(tm, 1); \
    tm += te; \
    t1 += OD_DCT_RSHIFT(tp, 1); \
    tp -= t1; \
    tu -= OD_DCT_RSHIFT(t6, 1); \
    t6 += tu; \
    ta -= OD_DCT_RSHIFT(td, 1); \
    td += ta; \
    tl += OD_DCT_RSHIFT(ti, 1); \
    ti -= tl; \
    t5 += OD_DCT_RSHIFT(tt, 1); \
    tt -= t5; \
    tq += OD_DCT_RSHIFT(t2, 1); \
    t2 -= tq; \
    \
    t8 -= tgh; \
    tg += t8; \
    tn += tfh; \
    tf -= tn; \
    t7 -= tvh; \
    tv += t7; \
    to -= t0h; \
    t0 += to; \
    tc += tbh; \
    tb -= tc; \
    tj += tkh; \
    tk -= tj; \
    ts += t4h; \
    t4 -= ts; \
    t3 += trh; \
    tr -= t3; \
    \
    tk = -tk; \
    \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    tc -= (tj*2485 + 4096) >> 13; \
    /* 18205/32768 ~= Sin[3*Pi/16] ~= 0.555570233019602 */ \
    tj += (tc*18205 + 16384) >> 15; \
    /* 2485/8192 ~= Tan[3*Pi/32] ~= 0.303346683607342 */ \
    tc -= (tj*2485 + 4096) >> 13; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    ts -= (t3*3227 + 16384) >> 15; \
    /* 6393/32768 ~= Sin[Pi/16] ~= 0.19509032201612825 */ \
    t3 += (ts*6393 + 16384) >> 15; \
    /* 3227/32768 ~= Tan[Pi/32] ~= 0.09849140335716425 */ \
    ts -= (t3*3227 + 16384) >> 15; \
    /* 17515/32768 ~= Tan[5*Pi/32] ~= 0.5345111359507916 */ \
    tk -= (tb*17515 + 16384) >> 15; \
    /* 13623/16384 ~= Sin[5*Pi/16] ~= 0.8314696123025452 */ \
    tb += (tk*13623 + 8192) >> 14; \
    /* 17515/32768 ~= Tan[5*Pi/32] ~= 0.5345111359507916 */ \
    tk -= (tb*17515 + 16384) >> 15; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.8206787908286602 */ \
    t4 -= (tr*6723 + 4096) >> 13; \
    /* 16069/16384 ~= Sin[7*Pi/16] ~= 0.9807852804032304 */ \
    tr += (t4*16069 + 8192) >> 14; \
    /* 6723/8192 ~= Tan[7*Pi/32] ~= 0.8206787908286602 */ \
    t4 -= (tr*6723 + 4096) >> 13; \
    \
    t4 = -t4; \
    \
    tp += tm; \
    tm -= OD_DCT_RSHIFT(tp, 1); \
    t9 -= t6; \
    t6 += OD_DCT_RSHIFT(t9, 1); \
    th -= t1; \
    t1 += OD_DCT_RSHIFT(th, 1); \
    tu -= te; \
    te += OD_DCT_RSHIFT(tu, 1); /* pass */ \
    t5 -= tl; \
    tl += OD_DCT_RSHIFT(t5, 1); \
    ta += tq; \
    tq -= OD_DCT_RSHIFT(ta, 1); \
    td += tt; \
    tt -= OD_DCT_RSHIFT(td, 1); \
    t2 -= ti; \
    ti += OD_DCT_RSHIFT(t2, 1); /* pass */ \
    t7 += t8; \
    t8 -= OD_DCT_RSHIFT(t7, 1); \
    tn -= to; \
    to += OD_DCT_RSHIFT(tn, 1); \
    tf -= tv; \
    tv += OD_DCT_RSHIFT(tf, 1); \
    t0 += tg; \
    tg -= OD_DCT_RSHIFT(t0, 1); /* pass */ \
    tj -= t3; \
    t3 += OD_DCT_RSHIFT(tj, 1); /* pass */ \
    ts -= tc; \
    tc += OD_DCT_RSHIFT(ts, 1); \
    t4 -= tb; \
    tb += OD_DCT_RSHIFT(t4, 1); /* pass */ \
    tk -= tr; \
    tr += OD_DCT_RSHIFT(tk, 1); \
    \
    t1 = -t1; \
    t3 = -t3; \
    t7 = -t7; \
    t8 = -t8; \
    tg = -tg; \
    tm = -tm; \
    to = -to; \
    \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    tm -= (t9*14341 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t9 += (tm*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    tm -= (t9*4161 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    tp -= (t6*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t6 += (tp*15137 + 8192) >> 14; \
    /* 28681/32768 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    tp -= (t6*28681 + 16384) >> 15; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    th += (te*19195 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    te += (th*11585 + 8192) >> 14; \
    /* 29957/32768 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    th -= (te*29957 + 16384) >> 15; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    tq -= (t5*14341 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t5 += (tq*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    tq -= (t5*4161 + 8192) >> 14; \
    /* 3259/8192 ~= 2*Tan[Pi/16] ~= 0.397824734759316 */ \
    ta -= (tl*3259 + 4096) >> 13; \
    /* 3135/16384 ~= Sin[Pi/8]/2 ~= 0.1913417161825449 */ \
    tl += (ta*3135 + 8192) >> 14; \
    /* 3259/8192 ~= 2*Tan[Pi/16] ~= 0.397824734759316 */ \
    ta -= (tl*3259 + 4096) >> 13; \
    /* 7489/8192 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    ti -= (td*7489 + 4096) >> 13; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    td += (ti*11585 + 8192) >> 14; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    ti += (td*19195 + 16384) >> 15; \
    /* 14341/16384 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    to -= (t7*14341 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t7 += (to*15137 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    to -= (t7*4161 + 8192) >> 14; \
    /* 4161/16384 ~= Tan[3*Pi/16] - Tan[Pi/8] ~= 0.253965075546204 */ \
    tn -= (t8*4161 + 8192) >> 14; \
    /* 15137/16384 ~= Sin[3*Pi/8] ~= 0.923879532511287 */ \
    t8 += (tn*15137 + 8192) >> 14; \
    /* 28681/32768 ~= Tan[3*Pi/16] + Tan[Pi/8]/2 ~= 0.875285419105846 */ \
    tn -= (t8*28681 + 16384) >> 15; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    tf += (tg*19195 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    tg += (tf*11585 + 8192) >> 14; \
    /* 29957/32768 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    tf -= (tg*29957 + 16384) >> 15; \
    /* -19195/32768 ~= Tan[Pi/8] - Tan[Pi/4] ~= -0.585786437626905 */ \
    tj += (tc*19195 + 16384) >> 15; \
    /* 11585/16384 ~= Sin[Pi/4] ~= 0.707106781186548 */ \
    tc += (tj*11585 + 8192) >> 14; \
    /* 29957/32768 ~= Tan[Pi/8] + Tan[Pi/4]/2 ~= 0.914213562373095 */ \
    tj -= (tc*29957 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    tk += (tb*13573 + 8192) >> 14; \
    /* 11585/32768 ~= Sin[Pi/4]/2 ~= 0.353553390593274 */ \
    tb -= (tk*11585 + 16384) >> 15; \
    /* 13573/16384 ~= 2*Tan[Pi/8] ~= 0.828427124746190 */ \
    tk += (tb*13573 + 8192) >> 14; \
    \
    tf = -tf; \
    \
  } \
  while (0)

#define OD_FDCT_64(u0, uw, ug, uM, u8, uE, uo, uU, u4, uA, uk, uQ, uc, uI, \
 us, uY, u2, uy, ui, uO, ua, uG, uq, uW, u6, uC, um, uS, ue, uK, uu, u_, u1, \
 ux, uh, uN, u9, uF, up, uV, u5, uB, ul, uR, ud, uJ, ut, uZ, u3, uz, uj, uP, \
 ub, uH, ur, uX, u7, uD, un, uT, uf, uL, uv, u) \
  /* Embedded 64-point orthonormal Type-II fDCT. */ \
  do { \
    int uwh; \
    int uxh; \
    int uyh; \
    int uzh; \
    int uAh; \
    int uBh; \
    int uCh; \
    int uDh; \
    int uEh; \
    int uFh; \
    int uGh; \
    int uHh; \
    int uIh; \
    int uJh; \
    int uKh; \
    int uLh; \
    int uMh; \
    int uNh; \
    int uOh; \
    int uPh; \
    int uQh; \
    int uRh; \
    int uSh; \
    int uTh; \
    int uUh; \
    int uVh; \
    int uWh; \
    int uXh; \
    int uYh; \
    int uZh; \
    int u_h; \
    int uh_; \
    u = u0 - u; \
    uh_ = OD_DCT_RSHIFT(u, 1); \
    u0 -= uh_; \
    u_ += u1; \
    u_h = OD_DCT_RSHIFT(u_, 1); \
    u1 = u_h - u1; \
    uZ = u2 - uZ; \
    uZh = OD_DCT_RSHIFT(uZ, 1); \
    u2 -= uZh; \
    uY += u3; \
    uYh = OD_DCT_RSHIFT(uY, 1); \
    u3 = uYh - u3; \
    uX = u4 - uX; \
    uXh = OD_DCT_RSHIFT(uX, 1); \
    u4 -= uXh; \
    uW += u5; \
    uWh = OD_DCT_RSHIFT(uW, 1); \
    u5 = uWh - u5; \
    uV = u6 - uV; \
    uVh = OD_DCT_RSHIFT(uV, 1); \
    u6 -= uVh; \
    uU += u7; \
    uUh = OD_DCT_RSHIFT(uU, 1); \
    u7 = uUh - u7; \
    uT = u8 - uT; \
    uTh = OD_DCT_RSHIFT(uT, 1); \
    u8 -= uTh; \
    uS += u9; \
    uSh = OD_DCT_RSHIFT(uS, 1); \
    u9 = uSh - u9; \
    uR = ua - uR; \
    uRh = OD_DCT_RSHIFT(uR, 1); \
    ua -= uRh; \
    uQ += ub; \
    uQh = OD_DCT_RSHIFT(uQ, 1); \
    ub = uQh - ub; \
    uP = uc - uP; \
    uPh = OD_DCT_RSHIFT(uP, 1); \
    uc -= uPh; \
    uO += ud; \
    uOh = OD_DCT_RSHIFT(uO, 1); \
    ud = uOh - ud; \
    uN = ue - uN; \
    uNh = OD_DCT_RSHIFT(uN, 1); \
    ue -= uNh; \
    uM += uf; \
    uMh = OD_DCT_RSHIFT(uM, 1); \
    uf = uMh - uf; \
    uL = ug - uL; \
    uLh = OD_DCT_RSHIFT(uL, 1); \
    ug -= uLh; \
    uK += uh; \
    uKh = OD_DCT_RSHIFT(uK, 1); \
    uh = uKh - uh; \
    uJ = ui - uJ; \
    uJh = OD_DCT_RSHIFT(uJ, 1); \
    ui -= uJh; \
    uI += uj; \
    uIh = OD_DCT_RSHIFT(uI, 1); \
    uj = uIh - uj; \
    uH = uk - uH; \
    uHh = OD_DCT_RSHIFT(uH, 1); \
    uk -= uHh; \
    uG += ul; \
    uGh = OD_DCT_RSHIFT(uG, 1); \
    ul = uGh - ul; \
    uF = um - uF; \
    uFh = OD_DCT_RSHIFT(uF, 1); \
    um -= uFh; \
    uE += un; \
    uEh = OD_DCT_RSHIFT(uE, 1); \
    un = uEh - un; \
    uD = uo - uD; \
    uDh = OD_DCT_RSHIFT(uD, 1); \
    uo -= uDh; \
    uC += up; \
    uCh = OD_DCT_RSHIFT(uC, 1); \
    up = uCh - up; \
    uB = uq - uB; \
    uBh = OD_DCT_RSHIFT(uB, 1); \
    uq -= uBh; \
    uA += ur; \
    uAh = OD_DCT_RSHIFT(uA, 1); \
    ur = uAh - ur; \
    uz = us - uz; \
    uzh = OD_DCT_RSHIFT(uz, 1); \
    us -= uzh; \
    uy += ut; \
    uyh = OD_DCT_RSHIFT(uy, 1); \
    ut = uyh - ut; \
    ux = uu - ux; \
    uxh = OD_DCT_RSHIFT(ux, 1); \
    uu -= uxh; \
    uw += uv; \
    uwh = OD_DCT_RSHIFT(uw, 1); \
    uv = uwh - uv; \
    OD_FDCT_32_ASYM(u0, uw, uwh, ug, uM, uMh, u8, uE, uEh, uo, uU, uUh, \
     u4, uA, uAh, uk, uQ, uQh, uc, uI, uIh, us, uY, uYh, u2, uy, uyh, \
     ui, uO, uOh, ua, uG, uGh, uq, uW, uWh, u6, uC, uCh, um, uS, uSh, \
     ue, uK, uKh, uu, u_, u_h); \
    OD_FDST_32_ASYM(u, uv, uL, uf, uT, un, uD, u7, uX, ur, uH, ub, uP, uj, \
     uz, u3, uZ, ut, uJ, ud, uR, ul, uB, u5, uV, up, uF, u9, uN, uh, ux, u1); \
  } \
  while (0)

#define OD_IDCT_64(u0, uw, ug, uM, u8, uE, uo, uU, u4, uA, uk, uQ, uc, uI, \
 us, uY, u2, uy, ui, uO, ua, uG, uq, uW, u6, uC, um, uS, ue, uK, uu, u_, u1, \
 ux, uh, uN, u9, uF, up, uV, u5, uB, ul, uR, ud, uJ, ut, uZ, u3, uz, uj, uP, \
 ub, uH, ur, uX, u7, uD, un, uT, uf, uL, uv, u) \
  /* Embedded 64-point orthonormal Type-II fDCT. */ \
  do { \
    int u1h; \
    int u3h; \
    int u5h; \
    int u7h; \
    int u9h; \
    int ubh; \
    int udh; \
    int ufh; \
    int uhh; \
    int ujh; \
    int ulh; \
    int unh; \
    int uph; \
    int urh; \
    int uth; \
    int uvh; \
    int uxh; \
    int uzh; \
    int uBh; \
    int uDh; \
    int uFh; \
    int uHh; \
    int uJh; \
    int uLh; \
    int uNh; \
    int uPh; \
    int uRh; \
    int uTh; \
    int uVh; \
    int uXh; \
    int uZh; \
    int uh_; \
    OD_IDST_32_ASYM(u, uL, uT, uD, uX, uH, uP, uz, uZ, uJ, uR, uB, uV, uF, \
     uN, ux, u_, uK, uS, uC, uW, uG, uO, uy, uY, uI, uQ, uA, uU, uE, uM, uw); \
    OD_IDCT_32_ASYM(u0, ug, u8, uo, u4, uk, uc, us, u2, ui, ua, uq, u6, um, \
     ue, uu, u1, u1h, uh, uhh, u9, u9h, up, uph, u5, u5h, ul, ulh, ud, udh, \
     ut, uth, u3, u3h, uj, ujh, ub, ubh, ur, urh, u7, u7h, un, unh, uf, ufh, \
     uv, uvh); \
    uh_ = OD_DCT_RSHIFT(u, 1); \
    u0 += uh_; \
    u = u0 - u; \
    u_ = u1h - u_; \
    u1 -= u_; \
    uZh = OD_DCT_RSHIFT(uZ, 1); \
    u2 += uZh; \
    uZ = u2 - uZ; \
    uY = u3h - uY; \
    u3 -= uY; \
    uXh = OD_DCT_RSHIFT(uX, 1); \
    u4 += uXh; \
    uX = u4 - uX; \
    uW = u5h - uW; \
    u5 -= uW; \
    uVh = OD_DCT_RSHIFT(uV, 1); \
    u6 += uVh; \
    uV = u6 - uV; \
    uU = u7h - uU; \
    u7 -= uU; \
    uTh = OD_DCT_RSHIFT(uT, 1); \
    u8 += uTh; \
    uT = u8 - uT; \
    uS = u9h - uS; \
    u9 -= uS; \
    uRh = OD_DCT_RSHIFT(uR, 1); \
    ua += uRh; \
    uR = ua - uR; \
    uQ = ubh - uQ; \
    ub -= uQ; \
    uPh = OD_DCT_RSHIFT(uP, 1); \
    uc += uPh; \
    uP = uc - uP; \
    uO = udh - uO; \
    ud -= uO; \
    uNh = OD_DCT_RSHIFT(uN, 1); \
    ue += uNh; \
    uN = ue - uN; \
    uM = ufh - uM; \
    uf -= uM; \
    uLh = OD_DCT_RSHIFT(uL, 1); \
    ug += uLh; \
    uL = ug - uL; \
    uK = uhh - uK; \
    uh -= uK; \
    uJh = OD_DCT_RSHIFT(uJ, 1); \
    ui += uJh; \
    uJ = ui - uJ; \
    uI = ujh - uI; \
    uj -= uI; \
    uHh = OD_DCT_RSHIFT(uH, 1); \
    uk += uHh; \
    uH = uk - uH; \
    uG = ulh - uG; \
    ul -= uG; \
    uFh = OD_DCT_RSHIFT(uF, 1); \
    um += uFh; \
    uF = um - uF; \
    uE = unh - uE; \
    un -= uE; \
    uDh = OD_DCT_RSHIFT(uD, 1); \
    uo += uDh; \
    uD = uo - uD; \
    uC = uph - uC; \
    up -= uC; \
    uBh = OD_DCT_RSHIFT(uB, 1); \
    uq += uBh; \
    uB = uq - uB; \
    uA = urh - uA; \
    ur -= uA; \
    uzh = OD_DCT_RSHIFT(uz, 1); \
    us += uzh; \
    uz = us - uz; \
    uy = uth - uy; \
    ut -= uy; \
    uxh = OD_DCT_RSHIFT(ux, 1); \
    uu += uxh; \
    ux = uu - ux; \
    uw = uvh - uw; \
    uv -= uw; \
  } while (0)

#if OD_USE_EMBEDDED_DCTS
void od_bin_fdct4(od_coeff y[4], const od_coeff *x, int xstride) {
  int t0;
  int t1;
  int t2;
  int t3;
  t0 = x[0*xstride];
  t2 = x[1*xstride];
  t1 = x[2*xstride];
  t3 = x[3*xstride];
  OD_FDCT_4(t0, t2, t1, t3);
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
}

void od_bin_idct4(od_coeff *x, int xstride, const od_coeff y[4]) {
  int t0;
  int t1;
  int t2;
  int t3;
  t0 = y[0];
  t2 = y[1];
  t1 = y[2];
  t3 = y[3];
  OD_IDCT_4(t0, t2, t1, t3);
  x[0*xstride] = t0;
  x[1*xstride] = t1;
  x[2*xstride] = t2;
  x[3*xstride] = t3;
}

void od_bin_fdct8(od_coeff y[8], const od_coeff *x, int xstride) {
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  t0 = x[0*xstride];
  t4 = x[1*xstride];
  t2 = x[2*xstride];
  t6 = x[3*xstride];
  t1 = x[4*xstride];
  t5 = x[5*xstride];
  t3 = x[6*xstride];
  t7 = x[7*xstride];
  OD_FDCT_8(t0, t4, t2, t6, t1, t5, t3, t7);
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
  y[4] = (od_coeff)t4;
  y[5] = (od_coeff)t5;
  y[6] = (od_coeff)t6;
  y[7] = (od_coeff)t7;
}

void od_bin_idct8(od_coeff *x, int xstride, const od_coeff y[8]) {
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  t0 = y[0];
  t4 = y[1];
  t2 = y[2];
  t6 = y[3];
  t1 = y[4];
  t5 = y[5];
  t3 = y[6];
  t7 = y[7];
  OD_IDCT_8(t0, t4, t2, t6, t1, t5, t3, t7);
  x[0*xstride] = (od_coeff)t0;
  x[1*xstride] = (od_coeff)t1;
  x[2*xstride] = (od_coeff)t2;
  x[3*xstride] = (od_coeff)t3;
  x[4*xstride] = (od_coeff)t4;
  x[5*xstride] = (od_coeff)t5;
  x[6*xstride] = (od_coeff)t6;
  x[7*xstride] = (od_coeff)t7;
}

void od_bin_fdct16(od_coeff y[16], const od_coeff *x, int xstride) {
  int s0;
  int s1;
  int s2;
  int s3;
  int s4;
  int s5;
  int s6;
  int s7;
  int s8;
  int s9;
  int sa;
  int sb;
  int sc;
  int sd;
  int se;
  int sf;
  s0 = x[0*xstride];
  s8 = x[1*xstride];
  s4 = x[2*xstride];
  sc = x[3*xstride];
  s2 = x[4*xstride];
  sa = x[5*xstride];
  s6 = x[6*xstride];
  se = x[7*xstride];
  s1 = x[8*xstride];
  s9 = x[9*xstride];
  s5 = x[10*xstride];
  sd = x[11*xstride];
  s3 = x[12*xstride];
  sb = x[13*xstride];
  s7 = x[14*xstride];
  sf = x[15*xstride];
  OD_FDCT_16(s0, s8, s4, sc, s2, sa, s6, se, s1, s9, s5, sd, s3, sb, s7, sf);
  y[0] = (od_coeff)s0;
  y[1] = (od_coeff)s1;
  y[2] = (od_coeff)s2;
  y[3] = (od_coeff)s3;
  y[4] = (od_coeff)s4;
  y[5] = (od_coeff)s5;
  y[6] = (od_coeff)s6;
  y[7] = (od_coeff)s7;
  y[8] = (od_coeff)s8;
  y[9] = (od_coeff)s9;
  y[10] = (od_coeff)sa;
  y[11] = (od_coeff)sb;
  y[12] = (od_coeff)sc;
  y[13] = (od_coeff)sd;
  y[14] = (od_coeff)se;
  y[15] = (od_coeff)sf;
}

void od_bin_idct16(od_coeff *x, int xstride, const od_coeff y[16]) {
  int s0;
  int s1;
  int s2;
  int s3;
  int s4;
  int s5;
  int s6;
  int s7;
  int s8;
  int s9;
  int sa;
  int sb;
  int sc;
  int sd;
  int se;
  int sf;
  s0 = y[0];
  s8 = y[1];
  s4 = y[2];
  sc = y[3];
  s2 = y[4];
  sa = y[5];
  s6 = y[6];
  se = y[7];
  s1 = y[8];
  s9 = y[9];
  s5 = y[10];
  sd = y[11];
  s3 = y[12];
  sb = y[13];
  s7 = y[14];
  sf = y[15];
  OD_IDCT_16(s0, s8, s4, sc, s2, sa, s6, se, s1, s9, s5, sd, s3, sb, s7, sf);
  x[0*xstride] = (od_coeff)s0;
  x[1*xstride] = (od_coeff)s1;
  x[2*xstride] = (od_coeff)s2;
  x[3*xstride] = (od_coeff)s3;
  x[4*xstride] = (od_coeff)s4;
  x[5*xstride] = (od_coeff)s5;
  x[6*xstride] = (od_coeff)s6;
  x[7*xstride] = (od_coeff)s7;
  x[8*xstride] = (od_coeff)s8;
  x[9*xstride] = (od_coeff)s9;
  x[10*xstride] = (od_coeff)sa;
  x[11*xstride] = (od_coeff)sb;
  x[12*xstride] = (od_coeff)sc;
  x[13*xstride] = (od_coeff)sd;
  x[14*xstride] = (od_coeff)se;
  x[15*xstride] = (od_coeff)sf;
}
#endif

void od_bin_fdct32(od_coeff y[32], const od_coeff *x, int xstride) {
  /*215 adds, 38 shifts, 87 "muls".*/
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  int t8;
  int t9;
  int ta;
  int tb;
  int tc;
  int td;
  int te;
  int tf;
  int tg;
  int th;
  int ti;
  int tj;
  int tk;
  int tl;
  int tm;
  int tn;
  int to;
  int tp;
  int tq;
  int tr;
  int ts;
  int tt;
  int tu;
  int tv;
  t0 = x[0*xstride];
  tg = x[1*xstride];
  t8 = x[2*xstride];
  to = x[3*xstride];
  t4 = x[4*xstride];
  tk = x[5*xstride];
  tc = x[6*xstride];
  ts = x[7*xstride];
  t2 = x[8*xstride];
  ti = x[9*xstride];
  ta = x[10*xstride];
  tq = x[11*xstride];
  t6 = x[12*xstride];
  tm = x[13*xstride];
  te = x[14*xstride];
  tu = x[15*xstride];
  t1 = x[16*xstride];
  th = x[17*xstride];
  t9 = x[18*xstride];
  tp = x[19*xstride];
  t5 = x[20*xstride];
  tl = x[21*xstride];
  td = x[22*xstride];
  tt = x[23*xstride];
  t3 = x[24*xstride];
  tj = x[25*xstride];
  tb = x[26*xstride];
  tr = x[27*xstride];
  t7 = x[28*xstride];
  tn = x[29*xstride];
  tf = x[30*xstride];
  tv = x[31*xstride];
  OD_FDCT_32(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, t6, tm, te, tu,
   t1, th, t9, tp, t5, tl, td, tt, t3, tj, tb, tr, t7, tn, tf, tv);
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
  y[4] = (od_coeff)t4;
  y[5] = (od_coeff)t5;
  y[6] = (od_coeff)t6;
  y[7] = (od_coeff)t7;
  y[8] = (od_coeff)t8;
  y[9] = (od_coeff)t9;
  y[10] = (od_coeff)ta;
  y[11] = (od_coeff)tb;
  y[12] = (od_coeff)tc;
  y[13] = (od_coeff)td;
  y[14] = (od_coeff)te;
  y[15] = (od_coeff)tf;
  y[16] = (od_coeff)tg;
  y[17] = (od_coeff)th;
  y[18] = (od_coeff)ti;
  y[19] = (od_coeff)tj;
  y[20] = (od_coeff)tk;
  y[21] = (od_coeff)tl;
  y[22] = (od_coeff)tm;
  y[23] = (od_coeff)tn;
  y[24] = (od_coeff)to;
  y[25] = (od_coeff)tp;
  y[26] = (od_coeff)tq;
  y[27] = (od_coeff)tr;
  y[28] = (od_coeff)ts;
  y[29] = (od_coeff)tt;
  y[30] = (od_coeff)tu;
  y[31] = (od_coeff)tv;
}

void od_bin_idct32(od_coeff *x, int xstride, const od_coeff y[32]) {
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  int t8;
  int t9;
  int ta;
  int tb;
  int tc;
  int td;
  int te;
  int tf;
  int tg;
  int th;
  int ti;
  int tj;
  int tk;
  int tl;
  int tm;
  int tn;
  int to;
  int tp;
  int tq;
  int tr;
  int ts;
  int tt;
  int tu;
  int tv;
  t0 = y[0];
  tg = y[1];
  t8 = y[2];
  to = y[3];
  t4 = y[4];
  tk = y[5];
  tc = y[6];
  ts = y[7];
  t2 = y[8];
  ti = y[9];
  ta = y[10];
  tq = y[11];
  t6 = y[12];
  tm = y[13];
  te = y[14];
  tu = y[15];
  t1 = y[16];
  th = y[17];
  t9 = y[18];
  tp = y[19];
  t5 = y[20];
  tl = y[21];
  td = y[22];
  tt = y[23];
  t3 = y[24];
  tj = y[25];
  tb = y[26];
  tr = y[27];
  t7 = y[28];
  tn = y[29];
  tf = y[30];
  tv = y[31];
  OD_IDCT_32(t0, tg, t8, to, t4, tk, tc, ts, t2, ti, ta, tq, t6, tm, te, tu,
   t1, th, t9, tp, t5, tl, td, tt, t3, tj, tb, tr, t7, tn, tf, tv);
  x[0*xstride] = (od_coeff)t0;
  x[1*xstride] = (od_coeff)t1;
  x[2*xstride] = (od_coeff)t2;
  x[3*xstride] = (od_coeff)t3;
  x[4*xstride] = (od_coeff)t4;
  x[5*xstride] = (od_coeff)t5;
  x[6*xstride] = (od_coeff)t6;
  x[7*xstride] = (od_coeff)t7;
  x[8*xstride] = (od_coeff)t8;
  x[9*xstride] = (od_coeff)t9;
  x[10*xstride] = (od_coeff)ta;
  x[11*xstride] = (od_coeff)tb;
  x[12*xstride] = (od_coeff)tc;
  x[13*xstride] = (od_coeff)td;
  x[14*xstride] = (od_coeff)te;
  x[15*xstride] = (od_coeff)tf;
  x[16*xstride] = (od_coeff)tg;
  x[17*xstride] = (od_coeff)th;
  x[18*xstride] = (od_coeff)ti;
  x[19*xstride] = (od_coeff)tj;
  x[20*xstride] = (od_coeff)tk;
  x[21*xstride] = (od_coeff)tl;
  x[22*xstride] = (od_coeff)tm;
  x[23*xstride] = (od_coeff)tn;
  x[24*xstride] = (od_coeff)to;
  x[25*xstride] = (od_coeff)tp;
  x[26*xstride] = (od_coeff)tq;
  x[27*xstride] = (od_coeff)tr;
  x[28*xstride] = (od_coeff)ts;
  x[29*xstride] = (od_coeff)tt;
  x[30*xstride] = (od_coeff)tu;
  x[31*xstride] = (od_coeff)tv;
}

void od_bin_fdct64(od_coeff y[64], const od_coeff *x, int xstride) {
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  int t8;
  int t9;
  int ta;
  int tb;
  int tc;
  int td;
  int te;
  int tf;
  int tg;
  int th;
  int ti;
  int tj;
  int tk;
  int tl;
  int tm;
  int tn;
  int to;
  int tp;
  int tq;
  int tr;
  int ts;
  int tt;
  int tu;
  int tv;
  int tw;
  int tx;
  int ty;
  int tz;
  int tA;
  int tB;
  int tC;
  int tD;
  int tE;
  int tF;
  int tG;
  int tH;
  int tI;
  int tJ;
  int tK;
  int tL;
  int tM;
  int tN;
  int tO;
  int tP;
  int tQ;
  int tR;
  int tS;
  int tT;
  int tU;
  int tV;
  int tW;
  int tX;
  int tY;
  int tZ;
  int t_;
  int t;
  t0 = x[0*xstride];
  tw = x[1*xstride];
  tg = x[2*xstride];
  tM = x[3*xstride];
  t8 = x[4*xstride];
  tE = x[5*xstride];
  to = x[6*xstride];
  tU = x[7*xstride];
  t4 = x[8*xstride];
  tA = x[9*xstride];
  tk = x[10*xstride];
  tQ = x[11*xstride];
  tc = x[12*xstride];
  tI = x[13*xstride];
  ts = x[14*xstride];
  tY = x[15*xstride];
  t2 = x[16*xstride];
  ty = x[17*xstride];
  ti = x[18*xstride];
  tO = x[19*xstride];
  ta = x[20*xstride];
  tG = x[21*xstride];
  tq = x[22*xstride];
  tW = x[23*xstride];
  t6 = x[24*xstride];
  tC = x[25*xstride];
  tm = x[26*xstride];
  tS = x[27*xstride];
  te = x[28*xstride];
  tK = x[29*xstride];
  tu = x[30*xstride];
  t_ = x[31*xstride];
  t1 = x[32*xstride];
  tx = x[33*xstride];
  th = x[34*xstride];
  tN = x[35*xstride];
  t9 = x[36*xstride];
  tF = x[37*xstride];
  tp = x[38*xstride];
  tV = x[39*xstride];
  t5 = x[40*xstride];
  tB = x[41*xstride];
  tl = x[42*xstride];
  tR = x[43*xstride];
  td = x[44*xstride];
  tJ = x[45*xstride];
  tt = x[46*xstride];
  tZ = x[47*xstride];
  t3 = x[48*xstride];
  tz = x[49*xstride];
  tj = x[50*xstride];
  tP = x[51*xstride];
  tb = x[52*xstride];
  tH = x[53*xstride];
  tr = x[54*xstride];
  tX = x[55*xstride];
  t7 = x[56*xstride];
  tD = x[57*xstride];
  tn = x[58*xstride];
  tT = x[59*xstride];
  tf = x[60*xstride];
  tL = x[61*xstride];
  tv = x[62*xstride];
  t = x[63*xstride];
  OD_FDCT_64(t0, tw, tg, tM, t8, tE, to, tU, t4, tA, tk, tQ, tc, tI, ts, tY,
   t2, ty, ti, tO, ta, tG, tq, tW, t6, tC, tm, tS, te, tK, tu, t_, t1, tx,
   th, tN, t9, tF, tp, tV, t5, tB, tl, tR, td, tJ, tt, tZ, t3, tz, tj, tP,
   tb, tH, tr, tX, t7, tD, tn, tT, tf, tL, tv, t);
  y[0] = (od_coeff)t0;
  y[1] = (od_coeff)t1;
  y[2] = (od_coeff)t2;
  y[3] = (od_coeff)t3;
  y[4] = (od_coeff)t4;
  y[5] = (od_coeff)t5;
  y[6] = (od_coeff)t6;
  y[7] = (od_coeff)t7;
  y[8] = (od_coeff)t8;
  y[9] = (od_coeff)t9;
  y[10] = (od_coeff)ta;
  y[11] = (od_coeff)tb;
  y[12] = (od_coeff)tc;
  y[13] = (od_coeff)td;
  y[14] = (od_coeff)te;
  y[15] = (od_coeff)tf;
  y[16] = (od_coeff)tg;
  y[17] = (od_coeff)th;
  y[18] = (od_coeff)ti;
  y[19] = (od_coeff)tj;
  y[20] = (od_coeff)tk;
  y[21] = (od_coeff)tl;
  y[22] = (od_coeff)tm;
  y[23] = (od_coeff)tn;
  y[24] = (od_coeff)to;
  y[25] = (od_coeff)tp;
  y[26] = (od_coeff)tq;
  y[27] = (od_coeff)tr;
  y[28] = (od_coeff)ts;
  y[29] = (od_coeff)tt;
  y[30] = (od_coeff)tu;
  y[31] = (od_coeff)tv;
  y[32] = (od_coeff)tw;
  y[33] = (od_coeff)tx;
  y[34] = (od_coeff)ty;
  y[35] = (od_coeff)tz;
  y[36] = (od_coeff)tA;
  y[37] = (od_coeff)tB;
  y[38] = (od_coeff)tC;
  y[39] = (od_coeff)tD;
  y[40] = (od_coeff)tE;
  y[41] = (od_coeff)tF;
  y[41] = (od_coeff)tF;
  y[42] = (od_coeff)tG;
  y[43] = (od_coeff)tH;
  y[44] = (od_coeff)tI;
  y[45] = (od_coeff)tJ;
  y[46] = (od_coeff)tK;
  y[47] = (od_coeff)tL;
  y[48] = (od_coeff)tM;
  y[49] = (od_coeff)tN;
  y[50] = (od_coeff)tO;
  y[51] = (od_coeff)tP;
  y[52] = (od_coeff)tQ;
  y[53] = (od_coeff)tR;
  y[54] = (od_coeff)tS;
  y[55] = (od_coeff)tT;
  y[56] = (od_coeff)tU;
  y[57] = (od_coeff)tV;
  y[58] = (od_coeff)tW;
  y[59] = (od_coeff)tX;
  y[60] = (od_coeff)tY;
  y[61] = (od_coeff)tZ;
  y[62] = (od_coeff)t_;
  y[63] = (od_coeff)t;
}

void od_bin_idct64(od_coeff *x, int xstride, const od_coeff y[64]) {
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t5;
  int t6;
  int t7;
  int t8;
  int t9;
  int ta;
  int tb;
  int tc;
  int td;
  int te;
  int tf;
  int tg;
  int th;
  int ti;
  int tj;
  int tk;
  int tl;
  int tm;
  int tn;
  int to;
  int tp;
  int tq;
  int tr;
  int ts;
  int tt;
  int tu;
  int tv;
  int tw;
  int tx;
  int ty;
  int tz;
  int tA;
  int tB;
  int tC;
  int tD;
  int tE;
  int tF;
  int tG;
  int tH;
  int tI;
  int tJ;
  int tK;
  int tL;
  int tM;
  int tN;
  int tO;
  int tP;
  int tQ;
  int tR;
  int tS;
  int tT;
  int tU;
  int tV;
  int tW;
  int tX;
  int tY;
  int tZ;
  int t_;
  int t;
  t0 = y[0];
  tw = y[1];
  tg = y[2];
  tM = y[3];
  t8 = y[4];
  tE = y[5];
  to = y[6];
  tU = y[7];
  t4 = y[8];
  tA = y[9];
  tk = y[10];
  tQ = y[11];
  tc = y[12];
  tI = y[13];
  ts = y[14];
  tY = y[15];
  t2 = y[16];
  ty = y[17];
  ti = y[18];
  tO = y[19];
  ta = y[20];
  tG = y[21];
  tq = y[22];
  tW = y[23];
  t6 = y[24];
  tC = y[25];
  tm = y[26];
  tS = y[27];
  te = y[28];
  tK = y[29];
  tu = y[30];
  t_ = y[31];
  t1 = y[32];
  tx = y[33];
  th = y[34];
  tN = y[35];
  t9 = y[36];
  tF = y[37];
  tp = y[38];
  tV = y[39];
  t5 = y[40];
  tB = y[41];
  tl = y[42];
  tR = y[43];
  td = y[44];
  tJ = y[45];
  tt = y[46];
  tZ = y[47];
  t3 = y[48];
  tz = y[49];
  tj = y[50];
  tP = y[51];
  tb = y[52];
  tH = y[53];
  tr = y[54];
  tX = y[55];
  t7 = y[56];
  tD = y[57];
  tn = y[58];
  tT = y[59];
  tf = y[60];
  tL = y[61];
  tv = y[62];
  t = y[63];
  OD_IDCT_64(t0, tw, tg, tM, t8, tE, to, tU, t4, tA, tk, tQ, tc, tI, ts, tY,
   t2, ty, ti, tO, ta, tG, tq, tW, t6, tC, tm, tS, te, tK, tu, t_, t1, tx,
   th, tN, t9, tF, tp, tV, t5, tB, tl, tR, td, tJ, tt, tZ, t3, tz, tj, tP,
   tb, tH, tr, tX, t7, tD, tn, tT, tf, tL, tv, t);
  x[0*xstride] = (od_coeff)t0;
  x[1*xstride] = (od_coeff)t1;
  x[2*xstride] = (od_coeff)t2;
  x[3*xstride] = (od_coeff)t3;
  x[4*xstride] = (od_coeff)t4;
  x[5*xstride] = (od_coeff)t5;
  x[6*xstride] = (od_coeff)t6;
  x[7*xstride] = (od_coeff)t7;
  x[8*xstride] = (od_coeff)t8;
  x[9*xstride] = (od_coeff)t9;
  x[10*xstride] = (od_coeff)ta;
  x[11*xstride] = (od_coeff)tb;
  x[12*xstride] = (od_coeff)tc;
  x[13*xstride] = (od_coeff)td;
  x[14*xstride] = (od_coeff)te;
  x[15*xstride] = (od_coeff)tf;
  x[16*xstride] = (od_coeff)tg;
  x[17*xstride] = (od_coeff)th;
  x[18*xstride] = (od_coeff)ti;
  x[19*xstride] = (od_coeff)tj;
  x[20*xstride] = (od_coeff)tk;
  x[21*xstride] = (od_coeff)tl;
  x[22*xstride] = (od_coeff)tm;
  x[23*xstride] = (od_coeff)tn;
  x[24*xstride] = (od_coeff)to;
  x[25*xstride] = (od_coeff)tp;
  x[26*xstride] = (od_coeff)tq;
  x[27*xstride] = (od_coeff)tr;
  x[28*xstride] = (od_coeff)ts;
  x[29*xstride] = (od_coeff)tt;
  x[30*xstride] = (od_coeff)tu;
  x[31*xstride] = (od_coeff)tv;
  x[32*xstride] = (od_coeff)tw;
  x[33*xstride] = (od_coeff)tx;
  x[34*xstride] = (od_coeff)ty;
  x[35*xstride] = (od_coeff)tz;
  x[36*xstride] = (od_coeff)tA;
  x[37*xstride] = (od_coeff)tB;
  x[38*xstride] = (od_coeff)tC;
  x[39*xstride] = (od_coeff)tD;
  x[40*xstride] = (od_coeff)tE;
  x[41*xstride] = (od_coeff)tF;
  x[41*xstride] = (od_coeff)tF;
  x[42*xstride] = (od_coeff)tG;
  x[43*xstride] = (od_coeff)tH;
  x[44*xstride] = (od_coeff)tI;
  x[45*xstride] = (od_coeff)tJ;
  x[46*xstride] = (od_coeff)tK;
  x[47*xstride] = (od_coeff)tL;
  x[48*xstride] = (od_coeff)tM;
  x[49*xstride] = (od_coeff)tN;
  x[50*xstride] = (od_coeff)tO;
  x[51*xstride] = (od_coeff)tP;
  x[52*xstride] = (od_coeff)tQ;
  x[53*xstride] = (od_coeff)tR;
  x[54*xstride] = (od_coeff)tS;
  x[55*xstride] = (od_coeff)tT;
  x[56*xstride] = (od_coeff)tU;
  x[57*xstride] = (od_coeff)tV;
  x[58*xstride] = (od_coeff)tW;
  x[59*xstride] = (od_coeff)tX;
  x[60*xstride] = (od_coeff)tY;
  x[61*xstride] = (od_coeff)tZ;
  x[62*xstride] = (od_coeff)t_;
  x[63*xstride] = (od_coeff)t;
}

void od_haar(od_coeff *y, int ystride,
  const od_coeff *x, int xstride, int ln) {
  int i;
  int j;
  int level;
  int tstride;
  int n;
  od_coeff tmp[OD_BSIZE_MAX*OD_BSIZE_MAX];
  n = 1 << ln;
  tstride = n;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      tmp[i*tstride + j] = x[i*xstride + j];
    }
  }
  for (level = 0; level < ln; level++) {
    int npairs;
    npairs = n >> level >> 1;
    for (i = 0; i < npairs; i++) {
      for (j = 0; j < npairs; j++) {
        od_coeff a;
        od_coeff b;
        od_coeff c;
        od_coeff d;
        a = tmp[2*i*tstride + 2*j];
        b = tmp[(2*i + 1)*tstride + 2*j];
        c = tmp[2*i*tstride + 2*j + 1];
        d = tmp[(2*i + 1)*tstride + 2*j + 1];
        OD_HAAR_KERNEL(a, b, c, d);
        tmp[i*tstride + j] = a;
        y[i*ystride + j + npairs] = b;
        y[(i + npairs)*ystride + j] = c;
        y[(i + npairs)*ystride + j + npairs] = d;
      }
    }
  }
  y[0] = tmp[0];
}

void od_haar_inv(od_coeff *x, int xstride,
 const od_coeff *y, int ystride, int ln) {
  int i;
  int j;
  int level;
  x[0] = y[0];
  for (level = ln - 1; level >= 0; level--) {
    int npairs;
    npairs = 1 << (ln - 1 - level);
    for (i = npairs - 1; i >= 0; i--) {
      for (j = npairs - 1; j >= 0; j--) {
        od_coeff a;
        od_coeff b;
        od_coeff c;
        od_coeff d;
        a = x[i*xstride + j];
        b = y[i*ystride + j + npairs];
        c = y[(i + npairs)*ystride + j];
        d = y[(i + npairs)*ystride + j + npairs];
        OD_HAAR_KERNEL(a, b, c, d);
        x[2*i*xstride + 2*j] = a;
        x[(2*i + 1)*xstride + 2*j] = b;
        x[2*i*xstride + 2*j + 1] = c;
        x[(2*i + 1)*xstride + 2*j + 1] = d;
      }
    }
  }
}

void od_bin_fdct32x32(od_coeff *y, int ystride,
 const od_coeff *x, int xstride) {
  od_coeff z[32*32];
  int i;
  for (i = 0; i < 32; i++) od_bin_fdct32(z + 32*i, x + i, xstride);
  for (i = 0; i < 32; i++) od_bin_fdct32(y + ystride*i, z + i, 32);
}

void od_bin_idct32x32(od_coeff *x, int xstride,
 const od_coeff *y, int ystride) {
  od_coeff z[32*32];
  int i;
  for (i = 0; i < 32; i++) od_bin_idct32(z + i, 32, y + ystride*i);
  for (i = 0; i < 32; i++) od_bin_idct32(x + i, xstride, z + 32*i);
}

void od_bin_fdct64x64(od_coeff *y, int ystride,
 const od_coeff *x, int xstride) {
  od_coeff z[64*64];
  int i;
  for (i = 0; i < 64; i++) od_bin_fdct64(z + 64*i, x + i, xstride);
  for (i = 0; i < 64; i++) od_bin_fdct64(y + ystride*i, z + i, 64);
}

void od_bin_idct64x64(od_coeff *x, int xstride,
 const od_coeff *y, int ystride) {
  od_coeff z[64*64];
  int i;
  for (i = 0; i < 64; i++) od_bin_idct64(z + i, 64, y + ystride*i);
  for (i = 0; i < 64; i++) od_bin_idct64(x + i, xstride, z + 64*i);
}

#if defined(OD_CHECKASM)
# include <stdio.h>

void od_dct_check(int bs, const od_coeff *ref, const od_coeff *x,
 int xstride) {
  int failed;
  int i;
  int j;
  int n;
  n = 4 << bs;
  failed = 0;
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      if (ref[i + j*n] != x[i + j*xstride]) {
        fprintf(stderr, "ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         ref[i + j*n], x[i + j*xstride], i, j);
        failed = 1;
      }
    }
  }
  if (failed) {
    fprintf(stderr, "od_bin %ix%i check failed.\n",
     n, n);
  }
  OD_ASSERT(!failed);
}
#endif

#if OD_DCT_TEST
/*Test code.*/
# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <string.h>

#if defined(OD_X86ASM)
# include "x86/cpu.h"
# include "x86/x86int.h"
#endif

#if defined(OD_ARMASM)
# include "arm/cpu.h"
# include "arm/armint.h"
#endif

/*The auto-correlation coefficent. 0.95 is a common value.*/
# define INPUT_AUTOCORR (0.95)
# define INPUT_AUTOCORR_2 (INPUT_AUTOCORR*INPUT_AUTOCORR)
# define INPUT_AUTOCORR_4 (INPUT_AUTOCORR_2*INPUT_AUTOCORR_2)
# define INPUT_AUTOCORR_8 (INPUT_AUTOCORR_4*INPUT_AUTOCORR_4)
# define INPUT_AUTOCORR_16 (INPUT_AUTOCORR_8*INPUT_AUTOCORR_8)
# define INPUT_AUTOCORR_32 (INPUT_AUTOCORR_16*INPUT_AUTOCORR_16)

/*An autocorrelation table.
  A common model for natural-image input is an AR-0 process with an
   autocorrelation coefficient of 0.95.
  This table contains various powers of 0.95, so that
   AUTOCORR[i-j+63]==(0.95)**abs(i-j), for i,j in [0..63].
  This makes it easy to compute element i,j of the covariance matrix.*/
static const double AUTOCORR[127] = {
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*
   INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*
   INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_2,
  INPUT_AUTOCORR,
  1,
  INPUT_AUTOCORR,
  INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*
   INPUT_AUTOCORR,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*
   INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_32*INPUT_AUTOCORR_16*INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*
   INPUT_AUTOCORR_2*INPUT_AUTOCORR
};

/*State information for the IEEE 1180 pseudo-random number generator.*/
static int ieee1180_rand_x;

/*Initializes the IEEE 1180 random number generator.*/
static void ieee1180_srand(int seed) {
  ieee1180_rand_x = seed;
}

/*Computes a random number between -l and h, inclusive, according to the
   specification in IEEE Std 1180-1990, "IEEE Standard Specifications for the
   Implementations of 8x8 Inverse Discrete Cosine Transform."*/
static int ieee1180_rand(int l, int h) {
  double x;
  int i;
  /*We do an unsigned multiply here to avoid signed overflow (which is
     undefined).
    The results should not change on a 2's-complement machine, since the low
     32-bits of the multiply are the same either way.*/
  ieee1180_rand_x = ieee1180_rand_x*1103515245U + 12345;
  i = ieee1180_rand_x & 0x7FFFFFFE;
  x = i/(double)0x7FFFFFFF*(l + h + 1);
  return (int)x - l;
}

static int od_exit_code = EXIT_SUCCESS;

static const char *ieee1180_meets(double val, double limit) {
  int meets;
  meets = fabs(val) <= limit;
  if (!meets) od_exit_code = EXIT_FAILURE;
  return meets ? "meets" : "FAILS";
}

/*The number of different input ranges.*/
# define IEEE1180_NRANGES (3)

/*The number of blocks of data to generate.*/
# define IEEE1180_NBLOCKS (10000)

static const int IEEE1180_L[IEEE1180_NRANGES] = { 256, 5, 300 };
static const int IEEE1180_H[IEEE1180_NRANGES] = { 255, 5, 300 };

/*TODO: This unscaling introduces noticeable bias, since all ties round up.*/
# define OD_COEFF_UNSCALE(x) \
 (((x) + (1 << OD_COEFF_SHIFT >> 1)) >> OD_COEFF_SHIFT)

typedef double od_dct_basis_row[OD_BSIZE_MAX];

/*The true forward 4-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const od_dct_basis_row OD_DCT4_BASIS[4] = {
  {
    0.5,                                 0.5,
    0.5,                                 0.5
  },
  {
     0.65328148243818826392832158671359,  0.27059805007309849219986160268319,
    -0.27059805007309849219986160268319, -0.65328148243818826392832158671359
  },
  {
     0.5,                                -0.5,
    -0.5,                                 0.5
  },
  {
    0.27059805007309849219986160268319, -0.65328148243818826392832158671359,
    0.65328148243818826392832158671359, -0.27059805007309849219986160268319
  },
};

/*The true forward 8-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const od_dct_basis_row OD_DCT8_BASIS[8] = {
  {
    0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
    0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
    0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
    0.35355339059327376220042218105242,  0.35355339059327376220042218105242
  },
  {
     0.49039264020161522456309111806712,  0.41573480615127261853939418880895,
     0.27778511650980111237141540697427,  0.097545161008064133924142434238511,
    -0.097545161008064133924142434238511,-0.27778511650980111237141540697427,
    -0.41573480615127261853939418880895, -0.49039264020161522456309111806712
  },
  {
     0.46193976625574337806409159469839,  0.19134171618254488586422999201520,
    -0.19134171618254488586422999201520, -0.46193976625574337806409159469839,
    -0.46193976625574337806409159469839, -0.19134171618254488586422999201520,
     0.19134171618254488586422999201520,  0.46193976625574337806409159469839
  },
  {
     0.41573480615127261853939418880895, -0.097545161008064133924142434238511,
    -0.49039264020161522456309111806712, -0.27778511650980111237141540697427,
     0.27778511650980111237141540697427,  0.49039264020161522456309111806712,
     0.097545161008064133924142434238511,-0.41573480615127261853939418880895
  },
  {
     0.35355339059327376220042218105242, -0.35355339059327376220042218105242,
    -0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242, -0.35355339059327376220042218105242,
    -0.35355339059327376220042218105242,  0.35355339059327376220042218105242
  },
  {
     0.27778511650980111237141540697427, -0.49039264020161522456309111806712,
     0.097545161008064133924142434238511, 0.41573480615127261853939418880895,
    -0.41573480615127261853939418880895, -0.097545161008064133924142434238511,
     0.49039264020161522456309111806712, -0.27778511650980111237141540697427
  },
  {
     0.19134171618254488586422999201520, -0.46193976625574337806409159469839,
     0.46193976625574337806409159469839, -0.19134171618254488586422999201520,
    -0.19134171618254488586422999201520,  0.46193976625574337806409159469839,
    -0.46193976625574337806409159469839,  0.19134171618254488586422999201520
  },
  {
     0.097545161008064133924142434238511,-0.27778511650980111237141540697427,
     0.41573480615127261853939418880895, -0.49039264020161522456309111806712,
     0.49039264020161522456309111806712, -0.41573480615127261853939418880895,
     0.27778511650980111237141540697427, -0.097545161008064133924142434238511
  }
};

/*The true forward 16-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const od_dct_basis_row OD_DCT16_BASIS[16] = {
  {
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25
  },
  {
     0.35185093438159561475651955574599,  0.33832950029358816956728612239141,
     0.31180625324666780814299756569942,  0.27330046675043937205975812924953,
     0.22429189658565907106266034730521,  0.16666391461943662432490137779072,
     0.10263113188058934528909230943878,  0.034654292299772865648933749133244,
    -0.034654292299772865648933749133244,-0.10263113188058934528909230943878,
    -0.16666391461943662432490137779072, -0.22429189658565907106266034730521,
    -0.27330046675043937205975812924953, -0.31180625324666780814299756569942,
    -0.33832950029358816956728612239141, -0.35185093438159561475651955574599
  },
  {
     0.34675996133053686545540479789161,  0.29396890060483967924361677615282,
     0.19642373959677554531947434191430,  0.068974844820735753083989390917343,
    -0.068974844820735753083989390917343,-0.19642373959677554531947434191430,
    -0.29396890060483967924361677615282, -0.34675996133053686545540479789161,
    -0.34675996133053686545540479789161, -0.29396890060483967924361677615282,
    -0.19642373959677554531947434191430, -0.068974844820735753083989390917343,
     0.068974844820735753083989390917343, 0.19642373959677554531947434191430,
     0.29396890060483967924361677615282,  0.34675996133053686545540479789161
  },
  {
     0.33832950029358816956728612239141,  0.22429189658565907106266034730521,
     0.034654292299772865648933749133244,-0.16666391461943662432490137779072,
    -0.31180625324666780814299756569942, -0.35185093438159561475651955574599,
    -0.27330046675043937205975812924953, -0.10263113188058934528909230943878,
     0.10263113188058934528909230943878,  0.27330046675043937205975812924953,
     0.35185093438159561475651955574599,  0.31180625324666780814299756569942,
     0.16666391461943662432490137779072, -0.034654292299772865648933749133244,
    -0.22429189658565907106266034730521, -0.33832950029358816956728612239141
  },
  {
     0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
     0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
     0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
     0.13529902503654924609993080134160,  0.32664074121909413196416079335680
  },
  {
     0.31180625324666780814299756569942,  0.034654292299772865648933749133244,
    -0.27330046675043937205975812924953, -0.33832950029358816956728612239141,
    -0.10263113188058934528909230943878,  0.22429189658565907106266034730521,
     0.35185093438159561475651955574599,  0.16666391461943662432490137779072,
    -0.16666391461943662432490137779072, -0.35185093438159561475651955574599,
    -0.22429189658565907106266034730521,  0.10263113188058934528909230943878,
     0.33832950029358816956728612239141,  0.27330046675043937205975812924953,
    -0.034654292299772865648933749133244,-0.31180625324666780814299756569942
  },
  {
     0.29396890060483967924361677615282, -0.068974844820735753083989390917343,
    -0.34675996133053686545540479789161, -0.19642373959677554531947434191430,
     0.19642373959677554531947434191430,  0.34675996133053686545540479789161,
     0.068974844820735753083989390917343,-0.29396890060483967924361677615282,
    -0.29396890060483967924361677615282,  0.068974844820735753083989390917343,
     0.34675996133053686545540479789161,  0.19642373959677554531947434191430,
    -0.19642373959677554531947434191430, -0.34675996133053686545540479789161,
    -0.068974844820735753083989390917343, 0.29396890060483967924361677615282
  },
  {
     0.27330046675043937205975812924953, -0.16666391461943662432490137779072,
    -0.33832950029358816956728612239141,  0.034654292299772865648933749133244,
     0.35185093438159561475651955574599,  0.10263113188058934528909230943878,
    -0.31180625324666780814299756569942, -0.22429189658565907106266034730521,
     0.22429189658565907106266034730521,  0.31180625324666780814299756569942,
    -0.10263113188058934528909230943878, -0.35185093438159561475651955574599,
    -0.034654292299772865648933749133244, 0.33832950029358816956728612239141,
     0.16666391461943662432490137779072, -0.27330046675043937205975812924953
  },
  {
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25
  },
  {
     0.22429189658565907106266034730521, -0.31180625324666780814299756569942,
    -0.10263113188058934528909230943878,  0.35185093438159561475651955574599,
    -0.034654292299772865648933749133244,-0.33832950029358816956728612239141,
     0.16666391461943662432490137779072,  0.27330046675043937205975812924953,
    -0.27330046675043937205975812924953, -0.16666391461943662432490137779072,
     0.33832950029358816956728612239141,  0.034654292299772865648933749133244,
    -0.35185093438159561475651955574599,  0.10263113188058934528909230943878,
     0.31180625324666780814299756569942, -0.22429189658565907106266034730521
  },
  {
     0.19642373959677554531947434191430, -0.34675996133053686545540479789161,
     0.068974844820735753083989390917343, 0.29396890060483967924361677615282,
    -0.29396890060483967924361677615282, -0.068974844820735753083989390917343,
     0.34675996133053686545540479789161, -0.19642373959677554531947434191430,
    -0.19642373959677554531947434191430,  0.34675996133053686545540479789161,
    -0.068974844820735753083989390917343,-0.29396890060483967924361677615282,
     0.29396890060483967924361677615282,  0.068974844820735753083989390917343,
    -0.34675996133053686545540479789161,  0.19642373959677554531947434191430
  },
  {
     0.16666391461943662432490137779072, -0.35185093438159561475651955574599,
     0.22429189658565907106266034730521,  0.10263113188058934528909230943878,
    -0.33832950029358816956728612239141,  0.27330046675043937205975812924953,
     0.034654292299772865648933749133244,-0.31180625324666780814299756569942,
     0.31180625324666780814299756569942, -0.034654292299772865648933749133244,
    -0.27330046675043937205975812924953,  0.33832950029358816956728612239141,
    -0.10263113188058934528909230943878, -0.22429189658565907106266034730521,
     0.35185093438159561475651955574599, -0.16666391461943662432490137779072
  },
  {
     0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
     0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
     0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
     0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680,  0.13529902503654924609993080134160
  },
  {
     0.10263113188058934528909230943878, -0.27330046675043937205975812924953,
     0.35185093438159561475651955574599, -0.31180625324666780814299756569942,
     0.16666391461943662432490137779072,  0.034654292299772865648933749133244,
    -0.22429189658565907106266034730521,  0.33832950029358816956728612239141,
    -0.33832950029358816956728612239141,  0.22429189658565907106266034730521,
    -0.034654292299772865648933749133244,-0.16666391461943662432490137779072,
     0.31180625324666780814299756569942, -0.35185093438159561475651955574599,
     0.27330046675043937205975812924953, -0.10263113188058934528909230943878
  },
  {
     0.068974844820735753083989390917343,-0.19642373959677554531947434191430,
     0.29396890060483967924361677615282, -0.34675996133053686545540479789161,
     0.34675996133053686545540479789161, -0.29396890060483967924361677615282,
     0.19642373959677554531947434191430, -0.068974844820735753083989390917343,
    -0.068974844820735753083989390917343, 0.19642373959677554531947434191430,
    -0.29396890060483967924361677615282,  0.34675996133053686545540479789161,
    -0.34675996133053686545540479789161,  0.29396890060483967924361677615282,
    -0.19642373959677554531947434191430,  0.068974844820735753083989390917343
  },
  {
     0.034654292299772865648933749133244,-0.10263113188058934528909230943878,
     0.16666391461943662432490137779072, -0.22429189658565907106266034730521,
     0.27330046675043937205975812924953, -0.31180625324666780814299756569942,
     0.33832950029358816956728612239141, -0.35185093438159561475651955574599,
     0.35185093438159561475651955574599, -0.33832950029358816956728612239141,
     0.31180625324666780814299756569942, -0.27330046675043937205975812924953,
     0.22429189658565907106266034730521, -0.16666391461943662432490137779072,
     0.10263113188058934528909230943878, -0.034654292299772865648933749133244
  }
};

/*The true forward 32-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const od_dct_basis_row OD_DCT32_BASIS[32] = {
  {
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621,  0.17677669529663688110021109052621
  },
  {
     0.24969886405129309817869290118978,  0.24729412749119524336291843450406,
     0.24250781329863599815099605182153,  0.23538601629575519460312735064988,
     0.22599732328086083289655007430763,  0.21443215250006801747556749607119,
     0.20080188287016122745166912824079,  0.18523778133873977279390422437379,
     0.16788973871175460015634421260686,  0.14892482612310833586675913220749,
     0.12852568604830543164842345974220,  0.10688877335757052358024171422220,
     0.084222463348055012672313303154787, 0.060745044975815972487068540519368,
     0.036682618613840437914712532411679, 0.012266918581854503563738744235671,
    -0.012266918581854503563738744235671,-0.036682618613840437914712532411679,
    -0.060745044975815972487068540519368,-0.084222463348055012672313303154787,
    -0.10688877335757052358024171422220, -0.12852568604830543164842345974220,
    -0.14892482612310833586675913220749, -0.16788973871175460015634421260686,
    -0.18523778133873977279390422437379, -0.20080188287016122745166912824079,
    -0.21443215250006801747556749607119, -0.22599732328086083289655007430763,
    -0.23538601629575519460312735064988, -0.24250781329863599815099605182153,
    -0.24729412749119524336291843450406, -0.24969886405129309817869290118978
  },
  {
     0.24879618166804922156120923827737,  0.23923508393305221623394947174507,
     0.22048031608708875742818921591510,  0.19325261334068424020272665243962,
     0.15859832104091137455379290330637,  0.11784918420649941213909690647631,
     0.072571169313615591909048093954349, 0.024504285082390150498548890972160,
    -0.024504285082390150498548890972160,-0.072571169313615591909048093954349,
    -0.11784918420649941213909690647631, -0.15859832104091137455379290330637,
    -0.19325261334068424020272665243962, -0.22048031608708875742818921591510,
    -0.23923508393305221623394947174507, -0.24879618166804922156120923827737,
    -0.24879618166804922156120923827737, -0.23923508393305221623394947174507,
    -0.22048031608708875742818921591510, -0.19325261334068424020272665243962,
    -0.15859832104091137455379290330637, -0.11784918420649941213909690647631,
    -0.072571169313615591909048093954349,-0.024504285082390150498548890972160,
     0.024504285082390150498548890972160, 0.072571169313615591909048093954349,
     0.11784918420649941213909690647631,  0.15859832104091137455379290330637,
     0.19325261334068424020272665243962,  0.22048031608708875742818921591510,
     0.23923508393305221623394947174507,  0.24879618166804922156120923827737
  },
  {
     0.24729412749119524336291843450406,  0.22599732328086083289655007430763,
     0.18523778133873977279390422437379,  0.12852568604830543164842345974220,
     0.060745044975815972487068540519368,-0.012266918581854503563738744235671,
    -0.084222463348055012672313303154787,-0.14892482612310833586675913220749,
    -0.20080188287016122745166912824079, -0.23538601629575519460312735064988,
    -0.24969886405129309817869290118978, -0.24250781329863599815099605182153,
    -0.21443215250006801747556749607119, -0.16788973871175460015634421260686,
    -0.10688877335757052358024171422220, -0.036682618613840437914712532411679,
     0.036682618613840437914712532411679, 0.10688877335757052358024171422220,
     0.16788973871175460015634421260686,  0.21443215250006801747556749607119,
     0.24250781329863599815099605182153,  0.24969886405129309817869290118978,
     0.23538601629575519460312735064988,  0.20080188287016122745166912824079,
     0.14892482612310833586675913220749,  0.084222463348055012672313303154787,
     0.012266918581854503563738744235671,-0.060745044975815972487068540519368,
    -0.12852568604830543164842345974220, -0.18523778133873977279390422437379,
    -0.22599732328086083289655007430763, -0.24729412749119524336291843450406
  },
  {
     0.24519632010080761228154555903356,  0.20786740307563630926969709440448,
     0.13889255825490055618570770348713,  0.048772580504032066962071217119256,
    -0.048772580504032066962071217119256,-0.13889255825490055618570770348713,
    -0.20786740307563630926969709440448, -0.24519632010080761228154555903356,
    -0.24519632010080761228154555903356, -0.20786740307563630926969709440448,
    -0.13889255825490055618570770348713, -0.048772580504032066962071217119256,
     0.048772580504032066962071217119256, 0.13889255825490055618570770348713,
     0.20786740307563630926969709440448,  0.24519632010080761228154555903356,
     0.24519632010080761228154555903356,  0.20786740307563630926969709440448,
     0.13889255825490055618570770348713,  0.048772580504032066962071217119256,
    -0.048772580504032066962071217119256,-0.13889255825490055618570770348713,
    -0.20786740307563630926969709440448, -0.24519632010080761228154555903356,
    -0.24519632010080761228154555903356, -0.20786740307563630926969709440448,
    -0.13889255825490055618570770348713, -0.048772580504032066962071217119256,
     0.048772580504032066962071217119256, 0.13889255825490055618570770348713,
     0.20786740307563630926969709440448,  0.24519632010080761228154555903356
  },
  {
     0.24250781329863599815099605182153,  0.18523778133873977279390422437379,
     0.084222463348055012672313303154787,-0.036682618613840437914712532411679,
    -0.14892482612310833586675913220749, -0.22599732328086083289655007430763,
    -0.24969886405129309817869290118978, -0.21443215250006801747556749607119,
    -0.12852568604830543164842345974220, -0.012266918581854503563738744235671,
     0.10688877335757052358024171422220,  0.20080188287016122745166912824079,
     0.24729412749119524336291843450406,  0.23538601629575519460312735064988,
     0.16788973871175460015634421260686,  0.060745044975815972487068540519368,
    -0.060745044975815972487068540519368,-0.16788973871175460015634421260686,
    -0.23538601629575519460312735064988, -0.24729412749119524336291843450406,
    -0.20080188287016122745166912824079, -0.10688877335757052358024171422220,
     0.012266918581854503563738744235671, 0.12852568604830543164842345974220,
     0.21443215250006801747556749607119,  0.24969886405129309817869290118978,
     0.22599732328086083289655007430763,  0.14892482612310833586675913220749,
     0.036682618613840437914712532411679,-0.084222463348055012672313303154787,
    -0.18523778133873977279390422437379, -0.24250781329863599815099605182153
  },
  {
     0.23923508393305221623394947174507,  0.15859832104091137455379290330637,
     0.024504285082390150498548890972160,-0.11784918420649941213909690647631,
    -0.22048031608708875742818921591510, -0.24879618166804922156120923827737,
    -0.19325261334068424020272665243962, -0.072571169313615591909048093954349,
     0.072571169313615591909048093954349, 0.19325261334068424020272665243962,
     0.24879618166804922156120923827737,  0.22048031608708875742818921591510,
     0.11784918420649941213909690647631, -0.024504285082390150498548890972160,
    -0.15859832104091137455379290330637, -0.23923508393305221623394947174507,
    -0.23923508393305221623394947174507, -0.15859832104091137455379290330637,
    -0.024504285082390150498548890972160, 0.11784918420649941213909690647631,
     0.22048031608708875742818921591510,  0.24879618166804922156120923827737,
     0.19325261334068424020272665243962,  0.072571169313615591909048093954349,
    -0.072571169313615591909048093954349,-0.19325261334068424020272665243962,
    -0.24879618166804922156120923827737, -0.22048031608708875742818921591510,
    -0.11784918420649941213909690647631,  0.024504285082390150498548890972160,
     0.15859832104091137455379290330637,  0.23923508393305221623394947174507
  },
  {
     0.23538601629575519460312735064988,  0.12852568604830543164842345974220,
    -0.036682618613840437914712532411679,-0.18523778133873977279390422437379,
    -0.24969886405129309817869290118978, -0.20080188287016122745166912824079,
    -0.060745044975815972487068540519368, 0.10688877335757052358024171422220,
     0.22599732328086083289655007430763,  0.24250781329863599815099605182153,
     0.14892482612310833586675913220749, -0.012266918581854503563738744235671,
    -0.16788973871175460015634421260686, -0.24729412749119524336291843450406,
    -0.21443215250006801747556749607119, -0.084222463348055012672313303154787,
     0.084222463348055012672313303154787, 0.21443215250006801747556749607119,
     0.24729412749119524336291843450406,  0.16788973871175460015634421260686,
     0.012266918581854503563738744235671,-0.14892482612310833586675913220749,
    -0.24250781329863599815099605182153, -0.22599732328086083289655007430763,
    -0.10688877335757052358024171422220,  0.060745044975815972487068540519368,
     0.20080188287016122745166912824079,  0.24969886405129309817869290118978,
     0.18523778133873977279390422437379,  0.036682618613840437914712532411679,
    -0.12852568604830543164842345974220, -0.23538601629575519460312735064988
  },
  {
     0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
     0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
     0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
     0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
     0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
     0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
     0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
     0.095670858091272442932114996007600, 0.23096988312782168903204579734920
  },
  {
     0.22599732328086083289655007430763,  0.060745044975815972487068540519368,
    -0.14892482612310833586675913220749, -0.24969886405129309817869290118978,
    -0.16788973871175460015634421260686,  0.036682618613840437914712532411679,
     0.21443215250006801747556749607119,  0.23538601629575519460312735064988,
     0.084222463348055012672313303154787,-0.12852568604830543164842345974220,
    -0.24729412749119524336291843450406, -0.18523778133873977279390422437379,
     0.012266918581854503563738744235671, 0.20080188287016122745166912824079,
     0.24250781329863599815099605182153,  0.10688877335757052358024171422220,
    -0.10688877335757052358024171422220, -0.24250781329863599815099605182153,
    -0.20080188287016122745166912824079, -0.012266918581854503563738744235671,
     0.18523778133873977279390422437379,  0.24729412749119524336291843450406,
     0.12852568604830543164842345974220, -0.084222463348055012672313303154787,
    -0.23538601629575519460312735064988, -0.21443215250006801747556749607119,
    -0.036682618613840437914712532411679, 0.16788973871175460015634421260686,
     0.24969886405129309817869290118978,  0.14892482612310833586675913220749,
    -0.060745044975815972487068540519368,-0.22599732328086083289655007430763
  },
  {
     0.22048031608708875742818921591510,  0.024504285082390150498548890972160,
    -0.19325261334068424020272665243962, -0.23923508393305221623394947174507,
    -0.072571169313615591909048093954349, 0.15859832104091137455379290330637,
     0.24879618166804922156120923827737,  0.11784918420649941213909690647631,
    -0.11784918420649941213909690647631, -0.24879618166804922156120923827737,
    -0.15859832104091137455379290330637,  0.072571169313615591909048093954349,
     0.23923508393305221623394947174507,  0.19325261334068424020272665243962,
    -0.024504285082390150498548890972160,-0.22048031608708875742818921591510,
    -0.22048031608708875742818921591510, -0.024504285082390150498548890972160,
     0.19325261334068424020272665243962,  0.23923508393305221623394947174507,
     0.072571169313615591909048093954349,-0.15859832104091137455379290330637,
    -0.24879618166804922156120923827737, -0.11784918420649941213909690647631,
     0.11784918420649941213909690647631,  0.24879618166804922156120923827737,
     0.15859832104091137455379290330637, -0.072571169313615591909048093954349,
    -0.23923508393305221623394947174507, -0.19325261334068424020272665243962,
     0.024504285082390150498548890972160, 0.22048031608708875742818921591510
  },
  {
     0.21443215250006801747556749607119, -0.012266918581854503563738744235671,
    -0.22599732328086083289655007430763, -0.20080188287016122745166912824079,
     0.036682618613840437914712532411679, 0.23538601629575519460312735064988,
     0.18523778133873977279390422437379, -0.060745044975815972487068540519368,
    -0.24250781329863599815099605182153, -0.16788973871175460015634421260686,
     0.084222463348055012672313303154787, 0.24729412749119524336291843450406,
     0.14892482612310833586675913220749, -0.10688877335757052358024171422220,
    -0.24969886405129309817869290118978, -0.12852568604830543164842345974220,
     0.12852568604830543164842345974220,  0.24969886405129309817869290118978,
     0.10688877335757052358024171422220, -0.14892482612310833586675913220749,
    -0.24729412749119524336291843450406, -0.084222463348055012672313303154787,
     0.16788973871175460015634421260686,  0.24250781329863599815099605182153,
     0.060745044975815972487068540519368,-0.18523778133873977279390422437379,
    -0.23538601629575519460312735064988, -0.036682618613840437914712532411679,
     0.20080188287016122745166912824079,  0.22599732328086083289655007430763,
     0.012266918581854503563738744235671,-0.21443215250006801747556749607119
  },
  {
     0.20786740307563630926969709440448, -0.048772580504032066962071217119256,
    -0.24519632010080761228154555903356, -0.13889255825490055618570770348713,
     0.13889255825490055618570770348713,  0.24519632010080761228154555903356,
     0.048772580504032066962071217119256,-0.20786740307563630926969709440448,
    -0.20786740307563630926969709440448,  0.048772580504032066962071217119256,
     0.24519632010080761228154555903356,  0.13889255825490055618570770348713,
    -0.13889255825490055618570770348713, -0.24519632010080761228154555903356,
    -0.048772580504032066962071217119256, 0.20786740307563630926969709440448,
     0.20786740307563630926969709440448, -0.048772580504032066962071217119256,
    -0.24519632010080761228154555903356, -0.13889255825490055618570770348713,
     0.13889255825490055618570770348713,  0.24519632010080761228154555903356,
     0.048772580504032066962071217119256,-0.20786740307563630926969709440448,
    -0.20786740307563630926969709440448,  0.048772580504032066962071217119256,
     0.24519632010080761228154555903356,  0.13889255825490055618570770348713,
    -0.13889255825490055618570770348713, -0.24519632010080761228154555903356,
    -0.048772580504032066962071217119256, 0.20786740307563630926969709440448
  },
  {
     0.20080188287016122745166912824079, -0.084222463348055012672313303154787,
    -0.24969886405129309817869290118978, -0.060745044975815972487068540519368,
     0.21443215250006801747556749607119,  0.18523778133873977279390422437379,
    -0.10688877335757052358024171422220, -0.24729412749119524336291843450406,
    -0.036682618613840437914712532411679, 0.22599732328086083289655007430763,
     0.16788973871175460015634421260686, -0.12852568604830543164842345974220,
    -0.24250781329863599815099605182153, -0.012266918581854503563738744235671,
     0.23538601629575519460312735064988,  0.14892482612310833586675913220749,
    -0.14892482612310833586675913220749, -0.23538601629575519460312735064988,
     0.012266918581854503563738744235671, 0.24250781329863599815099605182153,
     0.12852568604830543164842345974220, -0.16788973871175460015634421260686,
    -0.22599732328086083289655007430763,  0.036682618613840437914712532411679,
     0.24729412749119524336291843450406,  0.10688877335757052358024171422220,
    -0.18523778133873977279390422437379, -0.21443215250006801747556749607119,
     0.060745044975815972487068540519368, 0.24969886405129309817869290118978,
     0.084222463348055012672313303154787,-0.20080188287016122745166912824079
  },
  {
     0.19325261334068424020272665243962, -0.11784918420649941213909690647631,
    -0.23923508393305221623394947174507,  0.024504285082390150498548890972160,
     0.24879618166804922156120923827737,  0.072571169313615591909048093954349,
    -0.22048031608708875742818921591510, -0.15859832104091137455379290330637,
     0.15859832104091137455379290330637,  0.22048031608708875742818921591510,
    -0.072571169313615591909048093954349,-0.24879618166804922156120923827737,
    -0.024504285082390150498548890972160, 0.23923508393305221623394947174507,
     0.11784918420649941213909690647631, -0.19325261334068424020272665243962,
    -0.19325261334068424020272665243962,  0.11784918420649941213909690647631,
     0.23923508393305221623394947174507, -0.024504285082390150498548890972160,
    -0.24879618166804922156120923827737, -0.072571169313615591909048093954349,
     0.22048031608708875742818921591510,  0.15859832104091137455379290330637,
    -0.15859832104091137455379290330637, -0.22048031608708875742818921591510,
     0.072571169313615591909048093954349, 0.24879618166804922156120923827737,
     0.024504285082390150498548890972160,-0.23923508393305221623394947174507,
    -0.11784918420649941213909690647631,  0.19325261334068424020272665243962
  },
  {
     0.18523778133873977279390422437379, -0.14892482612310833586675913220749,
    -0.21443215250006801747556749607119,  0.10688877335757052358024171422220,
     0.23538601629575519460312735064988, -0.060745044975815972487068540519368,
    -0.24729412749119524336291843450406,  0.012266918581854503563738744235671,
     0.24969886405129309817869290118978,  0.036682618613840437914712532411679,
    -0.24250781329863599815099605182153, -0.084222463348055012672313303154787,
     0.22599732328086083289655007430763,  0.12852568604830543164842345974220,
    -0.20080188287016122745166912824079, -0.16788973871175460015634421260686,
     0.16788973871175460015634421260686,  0.20080188287016122745166912824079,
    -0.12852568604830543164842345974220, -0.22599732328086083289655007430763,
     0.084222463348055012672313303154787, 0.24250781329863599815099605182153,
    -0.036682618613840437914712532411679,-0.24969886405129309817869290118978,
    -0.012266918581854503563738744235671, 0.24729412749119524336291843450406,
     0.060745044975815972487068540519368,-0.23538601629575519460312735064988,
    -0.10688877335757052358024171422220,  0.21443215250006801747556749607119,
     0.14892482612310833586675913220749, -0.18523778133873977279390422437379
  },
  {
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621,
     0.17677669529663688110021109052621, -0.17677669529663688110021109052621,
    -0.17677669529663688110021109052621,  0.17677669529663688110021109052621
  },
  {
     0.16788973871175460015634421260686, -0.20080188287016122745166912824079,
    -0.12852568604830543164842345974220,  0.22599732328086083289655007430763,
     0.084222463348055012672313303154787,-0.24250781329863599815099605182153,
    -0.036682618613840437914712532411679, 0.24969886405129309817869290118978,
    -0.012266918581854503563738744235671,-0.24729412749119524336291843450406,
     0.060745044975815972487068540519368, 0.23538601629575519460312735064988,
    -0.10688877335757052358024171422220, -0.21443215250006801747556749607119,
     0.14892482612310833586675913220749,  0.18523778133873977279390422437379,
    -0.18523778133873977279390422437379, -0.14892482612310833586675913220749,
     0.21443215250006801747556749607119,  0.10688877335757052358024171422220,
    -0.23538601629575519460312735064988, -0.060745044975815972487068540519368,
     0.24729412749119524336291843450406,  0.012266918581854503563738744235671,
    -0.24969886405129309817869290118978,  0.036682618613840437914712532411679,
     0.24250781329863599815099605182153, -0.084222463348055012672313303154787,
    -0.22599732328086083289655007430763,  0.12852568604830543164842345974220,
     0.20080188287016122745166912824079, -0.16788973871175460015634421260686
  },
  {
     0.15859832104091137455379290330637, -0.22048031608708875742818921591510,
    -0.072571169313615591909048093954349, 0.24879618166804922156120923827737,
    -0.024504285082390150498548890972160,-0.23923508393305221623394947174507,
     0.11784918420649941213909690647631,  0.19325261334068424020272665243962,
    -0.19325261334068424020272665243962, -0.11784918420649941213909690647631,
     0.23923508393305221623394947174507,  0.024504285082390150498548890972160,
    -0.24879618166804922156120923827737,  0.072571169313615591909048093954349,
     0.22048031608708875742818921591510, -0.15859832104091137455379290330637,
    -0.15859832104091137455379290330637,  0.22048031608708875742818921591510,
     0.072571169313615591909048093954349,-0.24879618166804922156120923827737,
     0.024504285082390150498548890972160, 0.23923508393305221623394947174507,
    -0.11784918420649941213909690647631, -0.19325261334068424020272665243962,
     0.19325261334068424020272665243962,  0.11784918420649941213909690647631,
    -0.23923508393305221623394947174507, -0.024504285082390150498548890972160,
     0.24879618166804922156120923827737, -0.072571169313615591909048093954349,
    -0.22048031608708875742818921591510,  0.15859832104091137455379290330637
  },
  {
     0.14892482612310833586675913220749, -0.23538601629575519460312735064988,
    -0.012266918581854503563738744235671, 0.24250781329863599815099605182153,
    -0.12852568604830543164842345974220, -0.16788973871175460015634421260686,
     0.22599732328086083289655007430763,  0.036682618613840437914712532411679,
    -0.24729412749119524336291843450406,  0.10688877335757052358024171422220,
     0.18523778133873977279390422437379, -0.21443215250006801747556749607119,
    -0.060745044975815972487068540519368, 0.24969886405129309817869290118978,
    -0.084222463348055012672313303154787,-0.20080188287016122745166912824079,
     0.20080188287016122745166912824079,  0.084222463348055012672313303154787,
    -0.24969886405129309817869290118978,  0.060745044975815972487068540519368,
     0.21443215250006801747556749607119, -0.18523778133873977279390422437379,
    -0.10688877335757052358024171422220,  0.24729412749119524336291843450406,
    -0.036682618613840437914712532411679,-0.22599732328086083289655007430763,
     0.16788973871175460015634421260686,  0.12852568604830543164842345974220,
    -0.24250781329863599815099605182153,  0.012266918581854503563738744235671,
     0.23538601629575519460312735064988, -0.14892482612310833586675913220749
  },
  {
     0.13889255825490055618570770348713, -0.24519632010080761228154555903356,
     0.048772580504032066962071217119256, 0.20786740307563630926969709440448,
    -0.20786740307563630926969709440448, -0.048772580504032066962071217119256,
     0.24519632010080761228154555903356, -0.13889255825490055618570770348713,
    -0.13889255825490055618570770348713,  0.24519632010080761228154555903356,
    -0.048772580504032066962071217119256,-0.20786740307563630926969709440448,
     0.20786740307563630926969709440448,  0.048772580504032066962071217119256,
    -0.24519632010080761228154555903356,  0.13889255825490055618570770348713,
     0.13889255825490055618570770348713, -0.24519632010080761228154555903356,
     0.048772580504032066962071217119256, 0.20786740307563630926969709440448,
    -0.20786740307563630926969709440448, -0.048772580504032066962071217119256,
     0.24519632010080761228154555903356, -0.13889255825490055618570770348713,
    -0.13889255825490055618570770348713,  0.24519632010080761228154555903356,
    -0.048772580504032066962071217119256,-0.20786740307563630926969709440448,
     0.20786740307563630926969709440448,  0.048772580504032066962071217119256,
    -0.24519632010080761228154555903356,  0.13889255825490055618570770348713
  },
  {
     0.12852568604830543164842345974220, -0.24969886405129309817869290118978,
     0.10688877335757052358024171422220,  0.14892482612310833586675913220749,
    -0.24729412749119524336291843450406,  0.084222463348055012672313303154787,
     0.16788973871175460015634421260686, -0.24250781329863599815099605182153,
     0.060745044975815972487068540519368, 0.18523778133873977279390422437379,
    -0.23538601629575519460312735064988,  0.036682618613840437914712532411679,
     0.20080188287016122745166912824079, -0.22599732328086083289655007430763,
     0.012266918581854503563738744235671, 0.21443215250006801747556749607119,
    -0.21443215250006801747556749607119, -0.012266918581854503563738744235671,
     0.22599732328086083289655007430763, -0.20080188287016122745166912824079,
    -0.036682618613840437914712532411679, 0.23538601629575519460312735064988,
    -0.18523778133873977279390422437379, -0.060745044975815972487068540519368,
     0.24250781329863599815099605182153, -0.16788973871175460015634421260686,
    -0.084222463348055012672313303154787, 0.24729412749119524336291843450406,
    -0.14892482612310833586675913220749, -0.10688877335757052358024171422220,
     0.24969886405129309817869290118978, -0.12852568604830543164842345974220
  },
  {
     0.11784918420649941213909690647631, -0.24879618166804922156120923827737,
     0.15859832104091137455379290330637,  0.072571169313615591909048093954349,
    -0.23923508393305221623394947174507,  0.19325261334068424020272665243962,
     0.024504285082390150498548890972160,-0.22048031608708875742818921591510,
     0.22048031608708875742818921591510, -0.024504285082390150498548890972160,
    -0.19325261334068424020272665243962,  0.23923508393305221623394947174507,
    -0.072571169313615591909048093954349,-0.15859832104091137455379290330637,
     0.24879618166804922156120923827737, -0.11784918420649941213909690647631,
    -0.11784918420649941213909690647631,  0.24879618166804922156120923827737,
    -0.15859832104091137455379290330637, -0.072571169313615591909048093954349,
     0.23923508393305221623394947174507, -0.19325261334068424020272665243962,
    -0.024504285082390150498548890972160, 0.22048031608708875742818921591510,
    -0.22048031608708875742818921591510,  0.024504285082390150498548890972160,
     0.19325261334068424020272665243962, -0.23923508393305221623394947174507,
     0.072571169313615591909048093954349, 0.15859832104091137455379290330637,
    -0.24879618166804922156120923827737,  0.11784918420649941213909690647631
  },
  {
     0.10688877335757052358024171422220, -0.24250781329863599815099605182153,
     0.20080188287016122745166912824079, -0.012266918581854503563738744235671,
    -0.18523778133873977279390422437379,  0.24729412749119524336291843450406,
    -0.12852568604830543164842345974220, -0.084222463348055012672313303154787,
     0.23538601629575519460312735064988, -0.21443215250006801747556749607119,
     0.036682618613840437914712532411679, 0.16788973871175460015634421260686,
    -0.24969886405129309817869290118978,  0.14892482612310833586675913220749,
     0.060745044975815972487068540519368,-0.22599732328086083289655007430763,
     0.22599732328086083289655007430763, -0.060745044975815972487068540519368,
    -0.14892482612310833586675913220749,  0.24969886405129309817869290118978,
    -0.16788973871175460015634421260686, -0.036682618613840437914712532411679,
     0.21443215250006801747556749607119, -0.23538601629575519460312735064988,
     0.084222463348055012672313303154787, 0.12852568604830543164842345974220,
    -0.24729412749119524336291843450406,  0.18523778133873977279390422437379,
     0.012266918581854503563738744235671,-0.20080188287016122745166912824079,
     0.24250781329863599815099605182153, -0.10688877335757052358024171422220
  },
  {
     0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
     0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
     0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
     0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
     0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
     0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920,  0.095670858091272442932114996007600,
     0.095670858091272442932114996007600,-0.23096988312782168903204579734920,
     0.23096988312782168903204579734920, -0.095670858091272442932114996007600,
    -0.095670858091272442932114996007600, 0.23096988312782168903204579734920,
    -0.23096988312782168903204579734920,  0.095670858091272442932114996007600
  },
  {
     0.084222463348055012672313303154787,-0.21443215250006801747556749607119,
     0.24729412749119524336291843450406, -0.16788973871175460015634421260686,
     0.012266918581854503563738744235671, 0.14892482612310833586675913220749,
    -0.24250781329863599815099605182153,  0.22599732328086083289655007430763,
    -0.10688877335757052358024171422220, -0.060745044975815972487068540519368,
     0.20080188287016122745166912824079, -0.24969886405129309817869290118978,
     0.18523778133873977279390422437379, -0.036682618613840437914712532411679,
    -0.12852568604830543164842345974220,  0.23538601629575519460312735064988,
    -0.23538601629575519460312735064988,  0.12852568604830543164842345974220,
     0.036682618613840437914712532411679,-0.18523778133873977279390422437379,
     0.24969886405129309817869290118978, -0.20080188287016122745166912824079,
     0.060745044975815972487068540519368, 0.10688877335757052358024171422220,
    -0.22599732328086083289655007430763,  0.24250781329863599815099605182153,
    -0.14892482612310833586675913220749, -0.012266918581854503563738744235671,
     0.16788973871175460015634421260686, -0.24729412749119524336291843450406,
     0.21443215250006801747556749607119, -0.084222463348055012672313303154787
  },
  {
     0.072571169313615591909048093954349,-0.19325261334068424020272665243962,
     0.24879618166804922156120923827737, -0.22048031608708875742818921591510,
     0.11784918420649941213909690647631,  0.024504285082390150498548890972160,
    -0.15859832104091137455379290330637,  0.23923508393305221623394947174507,
    -0.23923508393305221623394947174507,  0.15859832104091137455379290330637,
    -0.024504285082390150498548890972160,-0.11784918420649941213909690647631,
     0.22048031608708875742818921591510, -0.24879618166804922156120923827737,
     0.19325261334068424020272665243962, -0.072571169313615591909048093954349,
    -0.072571169313615591909048093954349, 0.19325261334068424020272665243962,
    -0.24879618166804922156120923827737,  0.22048031608708875742818921591510,
    -0.11784918420649941213909690647631, -0.024504285082390150498548890972160,
     0.15859832104091137455379290330637, -0.23923508393305221623394947174507,
     0.23923508393305221623394947174507, -0.15859832104091137455379290330637,
     0.024504285082390150498548890972160, 0.11784918420649941213909690647631,
    -0.22048031608708875742818921591510,  0.24879618166804922156120923827737,
    -0.19325261334068424020272665243962,  0.072571169313615591909048093954349
  },
  {
     0.060745044975815972487068540519368,-0.16788973871175460015634421260686,
     0.23538601629575519460312735064988, -0.24729412749119524336291843450406,
     0.20080188287016122745166912824079, -0.10688877335757052358024171422220,
    -0.012266918581854503563738744235671, 0.12852568604830543164842345974220,
    -0.21443215250006801747556749607119,  0.24969886405129309817869290118978,
    -0.22599732328086083289655007430763,  0.14892482612310833586675913220749,
    -0.036682618613840437914712532411679,-0.084222463348055012672313303154787,
     0.18523778133873977279390422437379, -0.24250781329863599815099605182153,
     0.24250781329863599815099605182153, -0.18523778133873977279390422437379,
     0.084222463348055012672313303154787, 0.036682618613840437914712532411679,
    -0.14892482612310833586675913220749,  0.22599732328086083289655007430763,
    -0.24969886405129309817869290118978,  0.21443215250006801747556749607119,
    -0.12852568604830543164842345974220,  0.012266918581854503563738744235671,
     0.10688877335757052358024171422220, -0.20080188287016122745166912824079,
     0.24729412749119524336291843450406, -0.23538601629575519460312735064988,
     0.16788973871175460015634421260686, -0.060745044975815972487068540519368
  },
  {
     0.048772580504032066962071217119256,-0.13889255825490055618570770348713,
     0.20786740307563630926969709440448, -0.24519632010080761228154555903356,
     0.24519632010080761228154555903356, -0.20786740307563630926969709440448,
     0.13889255825490055618570770348713, -0.048772580504032066962071217119256,
    -0.048772580504032066962071217119256, 0.13889255825490055618570770348713,
    -0.20786740307563630926969709440448,  0.24519632010080761228154555903356,
    -0.24519632010080761228154555903356,  0.20786740307563630926969709440448,
    -0.13889255825490055618570770348713,  0.048772580504032066962071217119256,
     0.048772580504032066962071217119256,-0.13889255825490055618570770348713,
     0.20786740307563630926969709440448, -0.24519632010080761228154555903356,
     0.24519632010080761228154555903356, -0.20786740307563630926969709440448,
     0.13889255825490055618570770348713, -0.048772580504032066962071217119256,
    -0.048772580504032066962071217119256, 0.13889255825490055618570770348713,
    -0.20786740307563630926969709440448,  0.24519632010080761228154555903356,
    -0.24519632010080761228154555903356,  0.20786740307563630926969709440448,
    -0.13889255825490055618570770348713,  0.048772580504032066962071217119256
  },
  {
     0.036682618613840437914712532411679,-0.10688877335757052358024171422220,
     0.16788973871175460015634421260686, -0.21443215250006801747556749607119,
     0.24250781329863599815099605182153, -0.24969886405129309817869290118978,
     0.23538601629575519460312735064988, -0.20080188287016122745166912824079,
     0.14892482612310833586675913220749, -0.084222463348055012672313303154787,
     0.012266918581854503563738744235671, 0.060745044975815972487068540519368,
    -0.12852568604830543164842345974220,  0.18523778133873977279390422437379,
    -0.22599732328086083289655007430763,  0.24729412749119524336291843450406,
    -0.24729412749119524336291843450406,  0.22599732328086083289655007430763,
    -0.18523778133873977279390422437379,  0.12852568604830543164842345974220,
    -0.060745044975815972487068540519368,-0.012266918581854503563738744235671,
     0.084222463348055012672313303154787,-0.14892482612310833586675913220749,
     0.20080188287016122745166912824079, -0.23538601629575519460312735064988,
     0.24969886405129309817869290118978, -0.24250781329863599815099605182153,
     0.21443215250006801747556749607119, -0.16788973871175460015634421260686,
     0.10688877335757052358024171422220, -0.036682618613840437914712532411679
  },
  {
     0.024504285082390150498548890972160,-0.072571169313615591909048093954349,
     0.11784918420649941213909690647631, -0.15859832104091137455379290330637,
     0.19325261334068424020272665243962, -0.22048031608708875742818921591510,
     0.23923508393305221623394947174507, -0.24879618166804922156120923827737,
     0.24879618166804922156120923827737, -0.23923508393305221623394947174507,
     0.22048031608708875742818921591510, -0.19325261334068424020272665243962,
     0.15859832104091137455379290330637, -0.11784918420649941213909690647631,
     0.072571169313615591909048093954349,-0.024504285082390150498548890972160,
    -0.024504285082390150498548890972160, 0.072571169313615591909048093954349,
    -0.11784918420649941213909690647631,  0.15859832104091137455379290330637,
    -0.19325261334068424020272665243962,  0.22048031608708875742818921591510,
    -0.23923508393305221623394947174507,  0.24879618166804922156120923827737,
    -0.24879618166804922156120923827737,  0.23923508393305221623394947174507,
    -0.22048031608708875742818921591510,  0.19325261334068424020272665243962,
    -0.15859832104091137455379290330637,  0.11784918420649941213909690647631,
    -0.072571169313615591909048093954349, 0.024504285082390150498548890972160
  },
  {
     0.012266918581854503563738744235671,-0.036682618613840437914712532411679,
     0.060745044975815972487068540519368,-0.084222463348055012672313303154787,
     0.10688877335757052358024171422220, -0.12852568604830543164842345974220,
     0.14892482612310833586675913220749, -0.16788973871175460015634421260686,
     0.18523778133873977279390422437379, -0.20080188287016122745166912824079,
     0.21443215250006801747556749607119, -0.22599732328086083289655007430763,
     0.23538601629575519460312735064988, -0.24250781329863599815099605182153,
     0.24729412749119524336291843450406, -0.24969886405129309817869290118978,
     0.24969886405129309817869290118978, -0.24729412749119524336291843450406,
     0.24250781329863599815099605182153, -0.23538601629575519460312735064988,
     0.22599732328086083289655007430763, -0.21443215250006801747556749607119,
     0.20080188287016122745166912824079, -0.18523778133873977279390422437379,
     0.16788973871175460015634421260686, -0.14892482612310833586675913220749,
     0.12852568604830543164842345974220, -0.10688877335757052358024171422220,
     0.084222463348055012672313303154787,-0.060745044975815972487068540519368,
     0.036682618613840437914712532411679,-0.012266918581854503563738744235671
  }
};

/*The true forward 64-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const od_dct_basis_row OD_DCT64_BASIS[64] = {
  {
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125,
     0.125,                               0.125
  },
  {
     0.17672345346106673069472498415600,  0.17629771118253266395745838340776,
     0.17544725227590413843675978359671,  0.17417412557117862465636166724330,
     0.17248139814210053549033101371775,  0.17037314791731193194066488665533,
     0.16785445385626524004858887836470,  0.16493138371356506569293701686233,
     0.16161097942121587956084560348376,  0.15790124012399101798155186502293,
     0.15381110290879228162646885882803,  0.14935042127442479271101666423961,
     0.14452994139365530889220759879835,  0.13936127622474077458193396572651,
     0.13385687753479470495813284521973,  0.12803000590238956351843752239068,
     0.12189469877166149151050389002586,  0.11546573663487784907576643982592,
     0.10875860742493672495725784718877,  0.10178946920358000230129700949629,
     0.094575111235207343404082277365271, 0.087132913540067685045147230048925,
     0.079480805024268148817432153524406, 0.071637220287468842367069460402068,
     0.063621055212317597965816690486066, 0.055451621442613590037015233508646,
     0.047148599859865922722408214465630, 0.038731993170326232824561787558492,
     0.030222077716717299469868553789363, 0.021639354630747427065531392743213,
     0.013004500444088472936407629993685, 0.004338317276799999540542326666634,
    -0.004338317276799999540542326666634,-0.013004500444088472936407629993685,
    -0.021639354630747427065531392743213,-0.030222077716717299469868553789363,
    -0.038731993170326232824561787558492,-0.047148599859865922722408214465630,
    -0.055451621442613590037015233508646,-0.063621055212317597965816690486066,
    -0.071637220287468842367069460402068,-0.079480805024268148817432153524406,
    -0.087132913540067685045147230048925,-0.094575111235207343404082277365271,
    -0.10178946920358000230129700949629, -0.10875860742493672495725784718877,
    -0.11546573663487784907576643982592, -0.12189469877166149151050389002586,
    -0.12803000590238956351843752239068, -0.13385687753479470495813284521973,
    -0.13936127622474077458193396572651, -0.14452994139365530889220759879835,
    -0.14935042127442479271101666423961, -0.15381110290879228162646885882803,
    -0.15790124012399101798155186502293, -0.16161097942121587956084560348376,
    -0.16493138371356506569293701686233, -0.16785445385626524004858887836470,
    -0.17037314791731193194066488665533, -0.17248139814210053549033101371775,
    -0.17417412557117862465636166724330, -0.17544725227590413843675978359671,
    -0.17629771118253266395745838340776, -0.17672345346106673069472498415600
  },
  {
     0.17656376002524718647512421849032,  0.17486335449663478165921413022414,
     0.17147891927418672456199547790670,  0.16644304831921567823839589426492,
     0.15980423982190510363772032690233,  0.15162642913722598531903229617045,
     0.14198837305251784063881548345787,  0.13098289131657380087121582271272,
     0.11871597273471929730747707847705,  0.10530575443867740272410295104619,
     0.090881384161410012831963755651079, 0.075581776473850090965407023747544,
     0.059554274961645154658154180042717, 0.042953233225881292913572018164494,
     0.025938528373526445792454998016647, 0.008674021313492586318780005883468,
    -0.008674021313492586318780005883468,-0.025938528373526445792454998016647,
    -0.042953233225881292913572018164494,-0.059554274961645154658154180042717,
    -0.075581776473850090965407023747544,-0.090881384161410012831963755651079,
    -0.10530575443867740272410295104619, -0.11871597273471929730747707847705,
    -0.13098289131657380087121582271272, -0.14198837305251784063881548345787,
    -0.15162642913722598531903229617045, -0.15980423982190510363772032690233,
    -0.16644304831921567823839589426492, -0.17147891927418672456199547790670,
    -0.17486335449663478165921413022414, -0.17656376002524718647512421849032,
    -0.17656376002524718647512421849032, -0.17486335449663478165921413022414,
    -0.17147891927418672456199547790670, -0.16644304831921567823839589426492,
    -0.15980423982190510363772032690233, -0.15162642913722598531903229617045,
    -0.14198837305251784063881548345787, -0.13098289131657380087121582271272,
    -0.11871597273471929730747707847705, -0.10530575443867740272410295104619,
    -0.090881384161410012831963755651079,-0.075581776473850090965407023747544,
    -0.059554274961645154658154180042717,-0.042953233225881292913572018164494,
    -0.025938528373526445792454998016647,-0.008674021313492586318780005883468,
     0.008674021313492586318780005883468, 0.025938528373526445792454998016647,
     0.042953233225881292913572018164494, 0.059554274961645154658154180042717,
     0.075581776473850090965407023747544, 0.090881384161410012831963755651079,
     0.10530575443867740272410295104619,  0.11871597273471929730747707847705,
     0.13098289131657380087121582271272,  0.14198837305251784063881548345787,
     0.15162642913722598531903229617045,  0.15980423982190510363772032690233,
     0.16644304831921567823839589426492,  0.17147891927418672456199547790670,
     0.17486335449663478165921413022414,  0.17656376002524718647512421849032
  },
  {
     0.17629771118253266395745838340776,  0.17248139814210053549033101371775,
     0.16493138371356506569293701686233,  0.15381110290879228162646885882803,
     0.13936127622474077458193396572651,  0.12189469877166149151050389002586,
     0.10178946920358000230129700949629,  0.079480805024268148817432153524406,
     0.055451621442613590037015233508646, 0.030222077716717299469868553789363,
     0.004338317276799999540542326666634,-0.021639354630747427065531392743213,
    -0.047148599859865922722408214465630,-0.071637220287468842367069460402068,
    -0.094575111235207343404082277365271,-0.11546573663487784907576643982592,
    -0.13385687753479470495813284521973, -0.14935042127442479271101666423961,
    -0.16161097942121587956084560348376, -0.17037314791731193194066488665533,
    -0.17544725227590413843675978359671, -0.17672345346106673069472498415600,
    -0.17417412557117862465636166724330, -0.16785445385626524004858887836470,
    -0.15790124012399101798155186502293, -0.14452994139365530889220759879835,
    -0.12803000590238956351843752239068, -0.10875860742493672495725784718877,
    -0.087132913540067685045147230048925,-0.063621055212317597965816690486066,
    -0.038731993170326232824561787558492,-0.013004500444088472936407629993685,
     0.013004500444088472936407629993685, 0.038731993170326232824561787558492,
     0.063621055212317597965816690486066, 0.087132913540067685045147230048925,
     0.10875860742493672495725784718877,  0.12803000590238956351843752239068,
     0.14452994139365530889220759879835,  0.15790124012399101798155186502293,
     0.16785445385626524004858887836470,  0.17417412557117862465636166724330,
     0.17672345346106673069472498415600,  0.17544725227590413843675978359671,
     0.17037314791731193194066488665533,  0.16161097942121587956084560348376,
     0.14935042127442479271101666423961,  0.13385687753479470495813284521973,
     0.11546573663487784907576643982592,  0.094575111235207343404082277365271,
     0.071637220287468842367069460402068, 0.047148599859865922722408214465630,
     0.021639354630747427065531392743213,-0.004338317276799999540542326666634,
    -0.030222077716717299469868553789363,-0.055451621442613590037015233508646,
    -0.079480805024268148817432153524406,-0.10178946920358000230129700949629,
    -0.12189469877166149151050389002586, -0.13936127622474077458193396572651,
    -0.15381110290879228162646885882803, -0.16493138371356506569293701686233,
    -0.17248139814210053549033101371775, -0.17629771118253266395745838340776
  },
  {
     0.17592546719079780737825977787300,  0.16916475014679408478364306119571,
     0.15590312662333390407149878284971,  0.13665023337521968602987906462477,
     0.11214594829282953553133017365260,  0.083331957309718312162450688895359,
     0.051315565940294672644546154719392, 0.017327146149886432824466874566622,
    -0.017327146149886432824466874566622,-0.051315565940294672644546154719392,
    -0.083331957309718312162450688895359,-0.11214594829282953553133017365260,
    -0.13665023337521968602987906462477, -0.15590312662333390407149878284971,
    -0.16916475014679408478364306119571, -0.17592546719079780737825977787300,
    -0.17592546719079780737825977787300, -0.16916475014679408478364306119571,
    -0.15590312662333390407149878284971, -0.13665023337521968602987906462477,
    -0.11214594829282953553133017365260, -0.083331957309718312162450688895359,
    -0.051315565940294672644546154719392,-0.017327146149886432824466874566622,
     0.017327146149886432824466874566622, 0.051315565940294672644546154719392,
     0.083331957309718312162450688895359, 0.11214594829282953553133017365260,
     0.13665023337521968602987906462477,  0.15590312662333390407149878284971,
     0.16916475014679408478364306119571,  0.17592546719079780737825977787300,
     0.17592546719079780737825977787300,  0.16916475014679408478364306119571,
     0.15590312662333390407149878284971,  0.13665023337521968602987906462477,
     0.11214594829282953553133017365260,  0.083331957309718312162450688895359,
     0.051315565940294672644546154719392, 0.017327146149886432824466874566622,
    -0.017327146149886432824466874566622,-0.051315565940294672644546154719392,
    -0.083331957309718312162450688895359,-0.11214594829282953553133017365260,
    -0.13665023337521968602987906462477, -0.15590312662333390407149878284971,
    -0.16916475014679408478364306119571, -0.17592546719079780737825977787300,
    -0.17592546719079780737825977787300, -0.16916475014679408478364306119571,
    -0.15590312662333390407149878284971, -0.13665023337521968602987906462477,
    -0.11214594829282953553133017365260, -0.083331957309718312162450688895359,
    -0.051315565940294672644546154719392,-0.017327146149886432824466874566622,
     0.017327146149886432824466874566622, 0.051315565940294672644546154719392,
     0.083331957309718312162450688895359, 0.11214594829282953553133017365260,
     0.13665023337521968602987906462477,  0.15590312662333390407149878284971,
     0.16916475014679408478364306119571,  0.17592546719079780737825977787300
  },
  {
     0.17544725227590413843675978359671,  0.16493138371356506569293701686233,
     0.14452994139365530889220759879835,  0.11546573663487784907576643982592,
     0.079480805024268148817432153524406, 0.038731993170326232824561787558492,
    -0.004338317276799999540542326666634,-0.047148599859865922722408214465630,
    -0.087132913540067685045147230048925,-0.12189469877166149151050389002586,
    -0.14935042127442479271101666423961, -0.16785445385626524004858887836470,
    -0.17629771118253266395745838340776, -0.17417412557117862465636166724330,
    -0.16161097942121587956084560348376, -0.13936127622474077458193396572651,
    -0.10875860742493672495725784718877, -0.071637220287468842367069460402068,
    -0.030222077716717299469868553789363, 0.013004500444088472936407629993685,
     0.055451621442613590037015233508646, 0.094575111235207343404082277365271,
     0.12803000590238956351843752239068,  0.15381110290879228162646885882803,
     0.17037314791731193194066488665533,  0.17672345346106673069472498415600,
     0.17248139814210053549033101371775,  0.15790124012399101798155186502293,
     0.13385687753479470495813284521973,  0.10178946920358000230129700949629,
     0.063621055212317597965816690486066, 0.021639354630747427065531392743213,
    -0.021639354630747427065531392743213,-0.063621055212317597965816690486066,
    -0.10178946920358000230129700949629, -0.13385687753479470495813284521973,
    -0.15790124012399101798155186502293, -0.17248139814210053549033101371775,
    -0.17672345346106673069472498415600, -0.17037314791731193194066488665533,
    -0.15381110290879228162646885882803, -0.12803000590238956351843752239068,
    -0.094575111235207343404082277365271,-0.055451621442613590037015233508646,
    -0.013004500444088472936407629993685, 0.030222077716717299469868553789363,
     0.071637220287468842367069460402068, 0.10875860742493672495725784718877,
     0.13936127622474077458193396572651,  0.16161097942121587956084560348376,
     0.17417412557117862465636166724330,  0.17629771118253266395745838340776,
     0.16785445385626524004858887836470,  0.14935042127442479271101666423961,
     0.12189469877166149151050389002586,  0.087132913540067685045147230048925,
     0.047148599859865922722408214465630, 0.004338317276799999540542326666634,
    -0.038731993170326232824561787558492,-0.079480805024268148817432153524406,
    -0.11546573663487784907576643982592, -0.14452994139365530889220759879835,
    -0.16493138371356506569293701686233, -0.17544725227590413843675978359671
  },
  {
     0.17486335449663478165921413022414,  0.15980423982190510363772032690233,
     0.13098289131657380087121582271272,  0.090881384161410012831963755651079,
     0.042953233225881292913572018164494,-0.008674021313492586318780005883468,
    -0.059554274961645154658154180042717,-0.10530575443867740272410295104619,
    -0.14198837305251784063881548345787, -0.16644304831921567823839589426492,
    -0.17656376002524718647512421849032, -0.17147891927418672456199547790670,
    -0.15162642913722598531903229617045, -0.11871597273471929730747707847705,
    -0.075581776473850090965407023747544,-0.025938528373526445792454998016647,
     0.025938528373526445792454998016647, 0.075581776473850090965407023747544,
     0.11871597273471929730747707847705,  0.15162642913722598531903229617045,
     0.17147891927418672456199547790670,  0.17656376002524718647512421849032,
     0.16644304831921567823839589426492,  0.14198837305251784063881548345787,
     0.10530575443867740272410295104619,  0.059554274961645154658154180042717,
     0.008674021313492586318780005883468,-0.042953233225881292913572018164494,
    -0.090881384161410012831963755651079,-0.13098289131657380087121582271272,
    -0.15980423982190510363772032690233, -0.17486335449663478165921413022414,
    -0.17486335449663478165921413022414, -0.15980423982190510363772032690233,
    -0.13098289131657380087121582271272, -0.090881384161410012831963755651079,
    -0.042953233225881292913572018164494, 0.008674021313492586318780005883468,
     0.059554274961645154658154180042717, 0.10530575443867740272410295104619,
     0.14198837305251784063881548345787,  0.16644304831921567823839589426492,
     0.17656376002524718647512421849032,  0.17147891927418672456199547790670,
     0.15162642913722598531903229617045,  0.11871597273471929730747707847705,
     0.075581776473850090965407023747544, 0.025938528373526445792454998016647,
    -0.025938528373526445792454998016647,-0.075581776473850090965407023747544,
    -0.11871597273471929730747707847705, -0.15162642913722598531903229617045,
    -0.17147891927418672456199547790670, -0.17656376002524718647512421849032,
    -0.16644304831921567823839589426492, -0.14198837305251784063881548345787,
    -0.10530575443867740272410295104619, -0.059554274961645154658154180042717,
    -0.008674021313492586318780005883468, 0.042953233225881292913572018164494,
     0.090881384161410012831963755651079, 0.13098289131657380087121582271272,
     0.15980423982190510363772032690233,  0.17486335449663478165921413022414
  },
  {
     0.17417412557117862465636166724330,  0.15381110290879228162646885882803,
     0.11546573663487784907576643982592,  0.063621055212317597965816690486066,
     0.004338317276799999540542326666634,-0.055451621442613590037015233508646,
    -0.10875860742493672495725784718877, -0.14935042127442479271101666423961,
    -0.17248139814210053549033101371775, -0.17544725227590413843675978359671,
    -0.15790124012399101798155186502293, -0.12189469877166149151050389002586,
    -0.071637220287468842367069460402068,-0.013004500444088472936407629993685,
     0.047148599859865922722408214465630, 0.10178946920358000230129700949629,
     0.14452994139365530889220759879835,  0.17037314791731193194066488665533,
     0.17629771118253266395745838340776,  0.16161097942121587956084560348376,
     0.12803000590238956351843752239068,  0.079480805024268148817432153524406,
     0.021639354630747427065531392743213,-0.038731993170326232824561787558492,
    -0.094575111235207343404082277365271,-0.13936127622474077458193396572651,
    -0.16785445385626524004858887836470, -0.17672345346106673069472498415600,
    -0.16493138371356506569293701686233, -0.13385687753479470495813284521973,
    -0.087132913540067685045147230048925,-0.030222077716717299469868553789363,
     0.030222077716717299469868553789363, 0.087132913540067685045147230048925,
     0.13385687753479470495813284521973,  0.16493138371356506569293701686233,
     0.17672345346106673069472498415600,  0.16785445385626524004858887836470,
     0.13936127622474077458193396572651,  0.094575111235207343404082277365271,
     0.038731993170326232824561787558492,-0.021639354630747427065531392743213,
    -0.079480805024268148817432153524406,-0.12803000590238956351843752239068,
    -0.16161097942121587956084560348376, -0.17629771118253266395745838340776,
    -0.17037314791731193194066488665533, -0.14452994139365530889220759879835,
    -0.10178946920358000230129700949629, -0.047148599859865922722408214465630,
     0.013004500444088472936407629993685, 0.071637220287468842367069460402068,
     0.12189469877166149151050389002586,  0.15790124012399101798155186502293,
     0.17544725227590413843675978359671,  0.17248139814210053549033101371775,
     0.14935042127442479271101666423961,  0.10875860742493672495725784718877,
     0.055451621442613590037015233508646,-0.004338317276799999540542326666634,
    -0.063621055212317597965816690486066,-0.11546573663487784907576643982592,
    -0.15381110290879228162646885882803, -0.17417412557117862465636166724330
  },
  {
     0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
     0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
     0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
     0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
     0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
     0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
     0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
     0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
     0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
     0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
     0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
     0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
     0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
     0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
     0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
     0.14698445030241983962180838807641,  0.17337998066526843272770239894580
  },
  {
     0.17248139814210053549033101371775,  0.13936127622474077458193396572651,
     0.079480805024268148817432153524406, 0.004338317276799999540542326666634,
    -0.071637220287468842367069460402068,-0.13385687753479470495813284521973,
    -0.17037314791731193194066488665533, -0.17417412557117862465636166724330,
    -0.14452994139365530889220759879835, -0.087132913540067685045147230048925,
    -0.013004500444088472936407629993685, 0.063621055212317597965816690486066,
     0.12803000590238956351843752239068,  0.16785445385626524004858887836470,
     0.17544725227590413843675978359671,  0.14935042127442479271101666423961,
     0.094575111235207343404082277365271, 0.021639354630747427065531392743213,
    -0.055451621442613590037015233508646,-0.12189469877166149151050389002586,
    -0.16493138371356506569293701686233, -0.17629771118253266395745838340776,
    -0.15381110290879228162646885882803, -0.10178946920358000230129700949629,
    -0.030222077716717299469868553789363, 0.047148599859865922722408214465630,
     0.11546573663487784907576643982592,  0.16161097942121587956084560348376,
     0.17672345346106673069472498415600,  0.15790124012399101798155186502293,
     0.10875860742493672495725784718877,  0.038731993170326232824561787558492,
    -0.038731993170326232824561787558492,-0.10875860742493672495725784718877,
    -0.15790124012399101798155186502293, -0.17672345346106673069472498415600,
    -0.16161097942121587956084560348376, -0.11546573663487784907576643982592,
    -0.047148599859865922722408214465630, 0.030222077716717299469868553789363,
     0.10178946920358000230129700949629,  0.15381110290879228162646885882803,
     0.17629771118253266395745838340776,  0.16493138371356506569293701686233,
     0.12189469877166149151050389002586,  0.055451621442613590037015233508646,
    -0.021639354630747427065531392743213,-0.094575111235207343404082277365271,
    -0.14935042127442479271101666423961, -0.17544725227590413843675978359671,
    -0.16785445385626524004858887836470, -0.12803000590238956351843752239068,
    -0.063621055212317597965816690486066, 0.013004500444088472936407629993685,
     0.087132913540067685045147230048925, 0.14452994139365530889220759879835,
     0.17417412557117862465636166724330,  0.17037314791731193194066488665533,
     0.13385687753479470495813284521973,  0.071637220287468842367069460402068,
    -0.004338317276799999540542326666634,-0.079480805024268148817432153524406,
    -0.13936127622474077458193396572651, -0.17248139814210053549033101371775
  },
  {
     0.17147891927418672456199547790670,  0.13098289131657380087121582271272,
     0.059554274961645154658154180042717,-0.025938528373526445792454998016647,
    -0.10530575443867740272410295104619, -0.15980423982190510363772032690233,
    -0.17656376002524718647512421849032, -0.15162642913722598531903229617045,
    -0.090881384161410012831963755651079,-0.008674021313492586318780005883468,
     0.075581776473850090965407023747544, 0.14198837305251784063881548345787,
     0.17486335449663478165921413022414,  0.16644304831921567823839589426492,
     0.11871597273471929730747707847705,  0.042953233225881292913572018164494,
    -0.042953233225881292913572018164494,-0.11871597273471929730747707847705,
    -0.16644304831921567823839589426492, -0.17486335449663478165921413022414,
    -0.14198837305251784063881548345787, -0.075581776473850090965407023747544,
     0.008674021313492586318780005883468, 0.090881384161410012831963755651079,
     0.15162642913722598531903229617045,  0.17656376002524718647512421849032,
     0.15980423982190510363772032690233,  0.10530575443867740272410295104619,
     0.025938528373526445792454998016647,-0.059554274961645154658154180042717,
    -0.13098289131657380087121582271272, -0.17147891927418672456199547790670,
    -0.17147891927418672456199547790670, -0.13098289131657380087121582271272,
    -0.059554274961645154658154180042717, 0.025938528373526445792454998016647,
     0.10530575443867740272410295104619,  0.15980423982190510363772032690233,
     0.17656376002524718647512421849032,  0.15162642913722598531903229617045,
     0.090881384161410012831963755651079, 0.008674021313492586318780005883468,
    -0.075581776473850090965407023747544,-0.14198837305251784063881548345787,
    -0.17486335449663478165921413022414, -0.16644304831921567823839589426492,
    -0.11871597273471929730747707847705, -0.042953233225881292913572018164494,
     0.042953233225881292913572018164494, 0.11871597273471929730747707847705,
     0.16644304831921567823839589426492,  0.17486335449663478165921413022414,
     0.14198837305251784063881548345787,  0.075581776473850090965407023747544,
    -0.008674021313492586318780005883468,-0.090881384161410012831963755651079,
    -0.15162642913722598531903229617045, -0.17656376002524718647512421849032,
    -0.15980423982190510363772032690233, -0.10530575443867740272410295104619,
    -0.025938528373526445792454998016647, 0.059554274961645154658154180042717,
     0.13098289131657380087121582271272,  0.17147891927418672456199547790670
  },
  {
     0.17037314791731193194066488665533,  0.12189469877166149151050389002586,
     0.038731993170326232824561787558492,-0.055451621442613590037015233508646,
    -0.13385687753479470495813284521973, -0.17417412557117862465636166724330,
    -0.16493138371356506569293701686233, -0.10875860742493672495725784718877,
    -0.021639354630747427065531392743213, 0.071637220287468842367069460402068,
     0.14452994139365530889220759879835,  0.17629771118253266395745838340776,
     0.15790124012399101798155186502293,  0.094575111235207343404082277365271,
     0.004338317276799999540542326666634,-0.087132913540067685045147230048925,
    -0.15381110290879228162646885882803, -0.17672345346106673069472498415600,
    -0.14935042127442479271101666423961, -0.079480805024268148817432153524406,
     0.013004500444088472936407629993685, 0.10178946920358000230129700949629,
     0.16161097942121587956084560348376,  0.17544725227590413843675978359671,
     0.13936127622474077458193396572651,  0.063621055212317597965816690486066,
    -0.030222077716717299469868553789363,-0.11546573663487784907576643982592,
    -0.16785445385626524004858887836470, -0.17248139814210053549033101371775,
    -0.12803000590238956351843752239068, -0.047148599859865922722408214465630,
     0.047148599859865922722408214465630, 0.12803000590238956351843752239068,
     0.17248139814210053549033101371775,  0.16785445385626524004858887836470,
     0.11546573663487784907576643982592,  0.030222077716717299469868553789363,
    -0.063621055212317597965816690486066,-0.13936127622474077458193396572651,
    -0.17544725227590413843675978359671, -0.16161097942121587956084560348376,
    -0.10178946920358000230129700949629, -0.013004500444088472936407629993685,
     0.079480805024268148817432153524406, 0.14935042127442479271101666423961,
     0.17672345346106673069472498415600,  0.15381110290879228162646885882803,
     0.087132913540067685045147230048925,-0.004338317276799999540542326666634,
    -0.094575111235207343404082277365271,-0.15790124012399101798155186502293,
    -0.17629771118253266395745838340776, -0.14452994139365530889220759879835,
    -0.071637220287468842367069460402068, 0.021639354630747427065531392743213,
     0.10875860742493672495725784718877,  0.16493138371356506569293701686233,
     0.17417412557117862465636166724330,  0.13385687753479470495813284521973,
     0.055451621442613590037015233508646,-0.038731993170326232824561787558492,
    -0.12189469877166149151050389002586, -0.17037314791731193194066488665533
  },
  {
     0.16916475014679408478364306119571,  0.11214594829282953553133017365260,
     0.017327146149886432824466874566622,-0.083331957309718312162450688895359,
    -0.15590312662333390407149878284971, -0.17592546719079780737825977787300,
    -0.13665023337521968602987906462477, -0.051315565940294672644546154719392,
     0.051315565940294672644546154719392, 0.13665023337521968602987906462477,
     0.17592546719079780737825977787300,  0.15590312662333390407149878284971,
     0.083331957309718312162450688895359,-0.017327146149886432824466874566622,
    -0.11214594829282953553133017365260, -0.16916475014679408478364306119571,
    -0.16916475014679408478364306119571, -0.11214594829282953553133017365260,
    -0.017327146149886432824466874566622, 0.083331957309718312162450688895359,
     0.15590312662333390407149878284971,  0.17592546719079780737825977787300,
     0.13665023337521968602987906462477,  0.051315565940294672644546154719392,
    -0.051315565940294672644546154719392,-0.13665023337521968602987906462477,
    -0.17592546719079780737825977787300, -0.15590312662333390407149878284971,
    -0.083331957309718312162450688895359, 0.017327146149886432824466874566622,
     0.11214594829282953553133017365260,  0.16916475014679408478364306119571,
     0.16916475014679408478364306119571,  0.11214594829282953553133017365260,
     0.017327146149886432824466874566622,-0.083331957309718312162450688895359,
    -0.15590312662333390407149878284971, -0.17592546719079780737825977787300,
    -0.13665023337521968602987906462477, -0.051315565940294672644546154719392,
     0.051315565940294672644546154719392, 0.13665023337521968602987906462477,
     0.17592546719079780737825977787300,  0.15590312662333390407149878284971,
     0.083331957309718312162450688895359,-0.017327146149886432824466874566622,
    -0.11214594829282953553133017365260, -0.16916475014679408478364306119571,
    -0.16916475014679408478364306119571, -0.11214594829282953553133017365260,
    -0.017327146149886432824466874566622, 0.083331957309718312162450688895359,
     0.15590312662333390407149878284971,  0.17592546719079780737825977787300,
     0.13665023337521968602987906462477,  0.051315565940294672644546154719392,
    -0.051315565940294672644546154719392,-0.13665023337521968602987906462477,
    -0.17592546719079780737825977787300, -0.15590312662333390407149878284971,
    -0.083331957309718312162450688895359, 0.017327146149886432824466874566622,
     0.11214594829282953553133017365260,  0.16916475014679408478364306119571
  },
  {
     0.16785445385626524004858887836470,  0.10178946920358000230129700949629,
    -0.004338317276799999540542326666634,-0.10875860742493672495725784718877,
    -0.17037314791731193194066488665533, -0.16493138371356506569293701686233,
    -0.094575111235207343404082277365271, 0.013004500444088472936407629993685,
     0.11546573663487784907576643982592,  0.17248139814210053549033101371775,
     0.16161097942121587956084560348376,  0.087132913540067685045147230048925,
    -0.021639354630747427065531392743213,-0.12189469877166149151050389002586,
    -0.17417412557117862465636166724330, -0.15790124012399101798155186502293,
    -0.079480805024268148817432153524406, 0.030222077716717299469868553789363,
     0.12803000590238956351843752239068,  0.17544725227590413843675978359671,
     0.15381110290879228162646885882803,  0.071637220287468842367069460402068,
    -0.038731993170326232824561787558492,-0.13385687753479470495813284521973,
    -0.17629771118253266395745838340776, -0.14935042127442479271101666423961,
    -0.063621055212317597965816690486066, 0.047148599859865922722408214465630,
     0.13936127622474077458193396572651,  0.17672345346106673069472498415600,
     0.14452994139365530889220759879835,  0.055451621442613590037015233508646,
    -0.055451621442613590037015233508646,-0.14452994139365530889220759879835,
    -0.17672345346106673069472498415600, -0.13936127622474077458193396572651,
    -0.047148599859865922722408214465630, 0.063621055212317597965816690486066,
     0.14935042127442479271101666423961,  0.17629771118253266395745838340776,
     0.13385687753479470495813284521973,  0.038731993170326232824561787558492,
    -0.071637220287468842367069460402068,-0.15381110290879228162646885882803,
    -0.17544725227590413843675978359671, -0.12803000590238956351843752239068,
    -0.030222077716717299469868553789363, 0.079480805024268148817432153524406,
     0.15790124012399101798155186502293,  0.17417412557117862465636166724330,
     0.12189469877166149151050389002586,  0.021639354630747427065531392743213,
    -0.087132913540067685045147230048925,-0.16161097942121587956084560348376,
    -0.17248139814210053549033101371775, -0.11546573663487784907576643982592,
    -0.013004500444088472936407629993685, 0.094575111235207343404082277365271,
     0.16493138371356506569293701686233,  0.17037314791731193194066488665533,
     0.10875860742493672495725784718877,  0.004338317276799999540542326666634,
    -0.10178946920358000230129700949629, -0.16785445385626524004858887836470
  },
  {
     0.16644304831921567823839589426492,  0.090881384161410012831963755651079,
    -0.025938528373526445792454998016647,-0.13098289131657380087121582271272,
    -0.17656376002524718647512421849032, -0.14198837305251784063881548345787,
    -0.042953233225881292913572018164494, 0.075581776473850090965407023747544,
     0.15980423982190510363772032690233,  0.17147891927418672456199547790670,
     0.10530575443867740272410295104619, -0.008674021313492586318780005883468,
    -0.11871597273471929730747707847705, -0.17486335449663478165921413022414,
    -0.15162642913722598531903229617045, -0.059554274961645154658154180042717,
     0.059554274961645154658154180042717, 0.15162642913722598531903229617045,
     0.17486335449663478165921413022414,  0.11871597273471929730747707847705,
     0.008674021313492586318780005883468,-0.10530575443867740272410295104619,
    -0.17147891927418672456199547790670, -0.15980423982190510363772032690233,
    -0.075581776473850090965407023747544, 0.042953233225881292913572018164494,
     0.14198837305251784063881548345787,  0.17656376002524718647512421849032,
     0.13098289131657380087121582271272,  0.025938528373526445792454998016647,
    -0.090881384161410012831963755651079,-0.16644304831921567823839589426492,
    -0.16644304831921567823839589426492, -0.090881384161410012831963755651079,
     0.025938528373526445792454998016647, 0.13098289131657380087121582271272,
     0.17656376002524718647512421849032,  0.14198837305251784063881548345787,
     0.042953233225881292913572018164494,-0.075581776473850090965407023747544,
    -0.15980423982190510363772032690233, -0.17147891927418672456199547790670,
    -0.10530575443867740272410295104619,  0.008674021313492586318780005883468,
     0.11871597273471929730747707847705,  0.17486335449663478165921413022414,
     0.15162642913722598531903229617045,  0.059554274961645154658154180042717,
    -0.059554274961645154658154180042717,-0.15162642913722598531903229617045,
    -0.17486335449663478165921413022414, -0.11871597273471929730747707847705,
    -0.008674021313492586318780005883468, 0.10530575443867740272410295104619,
     0.17147891927418672456199547790670,  0.15980423982190510363772032690233,
     0.075581776473850090965407023747544,-0.042953233225881292913572018164494,
    -0.14198837305251784063881548345787, -0.17656376002524718647512421849032,
    -0.13098289131657380087121582271272, -0.025938528373526445792454998016647,
     0.090881384161410012831963755651079, 0.16644304831921567823839589426492
  },
  {
     0.16493138371356506569293701686233,  0.079480805024268148817432153524406,
    -0.047148599859865922722408214465630,-0.14935042127442479271101666423961,
    -0.17417412557117862465636166724330, -0.10875860742493672495725784718877,
     0.013004500444088472936407629993685, 0.12803000590238956351843752239068,
     0.17672345346106673069472498415600,  0.13385687753479470495813284521973,
     0.021639354630747427065531392743213,-0.10178946920358000230129700949629,
    -0.17248139814210053549033101371775, -0.15381110290879228162646885882803,
    -0.055451621442613590037015233508646, 0.071637220287468842367069460402068,
     0.16161097942121587956084560348376,  0.16785445385626524004858887836470,
     0.087132913540067685045147230048925,-0.038731993170326232824561787558492,
    -0.14452994139365530889220759879835, -0.17544725227590413843675978359671,
    -0.11546573663487784907576643982592,  0.004338317276799999540542326666634,
     0.12189469877166149151050389002586,  0.17629771118253266395745838340776,
     0.13936127622474077458193396572651,  0.030222077716717299469868553789363,
    -0.094575111235207343404082277365271,-0.17037314791731193194066488665533,
    -0.15790124012399101798155186502293, -0.063621055212317597965816690486066,
     0.063621055212317597965816690486066, 0.15790124012399101798155186502293,
     0.17037314791731193194066488665533,  0.094575111235207343404082277365271,
    -0.030222077716717299469868553789363,-0.13936127622474077458193396572651,
    -0.17629771118253266395745838340776, -0.12189469877166149151050389002586,
    -0.004338317276799999540542326666634, 0.11546573663487784907576643982592,
     0.17544725227590413843675978359671,  0.14452994139365530889220759879835,
     0.038731993170326232824561787558492,-0.087132913540067685045147230048925,
    -0.16785445385626524004858887836470, -0.16161097942121587956084560348376,
    -0.071637220287468842367069460402068, 0.055451621442613590037015233508646,
     0.15381110290879228162646885882803,  0.17248139814210053549033101371775,
     0.10178946920358000230129700949629, -0.021639354630747427065531392743213,
    -0.13385687753479470495813284521973, -0.17672345346106673069472498415600,
    -0.12803000590238956351843752239068, -0.013004500444088472936407629993685,
     0.10875860742493672495725784718877,  0.17417412557117862465636166724330,
     0.14935042127442479271101666423961,  0.047148599859865922722408214465630,
    -0.079480805024268148817432153524406,-0.16493138371356506569293701686233
  },
  {
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
     0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
     0.067649512518274623049965400670799, 0.16332037060954706598208039667840
  },
  {
     0.16161097942121587956084560348376,  0.055451621442613590037015233508646,
    -0.087132913540067685045147230048925,-0.17248139814210053549033101371775,
    -0.14452994139365530889220759879835, -0.021639354630747427065531392743213,
     0.11546573663487784907576643982592,  0.17672345346106673069472498415600,
     0.12189469877166149151050389002586, -0.013004500444088472936407629993685,
    -0.13936127622474077458193396572651, -0.17417412557117862465636166724330,
    -0.094575111235207343404082277365271, 0.047148599859865922722408214465630,
     0.15790124012399101798155186502293,  0.16493138371356506569293701686233,
     0.063621055212317597965816690486066,-0.079480805024268148817432153524406,
    -0.17037314791731193194066488665533, -0.14935042127442479271101666423961,
    -0.030222077716717299469868553789363, 0.10875860742493672495725784718877,
     0.17629771118253266395745838340776,  0.12803000590238956351843752239068,
    -0.004338317276799999540542326666634,-0.13385687753479470495813284521973,
    -0.17544725227590413843675978359671, -0.10178946920358000230129700949629,
     0.038731993170326232824561787558492, 0.15381110290879228162646885882803,
     0.16785445385626524004858887836470,  0.071637220287468842367069460402068,
    -0.071637220287468842367069460402068,-0.16785445385626524004858887836470,
    -0.15381110290879228162646885882803, -0.038731993170326232824561787558492,
     0.10178946920358000230129700949629,  0.17544725227590413843675978359671,
     0.13385687753479470495813284521973,  0.004338317276799999540542326666634,
    -0.12803000590238956351843752239068, -0.17629771118253266395745838340776,
    -0.10875860742493672495725784718877,  0.030222077716717299469868553789363,
     0.14935042127442479271101666423961,  0.17037314791731193194066488665533,
     0.079480805024268148817432153524406,-0.063621055212317597965816690486066,
    -0.16493138371356506569293701686233, -0.15790124012399101798155186502293,
    -0.047148599859865922722408214465630, 0.094575111235207343404082277365271,
     0.17417412557117862465636166724330,  0.13936127622474077458193396572651,
     0.013004500444088472936407629993685,-0.12189469877166149151050389002586,
    -0.17672345346106673069472498415600, -0.11546573663487784907576643982592,
     0.021639354630747427065531392743213, 0.14452994139365530889220759879835,
     0.17248139814210053549033101371775,  0.087132913540067685045147230048925,
    -0.055451621442613590037015233508646,-0.16161097942121587956084560348376
  },
  {
     0.15980423982190510363772032690233,  0.042953233225881292913572018164494,
    -0.10530575443867740272410295104619, -0.17656376002524718647512421849032,
    -0.11871597273471929730747707847705,  0.025938528373526445792454998016647,
     0.15162642913722598531903229617045,  0.16644304831921567823839589426492,
     0.059554274961645154658154180042717,-0.090881384161410012831963755651079,
    -0.17486335449663478165921413022414, -0.13098289131657380087121582271272,
     0.008674021313492586318780005883468, 0.14198837305251784063881548345787,
     0.17147891927418672456199547790670,  0.075581776473850090965407023747544,
    -0.075581776473850090965407023747544,-0.17147891927418672456199547790670,
    -0.14198837305251784063881548345787, -0.008674021313492586318780005883468,
     0.13098289131657380087121582271272,  0.17486335449663478165921413022414,
     0.090881384161410012831963755651079,-0.059554274961645154658154180042717,
    -0.16644304831921567823839589426492, -0.15162642913722598531903229617045,
    -0.025938528373526445792454998016647, 0.11871597273471929730747707847705,
     0.17656376002524718647512421849032,  0.10530575443867740272410295104619,
    -0.042953233225881292913572018164494,-0.15980423982190510363772032690233,
    -0.15980423982190510363772032690233, -0.042953233225881292913572018164494,
     0.10530575443867740272410295104619,  0.17656376002524718647512421849032,
     0.11871597273471929730747707847705, -0.025938528373526445792454998016647,
    -0.15162642913722598531903229617045, -0.16644304831921567823839589426492,
    -0.059554274961645154658154180042717, 0.090881384161410012831963755651079,
     0.17486335449663478165921413022414,  0.13098289131657380087121582271272,
    -0.008674021313492586318780005883468,-0.14198837305251784063881548345787,
    -0.17147891927418672456199547790670, -0.075581776473850090965407023747544,
     0.075581776473850090965407023747544, 0.17147891927418672456199547790670,
     0.14198837305251784063881548345787,  0.008674021313492586318780005883468,
    -0.13098289131657380087121582271272, -0.17486335449663478165921413022414,
    -0.090881384161410012831963755651079, 0.059554274961645154658154180042717,
     0.16644304831921567823839589426492,  0.15162642913722598531903229617045,
     0.025938528373526445792454998016647,-0.11871597273471929730747707847705,
    -0.17656376002524718647512421849032, -0.10530575443867740272410295104619,
     0.042953233225881292913572018164494, 0.15980423982190510363772032690233
  },
  {
     0.15790124012399101798155186502293,  0.030222077716717299469868553789363,
    -0.12189469877166149151050389002586, -0.17544725227590413843675978359671,
    -0.087132913540067685045147230048925, 0.071637220287468842367069460402068,
     0.17248139814210053549033101371775,  0.13385687753479470495813284521973,
    -0.013004500444088472936407629993685,-0.14935042127442479271101666423961,
    -0.16493138371356506569293701686233, -0.047148599859865922722408214465630,
     0.10875860742493672495725784718877,  0.17672345346106673069472498415600,
     0.10178946920358000230129700949629, -0.055451621442613590037015233508646,
    -0.16785445385626524004858887836470, -0.14452994139365530889220759879835,
    -0.004338317276799999540542326666634, 0.13936127622474077458193396572651,
     0.17037314791731193194066488665533,  0.063621055212317597965816690486066,
    -0.094575111235207343404082277365271,-0.17629771118253266395745838340776,
    -0.11546573663487784907576643982592,  0.038731993170326232824561787558492,
     0.16161097942121587956084560348376,  0.15381110290879228162646885882803,
     0.021639354630747427065531392743213,-0.12803000590238956351843752239068,
    -0.17417412557117862465636166724330, -0.079480805024268148817432153524406,
     0.079480805024268148817432153524406, 0.17417412557117862465636166724330,
     0.12803000590238956351843752239068, -0.021639354630747427065531392743213,
    -0.15381110290879228162646885882803, -0.16161097942121587956084560348376,
    -0.038731993170326232824561787558492, 0.11546573663487784907576643982592,
     0.17629771118253266395745838340776,  0.094575111235207343404082277365271,
    -0.063621055212317597965816690486066,-0.17037314791731193194066488665533,
    -0.13936127622474077458193396572651,  0.004338317276799999540542326666634,
     0.14452994139365530889220759879835,  0.16785445385626524004858887836470,
     0.055451621442613590037015233508646,-0.10178946920358000230129700949629,
    -0.17672345346106673069472498415600, -0.10875860742493672495725784718877,
     0.047148599859865922722408214465630, 0.16493138371356506569293701686233,
     0.14935042127442479271101666423961,  0.013004500444088472936407629993685,
    -0.13385687753479470495813284521973, -0.17248139814210053549033101371775,
    -0.071637220287468842367069460402068, 0.087132913540067685045147230048925,
     0.17544725227590413843675978359671,  0.12189469877166149151050389002586,
    -0.030222077716717299469868553789363,-0.15790124012399101798155186502293
  },
  {
     0.15590312662333390407149878284971,  0.017327146149886432824466874566622,
    -0.13665023337521968602987906462477, -0.16916475014679408478364306119571,
    -0.051315565940294672644546154719392, 0.11214594829282953553133017365260,
     0.17592546719079780737825977787300,  0.083331957309718312162450688895359,
    -0.083331957309718312162450688895359,-0.17592546719079780737825977787300,
    -0.11214594829282953553133017365260,  0.051315565940294672644546154719392,
     0.16916475014679408478364306119571,  0.13665023337521968602987906462477,
    -0.017327146149886432824466874566622,-0.15590312662333390407149878284971,
    -0.15590312662333390407149878284971, -0.017327146149886432824466874566622,
     0.13665023337521968602987906462477,  0.16916475014679408478364306119571,
     0.051315565940294672644546154719392,-0.11214594829282953553133017365260,
    -0.17592546719079780737825977787300, -0.083331957309718312162450688895359,
     0.083331957309718312162450688895359, 0.17592546719079780737825977787300,
     0.11214594829282953553133017365260, -0.051315565940294672644546154719392,
    -0.16916475014679408478364306119571, -0.13665023337521968602987906462477,
     0.017327146149886432824466874566622, 0.15590312662333390407149878284971,
     0.15590312662333390407149878284971,  0.017327146149886432824466874566622,
    -0.13665023337521968602987906462477, -0.16916475014679408478364306119571,
    -0.051315565940294672644546154719392, 0.11214594829282953553133017365260,
     0.17592546719079780737825977787300,  0.083331957309718312162450688895359,
    -0.083331957309718312162450688895359,-0.17592546719079780737825977787300,
    -0.11214594829282953553133017365260,  0.051315565940294672644546154719392,
     0.16916475014679408478364306119571,  0.13665023337521968602987906462477,
    -0.017327146149886432824466874566622,-0.15590312662333390407149878284971,
    -0.15590312662333390407149878284971, -0.017327146149886432824466874566622,
     0.13665023337521968602987906462477,  0.16916475014679408478364306119571,
     0.051315565940294672644546154719392,-0.11214594829282953553133017365260,
    -0.17592546719079780737825977787300, -0.083331957309718312162450688895359,
     0.083331957309718312162450688895359, 0.17592546719079780737825977787300,
     0.11214594829282953553133017365260, -0.051315565940294672644546154719392,
    -0.16916475014679408478364306119571, -0.13665023337521968602987906462477,
     0.017327146149886432824466874566622, 0.15590312662333390407149878284971
  },
  {
     0.15381110290879228162646885882803,  0.004338317276799999540542326666634,
    -0.14935042127442479271101666423961, -0.15790124012399101798155186502293,
    -0.013004500444088472936407629993685, 0.14452994139365530889220759879835,
     0.16161097942121587956084560348376,  0.021639354630747427065531392743213,
    -0.13936127622474077458193396572651, -0.16493138371356506569293701686233,
    -0.030222077716717299469868553789363, 0.13385687753479470495813284521973,
     0.16785445385626524004858887836470,  0.038731993170326232824561787558492,
    -0.12803000590238956351843752239068, -0.17037314791731193194066488665533,
    -0.047148599859865922722408214465630, 0.12189469877166149151050389002586,
     0.17248139814210053549033101371775,  0.055451621442613590037015233508646,
    -0.11546573663487784907576643982592, -0.17417412557117862465636166724330,
    -0.063621055212317597965816690486066, 0.10875860742493672495725784718877,
     0.17544725227590413843675978359671,  0.071637220287468842367069460402068,
    -0.10178946920358000230129700949629, -0.17629771118253266395745838340776,
    -0.079480805024268148817432153524406, 0.094575111235207343404082277365271,
     0.17672345346106673069472498415600,  0.087132913540067685045147230048925,
    -0.087132913540067685045147230048925,-0.17672345346106673069472498415600,
    -0.094575111235207343404082277365271, 0.079480805024268148817432153524406,
     0.17629771118253266395745838340776,  0.10178946920358000230129700949629,
    -0.071637220287468842367069460402068,-0.17544725227590413843675978359671,
    -0.10875860742493672495725784718877,  0.063621055212317597965816690486066,
     0.17417412557117862465636166724330,  0.11546573663487784907576643982592,
    -0.055451621442613590037015233508646,-0.17248139814210053549033101371775,
    -0.12189469877166149151050389002586,  0.047148599859865922722408214465630,
     0.17037314791731193194066488665533,  0.12803000590238956351843752239068,
    -0.038731993170326232824561787558492,-0.16785445385626524004858887836470,
    -0.13385687753479470495813284521973,  0.030222077716717299469868553789363,
     0.16493138371356506569293701686233,  0.13936127622474077458193396572651,
    -0.021639354630747427065531392743213,-0.16161097942121587956084560348376,
    -0.14452994139365530889220759879835,  0.013004500444088472936407629993685,
     0.15790124012399101798155186502293,  0.14935042127442479271101666423961,
    -0.004338317276799999540542326666634,-0.15381110290879228162646885882803
  },
  {
     0.15162642913722598531903229617045, -0.008674021313492586318780005883468,
    -0.15980423982190510363772032690233, -0.14198837305251784063881548345787,
     0.025938528373526445792454998016647, 0.16644304831921567823839589426492,
     0.13098289131657380087121582271272, -0.042953233225881292913572018164494,
    -0.17147891927418672456199547790670, -0.11871597273471929730747707847705,
     0.059554274961645154658154180042717, 0.17486335449663478165921413022414,
     0.10530575443867740272410295104619, -0.075581776473850090965407023747544,
    -0.17656376002524718647512421849032, -0.090881384161410012831963755651079,
     0.090881384161410012831963755651079, 0.17656376002524718647512421849032,
     0.075581776473850090965407023747544,-0.10530575443867740272410295104619,
    -0.17486335449663478165921413022414, -0.059554274961645154658154180042717,
     0.11871597273471929730747707847705,  0.17147891927418672456199547790670,
     0.042953233225881292913572018164494,-0.13098289131657380087121582271272,
    -0.16644304831921567823839589426492, -0.025938528373526445792454998016647,
     0.14198837305251784063881548345787,  0.15980423982190510363772032690233,
     0.008674021313492586318780005883468,-0.15162642913722598531903229617045,
    -0.15162642913722598531903229617045,  0.008674021313492586318780005883468,
     0.15980423982190510363772032690233,  0.14198837305251784063881548345787,
    -0.025938528373526445792454998016647,-0.16644304831921567823839589426492,
    -0.13098289131657380087121582271272,  0.042953233225881292913572018164494,
     0.17147891927418672456199547790670,  0.11871597273471929730747707847705,
    -0.059554274961645154658154180042717,-0.17486335449663478165921413022414,
    -0.10530575443867740272410295104619,  0.075581776473850090965407023747544,
     0.17656376002524718647512421849032,  0.090881384161410012831963755651079,
    -0.090881384161410012831963755651079,-0.17656376002524718647512421849032,
    -0.075581776473850090965407023747544, 0.10530575443867740272410295104619,
     0.17486335449663478165921413022414,  0.059554274961645154658154180042717,
    -0.11871597273471929730747707847705, -0.17147891927418672456199547790670,
    -0.042953233225881292913572018164494, 0.13098289131657380087121582271272,
     0.16644304831921567823839589426492,  0.025938528373526445792454998016647,
    -0.14198837305251784063881548345787, -0.15980423982190510363772032690233,
    -0.008674021313492586318780005883468, 0.15162642913722598531903229617045
  },
  {
     0.14935042127442479271101666423961, -0.021639354630747427065531392743213,
    -0.16785445385626524004858887836470, -0.12189469877166149151050389002586,
     0.063621055212317597965816690486066, 0.17629771118253266395745838340776,
     0.087132913540067685045147230048925,-0.10178946920358000230129700949629,
    -0.17417412557117862465636166724330, -0.047148599859865922722408214465630,
     0.13385687753479470495813284521973,  0.16161097942121587956084560348376,
     0.004338317276799999540542326666634,-0.15790124012399101798155186502293,
    -0.13936127622474077458193396572651,  0.038731993170326232824561787558492,
     0.17248139814210053549033101371775,  0.10875860742493672495725784718877,
    -0.079480805024268148817432153524406,-0.17672345346106673069472498415600,
    -0.071637220287468842367069460402068, 0.11546573663487784907576643982592,
     0.17037314791731193194066488665533,  0.030222077716717299469868553789363,
    -0.14452994139365530889220759879835, -0.15381110290879228162646885882803,
     0.013004500444088472936407629993685, 0.16493138371356506569293701686233,
     0.12803000590238956351843752239068, -0.055451621442613590037015233508646,
    -0.17544725227590413843675978359671, -0.094575111235207343404082277365271,
     0.094575111235207343404082277365271, 0.17544725227590413843675978359671,
     0.055451621442613590037015233508646,-0.12803000590238956351843752239068,
    -0.16493138371356506569293701686233, -0.013004500444088472936407629993685,
     0.15381110290879228162646885882803,  0.14452994139365530889220759879835,
    -0.030222077716717299469868553789363,-0.17037314791731193194066488665533,
    -0.11546573663487784907576643982592,  0.071637220287468842367069460402068,
     0.17672345346106673069472498415600,  0.079480805024268148817432153524406,
    -0.10875860742493672495725784718877, -0.17248139814210053549033101371775,
    -0.038731993170326232824561787558492, 0.13936127622474077458193396572651,
     0.15790124012399101798155186502293, -0.004338317276799999540542326666634,
    -0.16161097942121587956084560348376, -0.13385687753479470495813284521973,
     0.047148599859865922722408214465630, 0.17417412557117862465636166724330,
     0.10178946920358000230129700949629, -0.087132913540067685045147230048925,
    -0.17629771118253266395745838340776, -0.063621055212317597965816690486066,
     0.12189469877166149151050389002586,  0.16785445385626524004858887836470,
     0.021639354630747427065531392743213,-0.14935042127442479271101666423961
  },
  {
     0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
     0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
     0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
     0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
     0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
     0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
     0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
     0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
     0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
     0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
     0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
     0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
     0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
     0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
     0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
     0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672, 0.14698445030241983962180838807641
  },
  {
     0.14452994139365530889220759879835, -0.047148599859865922722408214465630,
    -0.17629771118253266395745838340776, -0.071637220287468842367069460402068,
     0.12803000590238956351843752239068,  0.15790124012399101798155186502293,
    -0.021639354630747427065531392743213,-0.17248139814210053549033101371775,
    -0.094575111235207343404082277365271, 0.10875860742493672495725784718877,
     0.16785445385626524004858887836470,  0.004338317276799999540542326666634,
    -0.16493138371356506569293701686233, -0.11546573663487784907576643982592,
     0.087132913540067685045147230048925, 0.17417412557117862465636166724330,
     0.030222077716717299469868553789363,-0.15381110290879228162646885882803,
    -0.13385687753479470495813284521973,  0.063621055212317597965816690486066,
     0.17672345346106673069472498415600,  0.055451621442613590037015233508646,
    -0.13936127622474077458193396572651, -0.14935042127442479271101666423961,
     0.038731993170326232824561787558492, 0.17544725227590413843675978359671,
     0.079480805024268148817432153524406,-0.12189469877166149151050389002586,
    -0.16161097942121587956084560348376,  0.013004500444088472936407629993685,
     0.17037314791731193194066488665533,  0.10178946920358000230129700949629,
    -0.10178946920358000230129700949629, -0.17037314791731193194066488665533,
    -0.013004500444088472936407629993685, 0.16161097942121587956084560348376,
     0.12189469877166149151050389002586, -0.079480805024268148817432153524406,
    -0.17544725227590413843675978359671, -0.038731993170326232824561787558492,
     0.14935042127442479271101666423961,  0.13936127622474077458193396572651,
    -0.055451621442613590037015233508646,-0.17672345346106673069472498415600,
    -0.063621055212317597965816690486066, 0.13385687753479470495813284521973,
     0.15381110290879228162646885882803, -0.030222077716717299469868553789363,
    -0.17417412557117862465636166724330, -0.087132913540067685045147230048925,
     0.11546573663487784907576643982592,  0.16493138371356506569293701686233,
    -0.004338317276799999540542326666634,-0.16785445385626524004858887836470,
    -0.10875860742493672495725784718877,  0.094575111235207343404082277365271,
     0.17248139814210053549033101371775,  0.021639354630747427065531392743213,
    -0.15790124012399101798155186502293, -0.12803000590238956351843752239068,
     0.071637220287468842367069460402068, 0.17629771118253266395745838340776,
     0.047148599859865922722408214465630,-0.14452994139365530889220759879835
  },
  {
     0.14198837305251784063881548345787, -0.059554274961645154658154180042717,
    -0.17656376002524718647512421849032, -0.042953233225881292913572018164494,
     0.15162642913722598531903229617045,  0.13098289131657380087121582271272,
    -0.075581776473850090965407023747544,-0.17486335449663478165921413022414,
    -0.025938528373526445792454998016647, 0.15980423982190510363772032690233,
     0.11871597273471929730747707847705, -0.090881384161410012831963755651079,
    -0.17147891927418672456199547790670, -0.008674021313492586318780005883468,
     0.16644304831921567823839589426492,  0.10530575443867740272410295104619,
    -0.10530575443867740272410295104619, -0.16644304831921567823839589426492,
     0.008674021313492586318780005883468, 0.17147891927418672456199547790670,
     0.090881384161410012831963755651079,-0.11871597273471929730747707847705,
    -0.15980423982190510363772032690233,  0.025938528373526445792454998016647,
     0.17486335449663478165921413022414,  0.075581776473850090965407023747544,
    -0.13098289131657380087121582271272, -0.15162642913722598531903229617045,
     0.042953233225881292913572018164494, 0.17656376002524718647512421849032,
     0.059554274961645154658154180042717,-0.14198837305251784063881548345787,
    -0.14198837305251784063881548345787,  0.059554274961645154658154180042717,
     0.17656376002524718647512421849032,  0.042953233225881292913572018164494,
    -0.15162642913722598531903229617045, -0.13098289131657380087121582271272,
     0.075581776473850090965407023747544, 0.17486335449663478165921413022414,
     0.025938528373526445792454998016647,-0.15980423982190510363772032690233,
    -0.11871597273471929730747707847705,  0.090881384161410012831963755651079,
     0.17147891927418672456199547790670,  0.008674021313492586318780005883468,
    -0.16644304831921567823839589426492, -0.10530575443867740272410295104619,
     0.10530575443867740272410295104619,  0.16644304831921567823839589426492,
    -0.008674021313492586318780005883468,-0.17147891927418672456199547790670,
    -0.090881384161410012831963755651079, 0.11871597273471929730747707847705,
     0.15980423982190510363772032690233, -0.025938528373526445792454998016647,
    -0.17486335449663478165921413022414, -0.075581776473850090965407023747544,
     0.13098289131657380087121582271272,  0.15162642913722598531903229617045,
    -0.042953233225881292913572018164494,-0.17656376002524718647512421849032,
    -0.059554274961645154658154180042717, 0.14198837305251784063881548345787
  },
  {
     0.13936127622474077458193396572651, -0.071637220287468842367069460402068,
    -0.17417412557117862465636166724330, -0.013004500444088472936407629993685,
     0.16785445385626524004858887836470,  0.094575111235207343404082277365271,
    -0.12189469877166149151050389002586, -0.15381110290879228162646885882803,
     0.047148599859865922722408214465630, 0.17672345346106673069472498415600,
     0.038731993170326232824561787558492,-0.15790124012399101798155186502293,
    -0.11546573663487784907576643982592,  0.10178946920358000230129700949629,
     0.16493138371356506569293701686233, -0.021639354630747427065531392743213,
    -0.17544725227590413843675978359671, -0.063621055212317597965816690486066,
     0.14452994139365530889220759879835,  0.13385687753479470495813284521973,
    -0.079480805024268148817432153524406,-0.17248139814210053549033101371775,
    -0.004338317276799999540542326666634, 0.17037314791731193194066488665533,
     0.087132913540067685045147230048925,-0.12803000590238956351843752239068,
    -0.14935042127442479271101666423961,  0.055451621442613590037015233508646,
     0.17629771118253266395745838340776,  0.030222077716717299469868553789363,
    -0.16161097942121587956084560348376, -0.10875860742493672495725784718877,
     0.10875860742493672495725784718877,  0.16161097942121587956084560348376,
    -0.030222077716717299469868553789363,-0.17629771118253266395745838340776,
    -0.055451621442613590037015233508646, 0.14935042127442479271101666423961,
     0.12803000590238956351843752239068, -0.087132913540067685045147230048925,
    -0.17037314791731193194066488665533,  0.004338317276799999540542326666634,
     0.17248139814210053549033101371775,  0.079480805024268148817432153524406,
    -0.13385687753479470495813284521973, -0.14452994139365530889220759879835,
     0.063621055212317597965816690486066, 0.17544725227590413843675978359671,
     0.021639354630747427065531392743213,-0.16493138371356506569293701686233,
    -0.10178946920358000230129700949629,  0.11546573663487784907576643982592,
     0.15790124012399101798155186502293, -0.038731993170326232824561787558492,
    -0.17672345346106673069472498415600, -0.047148599859865922722408214465630,
     0.15381110290879228162646885882803,  0.12189469877166149151050389002586,
    -0.094575111235207343404082277365271,-0.16785445385626524004858887836470,
     0.013004500444088472936407629993685, 0.17417412557117862465636166724330,
     0.071637220287468842367069460402068,-0.13936127622474077458193396572651
  },
  {
     0.13665023337521968602987906462477, -0.083331957309718312162450688895359,
    -0.16916475014679408478364306119571,  0.017327146149886432824466874566622,
     0.17592546719079780737825977787300,  0.051315565940294672644546154719392,
    -0.15590312662333390407149878284971, -0.11214594829282953553133017365260,
     0.11214594829282953553133017365260,  0.15590312662333390407149878284971,
    -0.051315565940294672644546154719392,-0.17592546719079780737825977787300,
    -0.017327146149886432824466874566622, 0.16916475014679408478364306119571,
     0.083331957309718312162450688895359,-0.13665023337521968602987906462477,
    -0.13665023337521968602987906462477,  0.083331957309718312162450688895359,
     0.16916475014679408478364306119571, -0.017327146149886432824466874566622,
    -0.17592546719079780737825977787300, -0.051315565940294672644546154719392,
     0.15590312662333390407149878284971,  0.11214594829282953553133017365260,
    -0.11214594829282953553133017365260, -0.15590312662333390407149878284971,
     0.051315565940294672644546154719392, 0.17592546719079780737825977787300,
     0.017327146149886432824466874566622,-0.16916475014679408478364306119571,
    -0.083331957309718312162450688895359, 0.13665023337521968602987906462477,
     0.13665023337521968602987906462477, -0.083331957309718312162450688895359,
    -0.16916475014679408478364306119571,  0.017327146149886432824466874566622,
     0.17592546719079780737825977787300,  0.051315565940294672644546154719392,
    -0.15590312662333390407149878284971, -0.11214594829282953553133017365260,
     0.11214594829282953553133017365260,  0.15590312662333390407149878284971,
    -0.051315565940294672644546154719392,-0.17592546719079780737825977787300,
    -0.017327146149886432824466874566622, 0.16916475014679408478364306119571,
     0.083331957309718312162450688895359,-0.13665023337521968602987906462477,
    -0.13665023337521968602987906462477,  0.083331957309718312162450688895359,
     0.16916475014679408478364306119571, -0.017327146149886432824466874566622,
    -0.17592546719079780737825977787300, -0.051315565940294672644546154719392,
     0.15590312662333390407149878284971,  0.11214594829282953553133017365260,
    -0.11214594829282953553133017365260, -0.15590312662333390407149878284971,
     0.051315565940294672644546154719392, 0.17592546719079780737825977787300,
     0.017327146149886432824466874566622,-0.16916475014679408478364306119571,
    -0.083331957309718312162450688895359, 0.13665023337521968602987906462477
  },
  {
     0.13385687753479470495813284521973, -0.094575111235207343404082277365271,
    -0.16161097942121587956084560348376,  0.047148599859865922722408214465630,
     0.17544725227590413843675978359671,  0.004338317276799999540542326666634,
    -0.17417412557117862465636166724330, -0.055451621442613590037015233508646,
     0.15790124012399101798155186502293,  0.10178946920358000230129700949629,
    -0.12803000590238956351843752239068, -0.13936127622474077458193396572651,
     0.087132913540067685045147230048925, 0.16493138371356506569293701686233,
    -0.038731993170326232824561787558492,-0.17629771118253266395745838340776,
    -0.013004500444088472936407629993685, 0.17248139814210053549033101371775,
     0.063621055212317597965816690486066,-0.15381110290879228162646885882803,
    -0.10875860742493672495725784718877,  0.12189469877166149151050389002586,
     0.14452994139365530889220759879835, -0.079480805024268148817432153524406,
    -0.16785445385626524004858887836470,  0.030222077716717299469868553789363,
     0.17672345346106673069472498415600,  0.021639354630747427065531392743213,
    -0.17037314791731193194066488665533, -0.071637220287468842367069460402068,
     0.14935042127442479271101666423961,  0.11546573663487784907576643982592,
    -0.11546573663487784907576643982592, -0.14935042127442479271101666423961,
     0.071637220287468842367069460402068, 0.17037314791731193194066488665533,
    -0.021639354630747427065531392743213,-0.17672345346106673069472498415600,
    -0.030222077716717299469868553789363, 0.16785445385626524004858887836470,
     0.079480805024268148817432153524406,-0.14452994139365530889220759879835,
    -0.12189469877166149151050389002586,  0.10875860742493672495725784718877,
     0.15381110290879228162646885882803, -0.063621055212317597965816690486066,
    -0.17248139814210053549033101371775,  0.013004500444088472936407629993685,
     0.17629771118253266395745838340776,  0.038731993170326232824561787558492,
    -0.16493138371356506569293701686233, -0.087132913540067685045147230048925,
     0.13936127622474077458193396572651,  0.12803000590238956351843752239068,
    -0.10178946920358000230129700949629, -0.15790124012399101798155186502293,
     0.055451621442613590037015233508646, 0.17417412557117862465636166724330,
    -0.004338317276799999540542326666634,-0.17544725227590413843675978359671,
    -0.047148599859865922722408214465630, 0.16161097942121587956084560348376,
     0.094575111235207343404082277365271,-0.13385687753479470495813284521973
  },
  {
     0.13098289131657380087121582271272, -0.10530575443867740272410295104619,
    -0.15162642913722598531903229617045,  0.075581776473850090965407023747544,
     0.16644304831921567823839589426492, -0.042953233225881292913572018164494,
    -0.17486335449663478165921413022414,  0.008674021313492586318780005883468,
     0.17656376002524718647512421849032,  0.025938528373526445792454998016647,
    -0.17147891927418672456199547790670, -0.059554274961645154658154180042717,
     0.15980423982190510363772032690233,  0.090881384161410012831963755651079,
    -0.14198837305251784063881548345787, -0.11871597273471929730747707847705,
     0.11871597273471929730747707847705,  0.14198837305251784063881548345787,
    -0.090881384161410012831963755651079,-0.15980423982190510363772032690233,
     0.059554274961645154658154180042717, 0.17147891927418672456199547790670,
    -0.025938528373526445792454998016647,-0.17656376002524718647512421849032,
    -0.008674021313492586318780005883468, 0.17486335449663478165921413022414,
     0.042953233225881292913572018164494,-0.16644304831921567823839589426492,
    -0.075581776473850090965407023747544, 0.15162642913722598531903229617045,
     0.10530575443867740272410295104619, -0.13098289131657380087121582271272,
    -0.13098289131657380087121582271272,  0.10530575443867740272410295104619,
     0.15162642913722598531903229617045, -0.075581776473850090965407023747544,
    -0.16644304831921567823839589426492,  0.042953233225881292913572018164494,
     0.17486335449663478165921413022414, -0.008674021313492586318780005883468,
    -0.17656376002524718647512421849032, -0.025938528373526445792454998016647,
     0.17147891927418672456199547790670,  0.059554274961645154658154180042717,
    -0.15980423982190510363772032690233, -0.090881384161410012831963755651079,
     0.14198837305251784063881548345787,  0.11871597273471929730747707847705,
    -0.11871597273471929730747707847705, -0.14198837305251784063881548345787,
     0.090881384161410012831963755651079, 0.15980423982190510363772032690233,
    -0.059554274961645154658154180042717,-0.17147891927418672456199547790670,
     0.025938528373526445792454998016647, 0.17656376002524718647512421849032,
     0.008674021313492586318780005883468,-0.17486335449663478165921413022414,
    -0.042953233225881292913572018164494, 0.16644304831921567823839589426492,
     0.075581776473850090965407023747544,-0.15162642913722598531903229617045,
    -0.10530575443867740272410295104619,  0.13098289131657380087121582271272
  },
  {
     0.12803000590238956351843752239068, -0.11546573663487784907576643982592,
    -0.13936127622474077458193396572651,  0.10178946920358000230129700949629,
     0.14935042127442479271101666423961, -0.087132913540067685045147230048925,
    -0.15790124012399101798155186502293,  0.071637220287468842367069460402068,
     0.16493138371356506569293701686233, -0.055451621442613590037015233508646,
    -0.17037314791731193194066488665533,  0.038731993170326232824561787558492,
     0.17417412557117862465636166724330, -0.021639354630747427065531392743213,
    -0.17629771118253266395745838340776,  0.004338317276799999540542326666634,
     0.17672345346106673069472498415600,  0.013004500444088472936407629993685,
    -0.17544725227590413843675978359671, -0.030222077716717299469868553789363,
     0.17248139814210053549033101371775,  0.047148599859865922722408214465630,
    -0.16785445385626524004858887836470, -0.063621055212317597965816690486066,
     0.16161097942121587956084560348376,  0.079480805024268148817432153524406,
    -0.15381110290879228162646885882803, -0.094575111235207343404082277365271,
     0.14452994139365530889220759879835,  0.10875860742493672495725784718877,
    -0.13385687753479470495813284521973, -0.12189469877166149151050389002586,
     0.12189469877166149151050389002586,  0.13385687753479470495813284521973,
    -0.10875860742493672495725784718877, -0.14452994139365530889220759879835,
     0.094575111235207343404082277365271, 0.15381110290879228162646885882803,
    -0.079480805024268148817432153524406,-0.16161097942121587956084560348376,
     0.063621055212317597965816690486066, 0.16785445385626524004858887836470,
    -0.047148599859865922722408214465630,-0.17248139814210053549033101371775,
     0.030222077716717299469868553789363, 0.17544725227590413843675978359671,
    -0.013004500444088472936407629993685,-0.17672345346106673069472498415600,
    -0.004338317276799999540542326666634, 0.17629771118253266395745838340776,
     0.021639354630747427065531392743213,-0.17417412557117862465636166724330,
    -0.038731993170326232824561787558492, 0.17037314791731193194066488665533,
     0.055451621442613590037015233508646,-0.16493138371356506569293701686233,
    -0.071637220287468842367069460402068, 0.15790124012399101798155186502293,
     0.087132913540067685045147230048925,-0.14935042127442479271101666423961,
    -0.10178946920358000230129700949629,  0.13936127622474077458193396572651,
     0.11546573663487784907576643982592, -0.12803000590238956351843752239068
  },
  {
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125,
     0.125,                              -0.125,
    -0.125,                               0.125
  },
  {
     0.12189469877166149151050389002586, -0.13385687753479470495813284521973,
    -0.10875860742493672495725784718877,  0.14452994139365530889220759879835,
     0.094575111235207343404082277365271,-0.15381110290879228162646885882803,
    -0.079480805024268148817432153524406, 0.16161097942121587956084560348376,
     0.063621055212317597965816690486066,-0.16785445385626524004858887836470,
    -0.047148599859865922722408214465630, 0.17248139814210053549033101371775,
     0.030222077716717299469868553789363,-0.17544725227590413843675978359671,
    -0.013004500444088472936407629993685, 0.17672345346106673069472498415600,
    -0.004338317276799999540542326666634,-0.17629771118253266395745838340776,
     0.021639354630747427065531392743213, 0.17417412557117862465636166724330,
    -0.038731993170326232824561787558492,-0.17037314791731193194066488665533,
     0.055451621442613590037015233508646, 0.16493138371356506569293701686233,
    -0.071637220287468842367069460402068,-0.15790124012399101798155186502293,
     0.087132913540067685045147230048925, 0.14935042127442479271101666423961,
    -0.10178946920358000230129700949629, -0.13936127622474077458193396572651,
     0.11546573663487784907576643982592,  0.12803000590238956351843752239068,
    -0.12803000590238956351843752239068, -0.11546573663487784907576643982592,
     0.13936127622474077458193396572651,  0.10178946920358000230129700949629,
    -0.14935042127442479271101666423961, -0.087132913540067685045147230048925,
     0.15790124012399101798155186502293,  0.071637220287468842367069460402068,
    -0.16493138371356506569293701686233, -0.055451621442613590037015233508646,
     0.17037314791731193194066488665533,  0.038731993170326232824561787558492,
    -0.17417412557117862465636166724330, -0.021639354630747427065531392743213,
     0.17629771118253266395745838340776,  0.004338317276799999540542326666634,
    -0.17672345346106673069472498415600,  0.013004500444088472936407629993685,
     0.17544725227590413843675978359671, -0.030222077716717299469868553789363,
    -0.17248139814210053549033101371775,  0.047148599859865922722408214465630,
     0.16785445385626524004858887836470, -0.063621055212317597965816690486066,
    -0.16161097942121587956084560348376,  0.079480805024268148817432153524406,
     0.15381110290879228162646885882803, -0.094575111235207343404082277365271,
    -0.14452994139365530889220759879835,  0.10875860742493672495725784718877,
     0.13385687753479470495813284521973, -0.12189469877166149151050389002586
  },
  {
     0.11871597273471929730747707847705, -0.14198837305251784063881548345787,
    -0.090881384161410012831963755651079, 0.15980423982190510363772032690233,
     0.059554274961645154658154180042717,-0.17147891927418672456199547790670,
    -0.025938528373526445792454998016647, 0.17656376002524718647512421849032,
    -0.008674021313492586318780005883468,-0.17486335449663478165921413022414,
     0.042953233225881292913572018164494, 0.16644304831921567823839589426492,
    -0.075581776473850090965407023747544,-0.15162642913722598531903229617045,
     0.10530575443867740272410295104619,  0.13098289131657380087121582271272,
    -0.13098289131657380087121582271272, -0.10530575443867740272410295104619,
     0.15162642913722598531903229617045,  0.075581776473850090965407023747544,
    -0.16644304831921567823839589426492, -0.042953233225881292913572018164494,
     0.17486335449663478165921413022414,  0.008674021313492586318780005883468,
    -0.17656376002524718647512421849032,  0.025938528373526445792454998016647,
     0.17147891927418672456199547790670, -0.059554274961645154658154180042717,
    -0.15980423982190510363772032690233,  0.090881384161410012831963755651079,
     0.14198837305251784063881548345787, -0.11871597273471929730747707847705,
    -0.11871597273471929730747707847705,  0.14198837305251784063881548345787,
     0.090881384161410012831963755651079,-0.15980423982190510363772032690233,
    -0.059554274961645154658154180042717, 0.17147891927418672456199547790670,
     0.025938528373526445792454998016647,-0.17656376002524718647512421849032,
     0.008674021313492586318780005883468, 0.17486335449663478165921413022414,
    -0.042953233225881292913572018164494,-0.16644304831921567823839589426492,
     0.075581776473850090965407023747544, 0.15162642913722598531903229617045,
    -0.10530575443867740272410295104619, -0.13098289131657380087121582271272,
     0.13098289131657380087121582271272,  0.10530575443867740272410295104619,
    -0.15162642913722598531903229617045, -0.075581776473850090965407023747544,
     0.16644304831921567823839589426492,  0.042953233225881292913572018164494,
    -0.17486335449663478165921413022414, -0.008674021313492586318780005883468,
     0.17656376002524718647512421849032, -0.025938528373526445792454998016647,
    -0.17147891927418672456199547790670,  0.059554274961645154658154180042717,
     0.15980423982190510363772032690233, -0.090881384161410012831963755651079,
    -0.14198837305251784063881548345787,  0.11871597273471929730747707847705
  },
  {
     0.11546573663487784907576643982592, -0.14935042127442479271101666423961,
    -0.071637220287468842367069460402068, 0.17037314791731193194066488665533,
     0.021639354630747427065531392743213,-0.17672345346106673069472498415600,
     0.030222077716717299469868553789363, 0.16785445385626524004858887836470,
    -0.079480805024268148817432153524406,-0.14452994139365530889220759879835,
     0.12189469877166149151050389002586,  0.10875860742493672495725784718877,
    -0.15381110290879228162646885882803, -0.063621055212317597965816690486066,
     0.17248139814210053549033101371775,  0.013004500444088472936407629993685,
    -0.17629771118253266395745838340776,  0.038731993170326232824561787558492,
     0.16493138371356506569293701686233, -0.087132913540067685045147230048925,
    -0.13936127622474077458193396572651,  0.12803000590238956351843752239068,
     0.10178946920358000230129700949629, -0.15790124012399101798155186502293,
    -0.055451621442613590037015233508646, 0.17417412557117862465636166724330,
     0.004338317276799999540542326666634,-0.17544725227590413843675978359671,
     0.047148599859865922722408214465630, 0.16161097942121587956084560348376,
    -0.094575111235207343404082277365271,-0.13385687753479470495813284521973,
     0.13385687753479470495813284521973,  0.094575111235207343404082277365271,
    -0.16161097942121587956084560348376, -0.047148599859865922722408214465630,
     0.17544725227590413843675978359671, -0.004338317276799999540542326666634,
    -0.17417412557117862465636166724330,  0.055451621442613590037015233508646,
     0.15790124012399101798155186502293, -0.10178946920358000230129700949629,
    -0.12803000590238956351843752239068,  0.13936127622474077458193396572651,
     0.087132913540067685045147230048925,-0.16493138371356506569293701686233,
    -0.038731993170326232824561787558492, 0.17629771118253266395745838340776,
    -0.013004500444088472936407629993685,-0.17248139814210053549033101371775,
     0.063621055212317597965816690486066, 0.15381110290879228162646885882803,
    -0.10875860742493672495725784718877, -0.12189469877166149151050389002586,
     0.14452994139365530889220759879835,  0.079480805024268148817432153524406,
    -0.16785445385626524004858887836470, -0.030222077716717299469868553789363,
     0.17672345346106673069472498415600, -0.021639354630747427065531392743213,
    -0.17037314791731193194066488665533,  0.071637220287468842367069460402068,
     0.14935042127442479271101666423961, -0.11546573663487784907576643982592
  },
  {
     0.11214594829282953553133017365260, -0.15590312662333390407149878284971,
    -0.051315565940294672644546154719392, 0.17592546719079780737825977787300,
    -0.017327146149886432824466874566622,-0.16916475014679408478364306119571,
     0.083331957309718312162450688895359, 0.13665023337521968602987906462477,
    -0.13665023337521968602987906462477, -0.083331957309718312162450688895359,
     0.16916475014679408478364306119571,  0.017327146149886432824466874566622,
    -0.17592546719079780737825977787300,  0.051315565940294672644546154719392,
     0.15590312662333390407149878284971, -0.11214594829282953553133017365260,
    -0.11214594829282953553133017365260,  0.15590312662333390407149878284971,
     0.051315565940294672644546154719392,-0.17592546719079780737825977787300,
     0.017327146149886432824466874566622, 0.16916475014679408478364306119571,
    -0.083331957309718312162450688895359,-0.13665023337521968602987906462477,
     0.13665023337521968602987906462477,  0.083331957309718312162450688895359,
    -0.16916475014679408478364306119571, -0.017327146149886432824466874566622,
     0.17592546719079780737825977787300, -0.051315565940294672644546154719392,
    -0.15590312662333390407149878284971,  0.11214594829282953553133017365260,
     0.11214594829282953553133017365260, -0.15590312662333390407149878284971,
    -0.051315565940294672644546154719392, 0.17592546719079780737825977787300,
    -0.017327146149886432824466874566622,-0.16916475014679408478364306119571,
     0.083331957309718312162450688895359, 0.13665023337521968602987906462477,
    -0.13665023337521968602987906462477, -0.083331957309718312162450688895359,
     0.16916475014679408478364306119571,  0.017327146149886432824466874566622,
    -0.17592546719079780737825977787300,  0.051315565940294672644546154719392,
     0.15590312662333390407149878284971, -0.11214594829282953553133017365260,
    -0.11214594829282953553133017365260,  0.15590312662333390407149878284971,
     0.051315565940294672644546154719392,-0.17592546719079780737825977787300,
     0.017327146149886432824466874566622, 0.16916475014679408478364306119571,
    -0.083331957309718312162450688895359,-0.13665023337521968602987906462477,
     0.13665023337521968602987906462477,  0.083331957309718312162450688895359,
    -0.16916475014679408478364306119571, -0.017327146149886432824466874566622,
     0.17592546719079780737825977787300, -0.051315565940294672644546154719392,
    -0.15590312662333390407149878284971,  0.11214594829282953553133017365260
  },
  {
     0.10875860742493672495725784718877, -0.16161097942121587956084560348376,
    -0.030222077716717299469868553789363, 0.17629771118253266395745838340776,
    -0.055451621442613590037015233508646,-0.14935042127442479271101666423961,
     0.12803000590238956351843752239068,  0.087132913540067685045147230048925,
    -0.17037314791731193194066488665533, -0.004338317276799999540542326666634,
     0.17248139814210053549033101371775, -0.079480805024268148817432153524406,
    -0.13385687753479470495813284521973,  0.14452994139365530889220759879835,
     0.063621055212317597965816690486066,-0.17544725227590413843675978359671,
     0.021639354630747427065531392743213, 0.16493138371356506569293701686233,
    -0.10178946920358000230129700949629, -0.11546573663487784907576643982592,
     0.15790124012399101798155186502293,  0.038731993170326232824561787558492,
    -0.17672345346106673069472498415600,  0.047148599859865922722408214465630,
     0.15381110290879228162646885882803, -0.12189469877166149151050389002586,
    -0.094575111235207343404082277365271, 0.16785445385626524004858887836470,
     0.013004500444088472936407629993685,-0.17417412557117862465636166724330,
     0.071637220287468842367069460402068, 0.13936127622474077458193396572651,
    -0.13936127622474077458193396572651, -0.071637220287468842367069460402068,
     0.17417412557117862465636166724330, -0.013004500444088472936407629993685,
    -0.16785445385626524004858887836470,  0.094575111235207343404082277365271,
     0.12189469877166149151050389002586, -0.15381110290879228162646885882803,
    -0.047148599859865922722408214465630, 0.17672345346106673069472498415600,
    -0.038731993170326232824561787558492,-0.15790124012399101798155186502293,
     0.11546573663487784907576643982592,  0.10178946920358000230129700949629,
    -0.16493138371356506569293701686233, -0.021639354630747427065531392743213,
     0.17544725227590413843675978359671, -0.063621055212317597965816690486066,
    -0.14452994139365530889220759879835,  0.13385687753479470495813284521973,
     0.079480805024268148817432153524406,-0.17248139814210053549033101371775,
     0.004338317276799999540542326666634, 0.17037314791731193194066488665533,
    -0.087132913540067685045147230048925,-0.12803000590238956351843752239068,
     0.14935042127442479271101666423961,  0.055451621442613590037015233508646,
    -0.17629771118253266395745838340776,  0.030222077716717299469868553789363,
     0.16161097942121587956084560348376, -0.10875860742493672495725784718877
  },
  {
     0.10530575443867740272410295104619, -0.16644304831921567823839589426492,
    -0.008674021313492586318780005883468, 0.17147891927418672456199547790670,
    -0.090881384161410012831963755651079,-0.11871597273471929730747707847705,
     0.15980423982190510363772032690233,  0.025938528373526445792454998016647,
    -0.17486335449663478165921413022414,  0.075581776473850090965407023747544,
     0.13098289131657380087121582271272, -0.15162642913722598531903229617045,
    -0.042953233225881292913572018164494, 0.17656376002524718647512421849032,
    -0.059554274961645154658154180042717,-0.14198837305251784063881548345787,
     0.14198837305251784063881548345787,  0.059554274961645154658154180042717,
    -0.17656376002524718647512421849032,  0.042953233225881292913572018164494,
     0.15162642913722598531903229617045, -0.13098289131657380087121582271272,
    -0.075581776473850090965407023747544, 0.17486335449663478165921413022414,
    -0.025938528373526445792454998016647,-0.15980423982190510363772032690233,
     0.11871597273471929730747707847705,  0.090881384161410012831963755651079,
    -0.17147891927418672456199547790670,  0.008674021313492586318780005883468,
     0.16644304831921567823839589426492, -0.10530575443867740272410295104619,
    -0.10530575443867740272410295104619,  0.16644304831921567823839589426492,
     0.008674021313492586318780005883468,-0.17147891927418672456199547790670,
     0.090881384161410012831963755651079, 0.11871597273471929730747707847705,
    -0.15980423982190510363772032690233, -0.025938528373526445792454998016647,
     0.17486335449663478165921413022414, -0.075581776473850090965407023747544,
    -0.13098289131657380087121582271272,  0.15162642913722598531903229617045,
     0.042953233225881292913572018164494,-0.17656376002524718647512421849032,
     0.059554274961645154658154180042717, 0.14198837305251784063881548345787,
    -0.14198837305251784063881548345787, -0.059554274961645154658154180042717,
     0.17656376002524718647512421849032, -0.042953233225881292913572018164494,
    -0.15162642913722598531903229617045,  0.13098289131657380087121582271272,
     0.075581776473850090965407023747544,-0.17486335449663478165921413022414,
     0.025938528373526445792454998016647, 0.15980423982190510363772032690233,
    -0.11871597273471929730747707847705, -0.090881384161410012831963755651079,
     0.17147891927418672456199547790670, -0.008674021313492586318780005883468,
    -0.16644304831921567823839589426492,  0.10530575443867740272410295104619
  },
  {
     0.10178946920358000230129700949629, -0.17037314791731193194066488665533,
     0.013004500444088472936407629993685, 0.16161097942121587956084560348376,
    -0.12189469877166149151050389002586, -0.079480805024268148817432153524406,
     0.17544725227590413843675978359671, -0.038731993170326232824561787558492,
    -0.14935042127442479271101666423961,  0.13936127622474077458193396572651,
     0.055451621442613590037015233508646,-0.17672345346106673069472498415600,
     0.063621055212317597965816690486066, 0.13385687753479470495813284521973,
    -0.15381110290879228162646885882803, -0.030222077716717299469868553789363,
     0.17417412557117862465636166724330, -0.087132913540067685045147230048925,
    -0.11546573663487784907576643982592,  0.16493138371356506569293701686233,
     0.004338317276799999540542326666634,-0.16785445385626524004858887836470,
     0.10875860742493672495725784718877,  0.094575111235207343404082277365271,
    -0.17248139814210053549033101371775,  0.021639354630747427065531392743213,
     0.15790124012399101798155186502293, -0.12803000590238956351843752239068,
    -0.071637220287468842367069460402068, 0.17629771118253266395745838340776,
    -0.047148599859865922722408214465630,-0.14452994139365530889220759879835,
     0.14452994139365530889220759879835,  0.047148599859865922722408214465630,
    -0.17629771118253266395745838340776,  0.071637220287468842367069460402068,
     0.12803000590238956351843752239068, -0.15790124012399101798155186502293,
    -0.021639354630747427065531392743213, 0.17248139814210053549033101371775,
    -0.094575111235207343404082277365271,-0.10875860742493672495725784718877,
     0.16785445385626524004858887836470, -0.004338317276799999540542326666634,
    -0.16493138371356506569293701686233,  0.11546573663487784907576643982592,
     0.087132913540067685045147230048925,-0.17417412557117862465636166724330,
     0.030222077716717299469868553789363, 0.15381110290879228162646885882803,
    -0.13385687753479470495813284521973, -0.063621055212317597965816690486066,
     0.17672345346106673069472498415600, -0.055451621442613590037015233508646,
    -0.13936127622474077458193396572651,  0.14935042127442479271101666423961,
     0.038731993170326232824561787558492,-0.17544725227590413843675978359671,
     0.079480805024268148817432153524406, 0.12189469877166149151050389002586,
    -0.16161097942121587956084560348376, -0.013004500444088472936407629993685,
     0.17037314791731193194066488665533, -0.10178946920358000230129700949629
  },
  {
     0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
     0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
     0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
     0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
     0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
     0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
     0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
     0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
     0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
     0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
     0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
     0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580,  0.098211869798387772659737170957152,
     0.098211869798387772659737170957152,-0.17337998066526843272770239894580,
     0.034487422410367876541994695458672, 0.14698445030241983962180838807641,
    -0.14698445030241983962180838807641, -0.034487422410367876541994695458672,
     0.17337998066526843272770239894580, -0.098211869798387772659737170957152,
    -0.098211869798387772659737170957152, 0.17337998066526843272770239894580,
    -0.034487422410367876541994695458672,-0.14698445030241983962180838807641,
     0.14698445030241983962180838807641,  0.034487422410367876541994695458672,
    -0.17337998066526843272770239894580,  0.098211869798387772659737170957152
  },
  {
     0.094575111235207343404082277365271,-0.17544725227590413843675978359671,
     0.055451621442613590037015233508646, 0.12803000590238956351843752239068,
    -0.16493138371356506569293701686233,  0.013004500444088472936407629993685,
     0.15381110290879228162646885882803, -0.14452994139365530889220759879835,
    -0.030222077716717299469868553789363, 0.17037314791731193194066488665533,
    -0.11546573663487784907576643982592, -0.071637220287468842367069460402068,
     0.17672345346106673069472498415600, -0.079480805024268148817432153524406,
    -0.10875860742493672495725784718877,  0.17248139814210053549033101371775,
    -0.038731993170326232824561787558492,-0.13936127622474077458193396572651,
     0.15790124012399101798155186502293,  0.004338317276799999540542326666634,
    -0.16161097942121587956084560348376,  0.13385687753479470495813284521973,
     0.047148599859865922722408214465630,-0.17417412557117862465636166724330,
     0.10178946920358000230129700949629,  0.087132913540067685045147230048925,
    -0.17629771118253266395745838340776,  0.063621055212317597965816690486066,
     0.12189469877166149151050389002586, -0.16785445385626524004858887836470,
     0.021639354630747427065531392743213, 0.14935042127442479271101666423961,
    -0.14935042127442479271101666423961, -0.021639354630747427065531392743213,
     0.16785445385626524004858887836470, -0.12189469877166149151050389002586,
    -0.063621055212317597965816690486066, 0.17629771118253266395745838340776,
    -0.087132913540067685045147230048925,-0.10178946920358000230129700949629,
     0.17417412557117862465636166724330, -0.047148599859865922722408214465630,
    -0.13385687753479470495813284521973,  0.16161097942121587956084560348376,
    -0.004338317276799999540542326666634,-0.15790124012399101798155186502293,
     0.13936127622474077458193396572651,  0.038731993170326232824561787558492,
    -0.17248139814210053549033101371775,  0.10875860742493672495725784718877,
     0.079480805024268148817432153524406,-0.17672345346106673069472498415600,
     0.071637220287468842367069460402068, 0.11546573663487784907576643982592,
    -0.17037314791731193194066488665533,  0.030222077716717299469868553789363,
     0.14452994139365530889220759879835, -0.15381110290879228162646885882803,
    -0.013004500444088472936407629993685, 0.16493138371356506569293701686233,
    -0.12803000590238956351843752239068, -0.055451621442613590037015233508646,
     0.17544725227590413843675978359671, -0.094575111235207343404082277365271
  },
  {
     0.090881384161410012831963755651079,-0.17656376002524718647512421849032,
     0.075581776473850090965407023747544, 0.10530575443867740272410295104619,
    -0.17486335449663478165921413022414,  0.059554274961645154658154180042717,
     0.11871597273471929730747707847705, -0.17147891927418672456199547790670,
     0.042953233225881292913572018164494, 0.13098289131657380087121582271272,
    -0.16644304831921567823839589426492,  0.025938528373526445792454998016647,
     0.14198837305251784063881548345787, -0.15980423982190510363772032690233,
     0.008674021313492586318780005883468, 0.15162642913722598531903229617045,
    -0.15162642913722598531903229617045, -0.008674021313492586318780005883468,
     0.15980423982190510363772032690233, -0.14198837305251784063881548345787,
    -0.025938528373526445792454998016647, 0.16644304831921567823839589426492,
    -0.13098289131657380087121582271272, -0.042953233225881292913572018164494,
     0.17147891927418672456199547790670, -0.11871597273471929730747707847705,
    -0.059554274961645154658154180042717, 0.17486335449663478165921413022414,
    -0.10530575443867740272410295104619, -0.075581776473850090965407023747544,
     0.17656376002524718647512421849032, -0.090881384161410012831963755651079,
    -0.090881384161410012831963755651079, 0.17656376002524718647512421849032,
    -0.075581776473850090965407023747544,-0.10530575443867740272410295104619,
     0.17486335449663478165921413022414, -0.059554274961645154658154180042717,
    -0.11871597273471929730747707847705,  0.17147891927418672456199547790670,
    -0.042953233225881292913572018164494,-0.13098289131657380087121582271272,
     0.16644304831921567823839589426492, -0.025938528373526445792454998016647,
    -0.14198837305251784063881548345787,  0.15980423982190510363772032690233,
    -0.008674021313492586318780005883468,-0.15162642913722598531903229617045,
     0.15162642913722598531903229617045,  0.008674021313492586318780005883468,
    -0.15980423982190510363772032690233,  0.14198837305251784063881548345787,
     0.025938528373526445792454998016647,-0.16644304831921567823839589426492,
     0.13098289131657380087121582271272,  0.042953233225881292913572018164494,
    -0.17147891927418672456199547790670,  0.11871597273471929730747707847705,
     0.059554274961645154658154180042717,-0.17486335449663478165921413022414,
     0.10530575443867740272410295104619,  0.075581776473850090965407023747544,
    -0.17656376002524718647512421849032,  0.090881384161410012831963755651079
  },
  {
     0.087132913540067685045147230048925,-0.17672345346106673069472498415600,
     0.094575111235207343404082277365271, 0.079480805024268148817432153524406,
    -0.17629771118253266395745838340776,  0.10178946920358000230129700949629,
     0.071637220287468842367069460402068,-0.17544725227590413843675978359671,
     0.10875860742493672495725784718877,  0.063621055212317597965816690486066,
    -0.17417412557117862465636166724330,  0.11546573663487784907576643982592,
     0.055451621442613590037015233508646,-0.17248139814210053549033101371775,
     0.12189469877166149151050389002586,  0.047148599859865922722408214465630,
    -0.17037314791731193194066488665533,  0.12803000590238956351843752239068,
     0.038731993170326232824561787558492,-0.16785445385626524004858887836470,
     0.13385687753479470495813284521973,  0.030222077716717299469868553789363,
    -0.16493138371356506569293701686233,  0.13936127622474077458193396572651,
     0.021639354630747427065531392743213,-0.16161097942121587956084560348376,
     0.14452994139365530889220759879835,  0.013004500444088472936407629993685,
    -0.15790124012399101798155186502293,  0.14935042127442479271101666423961,
     0.004338317276799999540542326666634,-0.15381110290879228162646885882803,
     0.15381110290879228162646885882803, -0.004338317276799999540542326666634,
    -0.14935042127442479271101666423961,  0.15790124012399101798155186502293,
    -0.013004500444088472936407629993685,-0.14452994139365530889220759879835,
     0.16161097942121587956084560348376, -0.021639354630747427065531392743213,
    -0.13936127622474077458193396572651,  0.16493138371356506569293701686233,
    -0.030222077716717299469868553789363,-0.13385687753479470495813284521973,
     0.16785445385626524004858887836470, -0.038731993170326232824561787558492,
    -0.12803000590238956351843752239068,  0.17037314791731193194066488665533,
    -0.047148599859865922722408214465630,-0.12189469877166149151050389002586,
     0.17248139814210053549033101371775, -0.055451621442613590037015233508646,
    -0.11546573663487784907576643982592,  0.17417412557117862465636166724330,
    -0.063621055212317597965816690486066,-0.10875860742493672495725784718877,
     0.17544725227590413843675978359671, -0.071637220287468842367069460402068,
    -0.10178946920358000230129700949629,  0.17629771118253266395745838340776,
    -0.079480805024268148817432153524406,-0.094575111235207343404082277365271,
     0.17672345346106673069472498415600, -0.087132913540067685045147230048925
  },
  {
     0.083331957309718312162450688895359,-0.17592546719079780737825977787300,
     0.11214594829282953553133017365260,  0.051315565940294672644546154719392,
    -0.16916475014679408478364306119571,  0.13665023337521968602987906462477,
     0.017327146149886432824466874566622,-0.15590312662333390407149878284971,
     0.15590312662333390407149878284971, -0.017327146149886432824466874566622,
    -0.13665023337521968602987906462477,  0.16916475014679408478364306119571,
    -0.051315565940294672644546154719392,-0.11214594829282953553133017365260,
     0.17592546719079780737825977787300, -0.083331957309718312162450688895359,
    -0.083331957309718312162450688895359, 0.17592546719079780737825977787300,
    -0.11214594829282953553133017365260, -0.051315565940294672644546154719392,
     0.16916475014679408478364306119571, -0.13665023337521968602987906462477,
    -0.017327146149886432824466874566622, 0.15590312662333390407149878284971,
    -0.15590312662333390407149878284971,  0.017327146149886432824466874566622,
     0.13665023337521968602987906462477, -0.16916475014679408478364306119571,
     0.051315565940294672644546154719392, 0.11214594829282953553133017365260,
    -0.17592546719079780737825977787300,  0.083331957309718312162450688895359,
     0.083331957309718312162450688895359,-0.17592546719079780737825977787300,
     0.11214594829282953553133017365260,  0.051315565940294672644546154719392,
    -0.16916475014679408478364306119571,  0.13665023337521968602987906462477,
     0.017327146149886432824466874566622,-0.15590312662333390407149878284971,
     0.15590312662333390407149878284971, -0.017327146149886432824466874566622,
    -0.13665023337521968602987906462477,  0.16916475014679408478364306119571,
    -0.051315565940294672644546154719392,-0.11214594829282953553133017365260,
     0.17592546719079780737825977787300, -0.083331957309718312162450688895359,
    -0.083331957309718312162450688895359, 0.17592546719079780737825977787300,
    -0.11214594829282953553133017365260, -0.051315565940294672644546154719392,
     0.16916475014679408478364306119571, -0.13665023337521968602987906462477,
    -0.017327146149886432824466874566622, 0.15590312662333390407149878284971,
    -0.15590312662333390407149878284971,  0.017327146149886432824466874566622,
     0.13665023337521968602987906462477, -0.16916475014679408478364306119571,
     0.051315565940294672644546154719392, 0.11214594829282953553133017365260,
    -0.17592546719079780737825977787300,  0.083331957309718312162450688895359
  },
  {
     0.079480805024268148817432153524406,-0.17417412557117862465636166724330,
     0.12803000590238956351843752239068,  0.021639354630747427065531392743213,
    -0.15381110290879228162646885882803,  0.16161097942121587956084560348376,
    -0.038731993170326232824561787558492,-0.11546573663487784907576643982592,
     0.17629771118253266395745838340776, -0.094575111235207343404082277365271,
    -0.063621055212317597965816690486066, 0.17037314791731193194066488665533,
    -0.13936127622474077458193396572651, -0.004338317276799999540542326666634,
     0.14452994139365530889220759879835, -0.16785445385626524004858887836470,
     0.055451621442613590037015233508646, 0.10178946920358000230129700949629,
    -0.17672345346106673069472498415600,  0.10875860742493672495725784718877,
     0.047148599859865922722408214465630,-0.16493138371356506569293701686233,
     0.14935042127442479271101666423961, -0.013004500444088472936407629993685,
    -0.13385687753479470495813284521973,  0.17248139814210053549033101371775,
    -0.071637220287468842367069460402068,-0.087132913540067685045147230048925,
     0.17544725227590413843675978359671, -0.12189469877166149151050389002586,
    -0.030222077716717299469868553789363, 0.15790124012399101798155186502293,
    -0.15790124012399101798155186502293,  0.030222077716717299469868553789363,
     0.12189469877166149151050389002586, -0.17544725227590413843675978359671,
     0.087132913540067685045147230048925, 0.071637220287468842367069460402068,
    -0.17248139814210053549033101371775,  0.13385687753479470495813284521973,
     0.013004500444088472936407629993685,-0.14935042127442479271101666423961,
     0.16493138371356506569293701686233, -0.047148599859865922722408214465630,
    -0.10875860742493672495725784718877,  0.17672345346106673069472498415600,
    -0.10178946920358000230129700949629, -0.055451621442613590037015233508646,
     0.16785445385626524004858887836470, -0.14452994139365530889220759879835,
     0.004338317276799999540542326666634, 0.13936127622474077458193396572651,
    -0.17037314791731193194066488665533,  0.063621055212317597965816690486066,
     0.094575111235207343404082277365271,-0.17629771118253266395745838340776,
     0.11546573663487784907576643982592,  0.038731993170326232824561787558492,
    -0.16161097942121587956084560348376,  0.15381110290879228162646885882803,
    -0.021639354630747427065531392743213,-0.12803000590238956351843752239068,
     0.17417412557117862465636166724330, -0.079480805024268148817432153524406
  },
  {
     0.075581776473850090965407023747544,-0.17147891927418672456199547790670,
     0.14198837305251784063881548345787, -0.008674021313492586318780005883468,
    -0.13098289131657380087121582271272,  0.17486335449663478165921413022414,
    -0.090881384161410012831963755651079,-0.059554274961645154658154180042717,
     0.16644304831921567823839589426492, -0.15162642913722598531903229617045,
     0.025938528373526445792454998016647, 0.11871597273471929730747707847705,
    -0.17656376002524718647512421849032,  0.10530575443867740272410295104619,
     0.042953233225881292913572018164494,-0.15980423982190510363772032690233,
     0.15980423982190510363772032690233, -0.042953233225881292913572018164494,
    -0.10530575443867740272410295104619,  0.17656376002524718647512421849032,
    -0.11871597273471929730747707847705, -0.025938528373526445792454998016647,
     0.15162642913722598531903229617045, -0.16644304831921567823839589426492,
     0.059554274961645154658154180042717, 0.090881384161410012831963755651079,
    -0.17486335449663478165921413022414,  0.13098289131657380087121582271272,
     0.008674021313492586318780005883468,-0.14198837305251784063881548345787,
     0.17147891927418672456199547790670, -0.075581776473850090965407023747544,
    -0.075581776473850090965407023747544, 0.17147891927418672456199547790670,
    -0.14198837305251784063881548345787,  0.008674021313492586318780005883468,
     0.13098289131657380087121582271272, -0.17486335449663478165921413022414,
     0.090881384161410012831963755651079, 0.059554274961645154658154180042717,
    -0.16644304831921567823839589426492,  0.15162642913722598531903229617045,
    -0.025938528373526445792454998016647,-0.11871597273471929730747707847705,
     0.17656376002524718647512421849032, -0.10530575443867740272410295104619,
    -0.042953233225881292913572018164494, 0.15980423982190510363772032690233,
    -0.15980423982190510363772032690233,  0.042953233225881292913572018164494,
     0.10530575443867740272410295104619, -0.17656376002524718647512421849032,
     0.11871597273471929730747707847705,  0.025938528373526445792454998016647,
    -0.15162642913722598531903229617045,  0.16644304831921567823839589426492,
    -0.059554274961645154658154180042717,-0.090881384161410012831963755651079,
     0.17486335449663478165921413022414, -0.13098289131657380087121582271272,
    -0.008674021313492586318780005883468, 0.14198837305251784063881548345787,
    -0.17147891927418672456199547790670,  0.075581776473850090965407023747544
  },
  {
     0.071637220287468842367069460402068,-0.16785445385626524004858887836470,
     0.15381110290879228162646885882803, -0.038731993170326232824561787558492,
    -0.10178946920358000230129700949629,  0.17544725227590413843675978359671,
    -0.13385687753479470495813284521973,  0.004338317276799999540542326666634,
     0.12803000590238956351843752239068, -0.17629771118253266395745838340776,
     0.10875860742493672495725784718877,  0.030222077716717299469868553789363,
    -0.14935042127442479271101666423961,  0.17037314791731193194066488665533,
    -0.079480805024268148817432153524406,-0.063621055212317597965816690486066,
     0.16493138371356506569293701686233, -0.15790124012399101798155186502293,
     0.047148599859865922722408214465630, 0.094575111235207343404082277365271,
    -0.17417412557117862465636166724330,  0.13936127622474077458193396572651,
    -0.013004500444088472936407629993685,-0.12189469877166149151050389002586,
     0.17672345346106673069472498415600, -0.11546573663487784907576643982592,
    -0.021639354630747427065531392743213, 0.14452994139365530889220759879835,
    -0.17248139814210053549033101371775,  0.087132913540067685045147230048925,
     0.055451621442613590037015233508646,-0.16161097942121587956084560348376,
     0.16161097942121587956084560348376, -0.055451621442613590037015233508646,
    -0.087132913540067685045147230048925, 0.17248139814210053549033101371775,
    -0.14452994139365530889220759879835,  0.021639354630747427065531392743213,
     0.11546573663487784907576643982592, -0.17672345346106673069472498415600,
     0.12189469877166149151050389002586,  0.013004500444088472936407629993685,
    -0.13936127622474077458193396572651,  0.17417412557117862465636166724330,
    -0.094575111235207343404082277365271,-0.047148599859865922722408214465630,
     0.15790124012399101798155186502293, -0.16493138371356506569293701686233,
     0.063621055212317597965816690486066, 0.079480805024268148817432153524406,
    -0.17037314791731193194066488665533,  0.14935042127442479271101666423961,
    -0.030222077716717299469868553789363,-0.10875860742493672495725784718877,
     0.17629771118253266395745838340776, -0.12803000590238956351843752239068,
    -0.004338317276799999540542326666634, 0.13385687753479470495813284521973,
    -0.17544725227590413843675978359671,  0.10178946920358000230129700949629,
     0.038731993170326232824561787558492,-0.15381110290879228162646885882803,
     0.16785445385626524004858887836470, -0.071637220287468842367069460402068
  },
  {
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799,
     0.067649512518274623049965400670799,-0.16332037060954706598208039667840,
     0.16332037060954706598208039667840, -0.067649512518274623049965400670799,
    -0.067649512518274623049965400670799, 0.16332037060954706598208039667840,
    -0.16332037060954706598208039667840,  0.067649512518274623049965400670799
  },
  {
     0.063621055212317597965816690486066,-0.15790124012399101798155186502293,
     0.17037314791731193194066488665533, -0.094575111235207343404082277365271,
    -0.030222077716717299469868553789363, 0.13936127622474077458193396572651,
    -0.17629771118253266395745838340776,  0.12189469877166149151050389002586,
    -0.004338317276799999540542326666634,-0.11546573663487784907576643982592,
     0.17544725227590413843675978359671, -0.14452994139365530889220759879835,
     0.038731993170326232824561787558492, 0.087132913540067685045147230048925,
    -0.16785445385626524004858887836470,  0.16161097942121587956084560348376,
    -0.071637220287468842367069460402068,-0.055451621442613590037015233508646,
     0.15381110290879228162646885882803, -0.17248139814210053549033101371775,
     0.10178946920358000230129700949629,  0.021639354630747427065531392743213,
    -0.13385687753479470495813284521973,  0.17672345346106673069472498415600,
    -0.12803000590238956351843752239068,  0.013004500444088472936407629993685,
     0.10875860742493672495725784718877, -0.17417412557117862465636166724330,
     0.14935042127442479271101666423961, -0.047148599859865922722408214465630,
    -0.079480805024268148817432153524406, 0.16493138371356506569293701686233,
    -0.16493138371356506569293701686233,  0.079480805024268148817432153524406,
     0.047148599859865922722408214465630,-0.14935042127442479271101666423961,
     0.17417412557117862465636166724330, -0.10875860742493672495725784718877,
    -0.013004500444088472936407629993685, 0.12803000590238956351843752239068,
    -0.17672345346106673069472498415600,  0.13385687753479470495813284521973,
    -0.021639354630747427065531392743213,-0.10178946920358000230129700949629,
     0.17248139814210053549033101371775, -0.15381110290879228162646885882803,
     0.055451621442613590037015233508646, 0.071637220287468842367069460402068,
    -0.16161097942121587956084560348376,  0.16785445385626524004858887836470,
    -0.087132913540067685045147230048925,-0.038731993170326232824561787558492,
     0.14452994139365530889220759879835, -0.17544725227590413843675978359671,
     0.11546573663487784907576643982592,  0.004338317276799999540542326666634,
    -0.12189469877166149151050389002586,  0.17629771118253266395745838340776,
    -0.13936127622474077458193396572651,  0.030222077716717299469868553789363,
     0.094575111235207343404082277365271,-0.17037314791731193194066488665533,
     0.15790124012399101798155186502293, -0.063621055212317597965816690486066
  },
  {
     0.059554274961645154658154180042717,-0.15162642913722598531903229617045,
     0.17486335449663478165921413022414, -0.11871597273471929730747707847705,
     0.008674021313492586318780005883468, 0.10530575443867740272410295104619,
    -0.17147891927418672456199547790670,  0.15980423982190510363772032690233,
    -0.075581776473850090965407023747544,-0.042953233225881292913572018164494,
     0.14198837305251784063881548345787, -0.17656376002524718647512421849032,
     0.13098289131657380087121582271272, -0.025938528373526445792454998016647,
    -0.090881384161410012831963755651079, 0.16644304831921567823839589426492,
    -0.16644304831921567823839589426492,  0.090881384161410012831963755651079,
     0.025938528373526445792454998016647,-0.13098289131657380087121582271272,
     0.17656376002524718647512421849032, -0.14198837305251784063881548345787,
     0.042953233225881292913572018164494, 0.075581776473850090965407023747544,
    -0.15980423982190510363772032690233,  0.17147891927418672456199547790670,
    -0.10530575443867740272410295104619, -0.008674021313492586318780005883468,
     0.11871597273471929730747707847705, -0.17486335449663478165921413022414,
     0.15162642913722598531903229617045, -0.059554274961645154658154180042717,
    -0.059554274961645154658154180042717, 0.15162642913722598531903229617045,
    -0.17486335449663478165921413022414,  0.11871597273471929730747707847705,
    -0.008674021313492586318780005883468,-0.10530575443867740272410295104619,
     0.17147891927418672456199547790670, -0.15980423982190510363772032690233,
     0.075581776473850090965407023747544, 0.042953233225881292913572018164494,
    -0.14198837305251784063881548345787,  0.17656376002524718647512421849032,
    -0.13098289131657380087121582271272,  0.025938528373526445792454998016647,
     0.090881384161410012831963755651079,-0.16644304831921567823839589426492,
     0.16644304831921567823839589426492, -0.090881384161410012831963755651079,
    -0.025938528373526445792454998016647, 0.13098289131657380087121582271272,
    -0.17656376002524718647512421849032,  0.14198837305251784063881548345787,
    -0.042953233225881292913572018164494,-0.075581776473850090965407023747544,
     0.15980423982190510363772032690233, -0.17147891927418672456199547790670,
     0.10530575443867740272410295104619,  0.008674021313492586318780005883468,
    -0.11871597273471929730747707847705,  0.17486335449663478165921413022414,
    -0.15162642913722598531903229617045,  0.059554274961645154658154180042717
  },
  {
     0.055451621442613590037015233508646,-0.14452994139365530889220759879835,
     0.17672345346106673069472498415600, -0.13936127622474077458193396572651,
     0.047148599859865922722408214465630, 0.063621055212317597965816690486066,
    -0.14935042127442479271101666423961,  0.17629771118253266395745838340776,
    -0.13385687753479470495813284521973,  0.038731993170326232824561787558492,
     0.071637220287468842367069460402068,-0.15381110290879228162646885882803,
     0.17544725227590413843675978359671, -0.12803000590238956351843752239068,
     0.030222077716717299469868553789363, 0.079480805024268148817432153524406,
    -0.15790124012399101798155186502293,  0.17417412557117862465636166724330,
    -0.12189469877166149151050389002586,  0.021639354630747427065531392743213,
     0.087132913540067685045147230048925,-0.16161097942121587956084560348376,
     0.17248139814210053549033101371775, -0.11546573663487784907576643982592,
     0.013004500444088472936407629993685, 0.094575111235207343404082277365271,
    -0.16493138371356506569293701686233,  0.17037314791731193194066488665533,
    -0.10875860742493672495725784718877,  0.004338317276799999540542326666634,
     0.10178946920358000230129700949629, -0.16785445385626524004858887836470,
     0.16785445385626524004858887836470, -0.10178946920358000230129700949629,
    -0.004338317276799999540542326666634, 0.10875860742493672495725784718877,
    -0.17037314791731193194066488665533,  0.16493138371356506569293701686233,
    -0.094575111235207343404082277365271,-0.013004500444088472936407629993685,
     0.11546573663487784907576643982592, -0.17248139814210053549033101371775,
     0.16161097942121587956084560348376, -0.087132913540067685045147230048925,
    -0.021639354630747427065531392743213, 0.12189469877166149151050389002586,
    -0.17417412557117862465636166724330,  0.15790124012399101798155186502293,
    -0.079480805024268148817432153524406,-0.030222077716717299469868553789363,
     0.12803000590238956351843752239068, -0.17544725227590413843675978359671,
     0.15381110290879228162646885882803, -0.071637220287468842367069460402068,
    -0.038731993170326232824561787558492, 0.13385687753479470495813284521973,
    -0.17629771118253266395745838340776,  0.14935042127442479271101666423961,
    -0.063621055212317597965816690486066,-0.047148599859865922722408214465630,
     0.13936127622474077458193396572651, -0.17672345346106673069472498415600,
     0.14452994139365530889220759879835, -0.055451621442613590037015233508646
  },
  {
     0.051315565940294672644546154719392,-0.13665023337521968602987906462477,
     0.17592546719079780737825977787300, -0.15590312662333390407149878284971,
     0.083331957309718312162450688895359, 0.017327146149886432824466874566622,
    -0.11214594829282953553133017365260,  0.16916475014679408478364306119571,
    -0.16916475014679408478364306119571,  0.11214594829282953553133017365260,
    -0.017327146149886432824466874566622,-0.083331957309718312162450688895359,
     0.15590312662333390407149878284971, -0.17592546719079780737825977787300,
     0.13665023337521968602987906462477, -0.051315565940294672644546154719392,
    -0.051315565940294672644546154719392, 0.13665023337521968602987906462477,
    -0.17592546719079780737825977787300,  0.15590312662333390407149878284971,
    -0.083331957309718312162450688895359,-0.017327146149886432824466874566622,
     0.11214594829282953553133017365260, -0.16916475014679408478364306119571,
     0.16916475014679408478364306119571, -0.11214594829282953553133017365260,
     0.017327146149886432824466874566622, 0.083331957309718312162450688895359,
    -0.15590312662333390407149878284971,  0.17592546719079780737825977787300,
    -0.13665023337521968602987906462477,  0.051315565940294672644546154719392,
     0.051315565940294672644546154719392,-0.13665023337521968602987906462477,
     0.17592546719079780737825977787300, -0.15590312662333390407149878284971,
     0.083331957309718312162450688895359, 0.017327146149886432824466874566622,
    -0.11214594829282953553133017365260,  0.16916475014679408478364306119571,
    -0.16916475014679408478364306119571,  0.11214594829282953553133017365260,
    -0.017327146149886432824466874566622,-0.083331957309718312162450688895359,
     0.15590312662333390407149878284971, -0.17592546719079780737825977787300,
     0.13665023337521968602987906462477, -0.051315565940294672644546154719392,
    -0.051315565940294672644546154719392, 0.13665023337521968602987906462477,
    -0.17592546719079780737825977787300,  0.15590312662333390407149878284971,
    -0.083331957309718312162450688895359,-0.017327146149886432824466874566622,
     0.11214594829282953553133017365260, -0.16916475014679408478364306119571,
     0.16916475014679408478364306119571, -0.11214594829282953553133017365260,
     0.017327146149886432824466874566622, 0.083331957309718312162450688895359,
    -0.15590312662333390407149878284971,  0.17592546719079780737825977787300,
    -0.13665023337521968602987906462477,  0.051315565940294672644546154719392
  },
  {
     0.047148599859865922722408214465630,-0.12803000590238956351843752239068,
     0.17248139814210053549033101371775, -0.16785445385626524004858887836470,
     0.11546573663487784907576643982592, -0.030222077716717299469868553789363,
    -0.063621055212317597965816690486066, 0.13936127622474077458193396572651,
    -0.17544725227590413843675978359671,  0.16161097942121587956084560348376,
    -0.10178946920358000230129700949629,  0.013004500444088472936407629993685,
     0.079480805024268148817432153524406,-0.14935042127442479271101666423961,
     0.17672345346106673069472498415600, -0.15381110290879228162646885882803,
     0.087132913540067685045147230048925, 0.004338317276799999540542326666634,
    -0.094575111235207343404082277365271, 0.15790124012399101798155186502293,
    -0.17629771118253266395745838340776,  0.14452994139365530889220759879835,
    -0.071637220287468842367069460402068,-0.021639354630747427065531392743213,
     0.10875860742493672495725784718877, -0.16493138371356506569293701686233,
     0.17417412557117862465636166724330, -0.13385687753479470495813284521973,
     0.055451621442613590037015233508646, 0.038731993170326232824561787558492,
    -0.12189469877166149151050389002586,  0.17037314791731193194066488665533,
    -0.17037314791731193194066488665533,  0.12189469877166149151050389002586,
    -0.038731993170326232824561787558492,-0.055451621442613590037015233508646,
     0.13385687753479470495813284521973, -0.17417412557117862465636166724330,
     0.16493138371356506569293701686233, -0.10875860742493672495725784718877,
     0.021639354630747427065531392743213, 0.071637220287468842367069460402068,
    -0.14452994139365530889220759879835,  0.17629771118253266395745838340776,
    -0.15790124012399101798155186502293,  0.094575111235207343404082277365271,
    -0.004338317276799999540542326666634,-0.087132913540067685045147230048925,
     0.15381110290879228162646885882803, -0.17672345346106673069472498415600,
     0.14935042127442479271101666423961, -0.079480805024268148817432153524406,
    -0.013004500444088472936407629993685, 0.10178946920358000230129700949629,
    -0.16161097942121587956084560348376,  0.17544725227590413843675978359671,
    -0.13936127622474077458193396572651,  0.063621055212317597965816690486066,
     0.030222077716717299469868553789363,-0.11546573663487784907576643982592,
     0.16785445385626524004858887836470, -0.17248139814210053549033101371775,
     0.12803000590238956351843752239068, -0.047148599859865922722408214465630
  },
  {
     0.042953233225881292913572018164494,-0.11871597273471929730747707847705,
     0.16644304831921567823839589426492, -0.17486335449663478165921413022414,
     0.14198837305251784063881548345787, -0.075581776473850090965407023747544,
    -0.008674021313492586318780005883468, 0.090881384161410012831963755651079,
    -0.15162642913722598531903229617045,  0.17656376002524718647512421849032,
    -0.15980423982190510363772032690233,  0.10530575443867740272410295104619,
    -0.025938528373526445792454998016647,-0.059554274961645154658154180042717,
     0.13098289131657380087121582271272, -0.17147891927418672456199547790670,
     0.17147891927418672456199547790670, -0.13098289131657380087121582271272,
     0.059554274961645154658154180042717, 0.025938528373526445792454998016647,
    -0.10530575443867740272410295104619,  0.15980423982190510363772032690233,
    -0.17656376002524718647512421849032,  0.15162642913722598531903229617045,
    -0.090881384161410012831963755651079, 0.008674021313492586318780005883468,
     0.075581776473850090965407023747544,-0.14198837305251784063881548345787,
     0.17486335449663478165921413022414, -0.16644304831921567823839589426492,
     0.11871597273471929730747707847705, -0.042953233225881292913572018164494,
    -0.042953233225881292913572018164494, 0.11871597273471929730747707847705,
    -0.16644304831921567823839589426492,  0.17486335449663478165921413022414,
    -0.14198837305251784063881548345787,  0.075581776473850090965407023747544,
     0.008674021313492586318780005883468,-0.090881384161410012831963755651079,
     0.15162642913722598531903229617045, -0.17656376002524718647512421849032,
     0.15980423982190510363772032690233, -0.10530575443867740272410295104619,
     0.025938528373526445792454998016647, 0.059554274961645154658154180042717,
    -0.13098289131657380087121582271272,  0.17147891927418672456199547790670,
    -0.17147891927418672456199547790670,  0.13098289131657380087121582271272,
    -0.059554274961645154658154180042717,-0.025938528373526445792454998016647,
     0.10530575443867740272410295104619, -0.15980423982190510363772032690233,
     0.17656376002524718647512421849032, -0.15162642913722598531903229617045,
     0.090881384161410012831963755651079,-0.008674021313492586318780005883468,
    -0.075581776473850090965407023747544, 0.14198837305251784063881548345787,
    -0.17486335449663478165921413022414,  0.16644304831921567823839589426492,
    -0.11871597273471929730747707847705,  0.042953233225881292913572018164494
  },
  {
     0.038731993170326232824561787558492,-0.10875860742493672495725784718877,
     0.15790124012399101798155186502293, -0.17672345346106673069472498415600,
     0.16161097942121587956084560348376, -0.11546573663487784907576643982592,
     0.047148599859865922722408214465630, 0.030222077716717299469868553789363,
    -0.10178946920358000230129700949629,  0.15381110290879228162646885882803,
    -0.17629771118253266395745838340776,  0.16493138371356506569293701686233,
    -0.12189469877166149151050389002586,  0.055451621442613590037015233508646,
     0.021639354630747427065531392743213,-0.094575111235207343404082277365271,
     0.14935042127442479271101666423961, -0.17544725227590413843675978359671,
     0.16785445385626524004858887836470, -0.12803000590238956351843752239068,
     0.063621055212317597965816690486066, 0.013004500444088472936407629993685,
    -0.087132913540067685045147230048925, 0.14452994139365530889220759879835,
    -0.17417412557117862465636166724330,  0.17037314791731193194066488665533,
    -0.13385687753479470495813284521973,  0.071637220287468842367069460402068,
     0.004338317276799999540542326666634,-0.079480805024268148817432153524406,
     0.13936127622474077458193396572651, -0.17248139814210053549033101371775,
     0.17248139814210053549033101371775, -0.13936127622474077458193396572651,
     0.079480805024268148817432153524406,-0.004338317276799999540542326666634,
    -0.071637220287468842367069460402068, 0.13385687753479470495813284521973,
    -0.17037314791731193194066488665533,  0.17417412557117862465636166724330,
    -0.14452994139365530889220759879835,  0.087132913540067685045147230048925,
    -0.013004500444088472936407629993685,-0.063621055212317597965816690486066,
     0.12803000590238956351843752239068, -0.16785445385626524004858887836470,
     0.17544725227590413843675978359671, -0.14935042127442479271101666423961,
     0.094575111235207343404082277365271,-0.021639354630747427065531392743213,
    -0.055451621442613590037015233508646, 0.12189469877166149151050389002586,
    -0.16493138371356506569293701686233,  0.17629771118253266395745838340776,
    -0.15381110290879228162646885882803,  0.10178946920358000230129700949629,
    -0.030222077716717299469868553789363,-0.047148599859865922722408214465630,
     0.11546573663487784907576643982592, -0.16161097942121587956084560348376,
     0.17672345346106673069472498415600, -0.15790124012399101798155186502293,
     0.10875860742493672495725784718877, -0.038731993170326232824561787558492
  },
  {
     0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
     0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
     0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
     0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
     0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
     0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
     0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
     0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
     0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
     0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
     0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
     0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152, 0.034487422410367876541994695458672,
     0.034487422410367876541994695458672,-0.098211869798387772659737170957152,
     0.14698445030241983962180838807641, -0.17337998066526843272770239894580,
     0.17337998066526843272770239894580, -0.14698445030241983962180838807641,
     0.098211869798387772659737170957152,-0.034487422410367876541994695458672,
    -0.034487422410367876541994695458672, 0.098211869798387772659737170957152,
    -0.14698445030241983962180838807641,  0.17337998066526843272770239894580,
    -0.17337998066526843272770239894580,  0.14698445030241983962180838807641,
    -0.098211869798387772659737170957152, 0.034487422410367876541994695458672
  },
  {
     0.030222077716717299469868553789363,-0.087132913540067685045147230048925,
     0.13385687753479470495813284521973, -0.16493138371356506569293701686233,
     0.17672345346106673069472498415600, -0.16785445385626524004858887836470,
     0.13936127622474077458193396572651, -0.094575111235207343404082277365271,
     0.038731993170326232824561787558492, 0.021639354630747427065531392743213,
    -0.079480805024268148817432153524406, 0.12803000590238956351843752239068,
    -0.16161097942121587956084560348376,  0.17629771118253266395745838340776,
    -0.17037314791731193194066488665533,  0.14452994139365530889220759879835,
    -0.10178946920358000230129700949629,  0.047148599859865922722408214465630,
     0.013004500444088472936407629993685,-0.071637220287468842367069460402068,
     0.12189469877166149151050389002586, -0.15790124012399101798155186502293,
     0.17544725227590413843675978359671, -0.17248139814210053549033101371775,
     0.14935042127442479271101666423961, -0.10875860742493672495725784718877,
     0.055451621442613590037015233508646, 0.004338317276799999540542326666634,
    -0.063621055212317597965816690486066, 0.11546573663487784907576643982592,
    -0.15381110290879228162646885882803,  0.17417412557117862465636166724330,
    -0.17417412557117862465636166724330,  0.15381110290879228162646885882803,
    -0.11546573663487784907576643982592,  0.063621055212317597965816690486066,
    -0.004338317276799999540542326666634,-0.055451621442613590037015233508646,
     0.10875860742493672495725784718877, -0.14935042127442479271101666423961,
     0.17248139814210053549033101371775, -0.17544725227590413843675978359671,
     0.15790124012399101798155186502293, -0.12189469877166149151050389002586,
     0.071637220287468842367069460402068,-0.013004500444088472936407629993685,
    -0.047148599859865922722408214465630, 0.10178946920358000230129700949629,
    -0.14452994139365530889220759879835,  0.17037314791731193194066488665533,
    -0.17629771118253266395745838340776,  0.16161097942121587956084560348376,
    -0.12803000590238956351843752239068,  0.079480805024268148817432153524406,
    -0.021639354630747427065531392743213,-0.038731993170326232824561787558492,
     0.094575111235207343404082277365271,-0.13936127622474077458193396572651,
     0.16785445385626524004858887836470, -0.17672345346106673069472498415600,
     0.16493138371356506569293701686233, -0.13385687753479470495813284521973,
     0.087132913540067685045147230048925,-0.030222077716717299469868553789363
  },
  {
     0.025938528373526445792454998016647,-0.075581776473850090965407023747544,
     0.11871597273471929730747707847705, -0.15162642913722598531903229617045,
     0.17147891927418672456199547790670, -0.17656376002524718647512421849032,
     0.16644304831921567823839589426492, -0.14198837305251784063881548345787,
     0.10530575443867740272410295104619, -0.059554274961645154658154180042717,
     0.008674021313492586318780005883468, 0.042953233225881292913572018164494,
    -0.090881384161410012831963755651079, 0.13098289131657380087121582271272,
    -0.15980423982190510363772032690233,  0.17486335449663478165921413022414,
    -0.17486335449663478165921413022414,  0.15980423982190510363772032690233,
    -0.13098289131657380087121582271272,  0.090881384161410012831963755651079,
    -0.042953233225881292913572018164494,-0.008674021313492586318780005883468,
     0.059554274961645154658154180042717,-0.10530575443867740272410295104619,
     0.14198837305251784063881548345787, -0.16644304831921567823839589426492,
     0.17656376002524718647512421849032, -0.17147891927418672456199547790670,
     0.15162642913722598531903229617045, -0.11871597273471929730747707847705,
     0.075581776473850090965407023747544,-0.025938528373526445792454998016647,
    -0.025938528373526445792454998016647, 0.075581776473850090965407023747544,
    -0.11871597273471929730747707847705,  0.15162642913722598531903229617045,
    -0.17147891927418672456199547790670,  0.17656376002524718647512421849032,
    -0.16644304831921567823839589426492,  0.14198837305251784063881548345787,
    -0.10530575443867740272410295104619,  0.059554274961645154658154180042717,
    -0.008674021313492586318780005883468,-0.042953233225881292913572018164494,
     0.090881384161410012831963755651079,-0.13098289131657380087121582271272,
     0.15980423982190510363772032690233, -0.17486335449663478165921413022414,
     0.17486335449663478165921413022414, -0.15980423982190510363772032690233,
     0.13098289131657380087121582271272, -0.090881384161410012831963755651079,
     0.042953233225881292913572018164494, 0.008674021313492586318780005883468,
    -0.059554274961645154658154180042717, 0.10530575443867740272410295104619,
    -0.14198837305251784063881548345787,  0.16644304831921567823839589426492,
    -0.17656376002524718647512421849032,  0.17147891927418672456199547790670,
    -0.15162642913722598531903229617045,  0.11871597273471929730747707847705,
    -0.075581776473850090965407023747544, 0.025938528373526445792454998016647
  },
  {
     0.021639354630747427065531392743213,-0.063621055212317597965816690486066,
     0.10178946920358000230129700949629, -0.13385687753479470495813284521973,
     0.15790124012399101798155186502293, -0.17248139814210053549033101371775,
     0.17672345346106673069472498415600, -0.17037314791731193194066488665533,
     0.15381110290879228162646885882803, -0.12803000590238956351843752239068,
     0.094575111235207343404082277365271,-0.055451621442613590037015233508646,
     0.013004500444088472936407629993685, 0.030222077716717299469868553789363,
    -0.071637220287468842367069460402068, 0.10875860742493672495725784718877,
    -0.13936127622474077458193396572651,  0.16161097942121587956084560348376,
    -0.17417412557117862465636166724330,  0.17629771118253266395745838340776,
    -0.16785445385626524004858887836470,  0.14935042127442479271101666423961,
    -0.12189469877166149151050389002586,  0.087132913540067685045147230048925,
    -0.047148599859865922722408214465630, 0.004338317276799999540542326666634,
     0.038731993170326232824561787558492,-0.079480805024268148817432153524406,
     0.11546573663487784907576643982592, -0.14452994139365530889220759879835,
     0.16493138371356506569293701686233, -0.17544725227590413843675978359671,
     0.17544725227590413843675978359671, -0.16493138371356506569293701686233,
     0.14452994139365530889220759879835, -0.11546573663487784907576643982592,
     0.079480805024268148817432153524406,-0.038731993170326232824561787558492,
    -0.004338317276799999540542326666634, 0.047148599859865922722408214465630,
    -0.087132913540067685045147230048925, 0.12189469877166149151050389002586,
    -0.14935042127442479271101666423961,  0.16785445385626524004858887836470,
    -0.17629771118253266395745838340776,  0.17417412557117862465636166724330,
    -0.16161097942121587956084560348376,  0.13936127622474077458193396572651,
    -0.10875860742493672495725784718877,  0.071637220287468842367069460402068,
    -0.030222077716717299469868553789363,-0.013004500444088472936407629993685,
     0.055451621442613590037015233508646,-0.094575111235207343404082277365271,
     0.12803000590238956351843752239068, -0.15381110290879228162646885882803,
     0.17037314791731193194066488665533, -0.17672345346106673069472498415600,
     0.17248139814210053549033101371775, -0.15790124012399101798155186502293,
     0.13385687753479470495813284521973, -0.10178946920358000230129700949629,
     0.063621055212317597965816690486066,-0.021639354630747427065531392743213
  },
  {
     0.017327146149886432824466874566622,-0.051315565940294672644546154719392,
     0.083331957309718312162450688895359,-0.11214594829282953553133017365260,
     0.13665023337521968602987906462477, -0.15590312662333390407149878284971,
     0.16916475014679408478364306119571, -0.17592546719079780737825977787300,
     0.17592546719079780737825977787300, -0.16916475014679408478364306119571,
     0.15590312662333390407149878284971, -0.13665023337521968602987906462477,
     0.11214594829282953553133017365260, -0.083331957309718312162450688895359,
     0.051315565940294672644546154719392,-0.017327146149886432824466874566622,
    -0.017327146149886432824466874566622, 0.051315565940294672644546154719392,
    -0.083331957309718312162450688895359, 0.11214594829282953553133017365260,
    -0.13665023337521968602987906462477,  0.15590312662333390407149878284971,
    -0.16916475014679408478364306119571,  0.17592546719079780737825977787300,
    -0.17592546719079780737825977787300,  0.16916475014679408478364306119571,
    -0.15590312662333390407149878284971,  0.13665023337521968602987906462477,
    -0.11214594829282953553133017365260,  0.083331957309718312162450688895359,
    -0.051315565940294672644546154719392, 0.017327146149886432824466874566622,
     0.017327146149886432824466874566622,-0.051315565940294672644546154719392,
     0.083331957309718312162450688895359,-0.11214594829282953553133017365260,
     0.13665023337521968602987906462477, -0.15590312662333390407149878284971,
     0.16916475014679408478364306119571, -0.17592546719079780737825977787300,
     0.17592546719079780737825977787300, -0.16916475014679408478364306119571,
     0.15590312662333390407149878284971, -0.13665023337521968602987906462477,
     0.11214594829282953553133017365260, -0.083331957309718312162450688895359,
     0.051315565940294672644546154719392,-0.017327146149886432824466874566622,
    -0.017327146149886432824466874566622, 0.051315565940294672644546154719392,
    -0.083331957309718312162450688895359, 0.11214594829282953553133017365260,
    -0.13665023337521968602987906462477,  0.15590312662333390407149878284971,
    -0.16916475014679408478364306119571,  0.17592546719079780737825977787300,
    -0.17592546719079780737825977787300,  0.16916475014679408478364306119571,
    -0.15590312662333390407149878284971,  0.13665023337521968602987906462477,
    -0.11214594829282953553133017365260,  0.083331957309718312162450688895359,
    -0.051315565940294672644546154719392, 0.017327146149886432824466874566622
  },
  {
     0.013004500444088472936407629993685,-0.038731993170326232824561787558492,
     0.063621055212317597965816690486066,-0.087132913540067685045147230048925,
     0.10875860742493672495725784718877, -0.12803000590238956351843752239068,
     0.14452994139365530889220759879835, -0.15790124012399101798155186502293,
     0.16785445385626524004858887836470, -0.17417412557117862465636166724330,
     0.17672345346106673069472498415600, -0.17544725227590413843675978359671,
     0.17037314791731193194066488665533, -0.16161097942121587956084560348376,
     0.14935042127442479271101666423961, -0.13385687753479470495813284521973,
     0.11546573663487784907576643982592, -0.094575111235207343404082277365271,
     0.071637220287468842367069460402068,-0.047148599859865922722408214465630,
     0.021639354630747427065531392743213, 0.004338317276799999540542326666634,
    -0.030222077716717299469868553789363, 0.055451621442613590037015233508646,
    -0.079480805024268148817432153524406, 0.10178946920358000230129700949629,
    -0.12189469877166149151050389002586,  0.13936127622474077458193396572651,
    -0.15381110290879228162646885882803,  0.16493138371356506569293701686233,
    -0.17248139814210053549033101371775,  0.17629771118253266395745838340776,
    -0.17629771118253266395745838340776,  0.17248139814210053549033101371775,
    -0.16493138371356506569293701686233,  0.15381110290879228162646885882803,
    -0.13936127622474077458193396572651,  0.12189469877166149151050389002586,
    -0.10178946920358000230129700949629,  0.079480805024268148817432153524406,
    -0.055451621442613590037015233508646, 0.030222077716717299469868553789363,
    -0.004338317276799999540542326666634,-0.021639354630747427065531392743213,
     0.047148599859865922722408214465630,-0.071637220287468842367069460402068,
     0.094575111235207343404082277365271,-0.11546573663487784907576643982592,
     0.13385687753479470495813284521973, -0.14935042127442479271101666423961,
     0.16161097942121587956084560348376, -0.17037314791731193194066488665533,
     0.17544725227590413843675978359671, -0.17672345346106673069472498415600,
     0.17417412557117862465636166724330, -0.16785445385626524004858887836470,
     0.15790124012399101798155186502293, -0.14452994139365530889220759879835,
     0.12803000590238956351843752239068, -0.10875860742493672495725784718877,
     0.087132913540067685045147230048925,-0.063621055212317597965816690486066,
     0.038731993170326232824561787558492,-0.013004500444088472936407629993685
  },
  {
     0.008674021313492586318780005883468,-0.025938528373526445792454998016647,
     0.042953233225881292913572018164494,-0.059554274961645154658154180042717,
     0.075581776473850090965407023747544,-0.090881384161410012831963755651079,
     0.10530575443867740272410295104619, -0.11871597273471929730747707847705,
     0.13098289131657380087121582271272, -0.14198837305251784063881548345787,
     0.15162642913722598531903229617045, -0.15980423982190510363772032690233,
     0.16644304831921567823839589426492, -0.17147891927418672456199547790670,
     0.17486335449663478165921413022414, -0.17656376002524718647512421849032,
     0.17656376002524718647512421849032, -0.17486335449663478165921413022414,
     0.17147891927418672456199547790670, -0.16644304831921567823839589426492,
     0.15980423982190510363772032690233, -0.15162642913722598531903229617045,
     0.14198837305251784063881548345787, -0.13098289131657380087121582271272,
     0.11871597273471929730747707847705, -0.10530575443867740272410295104619,
     0.090881384161410012831963755651079,-0.075581776473850090965407023747544,
     0.059554274961645154658154180042717,-0.042953233225881292913572018164494,
     0.025938528373526445792454998016647,-0.008674021313492586318780005883468,
    -0.008674021313492586318780005883468, 0.025938528373526445792454998016647,
    -0.042953233225881292913572018164494, 0.059554274961645154658154180042717,
    -0.075581776473850090965407023747544, 0.090881384161410012831963755651079,
    -0.10530575443867740272410295104619,  0.11871597273471929730747707847705,
    -0.13098289131657380087121582271272,  0.14198837305251784063881548345787,
    -0.15162642913722598531903229617045,  0.15980423982190510363772032690233,
    -0.16644304831921567823839589426492,  0.17147891927418672456199547790670,
    -0.17486335449663478165921413022414,  0.17656376002524718647512421849032,
    -0.17656376002524718647512421849032,  0.17486335449663478165921413022414,
    -0.17147891927418672456199547790670,  0.16644304831921567823839589426492,
    -0.15980423982190510363772032690233,  0.15162642913722598531903229617045,
    -0.14198837305251784063881548345787,  0.13098289131657380087121582271272,
    -0.11871597273471929730747707847705,  0.10530575443867740272410295104619,
    -0.090881384161410012831963755651079, 0.075581776473850090965407023747544,
    -0.059554274961645154658154180042717, 0.042953233225881292913572018164494,
    -0.025938528373526445792454998016647, 0.008674021313492586318780005883468
  },
  {
     0.004338317276799999540542326666634,-0.013004500444088472936407629993685,
     0.021639354630747427065531392743213,-0.030222077716717299469868553789363,
     0.038731993170326232824561787558492,-0.047148599859865922722408214465630,
     0.055451621442613590037015233508646,-0.063621055212317597965816690486066,
     0.071637220287468842367069460402068,-0.079480805024268148817432153524406,
     0.087132913540067685045147230048925,-0.094575111235207343404082277365271,
     0.10178946920358000230129700949629, -0.10875860742493672495725784718877,
     0.11546573663487784907576643982592, -0.12189469877166149151050389002586,
     0.12803000590238956351843752239068, -0.13385687753479470495813284521973,
     0.13936127622474077458193396572651, -0.14452994139365530889220759879835,
     0.14935042127442479271101666423961, -0.15381110290879228162646885882803,
     0.15790124012399101798155186502293, -0.16161097942121587956084560348376,
     0.16493138371356506569293701686233, -0.16785445385626524004858887836470,
     0.17037314791731193194066488665533, -0.17248139814210053549033101371775,
     0.17417412557117862465636166724330, -0.17544725227590413843675978359671,
     0.17629771118253266395745838340776, -0.17672345346106673069472498415600,
     0.17672345346106673069472498415600, -0.17629771118253266395745838340776,
     0.17544725227590413843675978359671, -0.17417412557117862465636166724330,
     0.17248139814210053549033101371775, -0.17037314791731193194066488665533,
     0.16785445385626524004858887836470, -0.16493138371356506569293701686233,
     0.16161097942121587956084560348376, -0.15790124012399101798155186502293,
     0.15381110290879228162646885882803, -0.14935042127442479271101666423961,
     0.14452994139365530889220759879835, -0.13936127622474077458193396572651,
     0.13385687753479470495813284521973, -0.12803000590238956351843752239068,
     0.12189469877166149151050389002586, -0.11546573663487784907576643982592,
     0.10875860742493672495725784718877, -0.10178946920358000230129700949629,
     0.094575111235207343404082277365271,-0.087132913540067685045147230048925,
     0.079480805024268148817432153524406,-0.071637220287468842367069460402068,
     0.063621055212317597965816690486066,-0.055451621442613590037015233508646,
     0.047148599859865922722408214465630,-0.038731993170326232824561787558492,
     0.030222077716717299469868553789363,-0.021639354630747427065531392743213,
     0.013004500444088472936407629993685,-0.004338317276799999540542326666634
  }
};

static const od_dct_basis_row *const OD_DCT_BASES[OD_NBSIZES] = {
  OD_DCT4_BASIS,
  OD_DCT8_BASIS,
  OD_DCT16_BASIS,
  OD_DCT32_BASIS,
  OD_DCT64_BASIS
};

static void idct(double x[], int bszi, const double y[]) {
  double t[OD_BSIZE_MAX];
  const od_dct_basis_row *basis;
  int i;
  int j;
  int n;
  basis = OD_DCT_BASES[bszi];
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (j = 0; j < n; j++) {
    t[j] = 0;
    for (i = 0; i < n; i++) t[j] += basis[i][j]*y[i];
  }
  /* workaround to prevent mis-compilation by clang-3.4 */
  OD_MOVE(x, t, n);
}

static void ieee1180_print_results(long sumerrs[OD_BSIZE_MAX][OD_BSIZE_MAX],
 long sumsqerrs[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int maxerr[OD_BSIZE_MAX][OD_BSIZE_MAX], int l, int h, int sign, int bszi) {
  double max;
  double total;
  int m;
  int i;
  int j;
  int n;
  printf("IEEE1180-1990 test results:\n");
  printf("Input range: [%i,%i]\n", -l, h);
  printf("Sign: %i\n", sign);
  printf("Iterations: %i\n\n", IEEE1180_NBLOCKS);
  printf("Peak absolute values of errors:\n");
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0, m = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (maxerr[i][j] > m) m = maxerr[i][j];
      printf("%4i", maxerr[i][j]);
    }
    printf("\n");
  }
  printf("Worst peak error = %i (%s spec limit 1)\n\n", m,
   ieee1180_meets((double)m, 1.0));
  printf("Mean square errors:\n");
  max = total = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      double err;
      err = sumsqerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += sumsqerrs[i][j];
      if (max < err) max = err;
    }
    printf("\n");
  }
  printf("Worst pmse = %.6f (%s spec limit 0.06)\n", max,
   ieee1180_meets(max, 0.06));
  total /= n*n*(double)IEEE1180_NBLOCKS;
  printf("Overall mse = %.6f (%s spec limit 0.02)\n\n", total,
   ieee1180_meets(total, 0.02));
  printf("Mean errors:\n");
  max = total = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      double err;
      err = sumerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += sumerrs[i][j];
      if (err < 0) err = -err;
      if (max < err) max = err;
    }
    printf("\n");
  }
  total /= n*n*(double)IEEE1180_NBLOCKS;
  printf("Worst mean error = %.6f (%s spec limit 0.015)\n", max,
   ieee1180_meets(max, 0.015));
  printf("Overall mean error = %.6f (%s spec limit 0.0015)\n\n", total,
   ieee1180_meets(total, 0.0015));
}

static const od_dct_func_2d *test_fdct_2d;
static const od_dct_func_2d *test_idct_2d;

static void ieee1180_test_block(long sumerrs[OD_BSIZE_MAX][OD_BSIZE_MAX],
 long sumsqerrs[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int maxerr[OD_BSIZE_MAX][OD_BSIZE_MAX], int l, int h, int sign, int bszi) {
  od_coeff block[OD_BSIZE_MAX][OD_BSIZE_MAX];
  od_coeff refcoefs[OD_BSIZE_MAX][OD_BSIZE_MAX];
  od_coeff refout[OD_BSIZE_MAX][OD_BSIZE_MAX];
  od_coeff testout[OD_BSIZE_MAX][OD_BSIZE_MAX];
  double floatcoefs[OD_BSIZE_MAX][OD_BSIZE_MAX];
  int global_maxerr;
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      block[i][j] = ieee1180_rand(l, h)*sign*OD_COEFF_SCALE;
    }
  }
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  (*test_fdct_2d[bszi])(refcoefs[0], OD_BSIZE_MAX, block[0], OD_BSIZE_MAX);
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      floatcoefs[i][j] = refcoefs[i][j]/(double)OD_COEFF_SCALE;
    }
  }
  for (i = 0; i < n; i++) idct(floatcoefs[i], bszi, floatcoefs[i]);
  for (j = 0; j < n; j++) {
    double x[OD_BSIZE_MAX];
    for (i = 0; i < n; i++) x[i] = floatcoefs[i][j];
    idct(x, bszi, x);
    for (i = 0; i < n; i++) floatcoefs[i][j] = x[i];
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      refout[i][j] = (od_coeff)floor(floatcoefs[i][j] + 0.5);
      refout[i][j] = OD_CLAMPI(-256, refout[i][j], 255);
    }
  }
  (*test_idct_2d[bszi])(testout[0], OD_BSIZE_MAX, refcoefs[0], OD_BSIZE_MAX);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (testout[i][j] != block[i][j]) {
        od_exit_code = EXIT_FAILURE;
        printf("Mismatch:\n");
        printf("in:\n");
        for (i = 0; i < n; i++) {
          printf("       ");
          for (j = 0; j < n; j++) printf(" %i", block[i][j]);
          printf("\n");
        }
        printf("xform:\n");
        for (i = 0; i < n; i++) {
          printf("       ");
          for (j = 0; j < n; j++) printf(" %i", refcoefs[i][j]);
          printf("\n");
        }
        printf("out:\n");
        for (i = 0; i < n; i++) {
          printf("       ");
          for (j = 0; j < n; j++) printf(" %i", testout[i][j]);
          printf("\n");
        }
        printf("\n");
        return;
      }
      testout[i][j] = OD_CLAMPI(-256, OD_COEFF_UNSCALE(testout[i][j]), 255);
    }
  }
  global_maxerr = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      int err;
      err = testout[i][j] - refout[i][j];
      sumerrs[i][j] += err;
      sumsqerrs[i][j] += err*err;
      if (err < 0) err = -err;
      maxerr[i][j] = OD_MAXI(maxerr[i][j], err);
      global_maxerr = OD_MAXI(global_maxerr, err);
    }
  }
  /*if (global_maxerr > 1) {
    int u;
    int v;
    printf("Excessive peak error: %i\n", global_maxerr);
    printf("Input:\n");
    for (u = 0; u < n; u++) {
      for (v = 0;v < n; v++) {
        printf("%5i", block[u][v]);
      }
      printf("\n");
    }
    printf("Forward transform coefficients:\n");
    for (u = 0; u < n; u++) {
      for (v = 0; v < n; v++) {
        printf("%5i", refcoefs[u][v]);
      }
      printf("\n");
    }
    printf("Reference inverse:\n");
    for (u = 0; u < n; u++) {
      for (v = 0; v < n; v++) {
        printf("%5i", refout[u][v]);
      }
      printf("\n");
    }
    printf("Integerized inverse:\n");
    for (u = 0; u < n; u++) {
      for (v = 0; v < n; v++) {
        printf("%5i", testout[u][v]);
      }
      printf("\n");
    }
  }*/
}

static void ieee1180_test(int bszi) {
  long sumerrs[OD_BSIZE_MAX][OD_BSIZE_MAX];
  long sumsqerrs[OD_BSIZE_MAX][OD_BSIZE_MAX];
  int maxerr[OD_BSIZE_MAX][OD_BSIZE_MAX];
  int i;
  int j;
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], 1, bszi);
    }
    ieee1180_print_results(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], 1, bszi);
  }
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], -1, bszi);
    }
    ieee1180_print_results(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], -1, bszi);
  }
}

static void print_basis(double basis[OD_BSIZE_MAX][OD_BSIZE_MAX], int bszi) {
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf("%8.5f%c", basis[i][j], j == n - 1 ? '\n' : ' ');
    }
  }
}

static void compute_fbasis(double basis[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int bszi) {
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    od_coeff x[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) x[j] = (i == j) << 8 ;
    (*OD_FDCT_1D[bszi])(x, x, 1);
    for (j = 0; j < n; j++) basis[j][i] = x[j]/256.0;
  }
}

static void compute_ibasis(double basis[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int bszi) {
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    od_coeff x[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) x[j] = (i == j) << 8;
    (*OD_IDCT_1D[bszi])(x, 1, x);
    for (j = 0; j < n; j++) basis[j][i] = x[j]/256.0;
  }
}

static void compute_ftrue_basis(double basis[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int bszi) {
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (j = 0; j < n; j++) {
    for (i = 0; i < n; i++) {
      basis[j][i] = sqrt(2.0/n)*cos((i + 0.5)*j*M_PI/n);
      if (j == 0) basis[j][i] *= M_SQRT1_2;
    }
  }
}

static void compute_itrue_basis(double basis[OD_BSIZE_MAX][OD_BSIZE_MAX],
 int bszi) {
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    double x[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) x[j] = (i == j);
    idct(x, bszi, x);
    for (j = 0; j < n; j++) basis[j][i] = x[j];
  }
}

static double compute_mse(double basis[OD_BSIZE_MAX][OD_BSIZE_MAX],
 double tbasis[OD_BSIZE_MAX][OD_BSIZE_MAX], int bszi) {
  double e[OD_BSIZE_MAX][OD_BSIZE_MAX];
  double ret;
  int i;
  int j;
  int k;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      e[i][j] = 0;
      for (k = 0; k < n; k++) {
        e[i][j] += (basis[i][k] - tbasis[i][k])*AUTOCORR[k - j + 63];
      }
    }
  }
  ret = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      ret += e[i][j]*(basis[i][j] - tbasis[i][j]);
    }
  }
  return ret/n;
}

static void check_bias(int bszi) {
  double rtacc[OD_BSIZE_MAX];
  double facc[OD_BSIZE_MAX];
  double q8acc[OD_BSIZE_MAX];
  double q7acc[OD_BSIZE_MAX];
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (j = 0; j < n; j++) q7acc[j] = q8acc[j] = facc[j] = rtacc[j] = 0;
  ieee1180_srand(1);
  for (i = 0; i < 10000000; i++) {
    od_coeff x[OD_BSIZE_MAX];
    od_coeff x2[OD_BSIZE_MAX];
    od_coeff y[OD_BSIZE_MAX];
    od_coeff y2[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) x[j] = ieee1180_rand(255, 255)*OD_COEFF_SCALE;
    (*OD_FDCT_1D[bszi])(y, x, 1);
    for (j = 0; j < n; j++) facc[j] += y[j];
    (*OD_IDCT_1D[bszi])(x2, 1, y);
    for (j = 0; j < n; j++) {
      if (x[j] != x2[j]) {
        od_exit_code = EXIT_FAILURE;
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < n; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < n; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < n; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
    for (j = 0; j < n; j++) x[j] = OD_COEFF_UNSCALE(x[j]);
    for (j = 0; j < n; j++) {
      y2[j] = y[j] + ieee1180_rand(1, 1)*OD_COEFF_SCALE;
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      rtacc[j] += x2[j] - x[j];
    }
    for (j = 0; j < n; j++) {
      y2[j] = (y[j] + (y[j] < 0 ? -4 : 4)*OD_COEFF_SCALE)/
       (8 << OD_COEFF_SHIFT)*(1 << (3 + OD_COEFF_SHIFT));
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      q8acc[j] += x2[j] - x[j];
    }
    for (j = 0; j < n; j++) {
      y2[j] = (y[j] + (y[j] < 0 ? -7 : 7)*OD_COEFF_SCALE/2)/
       (7 << OD_COEFF_SHIFT)*(7 << OD_COEFF_SHIFT);
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      q7acc[j] += x2[j] - x[j];
    }
  }
  printf("1-D Forward Bias:\n");
  for (j = 0; j < n; j++) {
    printf("% -18.15G%s", facc[j]/(i << OD_COEFF_SHIFT),
     (j & 3) == 3 ? "\n" : "  ");
  }
  printf("\n");
  printf("1-D Round-Trip Bias:\n");
  for (j = 0; j < n; j++) {
    printf("% -18.15G%s", rtacc[j]/i, (j & 3) == 3 ? "\n" : "  ");
  }
  printf("\n");
  printf("1-D Q=8 Bias:\n");
  for (j = 0; j < n; j++) {
    printf("% -18.15G%s", q8acc[j]/i, (j & 3) == 3 ? "\n" : "  ");
  }
  printf("\n");
  printf("1-D Q=7 Bias:\n");
  for (j = 0; j < n; j++) {
    printf("% -18.15G%s", q7acc[j]/i, (j & 3) == 3 ? "\n" : "  ");
  }
  printf("\n");
}

# if defined(OD_DCT_CHECK_OVERFLOW)

/*The total number of overflow checks we do.
  TODO: If we arranged the indices of these checks a bit better, we could save
   a bunch of work in dynamic_range().*/
#  define OD_DCT_NCHECK_OVERFLOW (388)
/*The overflow checks in the following ranges apply to the pre-/post- filters,
   and not the DCTs themselves.*/
#  define OD_FILTER_CHECK_OVERFLOW_MIN (51)
#  define OD_FILTER_CHECK_OVERFLOW_MAX (133)

od_coeff od_dct_check_min[OD_DCT_NCHECK_OVERFLOW];
od_coeff od_dct_check_max[OD_DCT_NCHECK_OVERFLOW];

/* Test up to 8-point filters to avoid a redesign if we want to turn them back
    on. */
#  define OD_MAX_FILT_SIZE (1)

static void od_bin_prefilter_columns_2d(
 od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2], int n, int f, int o) {
  od_coeff y[OD_BSIZE_MAX*2];
  int u;
  int v;
  /*Perform pre-filtering.*/
  for (v = 0; v < n*2; v++) {
    for (u = 0; u < n*2; u++) y[u] = x[u][v];
    (*OD_PRE_FILTER[f])(y + o, y + o);
    (*OD_PRE_FILTER[f])(y + n + o, y + n + o);
    for (u = 0; u < n*2; u++) x[u][v] = y[u];
  }
}

static void od_bin_fxform_2d(od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2],
 int bszi, int f) {
  od_coeff y[OD_BSIZE_MAX];
  int u;
  int v;
  int n;
  int fn;
  int o;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  fn = 1 << (OD_LOG_BSIZE0 + f);
  o = (n >> 1) - (fn >> 1);
  /*Perform pre-filtering.*/
  od_bin_prefilter_columns_2d(x, n, f, o);
  for (u = 0; u < n*2; u++) {
    (*OD_PRE_FILTER[f])(x[u] + o, x[u] + o);
    (*OD_PRE_FILTER[f])(x[u] + n + o, x[u] + n + o);
  }
  /*Perform DCT.*/
  for (u = n >> 1; u < n*3 >> 1; u++) {
    (*OD_FDCT_1D[bszi])(x[u] + (n >> 1), x[u] + (n >> 1), 1);
  }
  for (v = n >> 1; v < n*3 >> 1; v++) {
    for (u = n >> 1; u < n*3 >> 1; u++) y[u - (n >> 1)] = x[u][v];
    (*OD_FDCT_1D[bszi])(y, y, 1);
    for (u = n >> 1; u < n*3 >> 1; u++) x[u][v] = y[u - (n >> 1)];
  }
}

static void dynamic_range(int bszi) {
  int f;
  /*Test with all filters up to OD_MAX_FILT_SIZE.*/
  for (f = 0; f <= OD_MINI(bszi, OD_MAX_FILT_SIZE); f++) {
    static double
     basis2[OD_BSIZE_MAX][OD_BSIZE_MAX][OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
    static int8_t check_basis[OD_DCT_NCHECK_OVERFLOW][OD_BSIZE_MAX*2 + 1]
     [OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
    od_coeff old_check_min[OD_DCT_NCHECK_OVERFLOW];
    od_coeff old_check_max[OD_DCT_NCHECK_OVERFLOW];
    od_coeff min2[OD_BSIZE_MAX][OD_BSIZE_MAX];
    od_coeff max2[OD_BSIZE_MAX][OD_BSIZE_MAX];
    int i;
    int j;
    int k;
    int u;
    int v;
    int n;
    memcpy(old_check_min, od_dct_check_min, sizeof(od_dct_check_min));
    memcpy(old_check_max, od_dct_check_max, sizeof(od_dct_check_max));
    n = 1 << (OD_LOG_BSIZE0 + bszi);
    for (i = 0; i < n*2; i++) {
      for (j = 0; j < n*2; j++) {
        od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
        od_coeff y[OD_BSIZE_MAX];
        int fn;
        int o;
        /*Generate impulse.*/
        for (u = 0; u < n*2; u++) {
          for (v = 0; v < n*2; v++) {
            x[u][v] = (u == i && v == j) << (8 + OD_COEFF_SHIFT);
          }
        }
        memset(od_dct_check_min, 0, sizeof(od_dct_check_min));
        memset(od_dct_check_max, 0, sizeof(od_dct_check_max));
        fn = 1 << (OD_LOG_BSIZE0 + f);
        o = (n >> 1) - (fn >> 1);
        od_bin_prefilter_columns_2d(x, n, f, o);
        /*Retrieve basis elements for each row filter overflow check.
          Exactly one filter should have a non-zero response.*/
        for (k = OD_FILTER_CHECK_OVERFLOW_MIN;
         k < OD_FILTER_CHECK_OVERFLOW_MAX; k++) {
          check_basis[j][n*2][i][j] =
           (int8_t)(od_dct_check_min[k] < -od_dct_check_max[k] ? -1 : 1);
        }
        /*Prefilter rows.*/
        for (u = 0; u < n*2; u++) {
          memset(od_dct_check_min, 0, sizeof(od_dct_check_min));
          memset(od_dct_check_max, 0, sizeof(od_dct_check_max));
          (*OD_PRE_FILTER[f])(x[u] + o, x[u] + o);
          (*OD_PRE_FILTER[f])(x[u] + n + o, x[u] + n + o);
          /*Retrieve basis elements for each column filter overflow check.
            We compute a separate basis for each column to ensure that each
             check only has one non-zero input, which guarantees we get a
             basis from a consistent set of inputs (though in practice this
             probably does not matter).*/
          for (k = OD_FILTER_CHECK_OVERFLOW_MIN;
           k < OD_FILTER_CHECK_OVERFLOW_MAX; k++) {
            check_basis[k][u][i][j] =
             (int8_t)(od_dct_check_min[k] < -od_dct_check_max[k] ? -1 : 1);
          }
        }
        /*Transform rows.*/
        for (u = n >> 1; u < n*3 >> 1; u++) {
          memset(od_dct_check_min, 0, sizeof(od_dct_check_min));
          memset(od_dct_check_max, 0, sizeof(od_dct_check_max));
          (*OD_FDCT_1D[bszi])(x[u] + (n >> 1), x[u] + (n >> 1), 1);
          /*Retrieve basis elements for each row transform overflow check.*/
          for (k = 0; k < OD_FILTER_CHECK_OVERFLOW_MIN; k++) {
            check_basis[k][u][i][j] =
             (int8_t)(od_dct_check_min[k] < -od_dct_check_max[k] ? -1 : 1);
          }
          for (k = OD_FILTER_CHECK_OVERFLOW_MAX;
           k < OD_DCT_NCHECK_OVERFLOW; k++) {
            check_basis[k][u][i][j] =
             (int8_t)(od_dct_check_min[k] < -od_dct_check_max[k] ? -1 : 1);
          }
        }
        /*Transform columns.*/
        for (v = n >> 1; v < n*3 >> 1; v++) {
          memset(od_dct_check_min, 0, sizeof(od_dct_check_min));
          memset(od_dct_check_max, 0, sizeof(od_dct_check_max));
          for (u = n >> 1; u < n*3 >> 1; u++) y[u - (n >> 1)] = x[u][v];
          (*OD_FDCT_1D[bszi])(y, y, 1);
          for (u = n >> 1; u < n*3 >> 1; u++) x[u][v] = y[u - (n >> 1)];
          /*Retrieve basis elements for each column transform overflow
             check.*/
          for (k = 0; k < OD_FILTER_CHECK_OVERFLOW_MIN; k++) {
            check_basis[k][n + v - (n >> 1)][i][j] =
             (int8_t)(od_dct_check_min[k] < -od_dct_check_max[k] ? -1 : 1);
          }
          for (k = OD_FILTER_CHECK_OVERFLOW_MAX;
           k < OD_DCT_NCHECK_OVERFLOW; k++) {
            check_basis[k][n + v - (n >> 1)][i][j] =
             (int8_t)(od_dct_check_min[k] < -od_dct_check_max[k] ? -1 : 1);
          }
        }
        /*Retrieve output basis elements.*/
        for (u = 0; u < n; u++) {
          for (v = 0; v < n; v++) {
            basis2[u][v][i][j] =
             x[u + (n >> 1)][v + (n >> 1)]/(256.0*OD_COEFF_SCALE);
          }
        }
      }
    }
    /*We don't bother to update old_check_min/old_check_max while doing the
       transforms above, because the ones we're about to do should strictly
       exceed anything they encounter.*/
    memcpy(od_dct_check_min, old_check_min, sizeof(od_dct_check_min));
    memcpy(od_dct_check_max, old_check_max, sizeof(od_dct_check_max));
    /*Maximize/minimize each overflow check.*/
    for (k = 0; k < OD_DCT_NCHECK_OVERFLOW; k++) {
      int filter_k;
      filter_k = (k >= OD_FILTER_CHECK_OVERFLOW_MIN
       && k < OD_FILTER_CHECK_OVERFLOW_MAX);
      for (u = 0; u < 2*n + filter_k; u++) {
        od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
        for (i = 0; i < n*2; i++) {
          for (j = 0; j < n*2; j++) {
            x[i][j] =
             (check_basis[k][u][i][j] < 0 ? -255 : 255)*OD_COEFF_SCALE;
          }
        }
        od_bin_fxform_2d(x, bszi, f);
        for (i = 0; i < n*2; i++) {
          for (j = 0; j < n*2; j++) {
            x[i][j] =
             (check_basis[k][u][i][j] > 0 ? -255 : 255)*OD_COEFF_SCALE;
          }
        }
        od_bin_fxform_2d(x, bszi, f);
      }
    }
    for (u = 0; u < n; u++) {
      for (v = 0; v < n; v++) {
        od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
        for (i = 0; i < n*2; i++) {
          for (j = 0; j < n*2; j++) {
            x[i][j] = (basis2[u][v][i][j] < 0 ? -255 : 255)*OD_COEFF_SCALE;
          }
        }
        od_bin_fxform_2d(x, bszi, f);
        max2[u][v] = x[u + (n >> 1)][v + (n >> 1)];
        for (i = 0; i < n*2; i++) {
          for (j = 0; j < n*2; j++) {
            x[i][j] = (basis2[u][v][i][j] > 0 ? -255 : 255)*OD_COEFF_SCALE;
          }
        }
        od_bin_fxform_2d(x, bszi, f);
        min2[u][v] = x[u + (n >> 1)][v + (n >> 1)];
      }
    }
    printf("2-D scaled, prefiltered ranges (filter size = %i):\n",
     1 << (OD_LOG_BSIZE0 + f));
    for (u = 0; u < n; u++) {
      printf(" Min %2i:", u);
      for (v = 0; v < n; v++) printf(" %7i", min2[u][v]);
      printf("\n Max %2i:", u);
      for (v = 0; v < n; v++) printf(" %7i", max2[u][v]);
      printf("\n");
    }
  }
}

# endif

static void check_transform(int bszi) {
  od_coeff min[OD_BSIZE_MAX];
  od_coeff max[OD_BSIZE_MAX];
  double basis[OD_BSIZE_MAX][OD_BSIZE_MAX];
  double tbasis[OD_BSIZE_MAX][OD_BSIZE_MAX];
  int i;
  int j;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n; i++) {
    od_coeff x[OD_BSIZE_MAX];
    od_coeff y[OD_BSIZE_MAX];
    od_coeff x2[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) {
      x[j] = (i == j) << (8 + OD_COEFF_SHIFT);
    }
    OD_FDCT_1D[bszi](y, x, 1);
    OD_IDCT_1D[bszi](x2, 1, y);
    for (j = 0; j < n; j++) {
      basis[i][j] = y[j]/(256.0*OD_COEFF_SCALE);
    }
    for (j = 0; j < n; j++) {
      if (x[j] != x2[j]) {
        od_exit_code = EXIT_FAILURE;
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < n; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < n; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < n; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
  }
  for (i = 0; i < n; i++) {
    od_coeff x[OD_BSIZE_MAX];
    od_coeff y[OD_BSIZE_MAX];
    od_coeff x2[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) {
      x[j] = basis[j][i] < 0 ? -256 : 255;
    }
    OD_FDCT_1D[bszi](y, x, 1);
    OD_IDCT_1D[bszi](x2, 1, y);
    max[i] = y[i];
    for (j = 0; j < n; j++) {
      if (x[j] != x2[j]) {
        od_exit_code = EXIT_FAILURE;
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < n; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < n; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < n; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
    for (j = 0; j < n; j++) {
      x[j] = basis[j][i] > 0 ? -256 : 255;
    }
    OD_FDCT_1D[bszi](y, x, 1);
    OD_IDCT_1D[bszi](x2, 1, y);
    min[i] = y[i];
    for (j = 0; j < n; j++) {
      if (x[j] != x2[j]) {
        od_exit_code = EXIT_FAILURE;
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < n; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < n; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < n; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
  }
  printf("1-D ranges:\n");
  printf(" Min:");
  for (j = 0; j < n; j++) printf(" %5i", min[j]);
  printf("\n Max:");
  for (j = 0; j < n; j++) printf(" %5i", max[j]);
  printf("\n");
# if defined(OD_DCT_CHECK_OVERFLOW)
  dynamic_range(bszi);
# endif
  printf("od_bin_idct%i basis:\n", n);
  compute_ibasis(basis, bszi);
  print_basis(basis, bszi);
  printf("Scaled type-II iDCT basis:\n");
  compute_itrue_basis(tbasis, bszi);
  print_basis(tbasis, bszi);
  printf("\nod_bin_fdct%i basis:\n", n);
  compute_fbasis(basis, bszi);
  print_basis(basis, bszi);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis(tbasis, bszi);
  print_basis(tbasis, bszi);
  printf("MSE: %.32g\n\n", compute_mse(basis, tbasis, bszi));
  ieee1180_test(bszi);
  check_bias(bszi);
}

void run_test(void) {
  int bszi;
  for (bszi = 0; bszi <= OD_BLOCK_64X64; bszi++) check_transform(bszi);
}

int main(void) {
  test_fdct_2d = OD_FDCT_2D_C;
  test_idct_2d = OD_IDCT_2D_C;
  printf("Testing unoptimized DCT...\n");
  run_test();
#if defined(OD_X86ASM)
# if defined(OD_SSE2_INTRINSICS)
  if (od_cpu_flags_get() & OD_CPU_X86_SSE2) {
    static const od_dct_func_2d OD_FDCT_2D_SSE2[OD_NBSIZES + 1] = {
      od_bin_fdct4x4_sse2,
      od_bin_fdct8x8,
      od_bin_fdct16x16,
      od_bin_fdct32x32,
      od_bin_fdct64x64
    };
    static const od_dct_func_2d OD_IDCT_2D_SSE2[OD_NBSIZES + 1] = {
      od_bin_idct4x4_sse2,
      od_bin_idct8x8,
      od_bin_idct16x16,
      od_bin_idct32x32,
      od_bin_idct64x64
    };
    test_fdct_2d = OD_FDCT_2D_SSE2;
    test_idct_2d = OD_IDCT_2D_SSE2;
    printf("Testing SSE2 DCT...\n");
    run_test();
  }
# endif
# if defined(OD_SSE41_INTRINSICS)
  if (od_cpu_flags_get() & OD_CPU_X86_SSE4_1) {
    static const od_dct_func_2d OD_FDCT_2D_SSE41[OD_NBSIZES + 1] = {
      od_bin_fdct4x4_sse41,
      od_bin_fdct8x8,
      od_bin_fdct16x16,
      od_bin_fdct32x32,
      od_bin_fdct64x64
    };
    static const od_dct_func_2d OD_IDCT_2D_SSE41[OD_NBSIZES + 1] = {
      od_bin_idct4x4_sse41,
      od_bin_idct8x8,
      od_bin_idct16x16,
      od_bin_idct32x32,
      od_bin_idct64x64
    };
    test_fdct_2d = OD_FDCT_2D_SSE41;
    test_idct_2d = OD_IDCT_2D_SSE41;
    printf("Testing SSE4.1 DCT...\n");
    run_test();
  }
# endif
# if defined(OD_AVX2_INTRINSICS)
  if (od_cpu_flags_get() & OD_CPU_X86_AVX2) {
    static const od_dct_func_2d OD_FDCT_2D_AVX2[OD_NBSIZES + 1] = {
      od_bin_fdct4x4_sse41,
      od_bin_fdct8x8_avx2,
      od_bin_fdct16x16,
      od_bin_fdct32x32,
      od_bin_fdct64x64
    };
    static const od_dct_func_2d OD_IDCT_2D_AVX2[OD_NBSIZES + 1] = {
      od_bin_idct4x4_sse41,
      od_bin_idct8x8_avx2,
      od_bin_idct16x16,
      od_bin_idct32x32,
      od_bin_idct64x64
    };
    test_fdct_2d = OD_FDCT_2D_AVX2;
    test_idct_2d = OD_IDCT_2D_AVX2;
    printf("Testing AVX2 DCT...\n");
    run_test();
  }
# endif
#endif
#if defined(OD_ARMASM)
  if (od_cpu_flags_get() & OD_CPU_ARM_NEON) {
    static const od_dct_func_2d OD_FDCT_2D_NEON[OD_NBSIZES + 1] = {
      od_bin_fdct4x4,
      od_bin_fdct8x8,
      od_bin_fdct16x16,
      od_bin_fdct32x32,
      od_bin_fdct64x64
    };
    static const od_dct_func_2d OD_IDCT_2D_NEON[OD_NBSIZES + 1] = {
      od_bin_idct4x4,
      od_bin_idct8x8,
      od_bin_idct16x16,
      od_bin_idct32x32,
      od_bin_idct64x64
    };
    test_fdct_2d = OD_FDCT_2D_NEON;
    test_idct_2d = OD_IDCT_2D_NEON;
    run_test();
  }
#endif
  return od_exit_code;
}

#endif
