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

#define OD_FDCT_8(t0, t4, t2, t6, t1, t5, t3, t7) \
  /* Embedded 8-point orthonormal Type-II fDCT. */ \
  do { \
    int t4h; \
    int t6h; \
    int t7h; \
    t7 = t0 - t7; \
    t7h = OD_DCT_RSHIFT(t7, 1); \
    t0 -= t7h; \
    t4 += t3; \
    t4h = OD_DCT_RSHIFT(t4, 1); \
    t3 = t4h - t3; \
    t5 = t2 - t5; \
    t2 -= OD_DCT_RSHIFT(t5, 1); \
    t6 += t1; \
    t6h = OD_DCT_RSHIFT(t6, 1); \
    t1 = t6h - t1; \
    OD_FDCT_4_ASYM(t0, t4, t4h, t2, t6, t6h); \
    OD_FDST_4_ASYM(t7, t7h, t3, t5, t1); \
  } \
  while (0)

#define OD_IDCT_8(t0, t4, t2, t6, t1, t5, t3, t7) \
  /* Embedded 8-point orthonormal Type-II iDCT. */ \
  do { \
    int t1h_; \
    int t3h_; \
    int t7h_; \
    OD_IDST_4_ASYM(t7, t7h_, t5, t6, t4); \
    OD_IDCT_4_ASYM(t0, t2, t1, t1h_, t3, t3h_); \
    t4 = t3h_ - t4; \
    t3 -= t4; \
    t2 += OD_DCT_RSHIFT(t5, 1); \
    t5 = t2 - t5; \
    t6 = t1h_ - t6; \
    t1 -= t6; \
    t0 += t7h_; \
    t7 = t0 - t7; \
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

static const double OD_COS_TABLE64[] = {
  1.000000000000000, 0.999698818696204, 0.998795456205173, 0.997290456678690,
  0.995184726672198, 0.992479534598710, 0.989176509964781, 0.985277642388944,
  0.980785280403236, 0.975702130038531, 0.970031253194543, 0.963776065795439,
  0.956940335732209, 0.949528180593036, 0.941544065183021, 0.932992798834741,
  0.923879532511287, 0.914209755703533, 0.903989293123447, 0.893224301195517,
  0.881921264348357, 0.870086991108713, 0.857728610000273, 0.844853565249711,
  0.831469612302553, 0.817584813151584, 0.803207531480642, 0.788346427626608,
  0.773010453362743, 0.757208846506490, 0.740951125354958, 0.724247082951467,
  0.707106781186547, 0.689540544737065, 0.671558954847015, 0.653172842953774,
  0.634393284163666, 0.615231590580628, 0.595699304492432, 0.575808191417846,
  0.555570233019587, 0.534997619887108, 0.514102744193230, 0.492898192229781,
  0.471396736826001, 0.449611329654606, 0.427555093430283, 0.405241314004990,
  0.382683432365105, 0.359895036534988, 0.336889853392216, 0.313681740398892,
  0.290284677254455, 0.266712757474888, 0.242980179903264, 0.219101240156870,
  0.195090322016152, 0.170961888760305, 0.146730474455372, 0.122410675199216,
  0.098017140329573, 0.073564563599667, 0.049067674327418, 0.024541228522912,
  0.000000000000000,-0.024541228522909,-0.049067674327421,-0.073564563599662,
 -0.098017140329546,-0.122410675199216,-0.146730474455350,-0.170961888760287,
 -0.195090322016137,-0.219101240156846,-0.242980179903253,-0.266712757474898,
 -0.290284677254432,-0.313681740398898,-0.336889853392222,-0.359895036534979,
 -0.382683432365079,-0.405241314005000,-0.427555093430284,-0.449611329654584,
 -0.471396736825990,-0.492898192229771,-0.514102744193211,-0.534997619887093,
 -0.555570233019591,-0.575808191417836,-0.595699304492435,-0.615231590580616,
 -0.634393284163648,-0.653172842953771,-0.671558954847019,-0.689540544737068,
 -0.707106781186544,-0.724247082951464,-0.740951125354960,-0.757208846506483,
 -0.773010453362735,-0.788346427626606,-0.803207531480637,-0.817584813151590,
 -0.831469612302545,-0.844853565249693,-0.857728610000273,-0.870086991108708,
 -0.881921264348346,-0.893224301195518,-0.903989293123444,-0.914209755703528,
 -0.923879532511284,-0.932992798834737,-0.941544065183019,-0.949528180593039,
 -0.956940335732210,-0.963776065795440,-0.970031253194541,-0.975702130038530,
 -0.980785280403231,-0.985277642388942,-0.989176509964780,-0.992479534598710,
 -0.995184726672196,-0.997290456678690,-0.998795456205172,-0.999698818696204,
 -1.000000000000000,-0.999698818696204,-0.998795456205172,-0.997290456678690,
 -0.995184726672197,-0.992479534598710,-0.989176509964783,-0.985277642388942,
 -0.980785280403232,-0.975702130038528,-0.970031253194544,-0.963776065795441,
 -0.956940335732211,-0.949528180593039,-0.941544065183021,-0.932992798834739,
 -0.923879532511285,-0.914209755703527,-0.903989293123444,-0.893224301195516,
 -0.881921264348361,-0.870086991108712,-0.857728610000280,-0.844853565249711,
 -0.831469612302548,-0.817584813151587,-0.803207531480646,-0.788346427626606,
 -0.773010453362738,-0.757208846506486,-0.740951125354970,-0.724247082951469,
 -0.707106781186555,-0.689540544737071,-0.671558954847020,-0.653172842953785,
 -0.634393284163640,-0.615231590580628,-0.595699304492435,-0.575808191417845,
 -0.555570233019595,-0.534997619887102,-0.514102744193224,-0.492898192229784,
 -0.471396736825993,-0.449611329654609,-0.427555093430285,-0.405241314005004,
 -0.382683432365082,-0.359895036534990,-0.336889853392221,-0.313681740398902,
 -0.290284677254463,-0.266712757474901,-0.242980179903267,-0.219101240156870,
 -0.195090322016142,-0.170961888760302,-0.146730474455362,-0.122410675199236,
 -0.098017140329578,-0.073564563599667,-0.049067674327422,-0.024541228522914,
  0.000000000000000, 0.024541228522903, 0.049067674327414, 0.073564563599652,
  0.098017140329555, 0.122410675199211, 0.146730474455358, 0.170961888760289,
  0.195090322016119, 0.219101240156863, 0.242980179903246, 0.266712757474890,
  0.290284677254451, 0.313681740398885, 0.336889853392216, 0.359895036534983,
  0.382683432365075, 0.405241314004987, 0.427555093430278, 0.449611329654580,
  0.471396736825998, 0.492898192229774, 0.514102744193217, 0.534997619887098,
  0.555570233019582, 0.575808191417840, 0.595699304492440, 0.615231590580631,
  0.634393284163641, 0.653172842953777, 0.671558954847014, 0.689540544737053,
  0.707106781186547, 0.724247082951456, 0.740951125354955, 0.757208846506474,
  0.773010453362723, 0.788346427626613, 0.803207531480641, 0.817584813151593,
  0.831469612302535, 0.844853565249709, 0.857728610000262, 0.870086991108707,
  0.881921264348355, 0.893224301195506, 0.903989293123441, 0.914209755703528,
  0.923879532511285, 0.932992798834735, 0.941544065183018, 0.949528180593036,
  0.956940335732212, 0.963776065795440, 0.970031253194544, 0.975702130038526,
  0.980785280403229, 0.985277642388937, 0.989176509964782, 0.992479534598711,
  0.995184726672196, 0.997290456678690, 0.998795456205172, 0.999698818696204
};

void od_bin_fdct64(od_coeff y[64], const od_coeff *x, int xstride) {
  int i;
  double norm;
  norm = sqrt(2.0/64);
  for (i = 0; i < 64; i++) {
    int j;
    double sum;
    sum = 0;
    for (j = 0; j < 64; j++) {
      sum += x[j*xstride]*OD_COS_TABLE64[(i*(2*j + 1)) & 0xff];
    }
    y[i] = floor(0.5 + norm*sum*(i == 0 ? M_SQRT1_2 : 1));
  }
}

void od_bin_idct64(od_coeff *x, int xstride, const od_coeff y[64]) {
  int i;
  double norm;
  norm = sqrt(2.0/64);
  for (i = 0; i < 64; i++) {
    int j;
    double sum;
    sum = y[0]*M_SQRT1_2;
    for (j = 1; j < 64; j++) {
      sum += y[j]*OD_COS_TABLE64[(j*(2*i + 1)) & 0xff];
    }
    x[i*xstride] = floor(0.5 + norm*sum);
  }
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

# if defined(OD_DCT_CHECK_OVERFLOW)

int od_dct_check_min[388];
int od_dct_check_max[388];

#endif

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

static const od_dct_basis_row *const OD_DCT_BASES[OD_NBSIZES] = {
  OD_DCT4_BASIS,
  OD_DCT8_BASIS,
  OD_DCT16_BASIS,
  OD_DCT32_BASIS
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
      y2[j] = y[j] + (ieee1180_rand(1, 1)*OD_COEFF_SCALE);
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      rtacc[j] += x2[j] - x[j];
    }
    for (j = 0; j < n; j++) {
      y2[j] = (y[j] + ((y[j] < 0 ? -4 : 4)*OD_COEFF_SCALE))/
       (8 << OD_COEFF_SHIFT)*(1 << (3 + OD_COEFF_SHIFT));
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      q8acc[j] += x2[j] - x[j];
    }
    for (j = 0; j < n; j++) {
      y2[j] = (y[j] + (((y[j] < 0 ? -7 : 7)*OD_COEFF_SCALE)/2))/
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

static void od_bin_fxform_2d(od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2],
 int bszi) {
  od_coeff y[OD_BSIZE_MAX*2];
  int u;
  int v;
  int n;
  int f;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  /*Test with all filters up to OD_MAX_FILT_SIZE.*/
  for (f = 0; f <= OD_MINI(bszi, OD_MAX_FILT_SIZE); f++) {
    int fn;
    int o;
    fn = 1 << (OD_LOG_BSIZE0 + f);
    o = (n >> 1) - (fn >> 1);
    /*Perform pre-filtering.*/
    for (v = 0; v < n*2; v++) {
      for (u = 0; u < n*2; u++) y[u] = x[u][v];
      (*OD_PRE_FILTER[f])(y + o, y + o);
      (*OD_PRE_FILTER[f])(y + n + o, y + n + o);
      for (u = 0; u < n*2; u++) x[u][v] = y[u];
    }
    for (u = 0; u < n*2; u++) {
      (*OD_PRE_FILTER[f])(x[u] + o, x[u] + o);
      (*OD_PRE_FILTER[f])(x[u] + n + o, x[u] + n + o);
    }
    /*Perform DCT.*/
    for (u = n >> 1; u < n*3 >> 1; u++) {
      (*OD_FDCT_1D[bszi])(x[u] + (n >> 1), x[u] + (n >> 1), 1);
    }
    for (v = n >> 1; v < n*3 >> 1; v++) {
      for (u = n >> 1; u < n*3 >> 1; u++) y[u] = x[u][v];
      (*OD_FDCT_1D[bszi])(y + (n >> 1), y + (n >> 1), 1);
      for (u = n >> 1; u < n*3 >> 1; u++) x[u][v] = y[u];
    }
  }
}

static void dynamic_range(int bszi) {
  static double
   basis2[OD_BSIZE_MAX][OD_BSIZE_MAX][OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
  od_coeff min2[OD_BSIZE_MAX][OD_BSIZE_MAX];
  od_coeff max2[OD_BSIZE_MAX][OD_BSIZE_MAX];
  int i;
  int j;
  int u;
  int v;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  for (i = 0; i < n*2; i++) {
    for (j = 0; j < n*2; j++) {
      od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
      /*Generate impulse.*/
      for (u = 0; u < n*2; u++) {
        for (v = 0; v < n*2; v++) {
          x[u][v] = (u == i && v == j) << (8 + OD_COEFF_SHIFT);
        }
      }
      od_bin_fxform_2d(x, bszi);
      /*Retrieve basis elements.*/
      for (u = 0; u < n; u++) {
        for (v = 0; v < n; v++) {
          basis2[u][v][i][j] =
           x[u + (n >> 1)][v + (n >> 1)]/(256.0*OD_COEFF_SCALE);
        }
      }
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
      od_bin_fxform_2d(x, bszi);
      max2[u][v] = x[u + (n >> 1)][v + (n >> 1)];
      for (i = 0; i < n*2; i++) {
        for (j = 0; j < n*2; j++) {
          x[i][j] = (basis2[u][v][i][j] > 0 ? -255 : 255)*OD_COEFF_SCALE;
        }
      }
      od_bin_fxform_2d(x, bszi);
      min2[u][v] = x[u + (n >> 1)][v + (n >> 1)];
    }
  }
  printf("2-D scaled, prefiltered ranges:\n");
  for (u = 0; u < n; u++) {
    printf(" Min %2i:", u);
    for (v = 0; v < n; v++) printf(" %7i", min2[u][v]);
    printf("\n Max %2i:", u);
    for (v = 0; v < n; v++) printf(" %7i", max2[u][v]);
    printf("\n");
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
  for (bszi = 0; bszi <= OD_BLOCK_32X32; bszi++) check_transform(bszi);
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
      od_bin_fdct32x32
    };
    static const od_dct_func_2d OD_IDCT_2D_SSE2[OD_NBSIZES + 1] = {
      od_bin_idct4x4_sse2,
      od_bin_idct8x8,
      od_bin_idct16x16,
      od_bin_idct32x32
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
      od_bin_fdct32x32
    };
    static const od_dct_func_2d OD_IDCT_2D_SSE41[OD_NBSIZES + 1] = {
      od_bin_idct4x4_sse41,
      od_bin_idct8x8,
      od_bin_idct16x16,
      od_bin_idct32x32
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
      od_bin_fdct32x32
    };
    static const od_dct_func_2d OD_IDCT_2D_AVX2[OD_NBSIZES + 1] = {
      od_bin_idct4x4_sse41,
      od_bin_idct8x8_avx2,
      od_bin_idct16x16,
      od_bin_idct32x32
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
      od_bin_fdct32x32
    };
    static const od_dct_func_2d OD_IDCT_2D_NEON[OD_NBSIZES + 1] = {
      od_bin_idct4x4,
      od_bin_idct8x8,
      od_bin_idct16x16,
      od_bin_idct32x32
    };
    test_fdct_2d = OD_FDCT_2D_NEON;
    test_idct_2d = OD_IDCT_2D_NEON;
    run_test();
  }
#endif
  return od_exit_code;
}

#endif
