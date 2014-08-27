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

#include "dct.h"

/*Adjustments for the quantization step in Q15, a factor of 2^(1/4) per
 octave due to the 3/4 power in the quantizer.*/
const int OD_TRANS_QUANT_ADJ[3] = {
  32768, 27554, 23170
};

const unsigned char OD_ZIG4[16] = {
   0,  1,  5,  6,
   2,  4,  7, 12,
   3,  8, 11, 13,
   9, 10, 14, 15
};

const unsigned char OD_ZIG8[64] = {
   0,  1,  5,  6, 14, 15, 27, 28,
   2,  4,  7, 13, 16, 26, 29, 42,
   3,  8, 12, 17, 25, 30, 41, 43,
   9, 11, 18, 24, 31, 40, 44, 53,
  10, 19, 23, 32, 39, 45, 52, 54,
  20, 22, 33, 38, 46, 51, 55, 60,
  21, 34, 37, 47, 50, 56, 59, 61,
  35, 36, 48, 49, 57, 58, 62, 63
};

const unsigned char OD_ZIG16[256] = {
   0,   1,   5,   6,  14,  15,  27,  28,
  44,  45,  65,  66,  90,  91, 119, 120,
   2,   4,   7,  13,  16,  26,  29,  43,
  46,  64,  67,  89,  92, 118, 121, 150,
   3,   8,  12,  17,  25,  30,  42,  47,
  63,  68,  88,  93, 117, 122, 149, 151,
   9,  11,  18,  24,  31,  41,  48,  62,
  69,  87,  94, 116, 123, 148, 152, 177,
  10,  19,  23,  32,  40,  49,  61,  70,
  86,  95, 115, 124, 147, 153, 176, 178,
  20,  22,  33,  39,  50,  60,  71,  85,
  96, 114, 125, 146, 154, 175, 179, 200,
  21,  34,  38,  51,  59,  72,  84,  97,
 113, 126, 145, 155, 174, 180, 199, 201,
  35,  37,  52,  58,  73,  83,  98, 112,
 127, 144, 156, 173, 181, 198, 202, 219,
  36,  53,  57,  74,  82,  99, 111, 128,
 143, 157, 172, 182, 197, 203, 218, 220,
  54,  56,  75,  81, 100, 110, 129, 142,
 158, 171, 183, 196, 204, 217, 221, 234,
  55,  76,  80, 101, 109, 130, 141, 159,
 170, 184, 195, 205, 216, 222, 233, 235,
  77,  79, 102, 108, 131, 140, 160, 169,
 185, 194, 206, 215, 223, 232, 236, 245,
  78, 103, 107, 132, 139, 161, 168, 186,
 193, 207, 214, 224, 231, 237, 244, 246,
 104, 106, 133, 138, 162, 167, 187, 192,
 208, 213, 225, 230, 238, 243, 247, 252,
 105, 134, 137, 163, 166, 188, 191, 209,
 212, 226, 229, 239, 242, 248, 251, 253,
 135, 136, 164, 165, 189, 190, 210, 211,
 227, 228, 240, 241, 249, 250, 254, 255
};

const unsigned char *OD_DCT_ZIGS[OD_NBSIZES + 1] = {
  OD_ZIG4,
  OD_ZIG8,
  OD_ZIG16
};

/*Making function pointer tables at least one entry
   longer than needed makes it highly likely that an
   off-by-one will result in a null-pointer rather than
   another otherwise compatible function pointer.
  This can help avoid difficult to diagnose misbehavior.*/

const od_dct_func_2d OD_FDCT_2D[OD_NBSIZES + 1] = {
  od_bin_fdct4x4,
  od_bin_fdct8x8,
  od_bin_fdct16x16
};

const od_dct_func_2d OD_IDCT_2D[OD_NBSIZES + 1] = {
  od_bin_idct4x4,
  od_bin_idct8x8,
  od_bin_idct16x16
};

const od_fdct_func_1d OD_FDCT_1D[OD_NBSIZES + 1] = {
  od_bin_fdct4,
  od_bin_fdct8,
  od_bin_fdct16
};

const od_idct_func_1d OD_IDCT_1D[OD_NBSIZES + 1] = {
  od_bin_idct4,
  od_bin_idct8,
  od_bin_idct16
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

/*The auto-correlation coefficent. 0.95 is a common value.*/
# define INPUT_AUTOCORR (0.95)
# define INPUT_AUTOCORR_2 (INPUT_AUTOCORR*INPUT_AUTOCORR)
# define INPUT_AUTOCORR_4 (INPUT_AUTOCORR_2*INPUT_AUTOCORR_2)
# define INPUT_AUTOCORR_8 (INPUT_AUTOCORR_4*INPUT_AUTOCORR_4)

/*An autocorrelation table.
  A common model for natural-image input is an AR-0 process with an
   autocorrelation coefficient of 0.95.
  This table contains various powers of 0.95, so that
   AUTOCORR[i-j+15]==(0.95)**abs(i-j), for i,j in [0..15].
  This makes it easy to compute element i,j of the covariance matrix.*/
static const double AUTOCORR[31] = {
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
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR
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

static const od_dct_basis_row *const OD_DCT_BASES[OD_NBSIZES] = {
  OD_DCT4_BASIS,
  OD_DCT8_BASIS,
  OD_DCT16_BASIS
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
      block[i][j] = ieee1180_rand(l, h)*sign << OD_COEFF_SHIFT;
    }
  }
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  (*test_fdct_2d[bszi])(refcoefs[0], OD_BSIZE_MAX, block[0], OD_BSIZE_MAX);
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      floatcoefs[i][j] = refcoefs[i][j]/(double)(1 << OD_COEFF_SHIFT);
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
        e[i][j] += (basis[i][k] - tbasis[i][k])*AUTOCORR[k - j + 15];
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
    for (j = 0; j < n; j++) x[j] = ieee1180_rand(255, 255) << OD_COEFF_SHIFT;
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
      y2[j] = y[j] + (ieee1180_rand(1, 1) << OD_COEFF_SHIFT);
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      rtacc[j] += x2[j] - x[j];
    }
    for (j = 0; j < n; j++) {
      y2[j] = (y[j] + ((y[j] < 0 ? -4 : 4) << OD_COEFF_SHIFT))/
       (8 << OD_COEFF_SHIFT) << (3 + OD_COEFF_SHIFT);
    }
    (*OD_IDCT_1D[bszi])(x2, 1, y2);
    for (j = 0; j < n; j++) {
      x2[j] = OD_COEFF_UNSCALE(x2[j]);
      q8acc[j] += x2[j] - x[j];
    }
    for (j = 0; j < n; j++) {
      y2[j] = (y[j] + ((y[j] < 0 ? -7 : 7) << OD_COEFF_SHIFT >> 1))/
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

int od_dct_check_min[86];
int od_dct_check_max[86];

static void od_bin_fxform_2d(od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2],
 int bszi) {
  od_coeff y[OD_BSIZE_MAX*2];
  int u;
  int v;
  int n;
  n = 1 << (OD_LOG_BSIZE0 + bszi);
  /*Perform pre-filtering.*/
  for (v = 0; v < n*2; v++) {
    for (u = 0; u < n*2; u++) y[u] = x[u][v];
    (*OD_PRE_FILTER[bszi])(y, y);
    (*OD_PRE_FILTER[bszi])(y+n, y+n);
    for (u = 0; u < n*2; u++) x[u][v] = y[u];
  }
  for (u = 0; u < n*2; u++) {
    (*OD_PRE_FILTER[bszi])(x[u], x[u]);
    (*OD_PRE_FILTER[bszi])(x[u] + n, x[u] + n);
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
           x[u + (n >> 1)][v + (n >> 1)]/(256.0*(1 << OD_COEFF_SHIFT));
        }
      }
    }
  }
  for (u = 0; u < n; u++) {
    for (v = 0; v < n; v++) {
      od_coeff x[OD_BSIZE_MAX*2][OD_BSIZE_MAX*2];
      for (i = 0; i < n*2; i++) {
        for (j = 0; j < n*2; j++) {
          x[i][j] = (basis2[u][v][i][j] < 0 ? -255 : 255) << OD_COEFF_SHIFT;
        }
      }
      od_bin_fxform_2d(x, bszi);
      max2[u][v] = x[u + (n >> 1)][v + (n >> 1)];
      for (i = 0; i < n*2; i++) {
        for (j = 0; j < n*2; j++) {
          x[i][j] = (basis2[u][v][i][j] > 0 ? -255 : 255) << OD_COEFF_SHIFT;
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
  for (j = 0; j < n; j++) min[j] = max[j] = 0;
  for (i = 0; i < 1 << n; i++) {
    od_coeff x[OD_BSIZE_MAX];
    od_coeff y[OD_BSIZE_MAX];
    od_coeff x2[OD_BSIZE_MAX];
    for (j = 0; j < n; j++) x[j] = (i >> j & 1) ? 255 : -256;
    (*OD_FDCT_1D[bszi])(y, x, 1);
    (*OD_IDCT_1D[bszi])(x2, 1, y);
    for (j = 0; j < n; j++) {
      if (y[j] < min[j]) min[j] = y[j];
      else if (y[j] > max[j]) max[j] = y[j];
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
  for (bszi = 0; bszi < OD_NBSIZES; bszi++) check_transform(bszi);
}

int main(void) {
  test_fdct_2d = OD_FDCT_2D;
  test_idct_2d = OD_IDCT_2D;
  run_test();
  if (od_cpu_flags_get() & OD_CPU_X86_SSE2) {
    test_fdct_2d = OD_FDCT_2D_SSE2;
    test_idct_2d = OD_IDCT_2D_SSE2;
    run_test();
  }
  if (od_cpu_flags_get() & OD_CPU_X86_SSE4_1) {
    test_fdct_2d = OD_FDCT_2D_SSE4_1;
    test_idct_2d = OD_IDCT_2D_SSE4_1;
    run_test();
  }
  return od_exit_code;
}

#endif
