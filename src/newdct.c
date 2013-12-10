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

#ifdef HAVE_CONFIG_H
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

void od_bin_fdct4(od_coeff _y[4], const od_coeff *_x, int _xstride) {
  /*9 adds, 2 shifts, 3 "muls".*/
  int t0;
  int t1;
  int t2;
  int t2h;
  int t3;
  /*Initial permutation:*/
  t0 = *(_x + 0*_xstride);
  t2 = *(_x + 1*_xstride);
  t1 = *(_x + 2*_xstride);
  t3 = *(_x + 3*_xstride);
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
  /*23013/32768~=4*sin(\frac{\pi}{8})-2*tan(\frac{\pi}{8})~=
     0.70230660471416898931046248770220*/
  t3 -= (t1*23013+16384)>>15;
  /*21407/32768~=\sqrt{1/2}*cos(\frac{\pi}{8}))
     ~=0.65328148243818826392832158671359*/
  t1 += (t3*21407+16384)>>15;
  /*18293/16384~=4*sin(\frac{\pi}{8})-tan(\frac{\pi}{8})~=
     1.1165201670872640381121512119119*/
  t3 -= (t1*18293+8192)>>14;
  _y[0] = (od_coeff)t0;
  _y[1] = (od_coeff)t1;
  _y[2] = (od_coeff)t2;
  _y[3] = (od_coeff)t3;
}

void od_bin_idct4(od_coeff *_x, int _xstride, const od_coeff _y[4]) {
  int t0;
  int t1;
  int t2;
  int t2h;
  int t3;
  t0 = _y[0];
  t1 = _y[1];
  t2 = _y[2];
  t3 = _y[3];
  t3 += (t1*18293+8192)>>14;
  t1 -= (t3*21407+16384)>>15;
  t3 += (t1*23013+16384)>>15;
  t2 = t0 - t2;
  t2h = OD_DCT_RSHIFT(t2, 1);
  t0 -= t2h - OD_DCT_RSHIFT(t3, 1);
  t1 = t2h - t1;
  *(_x + 0*_xstride) = (od_coeff)t0;
  *(_x + 1*_xstride) = (od_coeff)(t2 - t1);
  *(_x + 2*_xstride) = (od_coeff)t1;
  *(_x + 3*_xstride) = (od_coeff)(t0 - t3);
}

void od_bin_fdct4x4(od_coeff *_y, int _ystride, const od_coeff *_x,
 int _xstride) {
  od_coeff z[4*4];
  int      i;
  for (i = 0; i < 4; i++) od_bin_fdct4(z + 4*i, _x + i, _xstride);
  for (i = 0; i < 4; i++) od_bin_fdct4(_y + _ystride*i, z + i, 4);
}

void od_bin_idct4x4(od_coeff *_x, int _xstride, const od_coeff *_y,
 int _ystride) {
  od_coeff z[4*4];
  int      i;
  for (i = 0; i < 4; i++) od_bin_idct4(z + i, 4, _y + _ystride*i);
  for (i = 0; i < 4; i++) od_bin_idct4(_x + i, _xstride, z + 4*i);
}

void od_bin_fdct8(od_coeff _y[8], const od_coeff *_x, int _xstride) {
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
  t0 = *(_x+0*_xstride);
  t4 = *(_x+1*_xstride);
  t2 = *(_x+2*_xstride);
  t6 = *(_x+3*_xstride);
  t7 = *(_x+4*_xstride);
  t3 = *(_x+5*_xstride);
  t5 = *(_x+6*_xstride);
  t1 = *(_x+7*_xstride);
  /*+1/-1 butterflies:*/
  t1 = t0-t1;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t0 -= t1h;
  t4 += t5;
  t4h = OD_DCT_RSHIFT(t4, 1);
  t5 -= t4h;
  t3 = t2-t3;
  t2 -= OD_DCT_RSHIFT(t3, 1);
  t6 += t7;
  t6h = OD_DCT_RSHIFT(t6, 1);
  t7 = t6h-t7;
  /*+ Embedded 4-point type-II DCT.*/
  t0 += t6h;
  t6 = t0-t6;
  t2 = t4h-t2;
  t4 = t2-t4;
  /*|-+ Embedded 2-point type-II DCT.*/
  /*13573/32768~=\sqrt{2}-1~=0.41421356237309504880168872420970*/
  t0 -= (t4*13573+16384)>>15;
  /*11585/16384~=\sqrt{\frac{1}{2}}~=0.70710678118654752440084436210485*/
  t4 += (t0*11585+8192)>>14;
  /*13573/32768~=\sqrt{2}-1~=0.41421356237309504880168872420970*/
  t0 -= (t4*13573+16384)>>15;
  /*|-+ Embedded 2-point type-IV DST.*/
  /*21895/32768~=\frac{1-cos(\frac{3\pi}{8})}{\sin(\frac{3\pi}{8})}~=
     0.66817863791929891999775768652308*/
  t6 -= (t2*21895+16384)>>15;
  /*15137/16384~=sin(\frac{3\pi}{8})~=0.92387953251128675612818318939679*/
  t2 += (t6*15137+8192)>>14;
  /*21895/32768~=\frac{1-cos(\frac{3\pi}{8})}{\sin(\frac{3\pi}{8})}~=
     0.66817863791929891999775768652308*/
  t6 -= (t2*21895+16384)>>15;
  /*+ Embedded 4-point type-IV DST.*/
  /*19195/32768~=2-\sqrt{2}~=0.58578643762690495119831127579030*/
  t3 += (t5*19195+16384)>>15;
  /*11585/16384~=\sqrt{\frac{1}{2}}~=0.70710678118654752440084436210485*/
  t5 += (t3*11585+8192)>>14;
  /*29957/32768~=\sqrt{2}-\frac{1}{2}~=0.91421356237309504880168872420970*/
  t3 -= (t5*29957+16384)>>15;
  t7 = OD_DCT_RSHIFT(t5, 1)-t7;
  t5 -= t7;
  t3 = t1h-t3;
  t1 -= t3;
  /*3227/32768~=\frac{1-cos(\frac{\pi}{16})}{sin(\frac{\pi}{16})}~=
     0.098491403357164253077197521291327*/
  t7 += (t1*3227+16384)>>15;
  /*6393/32768~=sin(\frac{\pi}{16})~=0.19509032201612826784828486847702*/
  t1 -= (t7*6393+16384)>>15;
  /*3227/32768~=\frac{1-cos(\frac{\pi}{16})}{sin(\frac{\pi}{16})}~=
     0.098491403357164253077197521291327*/
  t7 += (t1*3227+16384)>>15;
  /*2485/8192~=\frac{1-cos(\frac{3\pi}{16})}{sin(\frac{3\pi}{16})}~=
     0.30334668360734239167588394694130*/
  t5 += (t3*2485+4096)>>13;
  /*18205/32768~=sin(\frac{3\pi}{16})~=0.55557023301960222474283081394853*/
  t3 -= (t5*18205+16384)>>15;
  /*2485/8192~=\frac{1-cos(\frac{3\pi}{16})}{sin(\frac{3\pi}{16})}~=
     0.30334668360734239167588394694130*/
  t5 += (t3*2485+4096)>>13;
  _y[0] = (od_coeff)t0;
  _y[1] = (od_coeff)t1;
  _y[2] = (od_coeff)t2;
  _y[3] = (od_coeff)t3;
  _y[4] = (od_coeff)t4;
  _y[5] = (od_coeff)t5;
  _y[6] = (od_coeff)t6;
  _y[7] = (od_coeff)t7;
}

void od_bin_idct8(od_coeff *_x, int _xstride, const od_coeff _y[16]) {
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
  t0 = _y[0];
  t1 = _y[1];
  t2 = _y[2];
  t3 = _y[3];
  t4 = _y[4];
  t5 = _y[5];
  t6 = _y[6];
  t7 = _y[7];
  t5 -= (t3*2485+4096)>>13;
  t3 += (t5*18205+16384)>>15;
  t5 -= (t3*2485+4096)>>13;
  t7 -= (t1*3227+16384)>>15;
  t1 += (t7*6393+16384)>>15;
  t7 -= (t1*3227+16384)>>15;
  t1 += t3;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t3 = t1h-t3;
  t5 += t7;
  t7 = OD_DCT_RSHIFT(t5, 1)-t7;
  t3 += (t5*29957+16384)>>15;
  t5 -= (t3*11585+8192)>>14;
  t3 -= (t5*19195+16384)>>15;
  t6 += (t2*21895+16384)>>15;
  t2 -= (t6*15137+8192)>>14;
  t6 += (t2*21895+16384)>>15;
  t0 += (t4*13573+16384)>>15;
  t4 -= (t0*11585+8192)>>14;
  t0 += (t4*13573+16384)>>15;
  t4 = t2-t4;
  t4h = OD_DCT_RSHIFT(t4, 1);
  t2 = t4h-t2;
  t6 = t0-t6;
  t6h = OD_DCT_RSHIFT(t6, 1);
  t0 -= t6h;
  t7 = t6h-t7;
  t6 -= t7;
  t2 += OD_DCT_RSHIFT(t3, 1);
  t3 = t2-t3;
  t5 += t4h;
  t4 -= t5;
  t0 += t1h;
  t1 = t0-t1;
  *(_x+0*_xstride) = (od_coeff)t0;
  *(_x+1*_xstride) = (od_coeff)t4;
  *(_x+2*_xstride) = (od_coeff)t2;
  *(_x+3*_xstride) = (od_coeff)t6;
  *(_x+4*_xstride) = (od_coeff)t7;
  *(_x+5*_xstride) = (od_coeff)t3;
  *(_x+6*_xstride) = (od_coeff)t5;
  *(_x+7*_xstride) = (od_coeff)t1;
}

void od_bin_fdct8x8(od_coeff *_y, int _ystride, const od_coeff *_x,
 int _xstride) {
  od_coeff z[8*8];
  int      i;
  for (i = 0; i < 8; i++) od_bin_fdct8(z+8*i, _x+i, _xstride);
  for (i = 0; i < 8; i++) od_bin_fdct8(_y+_ystride*i, z+i, 8);
}

void od_bin_idct8x8(od_coeff *_x, int _xstride, const od_coeff *_y,
 int _ystride) {
  od_coeff z[8*8];
  int      i;
  for (i = 0; i < 8; i++) od_bin_idct8(z+i, 8, _y+_ystride*i);
  for (i = 0; i < 8; i++) od_bin_idct8(_x+i, _xstride, z+8*i);
}

void od_bin_fdct16(od_coeff _y[16], const od_coeff *_x, int _xstride) {
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
  t0 = *(_x+0*_xstride);
  t8 = *(_x+1*_xstride);
  t4 = *(_x+2*_xstride);
  tc = *(_x+3*_xstride);
  te = *(_x+4*_xstride);
  ta = *(_x+5*_xstride);
  t6 = *(_x+6*_xstride);
  t2 = *(_x+7*_xstride);
  t3 = *(_x+8*_xstride);
  td = *(_x+9*_xstride);
  t9 = *(_x+10*_xstride);
  tf = *(_x+11*_xstride);
  t1 = *(_x+12*_xstride);
  t7 = *(_x+13*_xstride);
  tb = *(_x+14*_xstride);
  t5 = *(_x+15*_xstride);
  /*+1/-1 butterflies:*/
  t5 = t0-t5;
  t8 += tb;
  t7 = t4-t7;
  tc += t1;
  tf = te-tf;
  ta += t9;
  td = t6-td;
  t2 += t3;
  t0 -= OD_DCT_RSHIFT(t5, 1);
  t8h = OD_DCT_RSHIFT(t8, 1);
  tb = t8h-tb;
  t4 -= OD_DCT_RSHIFT(t7, 1);
  tch = OD_DCT_RSHIFT(tc, 1);
  t1 = tch-t1;
  te -= OD_DCT_RSHIFT(tf, 1);
  tah = OD_DCT_RSHIFT(ta, 1);
  t9 = tah-t9;
  t6 -= OD_DCT_RSHIFT(td, 1);
  t2h = OD_DCT_RSHIFT(t2, 1);
  t3 = t2h-t3;
  /*+ Embedded 8-point type-II DCT.*/
  t0 += t2h;
  t6 = t8h-t6;
  t4 += tah;
  te = tch-te;
  t2 = t0-t2;
  t8 -= t6;
  ta = t4-ta;
  tc -= te;
  /*|-+ Embedded 4-point type-II DCT.*/
  tc = t0-tc;
  t8 += t4;
  t8h = OD_DCT_RSHIFT(t8, 1);
  t4 = t8h-t4;
  t0 -= OD_DCT_RSHIFT(tc, 1);
  /*|-|-+ Embedded 2-point type-II DCT.*/
  t0 += t8h;
  t8 = t0-t8;
  /*|-|-+ Embedded 2-point type-IV DST.*/
  /*32013/32768~=4*sin(\frac{\pi}{8})-2*tan(\frac{\pi}{8})
     ~=0.70230660471416898931046248770220*/
  tc -= (t4*23013+16384)>>15;
  /*21407~=\sqrt{1/2}*cos(\frac{\pi}{8}))~=0.65328148243818826392832158671359*/
  t4 += (tc*21407+16384)>>15;
  /*18293/16384~=4*sin(\frac{\pi}{8})-tan(\frac{\pi}{8})
     ~=1.1165201670872640381121512119119*/
  tc -= (t4*18293+8192)>>14;
  /*|-+ Embedded 4-point type-IV DST.*/
  /*13573/32768~=\sqrt{2}-1~=0.41421356237309504880168872420970*/
  t6 += (ta*13573+16384)>>15;
  /*11585/16384~=\sqrt{\frac{1}{2}}~=0.70710678118654752440084436210485*/
  ta -= (t6*11585+8192)>>14;
  /*13573/32768~=\sqrt{2}-1~=0.41421356237309504880168872420970*/
  t6 += (ta*13573+16384)>>15;
  ta += te;
  t2 += t6;
  te = OD_DCT_RSHIFT(ta, 1)-te;
  t6 = OD_DCT_RSHIFT(t2, 1)-t6;
  /*2775/2048~=\frac{\sqrt{2}-cos(\frac{\pi}{16})}{2sin(\frac{\pi}{16})}
     ~=1.1108400393486273201524536919723*/
  te += (t2*2275+1024)>>11;
  /*9041/32768~=\sqrt{2}sin(\frac{\pi}{16})
     ~=0.27589937928294301233595756366937*/
  t2 -= (te*9041+16384)>>15;
  /*2873/2048
     ~=\frac{cos(\frac{\pi}{16})-\sqrt{\frac{1}{2}}}{sin(\frac{\pi}{16})}
     ~=1.4028297067142967321050338435598*/
  te -= (t2*2873+1024)>>11;
  /*17185/32768~=\frac{\sqrt{2}-cos(\frac{3\pi}{16})}{2sin(\frac{3\pi}{16})}
     ~=0.52445569924008942966043945081053*/
  t6 -= (ta*17185+16384)>>15;
  /*12873/16384~=\sqrt{2}sin(\frac{3\pi}{16})
     ~=0.78569495838710218127789736765722*/
  ta += (t6*12873+8192)>>14;
  /*7335/32768
     ~=\frac{cos(\frac{3\pi}{16})-\sqrt{\frac{1}{2}}}{sin(\frac{3\pi}{16})}
     ~=0.22384718209265507914012811666071*/
  t6 += (ta*7335+16384)>>15;
  /*+ Embedded 8-point type-IV DST.*/
  /*1035/2048~=\frac{\sqrt{2}-cos(\frac{7\pi}{32})}{2sin(\frac{7\pi}{32})}
     ~=0.50536719493782972897642806316664*/
  t3 += (t5*1035+1024)>>11;
  /*14699/16384~=\sqrt{2}sin(\frac{7\pi}{32})
     ~=0.89716758634263628425064138922084*/
  t5 -= (t3*14699+8192)>>14;
  /*851/8192
    ~=\frac{cos(\frac{7\pi}{32})-\sqrt{\frac{1}{2}}}{sin(\frac{7\pi}{32}}
    ~=0.10388456785615844342131055214354*/
  t3 -= (t5*851+4096)>>13;
  /*17515/32768~=\frac{\sqrt{2}-cos(\frac{11\pi}{32})}{2sin(\frac{11\pi}{32})}
     ~=0.53452437516842143578098634302964*/
  tb += (td*17515+16384)>>15;
  /*40869/32768~=\sqrt{2}sin(\frac{11\pi}{32})
     ~=1.2472250129866712325719902627977*/
  td -= (tb*40869+16384)>>15;
  /*4379/16384
     ~=\frac{\sqrt{\frac{1}{2}}-cos(\frac{11\pi}{32})}{sin(\frac{11\pi}{32})}
     ~=0.26726880719302561523614336238196*/
  tb += (td*4379+8192)>>14;
  /*25809/32768~=\frac{\sqrt{2}-cos(\frac{3\pi}{32})}{2sin(\frac{3\pi}{32})}
     ~=0.78762894232967441973847776796517*/
  t9 += (t7*25809+16384)>>15;
  /*3363/8192~=\sqrt{2}sin(\frac{3\pi}{32})
     ~=0.41052452752235738115636923775513*/
  t7 -= (t9*3363+4096)>>13;
  /*14101/16384
     ~=\frac{cos(\frac{3\pi}{32})-\sqrt{\frac{1}{2}}}{sin(\frac{3\pi}{32})}
     ~=0.86065016213948579370059934044795*/
  t9 -= (t7*14101+8192)>>14;
  /*21669/32768~=\frac{\sqrt{2}-cos(\frac{15\pi}{32})}{2sin(\frac{15\pi}{32})}
     ~=0.66128246684651710406296283785232*/
  t1 += (tf*21669+16384)>>15;
  /*23059/16384~=\sqrt{2}sin(\frac{15\pi}{32})
     ~=1.4074037375263824590260782229840*/
  tf -= (t1*23059+8192)>>14;
  /*20055/32768
    ~=\frac{\sqrt{\frac{1}{2}}-cos(\frac{15\pi}{32})}{sin(\frac{15\pi}{32})}
    ~=0.61203676516793497752436407720666*/
  t1 += (tf*20055+16384)>>15;
  tf = t3-tf;
  td += t9;
  tfh = OD_DCT_RSHIFT(tf, 1);
  t3 -= tfh;
  tdh = OD_DCT_RSHIFT(td, 1);
  t9 = tdh-t9;
  t1 += t5;
  tb = t7-tb;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t5 = t1h-t5;
  tbh = OD_DCT_RSHIFT(tb, 1);
  t7 -= tbh;
  t3 += tbh;
  t5 = tdh-t5;
  t9 += tfh;
  t7 = t1h-t7;
  tb -= t3;
  td -= t5;
  tf = t9-tf;
  t1 -= t7;
  /*21895/32768~=\frac{1-cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})}
     ~=0.66817863791929891999775768652308*/
  t5 -= (tb*21895+16384)>>15;
  /*15137/16384~=sin(\frac{3\pi}{8})~=0.92387953251128675612818318939679*/
  tb += (t5*15137+8192)>>14;
  /*21895/32768~=\frac{1-cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})}
     ~=0.66817863791929891999775768652308*/
  t5 -= (tb*21895+16384)>>15;
  /*21895/32768~=\frac{1-cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})}
     ~=0.66817863791929891999775768652308*/
  td += (t3*21895+16384)>>15;
  /*15137/16384~=sin(\frac{3\pi}{8})~=0.92387953251128675612818318939679*/
  t3 -= (td*15137+8192)>>14;
  /*21895/32768~=\frac{1-cos(\frac{3\pi}{8})}{sin(\frac{3\pi}{8})}
     ~=0.66817863791929891999775768652308*/
  td += (t3*21895+16384)>>15;
  /*13573/32768~=\sqrt{2}-1~=0.41421356237309504880168872420970*/
  t1 -= (tf*13573+16384)>>15;
  /*11585/16384~=\sqrt{\frac{1}{2}}~=0.70710678118654752440084436210485*/
  tf += (t1*11585+8192)>>14;
  /*13573/32768~=\sqrt{2}-1~=0.41421356237309504880168872420970*/
  t1 -= (tf*13573+16384)>>15;
  _y[0] = (od_coeff)t0;
  _y[1] = (od_coeff)t1;
  _y[2] = (od_coeff)t2;
  _y[3] = (od_coeff)t3;
  _y[4] = (od_coeff)t4;
  _y[5] = (od_coeff)t5;
  _y[6] = (od_coeff)t6;
  _y[7] = (od_coeff)t7;
  _y[8] = (od_coeff)t8;
  _y[9] = (od_coeff)t9;
  _y[10] = (od_coeff)ta;
  _y[11] = (od_coeff)tb;
  _y[12] = (od_coeff)tc;
  _y[13] = (od_coeff)td;
  _y[14] = (od_coeff)te;
  _y[15] = (od_coeff)tf;
}

void od_bin_idct16(od_coeff *_x, int _xstride, const od_coeff _y[16]) {
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
  t0 = _y[0];
  t1 = _y[1];
  t2 = _y[2];
  t3 = _y[3];
  t4 = _y[4];
  t5 = _y[5];
  t6 = _y[6];
  t7 = _y[7];
  t8 = _y[8];
  t9 = _y[9];
  ta = _y[10];
  tb = _y[11];
  tc = _y[12];
  td = _y[13];
  te = _y[14];
  tf = _y[15];
  t1 += (tf*13573+16384)>>15;
  tf -= (t1*11585+8192)>>14;
  t1 += ((tf*13573+16384)>>15)+t7;
  td -= (t3*21895+16384)>>15;
  t3 += (td*15137+8192)>>14;
  t5 += (tb*21895+16384)>>15;
  tb -= (t5*15137+8192)>>14;
  t5 += (tb*21895+16384)>>15;
  td += t5-((t3*21895+16384)>>15);
  tf = t9-tf;
  tb += t3;
  tfh = OD_DCT_RSHIFT(tf, 1);
  t9 -= tfh;
  tbh = OD_DCT_RSHIFT(tb, 1);
  t3 += tfh-tbh;
  t1h = OD_DCT_RSHIFT(t1, 1);
  t7 = t1h-t7+tbh;
  tdh = OD_DCT_RSHIFT(td, 1);
  t5 += t1h-tdh;
  t9 = tdh-t9;
  td -= t9;
  tf = t3-tf;
  t1 -= t5+((tf*20055+16384)>>15);
  tf += (t1*23059+8192)>>14;
  t1 -= (tf*21669+16384)>>15;
  tb = t7-tb;
  t9 += (t7*14101+8192)>>14;
  t7 += (t9*3363+4096)>>13;
  t9 -= (t7*25809+16384)>>15;
  tb -= (td*4379+8192)>>14;
  td += (tb*40869+16384)>>15;
  tb -= (td*17515+16384)>>15;
  t3 += (t5*851+4096)>>13;
  t5 += (t3*14699+8192)>>14;
  t3 -= (t5*1035+1024)>>11;
  t6 -= (ta*7335+16384)>>15;
  ta -= (t6*12873+8192)>>14;
  te += (t2*2873+1024)>>11;
  t2 += (te*9041+16384)>>15;
  t6 = OD_DCT_RSHIFT(t2, 1)-t6-((ta*17185+16384)>>15);
  te = OD_DCT_RSHIFT(ta, 1)-te+((t2*2275+1024)>>11);
  t2 -= t6;
  ta -= te;
  t6 -= (ta*13573+16384)>>15;
  ta += (t6*11585+8192)>>14;
  t6 -= (ta*13573+16384)>>15;
  tc += (t4*18293+8192)>>14;
  t4 -= (tc*21407+16384)>>15;
  tc += (t4*23013+16384)>>15;
  t8 = t0-t8;
  t8h = OD_DCT_RSHIFT(t8, 1);
  t0 -= t8h-OD_DCT_RSHIFT(tc, 1);
  t4 = t8h-t4;
  t8 += t6-t4;
  tc = t0-tc+te;
  ta = t4-ta;
  t2 = t0-t2;
  tch = OD_DCT_RSHIFT(tc, 1);
  te = tch-te;
  tah = OD_DCT_RSHIFT(ta, 1);
  t4 -= tah;
  t8h = OD_DCT_RSHIFT(t8, 1);
  t6 = t8h-t6;
  t2h = OD_DCT_RSHIFT(t2, 1);
  t0 -= t2h;
  t3 = t2h-t3;
  t6 += OD_DCT_RSHIFT(td, 1);
  t9 = tah-t9;
  te += OD_DCT_RSHIFT(tf, 1);
  t1 = tch-t1;
  t4 += OD_DCT_RSHIFT(t7, 1);
  tb = t8h-tb;
  t0 += OD_DCT_RSHIFT(t5, 1);
  *(_x+0*_xstride) = (od_coeff)t0;
  *(_x+1*_xstride) = (od_coeff)(t8-tb);
  *(_x+2*_xstride) = (od_coeff)t4;
  *(_x+3*_xstride) = (od_coeff)(tc-t1);
  *(_x+4*_xstride) = (od_coeff)te;
  *(_x+5*_xstride) = (od_coeff)(ta-t9);
  *(_x+6*_xstride) = (od_coeff)t6;
  *(_x+7*_xstride) = (od_coeff)(t2-t3);
  *(_x+8*_xstride) = (od_coeff)t3;
  *(_x+9*_xstride) = (od_coeff)(t6-td);
  *(_x+10*_xstride) = (od_coeff)t9;
  *(_x+11*_xstride) = (od_coeff)(te-tf);
  *(_x+12*_xstride) = (od_coeff)t1;
  *(_x+13*_xstride) = (od_coeff)(t4-t7);
  *(_x+14*_xstride) = (od_coeff)tb;
  *(_x+15*_xstride) = (od_coeff)(t0-t5);
}

void od_bin_fdct16x16(od_coeff *_y, int _ystride,
 const od_coeff *_x, int _xstride) {
  od_coeff z[16*16];
  int      i;
  for (i = 0; i < 16; i++) od_bin_fdct16(z+16*i, _x+i, _xstride);
  for (i = 0; i < 16; i++) od_bin_fdct16(_y+_ystride*i, z+i, 16);
}

void od_bin_idct16x16(od_coeff *_x, int _xstride,
 const od_coeff *_y, int _ystride) {
  od_coeff z[16*16];
  int      i;
  for (i = 0; i < 16; i++) od_bin_idct16(z+i, 16, _y+_ystride*i);
  for (i = 0; i < 16; i++) od_bin_idct16(_x+i, _xstride, z+16*i);
}

#if OD_DCT_TEST
/*Test code.*/
# include <stdio.h>
# include <math.h>
# include <string.h>


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
static void ieee1180_srand(int _seed) {
  ieee1180_rand_x = _seed;
}

/*Computes a random number between -l and h, inclusive, according to the
   specification in IEEE Std 1180-1990, "IEEE Standard Specifications for the
   Implementations of 8x8 Inverse Discrete Cosine Transform."*/
static int ieee1180_rand(int _l, int _h) {
  double x;
  int    i;
  ieee1180_rand_x = (ieee1180_rand_x*1103515245)+12345;
  i = ieee1180_rand_x&0x7FFFFFFE;
  x = i/(double)0x7FFFFFFF*(_l+_h+1);
  return (int)x-_l;
}

static char *ieee1180_meets(double _val, double _limit) {
  return fabs(_val) <= _limit ? "meets" : "FAILS";
}

/*The number of different input ranges.*/
# define IEEE1180_NRANGES (3)

/*The number of blocks of data to generate.*/
# define IEEE1180_NBLOCKS (10000)

static const int IEEE1180_L[IEEE1180_NRANGES] = { 256, 5, 300 };
static const int IEEE1180_H[IEEE1180_NRANGES] = { 255, 5, 300 };


/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT4_ISCALE[4] = {
  1, 1, 1, 1
};

/*The true forward 4-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT4_BASIS[4][4] = {
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

void idct4(double _x[], const double _y[]) {
  double t[8];
  int    i;
  int    j;
  for (j = 0; j < 4; j++) {
    t[j] = 0;
    for (i = 0; i < 4; i++) t[j] += DCT4_BASIS[i][j]*_y[i];
  }
  for (j = 0; j < 4; j++) _x[j] = t[j];
}

void fdct4(double _x[], const double _y[]) {
  double t[8];
  int    i;
  int    j;
  for (j = 0; j < 4; j++) {
    t[j] = 0;
    for (i = 0; i < 4; i++) t[j] += DCT4_BASIS[j][i]*_y[i];
  }
  for (j = 0; j < 4; j++) _x[j] = t[j];
}

static void ieee1180_print_results4(long _sumerrs[4][4], long _sumsqerrs[4][4],
 int _maxerr[4][4], int _l, int _h, int _sign) {
  double max;
  double total;
  int    m;
  int    i;
  int    j;
  printf("IEEE1180-1990 test results:\n");
  printf("Input range: [%i,%i]\n", -_l, _h);
  printf("Sign: %i\n", _sign);
  printf("Iterations: %i\n\n", IEEE1180_NBLOCKS);
  printf("Peak absolute values of errors:\n");
  for (i = 0, m = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (_maxerr[i][j] > m) m = _maxerr[i][j];
      printf("%4i", _maxerr[i][j]);
    }
    printf("\n");
  }
  printf("Worst peak error = %i (%s spec limit 1)\n\n", m,
   ieee1180_meets((double)m, 1.0));
  printf("Mean square errors:\n");
  max = total = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      double err;
      err = _sumsqerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += _sumsqerrs[i][j];
      if (max < err) max = err;
    }
    printf("\n");
  }
  printf("Worst pmse = %.6f (%s spec limit 0.06)\n", max,
   ieee1180_meets(max, 0.06));
  total /= 4*4*(double)IEEE1180_NBLOCKS;
  printf("Overall mse = %.6f (%s spec limit 0.02)\n\n", total,
   ieee1180_meets(total, 0.02));
  printf("Mean errors:\n");
  max = total = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      double err;
      err = _sumerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += _sumerrs[i][j];
      if (err < 0) err = -err;
      if (max < err) max = err;
    }
    printf("\n");
  }
  total /= 4*4*(double)IEEE1180_NBLOCKS;
  printf("Worst mean error = %.6f (%s spec limit 0.015)\n", max,
   ieee1180_meets(max, 0.015));
  printf("Overall mean error = %.6f (%s spec limit 0.0015)\n\n", total,
   ieee1180_meets(total, 0.0015));
}

static void ieee1180_test_block4(long _sumerrs[4][4], long _sumsqerrs[4][4],
 int _maxerr[4][4], int _l, int _h, int _sign) {
  od_coeff block[4][4];
  od_coeff refcoefs[4][4];
  od_coeff refout[4][4];
  od_coeff testout[4][4];
  double   floatcoefs[4][4];
  int      maxerr;
  int      i;
  int      j;
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++) block[i][j] = ieee1180_rand(_l, _h)*_sign;
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  od_bin_fdct4x4(refcoefs[0], 4, block[0], 4);
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      /*Modification of IEEE1180: inputs to reference iDCT are scaled to match
         the scaling factors introduced by the forward integer transform.*/
      floatcoefs[i][j] = refcoefs[i][j]/(DCT4_ISCALE[i]*DCT4_ISCALE[j]);
    }
  }
  for (i = 0; i < 4; i++) idct4(floatcoefs[i], floatcoefs[i]);
  for (j = 0; j < 4; j++) {
    double x[4];
    for (i = 0; i < 4; i++) x[i] = floatcoefs[i][j];
    idct4(x, x);
    for (i = 0; i < 4; i++) floatcoefs[i][j] = x[i];
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      refout[i][j] = (od_coeff)(floatcoefs[i][j]+0.5);
      if (refout[i][j] > 255) refout[i][j] = 255;
      else if (refout[i][j] < -256) refout[i][j] = -256;
    }
  }
  od_bin_idct4x4(testout[0], 4, refcoefs[0], 4);
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      if (testout[i][j] != block[i][j]) {
        printf("Mismatch:\n");
        printf("in:\n");
        for (i = 0; i < 4; i++) {
          printf("       ");
          for (j = 0; j < 4; j++) printf(" %i", block[i][j]);
          printf("\n");
        }
        printf("xform:\n");
        for (i = 0; i < 4; i++) {
          printf("       ");
          for (j = 0; j < 4; j++) printf(" %i", refcoefs[i][j]);
          printf("\n");
        }
        printf("out:\n");
        for (i = 0; i < 4; i++) {
          printf("       ");
          for (j = 0; j < 4; j++) printf(" %i", testout[i][j]);
          printf("\n");
        }
        printf("\n");
        return;
      }
      if (testout[i][j] > 255) testout[i][j] = 255;
      else if (testout[i][j] < -256) testout[i][j] = -256;
    }
  }
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      int err;
      err = testout[i][j]-refout[i][j];
      _sumerrs[i][j] += err;
      _sumsqerrs[i][j] += err*err;
      if (err < 0) err = -err;
      if (_maxerr[i][j] < err) _maxerr[i][j] = err;
    }
  }
  for (i = 0, maxerr = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      int err;
      err = testout[i][j]-refout[i][j];
      if (err < 0) err = -err;
      if (err > maxerr) maxerr = err;
    }
  }
  /*if(maxerr>1){
    int u;
    int v;
    printf("Excessive peak error: %i\n",maxerr);
    printf("Input:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",block[u][v]);
      }
      printf("\n");
    }
    printf("Forward transform coefficients:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",refcoefs[u][v]);
      }
      printf("\n");
    }
    printf("Reference inverse:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",refout[u][v]);
      }
      printf("\n");
    }
    printf("Integerized inverse:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",testout[u][v]);
      }
      printf("\n");
    }
  }*/
}

static void ieee1180_test4(void) {
  long sumerrs[4][4];
  long sumsqerrs[4][4];
  int  maxerr[4][4];
  int  i;
  int  j;
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block4(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], 1);
    }
    ieee1180_print_results4(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], 1);
  }
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block4(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], -1);
    }
    ieee1180_print_results4(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], -1);
  }
}



static void print_basis4(double _basis[4][4]) {
  int i;
  int j;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      printf("%8.5lf%c", _basis[i][j], j == 3 ? '\n' : ' ');
    }
  }
}

static void compute_fbasis4(double _basis[4][4]) {
  int i;
  int j;
  for (i = 0; i < 4; i++) {
    od_coeff x[4];
    for (j = 0; j < 4; j++) x[j] = (i == j)<<8;
    od_bin_fdct4(x, x, 1);
    for (j = 0; j < 4; j++) _basis[j][i] = x[j]/256.0;
  }
}

static void compute_ftrue_basis4(double _basis[4][4]) {
  int i;
  int j;
  for (j = 0; j < 4; j++) {
    for (i = 0; i < 4; i++) {
      _basis[j][i] = sqrt(2.0/4)*cos((i+0.5)*j*M_PI/4)*DCT4_ISCALE[j];
      if (j == 0) _basis[j][i] *= M_SQRT1_2;
    }
  }
}

static double compute_mse4(double _basis[4][4], double _tbasis[4][4]) {
  double e[4][4];
  double ret;
  int    i;
  int    j;
  int    k;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      e[i][j] = 0;
      for (k = 0; k < 4; k++) {
        e[i][j] += (_basis[i][k]-_tbasis[i][k])*AUTOCORR[k-j+15];
      }
    }
  }
  ret = 0;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      ret += e[i][j]*(_basis[i][j]-_tbasis[i][j]);
    }
  }
  return ret/4;
}

static void check_bias4() {
  double rtacc[4];
  double facc[4];
  double q8acc[4];
  double q7acc[4];
  int    i;
  int    j;
  for (j = 0; j < 4; j++) q7acc[j] = q8acc[j] = facc[j] = rtacc[j] = 0;
  ieee1180_srand(1);
  for (i = 0; i < 10000000; i++) {
    od_coeff x[4];
    od_coeff x2[4];
    od_coeff y[4];
    od_coeff y2[4];
    for (j = 0; j < 4; j++) x[j] = ieee1180_rand(255, 255);
    od_bin_fdct4(y, x, 1);
    for (j = 0; j < 4; j++) facc[j] += y[j];
    od_bin_idct4(x2, 1, y);
    for (j = 0; j < 4; j++) {
      if (x[j] != x2[j]) {
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < 4; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < 4; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < 4; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
    for (j = 0; j < 4; j++) y2[j] = y[j]+ieee1180_rand(1, 1);
    od_bin_idct4(x2, 1, y2);
    for (j = 0; j < 4; j++) rtacc[j] += x2[j]-x[j];
    for (j = 0; j < 4; j++) y2[j] = y[j]/8<<3;
    od_bin_idct4(x2, 1, y2);
    for (j = 0; j < 4; j++) q8acc[j] += x2[j]-x[j];
    for (j = 0; j < 4; j++) y2[j] = y[j]/7*7;
    od_bin_idct4(x2, 1, y2);
    for (j = 0; j < 4; j++) q7acc[j] += x2[j]-x[j];
  }
  printf("1-D Forward Bias:\n");
  for (j = 0; j < 4; j++)
     printf("% -18.15G%s", facc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Round-Trip Bias:\n");
  for (j = 0; j < 4; j++)
     printf("% -18.15G%s", rtacc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Q=8 Bias:\n");
  for (j = 0; j < 4; j++)
     printf("% -18.15G%s", q8acc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Q=7 Bias:\n");
  for (j = 0; j < 4; j++)
     printf("% -18.15G%s", q7acc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
}

static void check4(void) {
  od_coeff min[4];
  od_coeff max[4];
  double   basis[4][4];
  double   tbasis[4][4];
  int      i;
  int      j;
  for (j = 0; j < 4; j++) min[j] = max[j] = 0;
  for (i = 0; i < 1<<4; i++) {
    od_coeff x[4];
    od_coeff y[4];
    od_coeff x2[4];
    for (j = 0; j < 4; j++) x[j] = (i>>j&1) ? 255 : -256;
    od_bin_fdct4(y, x, 1);
    od_bin_idct4(x2, 1, y);
    for (j = 0; j < 4; j++) {
      if (y[j] < min[j]) min[j] = y[j];
      else if (y[j] > max[j]) max[j] = y[j];
    }
    for (j = 0; j < 4; j++) {
      if (x[j] != x2[j]) {
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < 4; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < 4; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < 4; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
  }
  printf("Min:");
  for (j = 0; j < 4; j++) printf(" %5i", min[j]);
  printf("\nMax:");
  for (j = 0; j < 4; j++) printf(" %5i", max[j]);
  printf("\nod_bin_fdct4 basis:\n");
  compute_fbasis4(basis);
  print_basis4(basis);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis4(tbasis);
  print_basis4(tbasis);
  printf("MSE: %.32lg\n\n", compute_mse4(basis, tbasis));
  ieee1180_test4();
  check_bias4();
}


/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT8_ISCALE[8] = {
  1, 1, 1, 1, 1, 1, 1, 1
};

/*The true forward 8-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT8_BASIS[8][8] = {
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

void idct8(double _x[], const double _y[]) {
  double t[8];
  int    i;
  int    j;
  for (j = 0; j < 8; j++) {
    t[j] = 0;
    for (i = 0; i < 8; i++) t[j] += DCT8_BASIS[i][j]*_y[i];
  }
  for (j = 0; j < 8; j++) _x[j] = t[j];
}

static void ieee1180_print_results8(long _sumerrs[8][8], long _sumsqerrs[8][8],
 int _maxerr[8][8], int _l, int _h, int _sign) {
  double max;
  double total;
  int    m;
  int    i;
  int    j;
  printf("IEEE1180-1990 test results:\n");
  printf("Input range: [%i,%i]\n", -_l, _h);
  printf("Sign: %i\n", _sign);
  printf("Iterations: %i\n\n", IEEE1180_NBLOCKS);
  printf("Peak absolute values of errors:\n");
  for (i = 0, m = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      if (_maxerr[i][j] > m) m = _maxerr[i][j];
      printf("%4i", _maxerr[i][j]);
    }
    printf("\n");
  }
  printf("Worst peak error = %i (%s spec limit 1)\n\n", m,
   ieee1180_meets((double)m, 1.0));
  printf("Mean square errors:\n");
  max = total = 0;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      double err;
      err = _sumsqerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += _sumsqerrs[i][j];
      if (max < err) max = err;
    }
    printf("\n");
  }
  printf("Worst pmse = %.6f (%s spec limit 0.06)\n", max,
   ieee1180_meets(max, 0.06));
  total /= 8*8*(double)IEEE1180_NBLOCKS;
  printf("Overall mse = %.6f (%s spec limit 0.02)\n\n", total,
   ieee1180_meets(total, 0.02));
  printf("Mean errors:\n");
  max = total = 0;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      double err;
      err = _sumerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += _sumerrs[i][j];
      if (err < 0) err = -err;
      if (max < err) max = err;
    }
    printf("\n");
  }
  total /= 8*8*(double)IEEE1180_NBLOCKS;
  printf("Worst mean error = %.6f (%s spec limit 0.015)\n", max,
   ieee1180_meets(max, 0.015));
  printf("Overall mean error = %.6f (%s spec limit 0.0015)\n\n", total,
   ieee1180_meets(total, 0.0015));
}

static void ieee1180_test_block8(long _sumerrs[8][8], long _sumsqerrs[8][8],
 int _maxerr[8][8], int _l, int _h, int _sign) {
  od_coeff block[8][8];
  od_coeff refcoefs[8][8];
  od_coeff refout[8][8];
  od_coeff testout[8][8];
  double   floatcoefs[8][8];
  int      maxerr;
  int      i;
  int      j;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) block[i][j] = ieee1180_rand(_l, _h)*_sign;
  }
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  od_bin_fdct8x8(refcoefs[0], 8, block[0], 8);
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for (i = 0; i < 8; i++)
    for (j = 0; j < 8; j++) {
      /*Modification of IEEE1180: inputs to reference iDCT are scaled to match
         the scaling factors introduced by the forward integer transform.*/
      floatcoefs[i][j] = refcoefs[i][j]/(DCT8_ISCALE[i]*DCT8_ISCALE[j]);
    }
  for (i = 0; i < 8; i++) idct8(floatcoefs[i], floatcoefs[i]);
  for (j = 0; j < 8; j++) {
    double x[8];
    for (i = 0; i < 8; i++) x[i] = floatcoefs[i][j];
    idct8(x, x);
    for (i = 0; i < 8; i++) floatcoefs[i][j] = x[i];
  }
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      refout[i][j] = (od_coeff)(floatcoefs[i][j]+0.5);
      if (refout[i][j] > 255) refout[i][j] = 255;
      else if (refout[i][j] < -256) refout[i][j] = -256;
    }
  }
  od_bin_idct8x8(testout[0], 8, refcoefs[0], 8);
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      if (testout[i][j] != block[i][j]) {
        printf("Mismatch:\n");
        printf("in:\n");
        for (i = 0; i < 8; i++) {
          printf("       ");
          for (j = 0; j < 8; j++) printf(" %i", block[i][j]);
          printf("\n");
        }
        printf("xform:\n");
        for (i = 0; i < 8; i++) {
          printf("       ");
          for (j = 0; j < 8; j++) printf(" %i", refcoefs[i][j]);
          printf("\n");
        }
        printf("out:\n");
        for (i = 0; i < 8; i++) {
          printf("       ");
          for (j = 0; j < 8; j++) printf(" %i", testout[i][j]);
          printf("\n");
        }
        printf("\n");
        return;
      }
      if (testout[i][j] > 255) testout[i][j] = 255;
      else if (testout[i][j] < -256) testout[i][j] = -256;
    }
  }
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      int err;
      err = testout[i][j]-refout[i][j];
      _sumerrs[i][j] += err;
      _sumsqerrs[i][j] += err*err;
      if (err < 0) err = -err;
      if (_maxerr[i][j] < err) _maxerr[i][j] = err;
    }
  }
  for (i = 0, maxerr = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      int err;
      err = testout[i][j]-refout[i][j];
      if (err < 0) err = -err;
      if (err > maxerr) maxerr = err;
    }
  }
  /*if(maxerr>1){
    int u;
    int v;
    printf("Excessive peak error: %i\n",maxerr);
    printf("Input:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",block[u][v]);
      }
      printf("\n");
    }
    printf("Forward transform coefficients:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",refcoefs[u][v]);
      }
      printf("\n");
    }
    printf("Reference inverse:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",refout[u][v]);
      }
      printf("\n");
    }
    printf("Integerized inverse:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",testout[u][v]);
      }
      printf("\n");
    }
  }*/
}

static void ieee1180_test8(void) {
  long sumerrs[8][8];
  long sumsqerrs[8][8];
  int  maxerr[8][8];
  int  i;
  int  j;
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block8(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], 1);
    }
    ieee1180_print_results8(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], 1);
  }
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block8(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], -1);
    }
    ieee1180_print_results8(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], -1);
  }
}

static void print_basis8(double _basis[8][8]) {
  int i;
  int j;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      printf("%8.5lf%c", _basis[i][j], j == 8-1 ? '\n' : ' ');
    }
  }
}

static void compute_fbasis8(double _basis[8][8]) {
  int i;
  int j;
  for (i = 0; i < 8; i++) {
    od_coeff x[8];
    for (j = 0; j < 8; j++) x[j] = (i == j)*256;
    od_bin_fdct8(x, x, 1);
    for (j = 0; j < 8; j++) _basis[j][i] = x[j]/256.0;
  }
}

static void compute_ibasis8(double _basis[8][8]) {
  int i;
  int j;
  for (i = 0; i < 8; i++) {
    od_coeff x[8];
    for (j = 0; j < 8; j++) x[j] = (i == j)*256;
    od_bin_idct8(x, 1, x);
    for (j = 0; j < 8; j++) _basis[j][i] = x[j]/256.0;
  }
}

static void compute_ftrue_basis8(double _basis[8][8]) {
  int i;
  int j;
  for (j = 0; j < 8; j++) {
    for (i = 0; i < 8; i++) {
      _basis[j][i] = sqrt(2.0/8)*cos((i+0.5)*j*M_PI/8)*DCT8_ISCALE[j];
      if (j == 0) _basis[j][i] *= M_SQRT1_2;
    }
  }
}

static void compute_itrue_basis8(double _basis[8][8]) {
  int i;
  int j;
  for (i = 0; i < 8; i++) {
    double x[8];
    for (j = 0; j < 8; j++) x[j] = (i == j);
    idct8(x, x);
    for (j = 0; j < 8; j++) _basis[j][i] = x[j]/DCT8_ISCALE[i];
  }
}

static double compute_mse8(double _basis[8][8], double _tbasis[8][8]) {
  double e[8][8];
  double ret;
  int    i;
  int    j;
  int    k;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      e[i][j] = 0;
      for (k = 0; k < 8; k++) {
        e[i][j] += (_basis[i][k]-_tbasis[i][k])*AUTOCORR[k-j+15];
      }
    }
  }
  ret = 0;
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      ret += e[i][j]*(_basis[i][j]-_tbasis[i][j]);
    }
  }
  return ret/8;
}

static void check_bias8() {
  double rtacc[8];
  double facc[8];
  double q8acc[8];
  double q7acc[8];
  int    i;
  int    j;
  for (j = 0; j < 8; j++) q7acc[j] = q8acc[j] = facc[j] = rtacc[j] = 0;
  ieee1180_srand(1);
  for (i = 0; i < 10000000; i++) {
    od_coeff x[8];
    od_coeff x2[8];
    od_coeff y[8];
    od_coeff y2[8];
    for (j = 0; j < 8; j++) x[j] = ieee1180_rand(255, 255);
    od_bin_fdct8(y, x, 1);
    for (j = 0; j < 8; j++) facc[j] += y[j];
    od_bin_idct8(x2, 1, y);
    for (j = 0; j < 8; j++)
      if (x[j] != x2[j]) {
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < 8; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < 8; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < 8; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    for (j = 0; j < 8; j++) y2[j] = y[j]+ieee1180_rand(1, 1);
    od_bin_idct8(x2, 1, y2);
    for (j = 0; j < 8; j++) rtacc[j] += x2[j]-x[j];
    for (j = 0; j < 8; j++) y2[j] = y[j]/8<<3;
    od_bin_idct8(x2, 1, y2);
    for (j = 0; j < 8; j++) q8acc[j] += x2[j]-x[j];
    for (j = 0; j < 8; j++) y2[j] = y[j]/7*7;
    od_bin_idct8(x2, 1, y2);
    for (j = 0; j < 8; j++) q7acc[j] += x2[j]-x[j];
  }
  printf("1-D Forward Bias:\n");
  for (j = 0; j < 8; j++) printf("% -18.15G%s", facc[j]/i,
     (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Round-Trip Bias:\n");
  for (j = 0; j < 8; j++) printf("% -18.15G%s", rtacc[j]/i,
     (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Q=8 Bias:\n");
  for (j = 0; j < 8; j++) printf("% -18.15G%s", q8acc[j]/i,
     (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Q=7 Bias:\n");
  for (j = 0; j < 8; j++) printf("% -18.15G%s", q7acc[j]/i,
     (j&3) == 3 ? "\n" : "  ");
  printf("\n");
}

# if 0
static void bin_fxform_2d8(od_coeff _x[8*2][8*2]) {
  od_coeff y[8*2];
  int      u;
  int      v;
  /*Perform pre-filtering.*/
  for (u = 0; u < 8*2; u++) {
    od_pre_filter8(_x[u], _x[u]);
    od_pre_filter8(_x[u]+8, _x[u]+8);
  }
  for (v = 0; v < 8*2; v++) {
    for (u = 0; u < 8*2; u++) y[u] = _x[u][v];
    od_pre_filter8(y, y);
    od_pre_filter8(y+8, y+8);
    for (u = 0; u < 8*2; u++) _x[u][v] = y[u];
  }
  /*Perform DCT.*/
  for (u = 8/2; u < 8*3/2; u++) od_bin_fdct8(_x[u]+8/2, _x[u]+8/2, 1);
  for (v = 8/2; v < 8*3/2; v++) {
    for (u = 8/2; u < 8*3/2; u++) y[u] = _x[u][v];
    od_bin_fdct8(y+8/2, y+8/2, 1);
    for (u = 8/2; u < 8*3/2; u++) _x[u][v] = y[u];
  }
}

static void dynamic_range8(void) {
  double   basis2[8][8][8*2][8*2];
  od_coeff min2[8][8];
  od_coeff max2[8][8];
  int      i;
  int      j;
  int      u;
  int      v;
  for (i = 0; i < 8*2; i++) {
    for (j = 0; j < 8*2; j++) {
      od_coeff x[8*2][8*2];
      /*Generate impulse.*/
      for (u = 0; u < 8*2; u++) {
        for (v = 0; v < 8*2; v++) {
          x[u][v] = (u == i && v == j)*256;
        }
      }
      bin_fxform_2d8(x);
      /*Retrieve basis elements.*/
      for (u = 0; u < 8; u++) {
        for (v = 0; v < 8; v++) {
          basis2[u][v][i][j] = x[u+8/2][v+8/2]/256.0;
        }
      }
    }
  }
  for (u = 0; u < 8; u++) {
    for (v = 0; v < 8; v++) {
      od_coeff x[8*2][8*2];
      for (i = 0; i < 8*2; i++) {
        for (j = 0; j < 8*2; j++) {
          x[i][j] = basis2[u][v][i][j] < 0 ? -255 : 255;
        }
      }
      bin_fxform_2d8(x);
      max2[u][v] = x[u+8/2][v+8/2];
      for (i = 0; i < 8*2; i++) {
        for (j = 0; j < 8*2; j++) {
          x[i][j] = basis2[u][v][i][j] > 0 ? -255 : 255;
        }
      }
      bin_fxform_2d8(x);
      min2[u][v] = x[u+8/2][v+8/2];
    }
  }
  printf("2-D ranges:\n");
  for (u = 0; u < 8; u++) {
    printf("Min %2i:", u);
    for (v = 0; v < 8; v++) printf(" %6i", min2[u][v]);
    printf("\nMax %2i:", u);
    for (v = 0; v < 8; v++) printf(" %6i", max2[u][v]);
    printf("\n");
  }
}
# endif

static void check8(void) {
  od_coeff min[8];
  od_coeff max[8];
  double   basis[8][8];
  double   tbasis[8][8];
  int      i;
  int      j;
  /*dynamic_range8();*/
  for (j = 0; j < 8; j++) min[j] = max[j] = 0;
  for (i = 0; i < 1<<8; i++) {
    od_coeff x[8];
    od_coeff y[8];
    od_coeff x2[8];
    for (j = 0; j < 8; j++) x[j] = (i>>j&1) ? 255 : -256;
    od_bin_fdct8(y, x, 1);
    od_bin_idct8(x2, 1, y);
    for (j = 0; j < 8; j++) {
      if (y[j] < min[j]) min[j] = y[j];
      else if (y[j] > max[j]) max[j] = y[j];
    }
    for (j = 0; j < 8; j++) {
      if (x[j] != x2[j]) {
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < 8; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < 8; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < 8; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    }
  }
  printf("Min:");
  for (j = 0; j < 8; j++) printf(" %5i", min[j]);
  printf("\nMax:");
  for (j = 0; j < 8; j++) printf(" %5i", max[j]);
  printf("\nod_bin_idct8 basis:\n");
  compute_ibasis8(basis);
  print_basis8(basis);
  printf("Scaled type-II iDCT basis:\n");
  compute_itrue_basis8(tbasis);
  print_basis8(tbasis);
  printf("\nod_bin_fdct8 basis:\n");
  compute_fbasis8(basis);
  print_basis8(basis);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis8(tbasis);
  print_basis8(tbasis);
  printf("MSE: %.32lg\n\n", compute_mse8(basis, tbasis));
  ieee1180_test8();
  check_bias8();
}


/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT16_ISCALE[16] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

/*The true forward 16-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT16_BASIS[16][16] = {
  {
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25,
     0.25,                                0.25
  },{
     0.35185093438159561475651955574599,  0.33832950029358816956728612239141,
     0.31180625324666780814299756569942,  0.27330046675043937205975812924953,
     0.22429189658565907106266034730521,  0.16666391461943662432490137779072,
     0.10263113188058934528909230943878,  0.034654292299772865648933749133244,
    -0.034654292299772865648933749133244,-0.10263113188058934528909230943878,
    -0.16666391461943662432490137779072, -0.22429189658565907106266034730521,
    -0.27330046675043937205975812924953, -0.31180625324666780814299756569942,
    -0.33832950029358816956728612239141, -0.35185093438159561475651955574599
  },{
     0.34675996133053686545540479789161,  0.29396890060483967924361677615282,
     0.19642373959677554531947434191430,  0.068974844820735753083989390917343,
    -0.068974844820735753083989390917343,-0.19642373959677554531947434191430,
    -0.29396890060483967924361677615282, -0.34675996133053686545540479789161,
    -0.34675996133053686545540479789161, -0.29396890060483967924361677615282,
    -0.19642373959677554531947434191430, -0.068974844820735753083989390917343,
     0.068974844820735753083989390917343, 0.19642373959677554531947434191430,
     0.29396890060483967924361677615282,  0.34675996133053686545540479789161
  },{
     0.33832950029358816956728612239141,  0.22429189658565907106266034730521,
     0.034654292299772865648933749133244,-0.16666391461943662432490137779072,
    -0.31180625324666780814299756569942, -0.35185093438159561475651955574599,
    -0.27330046675043937205975812924953, -0.10263113188058934528909230943878,
     0.10263113188058934528909230943878,  0.27330046675043937205975812924953,
     0.35185093438159561475651955574599,  0.31180625324666780814299756569942,
     0.16666391461943662432490137779072, -0.034654292299772865648933749133244,
    -0.22429189658565907106266034730521, -0.33832950029358816956728612239141
  },{
     0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
     0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
     0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
     0.13529902503654924609993080134160,  0.32664074121909413196416079335680
  },{
     0.31180625324666780814299756569942,  0.034654292299772865648933749133244,
    -0.27330046675043937205975812924953, -0.33832950029358816956728612239141,
    -0.10263113188058934528909230943878,  0.22429189658565907106266034730521,
     0.35185093438159561475651955574599,  0.16666391461943662432490137779072,
    -0.16666391461943662432490137779072, -0.35185093438159561475651955574599,
    -0.22429189658565907106266034730521,  0.10263113188058934528909230943878,
     0.33832950029358816956728612239141,  0.27330046675043937205975812924953,
    -0.034654292299772865648933749133244,-0.31180625324666780814299756569942
  },{
     0.29396890060483967924361677615282, -0.068974844820735753083989390917343,
    -0.34675996133053686545540479789161, -0.19642373959677554531947434191430,
     0.19642373959677554531947434191430,  0.34675996133053686545540479789161,
     0.068974844820735753083989390917343,-0.29396890060483967924361677615282,
    -0.29396890060483967924361677615282,  0.068974844820735753083989390917343,
     0.34675996133053686545540479789161,  0.19642373959677554531947434191430,
    -0.19642373959677554531947434191430, -0.34675996133053686545540479789161,
    -0.068974844820735753083989390917343, 0.29396890060483967924361677615282
  },{
     0.27330046675043937205975812924953, -0.16666391461943662432490137779072,
    -0.33832950029358816956728612239141,  0.034654292299772865648933749133244,
     0.35185093438159561475651955574599,  0.10263113188058934528909230943878,
    -0.31180625324666780814299756569942, -0.22429189658565907106266034730521,
     0.22429189658565907106266034730521,  0.31180625324666780814299756569942,
    -0.10263113188058934528909230943878, -0.35185093438159561475651955574599,
    -0.034654292299772865648933749133244, 0.33832950029358816956728612239141,
     0.16666391461943662432490137779072, -0.27330046675043937205975812924953
  },{
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25,
     0.25,                               -0.25,
    -0.25,                                0.25
  },{
     0.22429189658565907106266034730521, -0.31180625324666780814299756569942,
    -0.10263113188058934528909230943878,  0.35185093438159561475651955574599,
    -0.034654292299772865648933749133244,-0.33832950029358816956728612239141,
     0.16666391461943662432490137779072,  0.27330046675043937205975812924953,
    -0.27330046675043937205975812924953, -0.16666391461943662432490137779072,
     0.33832950029358816956728612239141,  0.034654292299772865648933749133244,
    -0.35185093438159561475651955574599,  0.10263113188058934528909230943878,
     0.31180625324666780814299756569942, -0.22429189658565907106266034730521
  },{
     0.19642373959677554531947434191430, -0.34675996133053686545540479789161,
     0.068974844820735753083989390917343, 0.29396890060483967924361677615282,
    -0.29396890060483967924361677615282, -0.068974844820735753083989390917343,
     0.34675996133053686545540479789161, -0.19642373959677554531947434191430,
    -0.19642373959677554531947434191430,  0.34675996133053686545540479789161,
    -0.068974844820735753083989390917343,-0.29396890060483967924361677615282,
     0.29396890060483967924361677615282,  0.068974844820735753083989390917343,
    -0.34675996133053686545540479789161,  0.19642373959677554531947434191430
  },{
     0.16666391461943662432490137779072, -0.35185093438159561475651955574599,
     0.22429189658565907106266034730521,  0.10263113188058934528909230943878,
    -0.33832950029358816956728612239141,  0.27330046675043937205975812924953,
     0.034654292299772865648933749133244,-0.31180625324666780814299756569942,
     0.31180625324666780814299756569942, -0.034654292299772865648933749133244,
    -0.27330046675043937205975812924953,  0.33832950029358816956728612239141,
    -0.10263113188058934528909230943878, -0.22429189658565907106266034730521,
     0.35185093438159561475651955574599, -0.16666391461943662432490137779072
  },{
     0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
     0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680,  0.13529902503654924609993080134160,
     0.13529902503654924609993080134160, -0.32664074121909413196416079335680,
     0.32664074121909413196416079335680, -0.13529902503654924609993080134160,
    -0.13529902503654924609993080134160,  0.32664074121909413196416079335680,
    -0.32664074121909413196416079335680,  0.13529902503654924609993080134160
  },{
     0.10263113188058934528909230943878, -0.27330046675043937205975812924953,
     0.35185093438159561475651955574599, -0.31180625324666780814299756569942,
     0.16666391461943662432490137779072,  0.034654292299772865648933749133244,
    -0.22429189658565907106266034730521,  0.33832950029358816956728612239141,
    -0.33832950029358816956728612239141,  0.22429189658565907106266034730521,
    -0.034654292299772865648933749133244,-0.16666391461943662432490137779072,
     0.31180625324666780814299756569942, -0.35185093438159561475651955574599,
     0.27330046675043937205975812924953, -0.10263113188058934528909230943878
  },{
     0.068974844820735753083989390917343,-0.19642373959677554531947434191430,
     0.29396890060483967924361677615282, -0.34675996133053686545540479789161,
     0.34675996133053686545540479789161, -0.29396890060483967924361677615282,
     0.19642373959677554531947434191430, -0.068974844820735753083989390917343,
    -0.068974844820735753083989390917343, 0.19642373959677554531947434191430,
    -0.29396890060483967924361677615282,  0.34675996133053686545540479789161,
    -0.34675996133053686545540479789161,  0.29396890060483967924361677615282,
    -0.19642373959677554531947434191430,  0.068974844820735753083989390917343
  },{
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

void idct16(double _x[], const double _y[]) {
  double t[16];
  int    i;
  int    j;
  for (j = 0; j < 16; j++) {
    t[j] = 0;
    for (i = 0; i < 16; i++) t[j] += DCT16_BASIS[i][j]*_y[i];
  }
  for (j = 0; j < 16; j++) _x[j] = t[j];
}

static void ieee1180_print_results16(long _sumerrs[16][16],
 long _sumsqerrs[16][16], int _maxerr[16][16], int _l, int _h, int _sign) {
  double max;
  double total;
  int    m;
  int    i;
  int    j;
  printf("IEEE1180-1990 test results:\n");
  printf("Input range: [%i,%i]\n", -_l, _h);
  printf("Sign: %i\n", _sign);
  printf("Iterations: %i\n\n", IEEE1180_NBLOCKS);
  printf("Peak absolute values of errors:\n");
  for (i = 0, m = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      if (_maxerr[i][j] > m) m = _maxerr[i][j];
      printf("%4i", _maxerr[i][j]);
    }
    printf("\n");
  }
  printf("Worst peak error = %i (%s spec limit 1)\n\n", m,
   ieee1180_meets((double)m, 1.0));
  printf("Mean square errors:\n");
  max = total = 0;
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      double err;
      err = _sumsqerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += _sumsqerrs[i][j];
      if (max < err) max = err;
    }
    printf("\n");
  }
  printf("Worst pmse = %.6f (%s spec limit 0.06)\n", max,
   ieee1180_meets(max, 0.06));
  total /= 16*16*(double)IEEE1180_NBLOCKS;
  printf("Overall mse = %.6f (%s spec limit 0.02)\n\n", total,
   ieee1180_meets(total, 0.02));
  printf("Mean errors:\n");
  max = total = 0;
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      double err;
      err = _sumerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f", err);
      total += _sumerrs[i][j];
      if (err < 0) err = -err;
      if (max < err) max = err;
    }
    printf("\n");
  }
  total /= 16*16*(double)IEEE1180_NBLOCKS;
  printf("Worst mean error = %.6f (%s spec limit 0.015)\n", max,
   ieee1180_meets(max, 0.015));
  printf("Overall mean error = %.6f (%s spec limit 0.0015)\n\n", total,
   ieee1180_meets(total, 0.0015));
}

static void ieee1180_test_block16(long _sumerrs[16][16],
 long _sumsqerrs[16][16], int _maxerr[16][16], int _l, int _h, int _sign) {
  od_coeff block[16][16];
  od_coeff refcoefs[16][16];
  od_coeff refout[16][16];
  od_coeff testout[16][16];
  double   floatcoefs[16][16];
  int      maxerr;
  int      i;
  int      j;
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) block[i][j] = ieee1180_rand(_l, _h)*_sign;
  }
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  od_bin_fdct16x16(refcoefs[0], 16, block[0], 16);
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      /*Modification of IEEE1180: inputs to reference iDCT are scaled to match
         the scaling factors introduced by the forward integer transform.*/
      floatcoefs[i][j] = refcoefs[i][j]/(DCT16_ISCALE[i]*DCT16_ISCALE[j]);
    }
  }
  for (i = 0; i < 16; i++) idct16(floatcoefs[i], floatcoefs[i]);
  for (j = 0; j < 16; j++) {
    double x[16];
    for (i = 0; i < 16; i++) x[i] = floatcoefs[i][j];
    idct16(x, x);
    for (i = 0; i < 16; i++) floatcoefs[i][j] = x[i];
  }
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      refout[i][j] = (od_coeff)(floatcoefs[i][j]+0.5);
      if (refout[i][j] > 255) refout[i][j] = 255;
      else if (refout[i][j] < -256) refout[i][j] = -256;
    }
  }
  od_bin_idct16x16(testout[0], 16, refcoefs[0], 16);
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      if (testout[i][j] != block[i][j]) {
        printf("Mismatch:\n");
        printf("in:\n");
        for (i = 0; i < 16; i++) {
          printf("       ");
          for (j = 0; j < 16; j++) printf(" %i", block[i][j]);
          printf("\n");
        }
        printf("xform:\n");
        for (i = 0; i < 16; i++) {
          printf("       ");
          for (j = 0; j < 16; j++) printf(" %i", refcoefs[i][j]);
          printf("\n");
        }
        printf("out:\n");
        for (i = 0; i < 16; i++) {
          printf("       ");
          for (j = 0; j < 16; j++) printf(" %i", testout[i][j]);
          printf("\n");
        }
        printf("\n");
        return;
      }
      if (testout[i][j] > 255) testout[i][j] = 255;
      else if (testout[i][j] < -256) testout[i][j] = -256;
    }
  }
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      int err;
      err = testout[i][j]-refout[i][j];
      _sumerrs[i][j] += err;
      _sumsqerrs[i][j] += err*err;
      if (err < 0) err = -err;
      if (_maxerr[i][j] < err) _maxerr[i][j] = err;
    }
  }
  for (i = 0, maxerr = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      int err;
      err = testout[i][j]-refout[i][j];
      if (err < 0) err = -err;
      if (err > maxerr) maxerr = err;
    }
  }
  /*if(maxerr>1){
    int u;
    int v;
    printf("Excessive peak error: %i\n",maxerr);
    printf("Input:\n");
    for(u=0;u<16;u++){
      for(v=0;v<16;v++){
        printf("%5i",block[u][v]);
      }
      printf("\n");
    }
    printf("Forward transform coefficients:\n");
    for(u=0;u<16;u++){
      for(v=0;v<16;v++){
        printf("%5i",refcoefs[u][v]);
      }
      printf("\n");
    }
    printf("Reference inverse:\n");
    for(u=0;u<16;u++){
      for(v=0;v<16;v++){
        printf("%5i",refout[u][v]);
      }
      printf("\n");
    }
    printf("Integerized inverse:\n");
    for(u=0;u<16;u++){
      for(v=0;v<16;v++){
        printf("%5i",testout[u][v]);
      }
      printf("\n");
    }
  }*/
}

static void ieee1180_test16(void) {
  long sumerrs[16][16];
  long sumsqerrs[16][16];
  int  maxerr[16][16];
  int  i;
  int  j;
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block16(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], 1);
    }
    ieee1180_print_results16(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], 1);
  }
  ieee1180_srand(1);
  for (i = 0; i < IEEE1180_NRANGES; i++) {
    memset(sumerrs, 0, sizeof(sumerrs));
    memset(sumsqerrs, 0, sizeof(sumsqerrs));
    memset(maxerr, 0, sizeof(maxerr));
    for (j = 0; j < IEEE1180_NBLOCKS; j++) {
      ieee1180_test_block16(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
       IEEE1180_H[i], -1);
    }
    ieee1180_print_results16(sumerrs, sumsqerrs, maxerr, IEEE1180_L[i],
     IEEE1180_H[i], -1);
  }
}

static void print_basis16(double _basis[16][16]) {
  int i;
  int j;
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      printf("%8.5lf%c", _basis[i][j], j == 16-1 ? '\n' : ' ');
    }
  }
}

static void compute_fbasis16(double _basis[16][16]) {
  int i;
  int j;
  for (i = 0; i < 16; i++) {
    od_coeff x[16];
    for (j = 0; j < 16; j++) x[j] = (i == j)*256;
    od_bin_fdct16(x, x, 1);
    for (j = 0; j < 16; j++) _basis[j][i] = x[j]/256.0;
  }
}

static void compute_ibasis16(double _basis[16][16]) {
  int i;
  int j;
  for (i = 0; i < 16; i++) {
    od_coeff x[16];
    for (j = 0; j < 16; j++) x[j] = (i == j)*256;
    od_bin_idct16(x, 1, x);
    for (j = 0; j < 16; j++) _basis[j][i] = x[j]/256.0;
  }
}

static void compute_ftrue_basis16(double _basis[16][16]) {
  int i;
  int j;
  for (j = 0; j < 16; j++) {
    for (i = 0; i < 16; i++) {
      _basis[j][i] = sqrt(2.0/16)*cos((i+0.5)*j*M_PI/16)*DCT16_ISCALE[j];
      if (j == 0) _basis[j][i] *= M_SQRT1_2;
    }
  }
}

static void compute_itrue_basis16(double _basis[16][16]) {
  int i;
  int j;
  for (i = 0; i < 16; i++) {
    double x[16];
    for (j = 0; j < 16; j++) x[j] = (i == j);
    idct16(x, x);
    for (j = 0; j < 16; j++) _basis[j][i] = x[j]/DCT16_ISCALE[i];
  }
}

static double compute_mse16(double _basis[16][16], double _tbasis[16][16]) {
  double e[16][16];
  double ret;
  int    i;
  int    j;
  int    k;
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      e[i][j] = 0;
      for (k = 0; k < 16; k++) {
        e[i][j] += (_basis[i][k]-_tbasis[i][k])*AUTOCORR[k-j+15];
      }
    }
  }
  ret = 0;
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      ret += e[i][j]*(_basis[i][j]-_tbasis[i][j]);
    }
  }
  return ret/16;
}

static void check_bias16() {
  double rtacc[16];
  double facc[16];
  double q8acc[16];
  double q7acc[16];
  int    i;
  int    j;
  for (j = 0; j < 16; j++) q7acc[j] = q8acc[j] = facc[j] = rtacc[j] = 0;
  ieee1180_srand(1);
  for (i = 0; i < 10000000; i++) {
    od_coeff x[16];
    od_coeff x2[16];
    od_coeff y[16];
    od_coeff y2[16];
    for (j = 0; j < 16; j++) x[j] = ieee1180_rand(255, 255);
    od_bin_fdct16(y, x, 1);
    for (j = 0; j < 16; j++) facc[j] += y[j];
    od_bin_idct16(x2, 1, y);
    for (j = 0; j < 16; j++)
      if (x[j] != x2[j]) {
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < 16; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < 16; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < 16; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
    for (j = 0; j < 16; j++) y2[j] = y[j]+ieee1180_rand(1, 1);
    od_bin_idct16(x2, 1, y2);
    for (j = 0; j < 16; j++) rtacc[j] += x2[j]-x[j];
    for (j = 0; j < 16; j++) y2[j] = y[j]/8<<3;
    od_bin_idct16(x2, 1, y2);
    for (j = 0; j < 16; j++) q8acc[j] += x2[j]-x[j];
    for (j = 0; j < 16; j++) y2[j] = y[j]/7*7;
    od_bin_idct16(x2, 1, y2);
    for (j = 0; j < 16; j++) q7acc[j] += x2[j]-x[j];
  }
  printf("1-D Forward Bias:\n");
  for (j = 0; j < 16; j++)
    printf("% -18.15G%s", facc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Round-Trip Bias:\n");
  for (j = 0; j < 16; j++)
     printf("% -18.15G%s", rtacc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Q=8 Bias:\n");
  for (j = 0; j < 16; j++)
     printf("% -18.15G%s", q8acc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
  printf("1-D Q=7 Bias:\n");
  for (j = 0; j < 16; j++)
     printf("% -18.15G%s", q7acc[j]/i, (j&3) == 3 ? "\n" : "  ");
  printf("\n");
}

# if 0
static void bin_fxform_2d16(od_coeff _x[16*2][16*2]) {
  od_coeff y[16*2];
  int      u;
  int      v;
  /*Perform pre-filtering.*/
  for (u = 0; u < 16*2; u++) {
    od_pre_filter16(_x[u], _x[u]);
    od_pre_filter16(_x[u]+16, _x[u]+16);
  }
  for (v = 0; v < 16*2; v++) {
    for (u = 0; u < 16*2; u++) y[u] = _x[u][v];
    od_pre_filter16(y, y);
    od_pre_filter16(y+16, y+16);
    for (u = 0; u < 16*2; u++) _x[u][v] = y[u];
  }
  /*Perform DCT.*/
  for (u = 16/2; u < 16*3/2; u++) od_bin_fdct16(_x[u]+16/2, _x[u]+16/2, 1);
  for (v = 16/2; v < 16*3/2; v++) {
    for (u = 16/2; u < 16*3/2; u++) y[u] = _x[u][v];
    od_bin_fdct16(y+16/2, y+16/2, 1);
    for (u = 16/2; u < 16*3/2; u++) _x[u][v] = y[u];
  }
}

static void dynamic_range16(void) {
  double   basis2[16][16][16*2][16*2];
  od_coeff min2[16][16];
  od_coeff max2[16][16];
  int      i;
  int      j;
  int      u;
  int      v;
  for (i = 0; i < 16*2; i++) {
    for (j = 0; j < 16*2; j++) {
      od_coeff x[16*2][16*2];
      /*Generate impulse.*/
      for (u = 0; u < 16*2; u++) {
        for (v = 0; v < 16*2; v++) {
          x[u][v] = (u == i && v == j)*256;
        }
      }
      bin_fxform_2d16(x);
      /*Retrieve basis elements.*/
      for (u = 0; u < 16; u++) {
        for (v = 0; v < 16; v++) {
          basis2[u][v][i][j] = x[u+16/2][v+16/2]/256.0;
        }
      }
    }
  }
  for (u = 0; u < 16; u++) {
    for (v = 0; v < 16; v++) {
      od_coeff x[16*2][16*2];
      for (i = 0; i < 16*2; i++) {
        for (j = 0; j < 16*2; j++) {
          x[i][j] = basis2[u][v][i][j] < 0 ? -255 : 255;
        }
      }
      bin_fxform_2d16(x);
      max2[u][v] = x[u+16/2][v+16/2];
      for (i = 0; i < 16*2; i++) {
        for (j = 0; j < 16*2; j++) {
          x[i][j] = basis2[u][v][i][j] > 0 ? -255 : 255;
        }
      }
      bin_fxform_2d16(x);
      min2[u][v] = x[u+16/2][v+16/2];
    }
  }
  printf("2-D ranges:\n");
  for (u = 0; u < 16; u++) {
    printf("Min %2i:", u);
    for (v = 0; v < 16; v++) printf(" %6i", min2[u][v]);
    printf("\nMax %2i:", u);
    for (v = 0; v < 16; v++) printf(" %6i", max2[u][v]);
    printf("\n");
  }
}
# endif

static void check16(void) {
  od_coeff min[16];
  od_coeff max[16];
  double   basis[16][16];
  double   tbasis[16][16];
  int      i;
  int      j;
  /*dynamic_range16();*/
  for (j = 0; j < 16; j++) min[j] = max[j] = 0;
  for (i = 0; i < 1<<16; i++) {
    od_coeff x[16];
    od_coeff y[16];
    od_coeff x2[16];
    for (j = 0; j < 16; j++) x[j] = (i>>j&1) ? 255 : -256;
    od_bin_fdct16(y, x, 1);
    od_bin_idct16(x2, 1, y);
    for (j = 0; j < 16; j++) {
      if (y[j] < min[j]) min[j] = y[j];
      else if (y[j] > max[j]) max[j] = y[j];
    }
    for (j = 0; j < 16; j++)
      if (x[j] != x2[j]) {
        printf("Mismatch:\n");
        printf("in:    ");
        for (j = 0; j < 16; j++) printf(" %i", x[j]);
        printf("\nxform: ");
        for (j = 0; j < 16; j++) printf(" %i", y[j]);
        printf("\nout:   ");
        for (j = 0; j < 16; j++) printf(" %i", x2[j]);
        printf("\n\n");
        break;
      }
  }
  printf("Min:");
  for (j = 0; j < 16; j++) printf(" %5i", min[j]);
  printf("\nMax:");
  for (j = 0; j < 16; j++) printf(" %5i", max[j]);
  printf("\nod_bin_idct16 basis:\n");
  compute_ibasis16(basis);
  print_basis16(basis);
  printf("Scaled type-II iDCT basis:\n");
  compute_itrue_basis16(tbasis);
  print_basis16(tbasis);
  printf("\nod_bin_fdct16 basis:\n");
  compute_fbasis16(basis);
  print_basis16(basis);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis16(tbasis);
  print_basis16(tbasis);
  printf("MSE: %.32lg\n\n", compute_mse16(basis, tbasis));
  ieee1180_test16();
  check_bias16();
}


int main(void) {
  check4();
  check8();
  check16();
  return 0;
}
#endif
