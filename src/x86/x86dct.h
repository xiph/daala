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

/* This file should not be compiled as its own unit. Instead, it is included
   into avx2dct.c and sse2dct.c which provide implementations of the
   od_mm256_* functions. */

OD_SIMD_INLINE void od_mm256_overflow_check(od_m256i val, ogg_int32_t scale,
 ogg_int32_t offset, int idx) {
#if defined(OD_DCT_TEST) && defined(OD_DCT_CHECK_OVERFLOW)
  ogg_int32_t mem[8];
  int n;
  od_mm256_storeu_si256((od_m256i *)mem, val);
  for (n = 0; n < 4; n++) {
    OD_DCT_OVERFLOW_CHECK(mem[n], scale, offset, idx);
  }
#endif
  (void)val;
  (void)scale;
  (void)offset;
  (void)idx;
}

OD_SIMD_INLINE od_m256i od_mm256_unbiased_rshift32(od_m256i a, int b) {
  return od_mm256_srai_epi32(od_mm256_add_epi32(
   od_mm256_srli_epi32(a, 32 - b), a), b);
}

OD_SIMD_INLINE void od_mm256_transpose8(od_m256i *t0, od_m256i *t1,
 od_m256i *t2, od_m256i *t3, od_m256i *t4, od_m256i *t5, od_m256i *t6,
 od_m256i *t7) {
  od_m256i a;
  od_m256i b;
  od_m256i c;
  od_m256i d;
  od_m256i e;
  od_m256i f;
  od_m256i g;
  od_m256i h;
  od_m256i x;
  od_m256i y;
  a = od_mm256_unpacklo_epi32(*t0, *t1);
  b = od_mm256_unpacklo_epi32(*t2, *t3);
  c = od_mm256_unpackhi_epi32(*t0, *t1);
  d = od_mm256_unpackhi_epi32(*t2, *t3);
  e = od_mm256_unpacklo_epi32(*t4, *t5);
  f = od_mm256_unpacklo_epi32(*t6, *t7);
  g = od_mm256_unpackhi_epi32(*t4, *t5);
  h = od_mm256_unpackhi_epi32(*t6, *t7);
  x = od_mm256_unpacklo_epi64(a, b);
  y = od_mm256_unpacklo_epi64(e, f);
  *t0 = od_mm256_permute2x128_si256(x, y, 0 | (2 << 4));
  *t4 = od_mm256_permute2x128_si256(x, y, 1 | (3 << 4));
  x = od_mm256_unpackhi_epi64(a, b);
  y = od_mm256_unpackhi_epi64(e, f);
  *t1 = od_mm256_permute2x128_si256(x, y, 0 | (2 << 4));
  *t5 = od_mm256_permute2x128_si256(x, y, 1 | (3 << 4));
  x = od_mm256_unpacklo_epi64(c, d);
  y = od_mm256_unpacklo_epi64(g, h);
  *t2 = od_mm256_permute2x128_si256(x, y, 0 | (2 << 4));
  *t6 = od_mm256_permute2x128_si256(x, y, 1 | (3 << 4));
  x = od_mm256_unpackhi_epi64(c, d);
  y = od_mm256_unpackhi_epi64(g, h);
  *t3 = od_mm256_permute2x128_si256(x, y, 0 | (2 << 4));
  *t7 = od_mm256_permute2x128_si256(x, y, 1 | (3 << 4));
}

OD_SIMD_INLINE void load8(const od_coeff *x, int xstride,
 od_m256i *t0, od_m256i *t1, od_m256i *t2, od_m256i *t3,
 od_m256i *t4, od_m256i *t5, od_m256i *t6, od_m256i *t7) {
  *t0 = od_mm256_loadu_si256((const od_m256i *)x);
  *t1 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
  *t2 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
  *t3 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
  *t4 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
  *t5 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
  *t6 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
  *t7 = od_mm256_loadu_si256((const od_m256i *)(x += xstride));
}

OD_SIMD_INLINE void store8(od_coeff *x, int xstride,
 od_m256i t0, od_m256i t1, od_m256i t2, od_m256i t3,
 od_m256i t4, od_m256i t5, od_m256i t6, od_m256i t7) {
  od_mm256_storeu_si256((od_m256i *)x, t0);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t1);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t2);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t3);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t4);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t5);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t6);
  od_mm256_storeu_si256((od_m256i *)(x += xstride), t7);
}

#define OD_DCT_MUL_EPI32(a, scale, offset, shift) \
  od_mm256_srai_epi32(od_mm256_add_epi32(mul_epi32_256(a, scale), \
   od_mm256_set1_epi32(offset)), shift)
#define OD_DCT_MLS_EPI32(r, a, scale, offset, shift) \
  od_mm256_sub_epi32(r, OD_DCT_MUL_EPI32(a, scale, offset, shift))
#define OD_DCT_MLA_EPI32(r, a, scale, offset, shift) \
  od_mm256_add_epi32(r, OD_DCT_MUL_EPI32(a, scale, offset, shift))

OD_SIMD_INLINE void fdct8_kernel(od_m256i *x0, od_m256i *x1, od_m256i *x2,
 od_m256i *x3, od_m256i *x4, od_m256i *x5, od_m256i *x6, od_m256i *x7) {
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
  /*Initial permutation:*/
  od_m256i t0 = *x0;
  od_m256i t4 = *x1;
  od_m256i t2 = *x2;
  od_m256i t6 = *x3;
  od_m256i t7 = *x4;
  od_m256i t3 = *x5;
  od_m256i t5 = *x6;
  od_m256i t1 = *x7;
  od_m256i t1h;
  od_m256i t4h;
  od_m256i t6h;
  /*+1/-1 butterflies:*/
  t1 = od_mm256_sub_epi32(t0, t1);
  t1h = od_mm256_unbiased_rshift32(t1, 1);
  t0 = od_mm256_sub_epi32(t0, t1h);
  t4 = od_mm256_add_epi32(t4, t5);
  t4h = od_mm256_unbiased_rshift32(t4, 1);
  t5 = od_mm256_sub_epi32(t5, t4h);
  t3 = od_mm256_sub_epi32(t2, t3);
  t2 = od_mm256_sub_epi32(t2, od_mm256_unbiased_rshift32(t3, 1));
  t6 = od_mm256_add_epi32(t6, t7);
  t6h = od_mm256_unbiased_rshift32(t6, 1);
  t7 = od_mm256_sub_epi32(t6h, t7);
  /*+ Embedded 4-point type-II DCT.*/
  t0 = od_mm256_add_epi32(t0, t6h);
  t6 = od_mm256_sub_epi32(t0, t6);
  t2 = od_mm256_sub_epi32(t4h, t2);
  t4 = od_mm256_sub_epi32(t2, t4);
  /*|-+ Embedded 2-point type-II DCT.*/
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  od_mm256_overflow_check(t4, 13573, 16384, 3);
  t0 = OD_DCT_MLS_EPI32(t0, t4, 13573, 16384, 15);
  /*11585/16384 ~= \sqrt{\frac{1}{2}} ~= 0.70710678118654752440084436210485*/
  od_mm256_overflow_check(t0, 11585, 8192, 4);
  t4 = OD_DCT_MLA_EPI32(t4, t0, 11585, 8192, 14);
  /*13573/32768 ~= \sqrt{2} - 1 ~= 0.41421356237309504880168872420970*/
  od_mm256_overflow_check(t4, 13573, 16384, 5);
  t0 = OD_DCT_MLS_EPI32(t0, t4, 13573, 16384, 15);
  /*|-+ Embedded 2-point type-IV DST.*/
  /*21895/32768 ~= \frac{1 - cos(\frac{3\pi}{8})}{\sin(\frac{3\pi}{8})} ~=
     0.66817863791929891999775768652308*/
  od_mm256_overflow_check(t2, 21895, 16384, 6);
  t6 = OD_DCT_MLS_EPI32(t6, t2, 21895, 16384, 15);
  /*15137/16384~=sin(\frac{3\pi}{8})~=0.92387953251128675612818318939679*/
  od_mm256_overflow_check(t6, 15137, 8192, 7);
  t2 = OD_DCT_MLA_EPI32(t2, t6, 15137, 8192, 14);
  /*21895/32768 ~= \frac{1 - cos(\frac{3\pi}{8})}{\sin(\frac{3\pi}{8})}~=
     0.66817863791929891999775768652308*/
  od_mm256_overflow_check(t2, 21895, 16384, 8);
  t6 = OD_DCT_MLS_EPI32(t6, t2, 21895, 16384, 15);
  /*+ Embedded 4-point type-IV DST.*/
  /*19195/32768 ~= 2 - \sqrt{2} ~= 0.58578643762690495119831127579030*/
  od_mm256_overflow_check(t5, 19195, 16384, 9);
  t3 = OD_DCT_MLA_EPI32(t3, t5, 19195, 16384, 15);
  /*11585/16384 ~= \sqrt{\frac{1}{2}} ~= 0.70710678118654752440084436210485*/
  od_mm256_overflow_check(t3, 11585, 8192, 10);
  t5 = OD_DCT_MLA_EPI32(t5, t3, 11585, 8192, 14);
  /*7489/8192 ~= \sqrt{2}-\frac{1}{2} ~= 0.91421356237309504880168872420970*/
  od_mm256_overflow_check(t5, 7489, 4096, 11);
  t3 = OD_DCT_MLS_EPI32(t3, t5, 7489, 4096, 13);
  t7 = od_mm256_sub_epi32(od_mm256_unbiased_rshift32(t5, 1), t7);
  t5 = od_mm256_sub_epi32(t5, t7);
  t3 = od_mm256_sub_epi32(t1h, t3);
  t1 = od_mm256_sub_epi32(t1, t3);
  /*3227/32768 ~= \frac{1 - cos(\frac{\pi}{16})}{sin(\frac{\pi}{16})} ~=
     0.098491403357164253077197521291327*/
  od_mm256_overflow_check(t1, 3227, 16384, 12);
  t7 = OD_DCT_MLA_EPI32(t7, t1, 3227, 16384, 15);
  /*6393/32768 ~= sin(\frac{\pi}{16}) ~= 0.19509032201612826784828486847702*/
  od_mm256_overflow_check(t7, 6393, 16384, 13);
  t1 = OD_DCT_MLS_EPI32(t1, t7, 6393, 16384, 15);
  /*3227/32768 ~= \frac{1 - cos(\frac{\pi}{16})}{sin(\frac{\pi}{16})} ~=
     0.098491403357164253077197521291327*/
  od_mm256_overflow_check(t1, 3227, 16384, 14);
  t7 = OD_DCT_MLA_EPI32(t7, t1, 3227, 16384, 15);
  /*2485/8192 ~= \frac{1 - cos(\frac{3\pi}{16})}{sin(\frac{3\pi}{16})} ~=
     0.30334668360734239167588394694130*/
  od_mm256_overflow_check(t3, 2485, 4096, 15);
  t5 = OD_DCT_MLA_EPI32(t5, t3, 2485, 4096, 13);
  /*18205/32768 ~= sin(\frac{3\pi}{16}) ~= 0.55557023301960222474283081394853*/
  od_mm256_overflow_check(t5, 18205, 16384, 16);
  t3 = OD_DCT_MLS_EPI32(t3, t5, 18205, 16384, 15);
  /*2485/8192 ~= \frac{1 - cos(\frac{3\pi}{16})}{sin(\frac{3\pi}{16})} ~=
     0.30334668360734239167588394694130*/
  od_mm256_overflow_check(t3, 2485, 4096, 17);
  t5 = OD_DCT_MLA_EPI32(t5, t3, 2485, 4096, 13);
  od_mm256_transpose8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  *x0 = t0;
  *x1 = t1;
  *x2 = t2;
  *x3 = t3;
  *x4 = t4;
  *x5 = t5;
  *x6 = t6;
  *x7 = t7;
}

OD_SIMD_INLINE void idct8_kernel(od_m256i *y0, od_m256i *y1, od_m256i *y2,
 od_m256i *y3, od_m256i *y4, od_m256i *y5, od_m256i *y6, od_m256i *y7) {
  od_m256i t0 = *y0;
  od_m256i t1 = *y1;
  od_m256i t2 = *y2;
  od_m256i t3 = *y3;
  od_m256i t4 = *y4;
  od_m256i t5 = *y5;
  od_m256i t6 = *y6;
  od_m256i t7 = *y7;
  od_m256i t1h;
  od_m256i t4h;
  od_m256i t6h;
  od_mm256_transpose8(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  t5 = OD_DCT_MLS_EPI32(t5, t3, 2485, 4096, 13);
  t3 = OD_DCT_MLA_EPI32(t3, t5, 18205, 16384, 15);
  t5 = OD_DCT_MLS_EPI32(t5, t3, 2485, 4096, 13);
  t7 = OD_DCT_MLS_EPI32(t7, t1, 3227, 16384, 15);
  t1 = OD_DCT_MLA_EPI32(t1, t7, 6393, 16384, 15);
  t7 = OD_DCT_MLS_EPI32(t7, t1, 3227, 16384, 15);
  t1 = od_mm256_add_epi32(t1, t3);
  t1h = od_mm256_unbiased_rshift32(t1, 1);
  t3 = od_mm256_sub_epi32(t1h, t3);
  t5 = od_mm256_add_epi32(t5, t7);
  t7 = od_mm256_sub_epi32(od_mm256_unbiased_rshift32(t5, 1), t7);
  t3 = OD_DCT_MLA_EPI32(t3, t5, 7489, 4096, 13);
  t5 = OD_DCT_MLS_EPI32(t5, t3, 11585, 8192, 14);
  t3 = OD_DCT_MLS_EPI32(t3, t5, 19195, 16384, 15);
  t6 = OD_DCT_MLA_EPI32(t6, t2, 21895, 16384, 15);
  t2 = OD_DCT_MLS_EPI32(t2, t6, 15137, 8192, 14);
  t6 = OD_DCT_MLA_EPI32(t6, t2, 21895, 16384, 15);
  t0 = OD_DCT_MLA_EPI32(t0, t4, 13573, 16384, 15);
  t4 = OD_DCT_MLS_EPI32(t4, t0, 11585, 8192, 14);
  t0 = OD_DCT_MLA_EPI32(t0, t4, 13573, 16384, 15);
  t4 = od_mm256_sub_epi32(t2, t4);
  t4h = od_mm256_unbiased_rshift32(t4, 1);
  t2 = od_mm256_sub_epi32(t4h, t2);
  t6 = od_mm256_sub_epi32(t0, t6);
  t6h = od_mm256_unbiased_rshift32(t6, 1);
  t0 = od_mm256_sub_epi32(t0, t6h);
  t7 = od_mm256_sub_epi32(t6h, t7);
  t6 = od_mm256_sub_epi32(t6, t7);
  t2 = od_mm256_add_epi32(t2, od_mm256_unbiased_rshift32(t3, 1));
  t3 = od_mm256_sub_epi32(t2, t3);
  t5 = od_mm256_add_epi32(t5, t4h);
  t4 = od_mm256_sub_epi32(t4, t5);
  t0 = od_mm256_add_epi32(t0, t1h);
  t1 = od_mm256_sub_epi32(t0, t1);
  *y0 = t0;
  *y1 = t4;
  *y2 = t2;
  *y3 = t6;
  *y4 = t7;
  *y5 = t3;
  *y6 = t5;
  *y7 = t1;
}

void od_bin_fdct8x8_x86(od_coeff *y, int ystride,
 const od_coeff *x, int xstride) {
  od_m256i t0;
  od_m256i t1;
  od_m256i t2;
  od_m256i t3;
  od_m256i t4;
  od_m256i t5;
  od_m256i t6;
  od_m256i t7;
#if defined(OD_CHECKASM)
  od_coeff ref[8*8];
  od_bin_fdct8x8(ref, 8, x, xstride);
#endif
  load8(x, xstride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  fdct8_kernel(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  fdct8_kernel(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  store8(y, ystride, t0, t1, t2, t3, t4, t5, t6, t7);
#if defined(OD_CHECKASM)
  od_dct_check(1, ref, y, ystride);
#endif
}

void od_bin_idct8x8_x86(od_coeff *x, int xstride,
 const od_coeff *y, int ystride) {
  od_m256i t0;
  od_m256i t1;
  od_m256i t2;
  od_m256i t3;
  od_m256i t4;
  od_m256i t5;
  od_m256i t6;
  od_m256i t7;
#if defined(OD_CHECKASM)
  od_coeff ref[8*8];
  od_bin_idct8x8(ref, 8, y, ystride);
#endif
  load8(y, ystride, &t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  idct8_kernel(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  idct8_kernel(&t0, &t1, &t2, &t3, &t4, &t5, &t6, &t7);
  store8(x, xstride, t0, t1, t2, t3, t4, t5, t6, t7);
#if defined(OD_CHECKASM)
  od_dct_check(1, ref, x, xstride);
#endif
}
