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

#include <xmmintrin.h>

#include "../dct.h"

static __inline __m128i unbiased_rshift32(__m128i a, int b) {
  return _mm_srai_epi32(_mm_add_epi32(_mm_srli_epi32(a, 32-b), a), b);
}

static __inline __m128i mul_epi32(__m128i a, int b1) {
  __m128i b = _mm_set1_epi32(b1);
#ifdef __SSE4_1__
  return _mm_mullo_epi32(a, b);
#else
  __m128i lo = _mm_mul_epu32(a, b);
  __m128i hi = _mm_mul_epu32(_mm_srli_si128(a,4), _mm_srli_si128(b,4));
  return _mm_unpacklo_epi32(_mm_shuffle_epi32(lo, _MM_SHUFFLE(0,0,2,0)), _mm_shuffle_epi32(hi, _MM_SHUFFLE(0,0,2,0)));
#endif
}

static __inline void overflow_check(__m128i val, ogg_int32_t scale, ogg_int32_t offset, int idx) {
#if defined(OD_DCT_TEST) && defined(OD_DCT_CHECK_OVERFLOW)
  ogg_int32_t mem[4];
  int n;
  _mm_store_si128((__m128i *)mem, val);
  for (n = 0; n < 4; ++n) {
    OD_DCT_OVERFLOW_CHECK(mem[n], scale, offset, idx);
  }
#endif
  (void)val;
  (void)scale;
  (void)offset;
  (void)idx;
}

static __inline __m128i M(__m128i val, ogg_int32_t scale, ogg_int32_t offset, ogg_int32_t shift) {
  return _mm_srai_epi32(_mm_add_epi32(mul_epi32(val, scale), _mm_set1_epi32(offset)), shift);
}

static __inline void transpose4(__m128i *t0, __m128i *t1, __m128i *t2, __m128i *t3) {
  __m128i a = _mm_unpacklo_epi32(*t0, *t1);
  __m128i b = _mm_unpacklo_epi32(*t2, *t3);
  __m128i c = _mm_unpackhi_epi32(*t0, *t1);
  __m128i d = _mm_unpackhi_epi32(*t2, *t3);
  *t0 = _mm_unpacklo_epi64(a, b);
  *t1 = _mm_unpackhi_epi64(a, b);
  *t2 = _mm_unpacklo_epi64(c, d);
  *t3 = _mm_unpackhi_epi64(c, d);
}

static __inline void load4(const od_coeff *x, int xstride, __m128i *t0, __m128i *t1, __m128i *t2, __m128i *t3) {
  *t0 = _mm_load_si128((const __m128i *)x);
  *t1 = _mm_load_si128((const __m128i *)(x += xstride));
  *t2 = _mm_load_si128((const __m128i *)(x += xstride));
  *t3 = _mm_load_si128((const __m128i *)(x += xstride));
}

static __inline void store4(od_coeff *x, int xstride, __m128i t0, __m128i t1, __m128i t2, __m128i t3) {
  _mm_store_si128((__m128i *)x, t0);
  _mm_store_si128((__m128i *)(x += xstride), t1);
  _mm_store_si128((__m128i *)(x += xstride), t2);
  _mm_store_si128((__m128i *)(x += xstride), t3);
}

static __inline void fdct4_kernel(__m128i *x0, __m128i *x1, __m128i *x2, __m128i *x3) {
  /*9 adds, 2 shifts, 3 "muls".*/
  /*Initial permutation:*/
  __m128i t0 = *x0;
  __m128i t2 = *x1;
  __m128i t1 = *x2;
  __m128i t3 = *x3;
  __m128i t2h;
  /*+1/-1 butterflies:*/
  t3 = _mm_sub_epi32(t0, t3);
  t2 = _mm_add_epi32(t2, t1);
  t2h = unbiased_rshift32(t2, 1);
  t1 = _mm_sub_epi32(t2h, t1);
  t0 = _mm_sub_epi32(t0, unbiased_rshift32(t3, 1));
  /*+ Embedded 2-point type-II DCT.*/
  t0 = _mm_add_epi32(t0, t2h);
  t2 = _mm_sub_epi32(t0, t2);
  /*+ Embedded 2-point type-IV DST.*/
  /*23013/32768 ~= 4*sin(\frac{\pi}{8}) - 2*tan(\frac{\pi}{8}) ~=
     0.70230660471416898931046248770220*/
  overflow_check(t1, 23013, 16384, 0);
  t3 = _mm_sub_epi32(t3, M(t1, 23013, 16384, 15));
  /*21407/32768~=\sqrt{1/2}*cos(\frac{\pi}{8}))
     ~=0.65328148243818826392832158671359*/
  overflow_check(t3, 21407, 16384, 1);
  t1 = _mm_add_epi32(t1, M(t3, 21407, 16384, 15));
  /*18293/16384 ~= 4*sin(\frac{\pi}{8}) - tan(\frac{\pi}{8}) ~=
     1.1165201670872640381121512119119*/
  overflow_check(t1, 18293, 8192, 2);
  t3 = _mm_sub_epi32(t3, M(t1, 18293, 8192, 14));
  transpose4(&t0, &t1, &t2, &t3);
  *x0 = t0;
  *x1 = t1;
  *x2 = t2;
  *x3 = t3;
}

void od_bin_fdct4x4_sse2(od_coeff y[4], int ystride, const od_coeff *x, int xstride) {
  __m128i t0, t1, t2, t3;
  load4(x, xstride, &t0, &t1, &t2, &t3);
  fdct4_kernel(&t0, &t1, &t2, &t3);
  fdct4_kernel(&t0, &t1, &t2, &t3);
  store4(y, ystride, t0, t1, t2, t3);
}

static __inline void idct4_kernel(__m128i *y0, __m128i *y1, __m128i *y2, __m128i *y3) {
  __m128i t0 = *y0;
  __m128i t1 = *y1;
  __m128i t2 = *y2;
  __m128i t3 = *y3;
  __m128i t2h;
  transpose4(&t0, &t1, &t2, &t3);
  t3 = _mm_add_epi32(t3, M(t1, 18293, 8192, 14));
  t1 = _mm_sub_epi32(t1, M(t3, 21407, 16384, 15));
  t3 = _mm_add_epi32(t3, M(t1, 23013, 16384, 15));
  t2 = _mm_sub_epi32(t0, t2);
  t2h = unbiased_rshift32(t2, 1);
  t0 = _mm_sub_epi32(t0, _mm_sub_epi32(t2h, unbiased_rshift32(t3, 1)));
  t1 = _mm_sub_epi32(t2h, t1);
  *y0 = t0;
  *y1 = _mm_sub_epi32(t2, t1);
  *y2 = t1;
  *y3 = _mm_sub_epi32(t0, t3);
}

void od_bin_idct4x4_sse2(od_coeff *x, int xstride, const od_coeff *y, int ystride) {
  __m128i t0, t1, t2, t3;
  load4(y, ystride, &t0, &t1, &t2, &t3);
  idct4_kernel(&t0, &t1, &t2, &t3);
  idct4_kernel(&t0, &t1, &t2, &t3);
  store4(x, xstride, t0, t1, t2, t3);
}

const od_dct_func_2d OD_FDCT_2D_SSE2[OD_NBSIZES + 1] = {
  od_bin_fdct4x4_sse2,
  od_bin_fdct8x8,
  od_bin_fdct16x16
};

const od_dct_func_2d OD_IDCT_2D_SSE2[OD_NBSIZES + 1] = {
  od_bin_idct4x4_sse2,
  od_bin_idct8x8,
  od_bin_idct16x16
};
