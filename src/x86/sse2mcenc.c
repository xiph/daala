/*Daala video codec
Copyright (c) 2015 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
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
#include "config.h"
#endif

#include <stdio.h>
#include "x86enc.h"
#include "x86int.h"

#if defined(OD_X86ASM)
#include <immintrin.h>

OD_SIMD_INLINE void od_mc_butterfly_2x2_16x4(__m64 *t0, __m64 *t1,
 __m64 *t2, __m64 *t3) {
  __m64 a;
  __m64 b;
  __m64 c;
  __m64 d;
  /*a = t0 + t1, c = (t0 + t1) - (t1 + t1) = t0 - t1
    b = t2 + t3, d = (t2 + t3) - (t3 + t3) = t2 - t3*/
  a = _mm_add_pi16(*t0, *t1);
  c = _mm_add_pi16(*t1, *t1);
  c = _mm_sub_pi16(a, c);
  b = _mm_add_pi16(*t2, *t3);
  d = _mm_add_pi16(*t3, *t3);
  d = _mm_sub_pi16(b, d);
  *t0 = a;
  *t1 = b;
  *t2 = c;
  *t3 = d;
}

/*Transpose 4 vectors with 4 16-bit vectors.*/
OD_SIMD_INLINE void od_transpose16x4(__m64 *t0, __m64 *t1,
 __m64 *t2, __m64 *t3) {
  __m64 a;
  __m64 b;
  __m64 c;
  __m64 d;
  a = _mm_unpacklo_pi16(*t0, *t1);
  b = _mm_unpacklo_pi16(*t2, *t3);
  c = _mm_unpackhi_pi16(*t0, *t1);
  d = _mm_unpackhi_pi16(*t2, *t3);
  *t0 = _mm_unpacklo_pi32(a, b);
  *t1 = _mm_unpackhi_pi32(a, b);
  *t2 = _mm_unpacklo_pi32(c, d);
  *t3 = _mm_unpackhi_pi32(c, d);
}

OD_SIMD_INLINE __m64 od_load_convert_subtract_x4(const unsigned char *src_p,
 const unsigned char *ref_p) {
  __m64 src_vec;
  __m64 ref_vec;
  src_vec = _mm_cvtsi32_si64(*((uint32_t *)src_p));
  ref_vec = _mm_cvtsi32_si64(*((uint32_t *)ref_p));
  src_vec = _mm_unpacklo_pi8(src_vec, ref_vec);
  ref_vec = _mm_unpacklo_pi8(ref_vec, ref_vec);
  return _mm_sub_pi16(src_vec, ref_vec);
}

int32_t od_mc_compute_satd8_4x4_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t satd;
  __m64 sums;
  __m64 a;
  __m64 b;
  __m64 c;
  __m64 d;
  a = od_load_convert_subtract_x4(src + 0*systride, ref + 0*rystride);
  b = od_load_convert_subtract_x4(src + 1*systride, ref + 1*rystride);
  c = od_load_convert_subtract_x4(src + 2*systride, ref + 2*rystride);
  d = od_load_convert_subtract_x4(src + 3*systride, ref + 3*rystride);
  /*Vertical 1D transform.*/
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  od_transpose16x4(&a, &b, &c, &d);
  /*Horizontal 1D transform.*/
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  /*Use the fact that (abs(a+b)+abs(a-b))/2=max(abs(a),abs(b)) to merge the
     final butterfly stage with the calculating the absolute values and the
     first stage of accumulation.
    Calculates (abs(a+b)+abs(a-b))/2-0x7FFF.
    An offset must be added to the final sum before rounding to account for
     subtracting 0x7FFF.*/
  a = _mm_sub_pi16(_mm_max_pi16(a, b), _mm_adds_pi16(_mm_add_pi16(a, b),
   _mm_set1_pi16(0x7FFF)));
  c = _mm_sub_pi16(_mm_max_pi16(c, d), _mm_adds_pi16(_mm_add_pi16(c, d),
   _mm_set1_pi16(0x7FFF)));
  /*Take the sum of all the absolute values.*/
  sums = _mm_add_pi16(a, c);
  /*Sum the elements of the vector.*/
  sums = _mm_add_pi16(sums, _mm_shuffle_pi16(sums, _MM_SHUFFLE(0, 1, 2, 3)));
  sums = _mm_add_pi16(sums, _mm_shuffle_pi16(sums, _MM_SHUFFLE(2, 3, 0, 1)));
  sums = _mm_unpacklo_pi16(sums, _mm_setzero_si64());
  satd = _mm_cvtsi64_si32(sums);
  /*Subtract the offset (8) and round.*/
  satd = (satd + 1 - 8) >> 1;
#if defined(OD_CHECKASM)
  {
    int32_t c_satd;
    c_satd = od_mc_compute_satd8_4x4_c(src, systride, ref, rystride);
    if (satd != c_satd) {
      fprintf(stderr, "od_mc_compute_satd %ix%i check failed: %i!=%i\n",
       4, 4, satd, c_satd);
    }
  }
#endif
  return satd;
}

OD_SIMD_INLINE void od_mc_butterfly_2x2_16x8(__m128i *t0, __m128i *t1,
 __m128i *t2, __m128i *t3) {
  __m128i a;
  __m128i b;
  __m128i c;
  __m128i d;
  /*a = t0 + t1, c = (t0 + t1) - (t1 + t1) = t0 - t1
    b = t2 + t3, d = (t2 + t3) - (t3 + t3) = t2 - t3*/
  a = _mm_add_epi16(*t0, *t1);
  c = _mm_add_epi16(*t1, *t1);
  c = _mm_sub_epi16(a, c);
  b = _mm_add_epi16(*t2, *t3);
  d = _mm_add_epi16(*t3, *t3);
  d = _mm_sub_epi16(b, d);
  *t0 = a;
  *t1 = b;
  *t2 = c;
  *t3 = d;
}

/*Transpose 8 vectors with 8 16-bit values.*/
OD_SIMD_INLINE void od_transpose16x8(__m128i *t0, __m128i *t1,
 __m128i *t2, __m128i *t3,  __m128i *t4, __m128i *t5,
 __m128i *t6, __m128i *t7) {
  __m128i a0;
  __m128i b0;
  __m128i c0;
  __m128i d0;
  __m128i e0;
  __m128i f0;
  __m128i g0;
  __m128i h0;
  __m128i a1;
  __m128i b1;
  __m128i c1;
  __m128i d1;
  __m128i e1;
  __m128i f1;
  __m128i g1;
  __m128i h1;
  /*00112233*/
  a0 = _mm_unpacklo_epi16(*t0, *t1);
  b0 = _mm_unpacklo_epi16(*t2, *t3);
  c0 = _mm_unpacklo_epi16(*t4, *t5);
  d0 = _mm_unpacklo_epi16(*t6, *t7);
  /*44556677*/
  e0 = _mm_unpackhi_epi16(*t0, *t1);
  f0 = _mm_unpackhi_epi16(*t2, *t3);
  g0 = _mm_unpackhi_epi16(*t4, *t5);
  h0 = _mm_unpackhi_epi16(*t6, *t7);
  /*00001111*/
  a1 = _mm_unpacklo_epi32(a0, b0);
  b1 = _mm_unpacklo_epi32(c0, d0);
  /*22223333*/
  c1 = _mm_unpackhi_epi32(a0, b0);
  d1 = _mm_unpackhi_epi32(c0, d0);
  /*44445555*/
  e1 = _mm_unpacklo_epi32(e0, f0);
  f1 = _mm_unpacklo_epi32(g0, h0);
  /*66667777*/
  g1 = _mm_unpackhi_epi32(e0, f0);
  h1 = _mm_unpackhi_epi32(g0, h0);
  *t0 = _mm_unpacklo_epi64(a1, b1);
  *t1 = _mm_unpackhi_epi64(a1, b1);
  *t2 = _mm_unpacklo_epi64(c1, d1);
  *t3 = _mm_unpackhi_epi64(c1, d1);
  *t4 = _mm_unpacklo_epi64(e1, f1);
  *t5 = _mm_unpackhi_epi64(e1, f1);
  *t6 = _mm_unpacklo_epi64(g1, h1);
  *t7 = _mm_unpackhi_epi64(g1, h1);
}

OD_SIMD_INLINE __m128i od_load_convert_subtract_x8(const unsigned char *src_p,
 const unsigned char *ref_p) {
  __m128i src_vec;
  __m128i ref_vec;
  src_vec = _mm_loadl_epi64((__m128i *)src_p);
  ref_vec = _mm_loadl_epi64((__m128i *)ref_p);
  src_vec = _mm_unpacklo_epi8(src_vec, ref_vec);
  ref_vec = _mm_unpacklo_epi8(ref_vec, ref_vec);
  return _mm_sub_epi16(src_vec, ref_vec);
}

OD_SIMD_INLINE int32_t od_mc_compute_satd8_8x8_part(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t satd;
  __m128i sums;
  __m128i a;
  __m128i b;
  __m128i c;
  __m128i d;
  __m128i e;
  __m128i f;
  __m128i g;
  __m128i h;
  a = od_load_convert_subtract_x8(src + 0*systride, ref + 0*rystride);
  b = od_load_convert_subtract_x8(src + 1*systride, ref + 1*rystride);
  c = od_load_convert_subtract_x8(src + 2*systride, ref + 2*rystride);
  d = od_load_convert_subtract_x8(src + 3*systride, ref + 3*rystride);
  e = od_load_convert_subtract_x8(src + 4*systride, ref + 4*rystride);
  f = od_load_convert_subtract_x8(src + 5*systride, ref + 5*rystride);
  g = od_load_convert_subtract_x8(src + 6*systride, ref + 6*rystride);
  h = od_load_convert_subtract_x8(src + 7*systride, ref + 7*rystride);
  /*Vertical 1D transform.*/
  od_mc_butterfly_2x2_16x8(&a, &b, &c, &d);
  od_mc_butterfly_2x2_16x8(&e, &f, &g, &h);
  od_mc_butterfly_2x2_16x8(&a, &b, &e, &f);
  od_mc_butterfly_2x2_16x8(&c, &d, &g, &h);
  od_mc_butterfly_2x2_16x8(&a, &b, &e, &f);
  od_mc_butterfly_2x2_16x8(&c, &d, &g, &h);
  od_transpose16x8(&a, &c, &b, &d, &e, &g, &f, &h);
  /*Horizontal 1D transform.*/
  od_mc_butterfly_2x2_16x8(&a, &b, &c, &d);
  od_mc_butterfly_2x2_16x8(&e, &f, &g, &h);
  od_mc_butterfly_2x2_16x8(&a, &b, &e, &f);
  od_mc_butterfly_2x2_16x8(&c, &d, &g, &h);
  /*Use the fact that (abs(a+b)+abs(a-b))/2=max(abs(a),abs(b)) to merge the
     final butterfly stage with the calculating the absolute values and the
     first stage of accumulation.
    Calculates (abs(a+b)+abs(a-b))/2-0x7FFF.
    An offset must be added to the final sum before rounding to account for
     subtracting 0x7FFF.*/
  a = _mm_sub_epi16(_mm_max_epi16(a, b), _mm_adds_epi16(_mm_add_epi16(a, b),
   _mm_set1_epi16(0x7FFF)));
  e = _mm_sub_epi16(_mm_max_epi16(e, f), _mm_adds_epi16(_mm_add_epi16(e, f),
   _mm_set1_epi16(0x7FFF)));
  c = _mm_sub_epi16(_mm_max_epi16(c, d), _mm_adds_epi16(_mm_add_epi16(c, d),
   _mm_set1_epi16(0x7FFF)));
  g = _mm_sub_epi16(_mm_max_epi16(g, h), _mm_adds_epi16(_mm_add_epi16(g, h),
   _mm_set1_epi16(0x7FFF)));
  a = _mm_add_epi16(a, e);
  c = _mm_add_epi16(c, g);
  /*Convert to 32-bit unsigned integers and sum horizontally using madd to
    avoid overflowing 16-bit unsigned integers.
    The naively calculated max values of a and c are
     ((8 rows * 8 cols * 256) * 2 / 2 + 1 (offset)) * 2 or 0x8002.
    The actual max is lower so it is safe to use _mm_madd_epi16.*/
  a = _mm_madd_epi16(a, _mm_set1_epi16(1));
  c = _mm_madd_epi16(c, _mm_set1_epi16(1));
  sums = _mm_add_epi32(a, c);
  /*Sum the elements of the vector.*/
  sums = _mm_add_epi32(sums, _mm_shuffle_epi32(sums, _MM_SHUFFLE(0, 1, 2, 3)));
  sums = _mm_add_epi32(sums, _mm_shuffle_epi32(sums, _MM_SHUFFLE(2, 3, 0, 1)));
  satd = _mm_cvtsi128_si32(sums);
  /*Subtract the offset (32) and round.*/
  satd = (satd + 2 - 32) >> 2;
#if defined(OD_CHECKASM)
  {
    int32_t c_satd;
    c_satd = od_mc_compute_satd8_8x8_c(src, systride, ref, rystride);
    if (satd != c_satd) {
      fprintf(stderr, "od_mc_compute_satd %ix%i check failed: %i!=%i\n",
      8, 8, satd, c_satd);
    }
  }
#endif
  return satd;
}

/*Perform SATD on 8x8 blocks within src and ref then sum the results of
   each one.*/
OD_SIMD_INLINE int32_t od_mc_compute_sum_8x8_satd8(int ln,
 const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int n;
  int i;
  int j;
  int32_t satd;
  n = 1 << ln;
  OD_ASSERT(n >= 8);
  satd = 0;
  for (i = 0; i < n; i += 8) {
    for (j = 0; j < n; j += 8) {
      satd += od_mc_compute_satd8_8x8_part(src + i*systride + j, systride,
       ref + i*rystride + j, rystride);
    }
  }
  return satd;
}

int32_t od_mc_compute_satd8_8x8_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  return od_mc_compute_sum_8x8_satd8(3, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_16x16_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  return od_mc_compute_sum_8x8_satd8(4, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_32x32_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  return od_mc_compute_sum_8x8_satd8(5, src, systride, ref, rystride);
}

int32_t od_mc_compute_satd8_64x64_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  return od_mc_compute_sum_8x8_satd8(6, src, systride, ref, rystride);
}

#endif
