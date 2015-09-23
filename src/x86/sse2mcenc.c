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

OD_SIMD_INLINE void od_mc_butterfly_2x2_16x8(__m128i *t0, __m128i *t1,
 __m128i *t2, __m128i *t3) {
  __m128i in0 = *t0;
  __m128i in1 = *t1;
  __m128i in2 = *t2;
  __m128i in3 = *t3;
  *t0 = _mm_add_epi16(in0, in1);
  *t1 = _mm_add_epi16(in2, in3);
  *t2 = _mm_sub_epi16(in0, in1);
  *t3 = _mm_sub_epi16(in2, in3);
}

OD_SIMD_INLINE void od_mc_butterfly_2x2_32x4(__m128i *t0, __m128i *t1,
 __m128i *t2, __m128i *t3) {
  __m128i in0 = *t0;
  __m128i in1 = *t1;
  __m128i in2 = *t2;
  __m128i in3 = *t3;
  *t0 = _mm_add_epi32(in0, in1);
  *t1 = _mm_add_epi32(in2, in3);
  *t2 = _mm_sub_epi32(in0, in1);
  *t3 = _mm_sub_epi32(in2, in3);
}

OD_SIMD_INLINE void od_mc_butterfly_2x2_16x4(__m64 *t0, __m64 *t1,
 __m64 *t2, __m64 *t3) {
  __m64 in0 = *t0;
  __m64 in1 = *t1;
  __m64 in2 = *t2;
  __m64 in3 = *t3;
  *t0 = _mm_add_pi16(in0, in1);
  *t1 = _mm_add_pi16(in2, in3);
  *t2 = _mm_sub_pi16(in0, in1);
  *t3 = _mm_sub_pi16(in2, in3);
}

/*Convert 4 unsigned 8 bit values to 4 unsigned 16 bit.
  Both input and output are 64 bit (mmx) vectors.*/
OD_SIMD_INLINE __m64 od_cvtu8x4_u16x4(__m64 in) {
  return _mm_unpacklo_pi8(in, _mm_setzero_si64());
}

/*Convert 8 unsigned 8 bit values to 8 unsigned 16 bit.
  Both input and output are 128 bit (sse) integer vectors.*/
OD_SIMD_INLINE __m128i od_cvtu8x8_u16x8(__m128i in) {
  return _mm_unpacklo_epi8(in, _mm_setzero_si128());
}

/*Convert 16 unsigned integers to */
OD_SIMD_INLINE void od_cvtu8x16_2xu16x8(__m128i *outlo, __m128i *outhi,
 __m128i in)
{
  *outlo = _mm_unpacklo_epi8(in, _mm_setzero_si128());
  *outhi = _mm_unpackhi_epi8(in, _mm_setzero_si128());
}

/*Corresponds to _mm_cvtepi16_epi32 (sse4.1).*/
OD_SIMD_INLINE __m128i od_cvti16x4_i32x4(__m128i in) {
  /*_mm_cmplt_epi16 returns */
  __m128i extend = _mm_cmplt_epi16(in, _mm_setzero_si128());
  return _mm_unpacklo_epi16(in, extend);
}

/*Corresponds to _mm_abs_pi16 (ssse3).*/
OD_SIMD_INLINE __m64 od_abs16x4(__m64 in) {
  __m64 mask;
  mask = _mm_srai_pi16(in, 15);
  return _mm_xor_si64(_mm_add_pi16(in, mask), mask);
}

/*Corresponds to _mm_abs_epi32 (ssse3).*/
OD_SIMD_INLINE __m128i od_abs32x4(__m128i in) {
  __m128i mask;
  mask = _mm_srai_epi32(in, 31);
  return _mm_xor_si128(_mm_add_epi32(in, mask), mask);
}

/*Transpose 8 vectors with 8 16 bit vectors each.*/
OD_SIMD_INLINE void od_transpose16x8(__m128i *t0, __m128i *t1,
 __m128i *t2, __m128i *t3,  __m128i *t4, __m128i *t5, __m128i *t6, __m128i *t7) {
  __m128i a1;
  __m128i b1;
  __m128i c1;
  __m128i d1;
  __m128i e1;
  __m128i f1;
  __m128i g1;
  __m128i h1;
  /*00112233*/
  __m128i a0 = _mm_unpacklo_epi16(*t0, *t1);
  __m128i b0 = _mm_unpacklo_epi16(*t2, *t3);
  __m128i c0 = _mm_unpacklo_epi16(*t4, *t5);
  __m128i d0 = _mm_unpacklo_epi16(*t6, *t7);
  /*44556677*/
  __m128i e0 = _mm_unpackhi_epi16(*t0, *t1);
  __m128i f0 = _mm_unpackhi_epi16(*t2, *t3);
  __m128i g0 = _mm_unpackhi_epi16(*t4, *t5);
  __m128i h0 = _mm_unpackhi_epi16(*t6, *t7);
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

/*Transpose 4 vectors with 4 16 bit vectors each.*/
OD_SIMD_INLINE void od_transpose16x4(__m64 *t0, __m64 *t1,
 __m64 *t2, __m64 *t3) {
  __m64 a = _mm_unpacklo_pi16(*t0, *t1);
  __m64 b = _mm_unpacklo_pi16(*t2, *t3);
  __m64 c = _mm_unpackhi_pi16(*t0, *t1);
  __m64 d = _mm_unpackhi_pi16(*t2, *t3);
  *t0 = _mm_unpacklo_pi32(a, b);
  *t1 = _mm_unpackhi_pi32(a, b);
  *t2 = _mm_unpacklo_pi32(c, d);
  *t3 = _mm_unpackhi_pi32(c, d);
}

/*Get the sum of length elements in src.*/
OD_SIMD_INLINE int od_sum_array_32i(const int32_t *src, const int length)
{
  int sum;
  __m128i sums;
  const int32_t *src_p;
  int i;
  /*Length should be a power of 2 and greater than or equal to 4.*/
  OD_ASSERT(!((length - 1) & length));
  OD_ASSERT(a >= 4);
  sum = 0;
  sums = _mm_setzero_si128();
  src_p = src;
  /*Vertically sum the elements of the array into one vector.*/
  for(i = 0; i < length; i+=4) {
    __m128i tmp = _mm_load_si128((__m128i *)src_p);
    tmp = od_abs32x4(tmp);
    sums = _mm_add_epi32(sums, tmp);
    src_p += 4;
  }
  /*Sum the elements of the vector.*/
  sums = _mm_add_epi32(sums, _mm_srli_si128(sums, 8));
  sums = _mm_add_epi32(sums, _mm_srli_si128(sums, 4));
  sum = _mm_cvtsi128_si32(sums);
  return sum;
}

int od_mc_compute_satd_4x4_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int satd;
  __m64 sums;
  __m64 diff[4];
  __m64 *diff_p;
  const unsigned char *src_p;
  const unsigned char *ref_p;
  __m64 a;
  __m64 b;
  __m64 c;
  __m64 d;
  int i;
  diff_p = diff;
  src_p = src;
  ref_p = ref;
  for(i = 0; i < 4; i++) {
    __m64 src_vec;
    __m64 ref_vec;
    src_vec = _mm_cvtsi32_si64(*((uint32_t *)src_p));
    src_vec = od_cvtu8x4_u16x4(src_vec);
    ref_vec = _mm_cvtsi32_si64(*((uint32_t *)ref_p));
    ref_vec = od_cvtu8x4_u16x4(ref_vec);
    *diff_p = _mm_sub_pi16(src_vec, ref_vec);
    diff_p++;
    src_p += systride;
    ref_p += rystride;
  }
  a = diff[0];
  b = diff[1];
  c = diff[2];
  d = diff[3];
  /*Horizontal 1D transform.*/
  od_transpose16x4(&a, &b, &c, &d);
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  od_transpose16x4(&a, &b, &c, &d);
  /*Vertical 1D transform.*/
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  od_mc_butterfly_2x2_16x4(&a, &b, &c, &d);
  /*Take the sum of all the absolute values.*/
  a = od_abs16x4(a);
  b = od_abs16x4(b);
  c = od_abs16x4(c);
  d = od_abs16x4(d);
  a = _mm_add_pi16(a, b);
  c = _mm_add_pi16(c, d);
  sums = _mm_add_pi16(a, c);
  /*Sum the elements of the vector.*/
  sums = _mm_add_pi16(sums, _mm_srli_si64(sums, 4*8));
  sums = _mm_unpacklo_pi16(sums, _mm_setzero_si64());
  sums = _mm_add_pi32(sums, _mm_srli_si64(sums, 4*8));
  satd = _mm_cvtsi64_si32(sums);
  satd = (satd + 2) >> 2;
#if defined(OD_CHECKASM)
  {
    int c_satd = od_mc_compute_satd_4x4_c(src, systride, ref, rystride);
    if (satd != c_satd) {
      fprintf(stderr, "od_mc_compute_satd %ix%i check failed: %i!=%i\n",
       4, 4, satd, c_satd);
    }
  }
#endif
  return satd;
}

static void od_mc_compute_satd_8x8_hor_sse2(int16_t *dest, int dystride,
 const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  __m128i *diff_p;
  int16_t diff[8*8];
  int16_t *dest_p;
  const unsigned char *src_p;
  const unsigned char *ref_p;
  int i;
  diff_p = (__m128i *)diff;
  src_p = src;
  ref_p = ref;
  for (i = 0; i < 8; i++) {
    __m128i src_vec;
    __m128i ref_vec;
    src_vec = _mm_loadl_epi64((__m128i *)src_p);
    src_vec = od_cvtu8x8_u16x8(src_vec);
    ref_vec = _mm_loadl_epi64((__m128i *)ref_p);
    ref_vec = od_cvtu8x8_u16x8(ref_vec);
    *diff_p = _mm_sub_epi16(src_vec, ref_vec);
    diff_p++;
    src_p += systride;
    ref_p += rystride;
  }
  /*Horizontal 1D transform.*/
  dest_p = dest;
  diff_p = (__m128i *)diff;
  for (i = 0; i < 8; i += 8) {
    __m128i a;
    __m128i b;
    __m128i c;
    __m128i d;
    __m128i e;
    __m128i f;
    __m128i g;
    __m128i h;
    a = diff_p[0];
    b = diff_p[1];
    c = diff_p[2];
    d = diff_p[3];
    e = diff_p[4];
    f = diff_p[5];
    g = diff_p[6];
    h = diff_p[7];
    od_transpose16x8(&a, &b, &c, &d, &e, &f, &g, &h);
    od_mc_butterfly_2x2_16x8(&a, &b, &c, &d);
    od_mc_butterfly_2x2_16x8(&e, &f, &g, &h);
    od_mc_butterfly_2x2_16x8(&a, &b, &e, &f);
    od_mc_butterfly_2x2_16x8(&c, &d, &g, &h);
    od_mc_butterfly_2x2_16x8(&a, &b, &e, &f);
    od_mc_butterfly_2x2_16x8(&c, &d, &g, &h);
    od_transpose16x8(&a, &c, &b, &d, &e, &g, &f, &h);
    _mm_store_si128((__m128i *)(dest_p + 0*dystride), a);
    _mm_store_si128((__m128i *)(dest_p + 1*dystride), c);
    _mm_store_si128((__m128i *)(dest_p + 2*dystride), b);
    _mm_store_si128((__m128i *)(dest_p + 3*dystride), d);
    _mm_store_si128((__m128i *)(dest_p + 4*dystride), e);
    _mm_store_si128((__m128i *)(dest_p + 5*dystride), g);
    _mm_store_si128((__m128i *)(dest_p + 6*dystride), f);
    _mm_store_si128((__m128i *)(dest_p + 7*dystride), h);
    dest_p += 8*dystride;
    diff_p += 8;
  }
}

static void od_mc_compute_satd_8x8_ver_sse2(int32_t *dest, int dystride,
  int16_t *src, int systride) {
  int32_t *dest_p;
  int16_t *src_p;
  int i;
  src_p = src;
  dest_p = dest;
  /*Vertical 1D transform.*/
  for (i = 0; i < 8; i+=4) {
    __m128i a;
    __m128i b;
    __m128i c;
    __m128i d;
    __m128i e;
    __m128i f;
    __m128i g;
    __m128i h;
    a = _mm_loadl_epi64((__m128i *)(src_p + 0*systride));
    a = od_cvti16x4_i32x4(a);
    b = _mm_loadl_epi64((__m128i *)(src_p + 1*systride));
    b = od_cvti16x4_i32x4(b);
    c = _mm_loadl_epi64((__m128i *)(src_p + 2*systride));
    c = od_cvti16x4_i32x4(c);
    d = _mm_loadl_epi64((__m128i *)(src_p + 3*systride));
    d = od_cvti16x4_i32x4(d);
    e = _mm_loadl_epi64((__m128i *)(src_p + 4*systride));
    e = od_cvti16x4_i32x4(e);
    f = _mm_loadl_epi64((__m128i *)(src_p + 5*systride));
    f = od_cvti16x4_i32x4(f);
    g = _mm_loadl_epi64((__m128i *)(src_p + 6*systride));
    g = od_cvti16x4_i32x4(g);
    h = _mm_loadl_epi64((__m128i *)(src_p + 7*systride));
    h = od_cvti16x4_i32x4(h);
    od_mc_butterfly_2x2_32x4(&a, &b, &c, &d);
    od_mc_butterfly_2x2_32x4(&e, &f, &g, &h);
    od_mc_butterfly_2x2_32x4(&a, &b, &e, &f);
    od_mc_butterfly_2x2_32x4(&c, &d, &g, &h);
    od_mc_butterfly_2x2_32x4(&a, &b, &e, &f);
    od_mc_butterfly_2x2_32x4(&c, &d, &g, &h);
    _mm_store_si128((__m128i *)(dest_p + 0*dystride), a);
    _mm_store_si128((__m128i *)(dest_p + 1*dystride), c);
    _mm_store_si128((__m128i *)(dest_p + 2*dystride), b);
    _mm_store_si128((__m128i *)(dest_p + 3*dystride), d);
    _mm_store_si128((__m128i *)(dest_p + 4*dystride), e);
    _mm_store_si128((__m128i *)(dest_p + 5*dystride), f);
    _mm_store_si128((__m128i *)(dest_p + 6*dystride), g);
    _mm_store_si128((__m128i *)(dest_p + 7*dystride), h);
    src_p+=4;
    dest_p+=4;
  }
}

int od_mc_compute_satd_8x8_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int satd;
  int16_t buff1[8*8];
  int32_t buff2[8*8];
  od_mc_compute_satd_8x8_hor_sse2(buff1, 8, src, systride, ref, rystride);
  od_mc_compute_satd_8x8_ver_sse2(buff2, 8, buff1, 8);
  satd = od_sum_array_32i(buff2, 8*8);
  satd = (satd + 4) >> 3;
#if defined(OD_CHECKASM)
  {
    int c_satd = od_mc_compute_satd_8x8_c(src, systride, ref, rystride);
    if (satd != c_satd) {
      fprintf(stderr, "od_mc_compute_satd %ix%i check failed: %i!=%i\n",
      8, 8, satd, c_satd);
    }
  }
#endif
  return satd;
}

static void od_mc_compute_satd_16x16_hor_sse2(int16_t *dest, int dystride,
 const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int16_t buff[16*16];
  int blk_size;
  int i;
  int j;
  int row_ptr;
  int row_ptr2;
  int dest_row_ptr;
  int dest_row_ptr2;
  blk_size = 16;
  od_mc_compute_satd_8x8_hor_sse2(buff, blk_size, src, systride,
   ref, rystride);
  od_mc_compute_satd_8x8_hor_sse2(buff + 8, blk_size,
   src + 8, systride, ref + 8, rystride);
  od_mc_compute_satd_8x8_hor_sse2(buff + blk_size*8, blk_size,
   src + systride*8, systride, ref + rystride*8, rystride);
  od_mc_compute_satd_8x8_hor_sse2(buff + blk_size*8 + 8, blk_size,
   src + systride*8 + 8, systride, ref + rystride*8 + 8, rystride);
  for (j = 0; j < 8; j++) {
    row_ptr = j*blk_size;
    row_ptr2 = (j + 8)*blk_size;
    dest_row_ptr = j*dystride;
    dest_row_ptr2 = (j + 8)*dystride;
    for (i = 0; i < 8; i += 8) {
      __m128i a0;
      __m128i a1;
      __m128i b0;
      __m128i b1;
      __m128i c0;
      __m128i c1;
      __m128i d0;
      __m128i d1;
      a0 = _mm_load_si128((__m128i *)(buff + row_ptr + i));
      b0 = _mm_load_si128((__m128i *)(buff + row_ptr + i + 8));
      c0 = _mm_load_si128((__m128i *)(buff + row_ptr2 + i));
      d0 = _mm_load_si128((__m128i *)(buff + row_ptr2 + i + 8));
      a1 = _mm_add_epi16(a0, b0);
      b1 = _mm_sub_epi16(a0, b0);
      c1 = _mm_add_epi16(c0, d0);
      d1 = _mm_sub_epi16(c0, d0);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr + i), a1);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr + i + 8), b1);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr2 + i), c1);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr2 + i + 8), d1);
    }
  }
}

static void od_mc_compute_satd_16x16_ver_sse2(int32_t *dest, int dystride,
 int16_t *src, int systride) {
  int32_t buff[16*16];
  int blk_size;
  int i;
  int j;
  int row_ptr;
  int row_ptr2;
  int dest_row_ptr;
  int dest_row_ptr2;
  blk_size = 16;
  od_mc_compute_satd_8x8_ver_sse2(buff, blk_size, src, systride);
  od_mc_compute_satd_8x8_ver_sse2(buff + 8, blk_size, src + 8, systride);
  od_mc_compute_satd_8x8_ver_sse2(buff + blk_size*8, blk_size,
   src + systride*8, systride);
  od_mc_compute_satd_8x8_ver_sse2(buff + blk_size*8 + 8, blk_size,
   src + systride*8 + 8, systride);
  for (j = 0; j < 8; j++) {
    row_ptr = j*blk_size;
    row_ptr2 = (j + 8)*blk_size;
    dest_row_ptr = j*dystride;
    dest_row_ptr2 = (j + 8)*dystride;
    for (i = 0; i < 8; i += 4) {
      __m128i a0;
      __m128i a1;
      __m128i b0;
      __m128i b1;
      __m128i c0;
      __m128i c1;
      __m128i d0;
      __m128i d1;
      a0 = _mm_load_si128((__m128i *)(buff + row_ptr + i));
      b0 = _mm_load_si128((__m128i *)(buff + row_ptr + i + 8));
      c0 = _mm_load_si128((__m128i *)(buff + row_ptr2 + i));
      d0 = _mm_load_si128((__m128i *)(buff + row_ptr2 + i + 8));
      a1 = _mm_add_epi32(a0, c0);
      b1 = _mm_sub_epi32(a0, c0);
      c1 = _mm_add_epi32(b0, d0);
      d1 = _mm_sub_epi32(b0, d0);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr + i), a1);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr2 + i), b1);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr + i + 8), c1);
      _mm_store_si128((__m128i *)(dest + dest_row_ptr2 + i + 8), d1);
    }
  }
}

int od_mc_compute_satd_16x16_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int satd;
  int16_t buff1[16*16];
  int32_t buff2[16*16];
  od_mc_compute_satd_16x16_hor_sse2(buff1, 16, src, systride, ref, rystride);
  od_mc_compute_satd_16x16_ver_sse2(buff2, 16, buff1, 16);
  satd = od_sum_array_32i(buff2, 16*16);
  satd = (satd + 8) >> 4;
#if defined(OD_CHECKASM)
  {
    int c_satd = od_mc_compute_satd_16x16_c(src, systride, ref, rystride);
    if (satd != c_satd) {
      fprintf(stderr, "od_mc_compute_satd %ix%i check failed: %i!=%i\n",
      16, 16, satd, c_satd);
    }
  }
#endif
  return satd;
}

int od_mc_compute_satd_32x32_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int rystride) {
  int32_t satd;
  int16_t buff[32*32];
  int16_t buff2[32*32];
  int32_t buff3[32*32];
  int32_t buff4[32*32];
  int blk_size;
  int i;
  int j;
  int row_ptr;
  int row_ptr2;
  blk_size = 32;
  /*Horizontal 1D transform.*/
  od_mc_compute_satd_16x16_hor_sse2(buff, blk_size, src, systride,
   ref, rystride);
  od_mc_compute_satd_16x16_hor_sse2(buff + 16, blk_size,
   src + 16, systride, ref + 16, rystride);
  od_mc_compute_satd_16x16_hor_sse2(buff + blk_size*16, blk_size,
   src + systride*16, systride, ref + rystride*16, rystride);
  od_mc_compute_satd_16x16_hor_sse2(buff + blk_size*16 + 16, blk_size,
   src + systride*16 + 16, systride, ref + rystride*16 + 16, rystride);
  for (j = 0; j < 16; j++) {
    row_ptr = j*blk_size;
    row_ptr2 = (j + 16)*blk_size;
    for (i = 0; i < 16; i+=8) {
      __m128i a0;
      __m128i a1;
      __m128i b0;
      __m128i b1;
      __m128i c0;
      __m128i c1;
      __m128i d0;
      __m128i d1;
      a0 = _mm_load_si128((__m128i *)(buff + row_ptr + i));
      b0 = _mm_load_si128((__m128i *)(buff + row_ptr + i + 16));
      c0 = _mm_load_si128((__m128i *)(buff + row_ptr2 + i));
      d0 = _mm_load_si128((__m128i *)(buff + row_ptr2 + i + 16));
      a1 = _mm_add_epi16(a0, b0);
      b1 = _mm_sub_epi16(a0, b0);
      c1 = _mm_add_epi16(c0, d0);
      d1 = _mm_sub_epi16(c0, d0);
      _mm_store_si128((__m128i *)(buff2 + row_ptr + i), a1);
      _mm_store_si128((__m128i *)(buff2 + row_ptr + i + 16), b1);
      _mm_store_si128((__m128i *)(buff2 + row_ptr2 + i), c1);
      _mm_store_si128((__m128i *)(buff2 + row_ptr2 + i + 16), d1);
    }
  }
  /*Vertical 1D transform.*/
  od_mc_compute_satd_16x16_ver_sse2(buff3, blk_size, buff2, blk_size);
  od_mc_compute_satd_16x16_ver_sse2(buff3 + 16, blk_size,
   buff2 + 16, blk_size);
  od_mc_compute_satd_16x16_ver_sse2(buff3 + blk_size*16, blk_size,
   buff2 + blk_size*16, blk_size);
  od_mc_compute_satd_16x16_ver_sse2(buff3 + blk_size*16 + 16, blk_size,
   buff2 + blk_size*16 + 16, blk_size);
  for (j = 0; j < 16; j++) {
    row_ptr = j*blk_size;
    row_ptr2 = (j + 16)*blk_size;
    for (i = 0; i < 16; i+=4) {
      __m128i a0;
      __m128i a1;
      __m128i b0;
      __m128i b1;
      __m128i c0;
      __m128i c1;
      __m128i d0;
      __m128i d1;
      a0 = _mm_load_si128((__m128i *)(buff3 + row_ptr + i));
      b0 = _mm_load_si128((__m128i *)(buff3 + row_ptr + i + 16));
      c0 = _mm_load_si128((__m128i *)(buff3 + row_ptr2 + i));
      d0 = _mm_load_si128((__m128i *)(buff3 + row_ptr2 + i + 16));
      a1 = _mm_add_epi32(a0, c0);
      b1 = _mm_sub_epi32(a0, c0);
      c1 = _mm_add_epi32(b0, d0);
      d1 = _mm_sub_epi32(b0, d0);
      _mm_store_si128((__m128i *)(buff4 + row_ptr + i), a1);
      _mm_store_si128((__m128i *)(buff4 + row_ptr2 + i), b1);
      _mm_store_si128((__m128i *)(buff4 + row_ptr + i + 16), c1);
      _mm_store_si128((__m128i *)(buff4 + row_ptr2 + i + 16), d1);
    }
  }
  satd = od_sum_array_32i(buff4, 32*32);
  satd = (satd + 16) >> 5;
#if defined(OD_CHECKASM)
  {
    int c_satd = od_mc_compute_satd_32x32_c(src, systride, ref, rystride);
    if (satd != c_satd) {
      fprintf(stderr, "od_mc_compute_satd %ix%i check failed: %i!=%i\n",
      32, 32, satd, c_satd);
    }
  }
#endif
  return satd;
}
#endif
