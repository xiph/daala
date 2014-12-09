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

#include <immintrin.h>
#include "../dct.h"
#include "x86int.h"

typedef __m256i od_m256i;

OD_SIMD_INLINE od_m256i od_mm256_srai_epi32(od_m256i a, int c) {
  return _mm256_srai_epi32(a, c);
}

OD_SIMD_INLINE od_m256i od_mm256_srli_epi32(od_m256i a, int c) {
  return _mm256_srli_epi32(a, c);
}

OD_SIMD_INLINE od_m256i od_mm256_add_epi32(od_m256i a, od_m256i b) {
  return _mm256_add_epi32(a, b);
}

OD_SIMD_INLINE od_m256i od_mm256_sub_epi32(od_m256i a, od_m256i b) {
  return _mm256_sub_epi32(a, b);
}

OD_SIMD_INLINE od_m256i od_mm256_set1_epi32(int c) {
  return _mm256_set1_epi32(c);
}

OD_SIMD_INLINE od_m256i od_mm256_unpacklo_epi32(od_m256i a, od_m256i b) {
  return _mm256_unpacklo_epi32(a, b);
}

OD_SIMD_INLINE od_m256i od_mm256_unpackhi_epi32(od_m256i a, od_m256i b) {
  return _mm256_unpackhi_epi32(a, b);
}

OD_SIMD_INLINE od_m256i od_mm256_unpacklo_epi64(od_m256i a, od_m256i b) {
  return _mm256_unpacklo_epi64(a, b);
}

OD_SIMD_INLINE od_m256i od_mm256_unpackhi_epi64(od_m256i a, od_m256i b) {
  return _mm256_unpackhi_epi64(a, b);
}

OD_SIMD_INLINE od_m256i od_mm256_permute2x128_si256(od_m256i a,
 od_m256i b, const int c) {
  return _mm256_permute2x128_si256(a, b, c);
}

OD_SIMD_INLINE od_m256i od_mm256_loadu_si256(const od_m256i *ptr) {
  return _mm256_loadu_si256(ptr);
}

OD_SIMD_INLINE void od_mm256_storeu_si256(od_m256i *ptr, od_m256i a) {
  _mm256_storeu_si256(ptr, a);
}

OD_SIMD_INLINE od_m256i mul_epi32_256(od_m256i a, int b1) {
  return _mm256_mullo_epi32(a, _mm256_set1_epi32(b1));
}

#define od_bin_fdct8x8_x86 od_bin_fdct8x8_avx2
#define od_bin_idct8x8_x86 od_bin_idct8x8_avx2

#include "x86dct.h"
