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

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "x86int.h"
#include "cpu.h"

#include <stdio.h>

#if defined(OD_X86ASM)
#include <emmintrin.h>
#include "../filter.h"

#if defined(OD_CHECKASM)
void od_filter_dering_direction_check(int16_t *y, int ystride, int16_t *in,
 int ln, int threshold, int dir) {
  int16_t dst[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int i;
  int j;
  int failed;
  od_filter_dering_direction_c(dst, OD_BSIZE_MAX, in, ln,
   threshold, dir);
  failed = 0;
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      if (dst[i*OD_BSIZE_MAX + j] != y[i*ystride + j]) {
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         dst[i*OD_BSIZE_MAX + j],y[i*ystride + j],i,j);
        failed = 1;
      }
    }
  }
  if (failed) {
    fprintf(stderr, "od_filter_dering_direction check failed.\n");
  }
}
#endif

typedef void (*od_filter_dering_direction_fixed_func)(int16_t *y, int ystride,
 int16_t *in, int threshold, int dir);

/*Corresponds to _mm_abs_epi16 (ssse3).*/
OD_SIMD_INLINE __m128i _mm_abs_epi16(__m128i in) {
  __m128i mask;
  mask = _mm_cmpgt_epi16(_mm_setzero_si128(), in);
  return _mm_sub_epi16(_mm_xor_si128(in, mask), mask);
}

void od_filter_dering_direction_4x4(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  int i;
  int k;
  static const int taps[3] = {3, 2, 2};
  __m128i sum;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i res;
  for (i = 0; i < 8; i++) {
    sum = _mm_set1_epi16(0);
    row = _mm_loadl_epi64((__m128i*)&in[i*OD_FILT_BSTRIDE]);
    for (k = 0; k < 3; k++) {
      /* p = in[i*OD_FILT_BSTRIDE + offset] - row */;
      p = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)&in[i*OD_FILT_BSTRIDE +
       direction_offsets_table[dir][k]]), row);
      /* if (abs(p) < threshold) sum += taps[k]*p; */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), _mm_set1_epi16(threshold));
      p = _mm_mullo_epi16(p, _mm_set1_epi16(taps[k]));
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);

      /* p = in[i*OD_FILT_BSTRIDE + offset] - row */;
      p = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)&in[i*OD_FILT_BSTRIDE -
       direction_offsets_table[dir][k]]), row);
      /* if (abs(p) < threshold) sum += taps[k]*p1; */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), _mm_set1_epi16(threshold));
      p = _mm_mullo_epi16(p, _mm_set1_epi16(taps[k]));
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);
    }
    /* res = row + ((sum + 8) >> 4) */
    res = _mm_add_epi16(sum, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    res = _mm_add_epi16(row, res);
    _mm_storel_epi64((__m128i*)&y[i*ystride], res);
  }
}

void od_filter_dering_direction_8x8(int16_t *y, int ystride, int16_t *in,
 int threshold, int dir) {
  int i;
  int k;
  static const int taps[3] = {3, 2, 2};
  __m128i sum;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i res;
  for (i = 0; i < 8; i++) {
    sum = _mm_set1_epi16(0);
    row = _mm_loadu_si128((__m128i*)&in[i*OD_FILT_BSTRIDE]);
    for (k = 0; k < 3; k++) {
      /* p = in[i*OD_FILT_BSTRIDE + offset] - row */;
      p = _mm_sub_epi16(_mm_loadu_si128((__m128i*)&in[i*OD_FILT_BSTRIDE +
       direction_offsets_table[dir][k]]), row);
      /* if (abs(p) < threshold) sum += taps[k]*p; */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), _mm_set1_epi16(threshold));
      p = _mm_mullo_epi16(p, _mm_set1_epi16(taps[k]));
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);

      /* p = in[i*OD_FILT_BSTRIDE + offset] - row */;
      p = _mm_sub_epi16(_mm_loadu_si128((__m128i*)&in[i*OD_FILT_BSTRIDE -
       direction_offsets_table[dir][k]]), row);
      /* if (abs(p) < threshold) sum += taps[k]*p1; */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), _mm_set1_epi16(threshold));
      p = _mm_mullo_epi16(p, _mm_set1_epi16(taps[k]));
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);
    }
    /* res = row + ((sum + 8) >> 4) */
    res = _mm_add_epi16(sum, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    res = _mm_add_epi16(row, res);
    _mm_storeu_si128((__m128i*)&y[i*ystride], res);
  }
}

void od_filter_dering_direction_sse2(int16_t *y, int ystride, int16_t *in,
 int ln, int threshold, int dir) {
  static const od_filter_dering_direction_fixed_func VTBL[4] = {
    NULL,
    NULL,
    od_filter_dering_direction_4x4,
    od_filter_dering_direction_8x8,
  };
  OD_ASSERT(ln >= 2 && ln < 4);
  VTBL[ln](y, ystride, in, threshold, dir);
#if defined(OD_CHECKASM)
  od_filter_dering_direction_check(y, ystride, in, ln, threshold,
   dir);
#endif
}

#if defined(OD_CHECKASM)
void od_filter_dering_orthogonal_check(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int ln, int threshold, int dir) {
  int16_t dst[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int i;
  int j;
  int failed;
  od_filter_dering_orthogonal_c(dst, OD_BSIZE_MAX, in, x, xstride, ln,
   threshold, dir);
  failed = 0;
  for (i = 0; i < 1 << ln; i++) {
    for (j = 0; j < 1 << ln; j++) {
      if (dst[i*OD_BSIZE_MAX + j] != y[i*ystride + j]) {
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         dst[i*OD_BSIZE_MAX + j],y[i*ystride + j],i,j);
        failed = 1;
      }
    }
  }
  if (failed) {
    fprintf(stderr, "od_filter_dering_orthogonal_check check failed.\n");
  }
}
#endif

void od_filter_dering_orthogonal_4x4(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int threshold, int dir) {
  int i;
  int k;
  int offset;
  __m128i res;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i sum;
  __m128i athresh;
  if (dir <= 4) offset = OD_FILT_BSTRIDE;
  else offset = 1;
  for (i = 0; i < 8; i++) {
    sum = _mm_set1_epi16(0);
    row = _mm_loadl_epi64((__m128i*)&in[i*OD_FILT_BSTRIDE]);
    /* athresh = OD_MINI(threshold, threshold/3
       + abs(in[i*OD_FILT_BSTRIDE] - x[i*xstride])) */
    athresh = _mm_min_epi16(_mm_set1_epi16(threshold),
     _mm_add_epi16(_mm_set1_epi16(threshold/3),
     _mm_abs_epi16(_mm_sub_epi16(row,
     _mm_loadl_epi64((__m128i*)&x[i*xstride])))));

    for (k = 1; k <= 2; k++) {
      /* p = in[i*OD_FILT_BSTRIDE + k*offset] - row */
      p = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)&in[i*OD_FILT_BSTRIDE +
       k*offset]), row);
      /* if (abs(p) < athresh) sum += p */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), athresh);
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);

      /* p = in[i*OD_FILT_BSTRIDE - k*offset] - row */
      p = _mm_sub_epi16(_mm_loadl_epi64((__m128i*)&in[i*OD_FILT_BSTRIDE -
       k*offset]), row);
      /* if (abs(p) < athresh) sum += p */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), athresh);
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);
    }

    /* row + ((3*sum + 8) >> 4) */
    res = _mm_mullo_epi16(sum, _mm_set1_epi16(3));
    res = _mm_add_epi16(res, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    res = _mm_add_epi16(res, row);
    _mm_storel_epi64((__m128i*)&y[i*ystride], res);
  }
}

void od_filter_dering_orthogonal_8x8(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int threshold, int dir) {
  int i;
  int k;
  int offset;
  __m128i res;
  __m128i p;
  __m128i cmp;
  __m128i row;
  __m128i sum;
  __m128i athresh;
  if (dir <= 4) offset = OD_FILT_BSTRIDE;
  else offset = 1;
  for (i = 0; i < 8; i++) {
    sum = _mm_set1_epi16(0);
    row = _mm_loadu_si128((__m128i*)&in[i*OD_FILT_BSTRIDE]);
    /* athresh = OD_MINI(threshold, threshold/3
        + abs(in[i*OD_FILT_BSTRIDE] - x[i*xstride])) */
    athresh = _mm_min_epi16(_mm_set1_epi16(threshold),
     _mm_add_epi16(_mm_set1_epi16(threshold/3),
     _mm_abs_epi16(_mm_sub_epi16(row,
     _mm_loadu_si128((__m128i*)&x[i*xstride])))));

    for (k = 1; k <= 2; k++) {
      /* p = in[i*OD_FILT_BSTRIDE + k*offset] - row */
      p = _mm_sub_epi16(_mm_loadu_si128((__m128i*)&in[i*OD_FILT_BSTRIDE +
       k*offset]), row);
      /* if (abs(p) < athresh) sum += p */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), athresh);
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);

      /* p = in[i*OD_FILT_BSTRIDE - k*offset] - row */
      p = _mm_sub_epi16(_mm_loadu_si128((__m128i*)&in[i*OD_FILT_BSTRIDE -
       k*offset]), row);
      /* if (abs(p) < athresh) sum += p */
      cmp = _mm_cmplt_epi16(_mm_abs_epi16(p), athresh);
      p = _mm_and_si128(p, cmp);
      sum = _mm_add_epi16(sum, p);
    }

    /* row + ((3*sum + 8) >> 4) */
    res = _mm_mullo_epi16(sum, _mm_set1_epi16(3));
    res = _mm_add_epi16(res, _mm_set1_epi16(8));
    res = _mm_srai_epi16(res, 4);
    res = _mm_add_epi16(res, row);
    _mm_storeu_si128((__m128i*)&y[i*ystride], res);
  }
}

typedef void (*od_filter_dering_orthogonal_fixed_func)(int16_t *y,
 int ystride, int16_t *in, int16_t *x, int xstride, int threshold, int dir);

void od_filter_dering_orthogonal_sse2(int16_t *y, int ystride, int16_t *in,
 int16_t *x, int xstride, int ln, int threshold, int dir) {
  static const od_filter_dering_orthogonal_fixed_func VTBL[4] = {
    NULL,
    NULL,
    od_filter_dering_orthogonal_4x4,
    od_filter_dering_orthogonal_8x8,
  };
  OD_ASSERT(ln >= 2 && ln < 4);
  VTBL[ln](y, ystride, in, x, xstride, threshold, dir);
#if defined(OD_CHECKASM)
  od_filter_dering_orthogonal_check(y, ystride, in, x, xstride, ln, threshold,
   dir);
#endif
}

#endif
