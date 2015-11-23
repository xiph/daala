/*Daala video codec
Copyright (c) 2006-2015 Daala project contributors.  All rights reserved.

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

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "x86int.h"
#include "cpu.h"
#include "../mc.h"
#include "../util.h"

/*Block copy functions. Copying of overlapping regions has undefined
   behavior. Only 16x16, 32x32, 64x64 show a SIMD improvement.*/

#if defined(OD_X86ASM)
#include <xmmintrin.h>

#define OD_IM_LOAD_1(_rega, _regb, _regc, _regd) \
  "#OD_IM_LOAD_1\n\t" \
  "movdqu (%[src]), " _rega "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regb "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regc "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regd "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \

#define OD_IM_STORE_1(_rega, _regb, _regc, _regd) \
  "#OD_IM_STORE_1\n\t" \
  "movdqu " _rega ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst] \n\t" \
  "movdqu " _regb ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst] \n\t" \
  "movdqu " _regc ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \
  "movdqu " _regd ",(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \

void od_copy_16x16_8_sse2(unsigned char *_dst, int _dstride,
  const unsigned char *_src, int _sstride) {
  __asm__ __volatile__(
    OD_IM_LOAD_1 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_LOAD_1 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    OD_IM_STORE_1("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_STORE_1("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    OD_IM_LOAD_1 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_LOAD_1 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    OD_IM_STORE_1("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
    OD_IM_STORE_1("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
    :[dst]"+r"(_dst),[src]"+r"(_src)
    :[dstride]"r"((ptrdiff_t)_dstride),[sstride]"r"((ptrdiff_t)_sstride)
    :"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"
  );
}

#define OD_IM_LOAD_2(_rega, _regb, _regc, _regd) \
  "#OD_IM_LOAD_2\n\t" \
  "movdqu (%[src]), " _rega "\n\t" \
  "movdqu 16(%[src]), " _regb "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \
  "movdqu (%[src]), " _regc "\n\t" \
  "movdqu 16(%[src]), " _regd "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \

#define OD_IM_STORE_2(_rega, _regb, _regc, _regd) \
  "#OD_IM_STORE_1\n\t" \
  "movdqu " _rega ",(%[dst]) \n\t" \
  "movdqu " _regb ",16(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst] \n\t" \
  "movdqu " _regc ",(%[dst]) \n\t" \
  "movdqu " _regd ",16(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \

void od_copy_32x32_8_sse2(unsigned char *_dst, int _dstride,
  const unsigned char *_src, int _sstride) {
  int i;
  for (i = 0; i < 4; i++) {
    __asm__ __volatile__(
      OD_IM_LOAD_2 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_2 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_2("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_2("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_LOAD_2 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_2 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_2("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_2("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      :[dst]"+r"(_dst),[src]"+r"(_src)
      :[dstride]"r"((ptrdiff_t)_dstride),[sstride]"r"((ptrdiff_t)_sstride)
      :"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"
    );
  }
}

#define OD_IM_LOAD_4(_rega, _regb, _regc, _regd) \
  "#OD_IM_LOAD_4\n\t" \
  "movdqu (%[src]), " _rega "\n\t" \
  "movdqu 16(%[src]), " _regb "\n\t" \
  "movdqu 32(%[src]), " _regc "\n\t" \
  "movdqu 48(%[src]), " _regd "\n\t" \
  "lea (%[src],%[sstride]),%[src] \n\t" \

#define OD_IM_STORE_4(_rega, _regb, _regc, _regd) \
  "#OD_IM_STORE_4\n\t" \
  "movdqu " _rega ",(%[dst]) \n\t" \
  "movdqu " _regb ",16(%[dst]) \n\t" \
  "movdqu " _regc ",32(%[dst]) \n\t" \
  "movdqu " _regd ",48(%[dst]) \n\t" \
  "lea (%[dst],%[dstride]),%[dst]\n\t" \

void od_copy_64x64_8_sse2(unsigned char *_dst, int _dstride,
  const unsigned char *_src, int _sstride) {
  int i;
  for (i = 0; i < 16; i++) {
    __asm__ __volatile__(
      OD_IM_LOAD_4 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_4 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_4("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_4("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_LOAD_4 ("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_LOAD_4 ("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      OD_IM_STORE_4("%%xmm0", "%%xmm1", "%%xmm2", "%%xmm3")
      OD_IM_STORE_4("%%xmm4", "%%xmm5", "%%xmm6", "%%xmm7")
      :[dst]"+r"(_dst),[src]"+r"(_src)
      :[dstride]"r"((ptrdiff_t)_dstride),[sstride]"r"((ptrdiff_t)_sstride)
      :"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7"
    );
  }
}

#if defined(OD_CHECKASM)
void od_coeff_to_ref_buf_check(od_state *state,
 unsigned char *dst, int dst_xstride, int dst_ystride,
 od_coeff *src, int src_ystride, int lossless_p,
 int w, int h) {
  int x;
  int y;
  unsigned char *tmp;
  unsigned char *_tmp;
  _tmp = tmp = malloc(h * dst_ystride * dst_xstride);
  od_coeff_to_ref_buf_c(state, tmp, dst_xstride, dst_ystride,
   src, src_ystride, lossless_p, w, h);
  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
      if (dst_xstride == 1 ? tmp[x] != dst[x] :
       ((int16_t *)tmp)[x] != ((int16_t *)dst)[x]) {
        fprintf(stderr, "od_coeff_to_ref_buf_sse2 check failed.\n");
        OD_ASSERT(0);
        free(_tmp);
        return;
      }
    }
    tmp += dst_ystride;
    dst += dst_ystride;
  }
  free(_tmp);
}
#endif

void od_coeff_to_ref_buf_sse2(od_state *state,
 unsigned char *dst, int dst_xstride, int dst_ystride,
 od_coeff *src, int src_ystride, int lossless_p,
 int w, int h) {
  int x;
  int y;
  int coeff_shift;
  __m128i a;
  __m128i b;
  __m128i c;
  __m128i round;
  __m128i offset;
#if defined(OD_CHECKASM)
  unsigned char *_dst;
  od_coeff *_src;
  _dst = dst;
  _src = src;
#endif
  /*Use SSE2 version only when w is a multiple of 16.*/
  if (w & 7) {
    od_coeff_to_ref_buf_c(state, dst, dst_xstride, dst_ystride,
     src, src_ystride, lossless_p, w, h);
    return;
  }
  OD_ASSERT(sizeof(od_coeff) == 4);
  if (dst_xstride == 1) {
    /*The references are running at 8 bits.
      The transforms and coefficients may be operating at 8 bits (during
       lossless coding) or 8 + OD_COEFF_SHIFT bits (during lossy coding).*/
    coeff_shift = lossless_p ?
     (state->info.bitdepth_mode - OD_BITDEPTH_MODE_8)*2 :
     OD_COEFF_SHIFT;
    round = _mm_set1_epi16(1 << coeff_shift >> 1);
    offset = _mm_set1_epi16(128);
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x += 16) {
        /*Load 4 32-bit vectors and pack them into 2 16-bit vectors.
          The saturating 32-bit to 16-bit pack operation doesn't quite
           match the C version which only saturates after shifting and
           adding the offset.*/
        a = _mm_packs_epi32(
         _mm_loadu_si128((__m128i*)(src + x)),
         _mm_loadu_si128((__m128i*)(src + x + 4)));
        b = _mm_packs_epi32(
         _mm_loadu_si128((__m128i*)(src + x + 8)),
         _mm_loadu_si128((__m128i*)(src + x + 12)));
        /*Round, shift and offset.*/
        a = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(a, round), coeff_shift),
         offset);
        b = _mm_add_epi16(_mm_srai_epi16(_mm_add_epi16(b, round), coeff_shift),
         offset);
        /*Saturate, pack and store.*/
        c = _mm_packus_epi16(a, b);
        _mm_storeu_si128((__m128i*)(dst + x), c);
      }
      dst += dst_ystride;
      src += src_ystride;
    }
  } else {
    /*The references are running at greater than 8 bits, implying FPR.
      An FPR reference must run at full depth (8 + OD_COEFF_SHIFT).
      The transforms and coefficients may be operating at any supported
       bit depth (8, 10 or 12) */
    coeff_shift = lossless_p
     ? OD_COEFF_SHIFT - (state->info.bitdepth_mode - OD_BITDEPTH_MODE_8)*2
     : 0;
    round = _mm_set1_epi16(1 << 8 + OD_COEFF_SHIFT >> 1);
    for (y = 0; y < h; y++) {
      for (x = 0; x < w; x += 16) {
        a = _mm_packs_epi32(
         _mm_loadu_si128((__m128i*)(src + x)),
         _mm_loadu_si128((__m128i*)(src + x + 4)));
        b = _mm_packs_epi32(
         _mm_loadu_si128((__m128i*)(src + x + 8)),
         _mm_loadu_si128((__m128i*)(src + x + 12)));
        /*Round and shift.*/
        a = _mm_add_epi16(_mm_srai_epi16(a, coeff_shift), round);
        b = _mm_add_epi16(_mm_srai_epi16(b, coeff_shift), round);
        /*Store.*/
        _mm_storeu_si128((__m128i*)(dst + x*sizeof(int16_t)), a);
        _mm_storeu_si128((__m128i*)(dst + x*sizeof(int16_t) + 16), b);
      }
      dst += dst_ystride;
      src += src_ystride;
    }
  }
#if defined(OD_CHECKASM)
  od_coeff_to_ref_buf_check(state, _dst, dst_xstride, dst_ystride,
   _src, src_ystride, lossless_p, w, h);
#endif
}
#endif
