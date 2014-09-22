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

#include <arm_neon.h>
#include "../dct.h"
#include "armint.h"

#define OD_UNBIASED_RSHIFT(a, b) \
  vshrq_n_s32(vaddq_s32(vshlq_n_s32(a, 32 - b), a), b)

OD_SIMD_INLINE void od_overflow_check_epi32(int32x4_t val, ogg_int32_t scale,
 ogg_int32_t offset, int idx) {
#if defined(OD_DCT_TEST) && defined(OD_DCT_CHECK_OVERFLOW)
  ogg_int32_t mem[4];
  int n;
  _mm_store_si128((int32x4_t *)mem, val);
  for (n = 0; n < 4; n++) {
    OD_DCT_OVERFLOW_CHECK(mem[n], scale, offset, idx);
  }
#endif
  (void)val;
  (void)scale;
  (void)offset;
  (void)idx;
}

#define OD_DCT_MUL(val, scale, offset, shift) \
  vshrq_n_s32(vmlaq_s32(vdupq_n_s32(offset), val, vdupq_n_s32(scale)), shift)

OD_SIMD_INLINE uint64x2x2_t od_vswpq_u64(uint64x2_t a, uint64x2_t b) {
  uint64x2x2_t x;
  x.val[0] = vcombine_u64(vget_low_u64(a), vget_low_u64(b));
  x.val[1] = vcombine_u64(vget_high_u64(a), vget_high_u64(b));
  return x;
}

OD_SIMD_INLINE void od_transpose4(int32x4_t *t0, int32x4_t *t1,
 int32x4_t *t2, int32x4_t *t3) {
  uint64x2x2_t a, b;
  int32x4x2_t x;
  a = od_vswpq_u64(vreinterpretq_u64_s32(*t0), vreinterpretq_u64_s32(*t2));
  b = od_vswpq_u64(vreinterpretq_u64_s32(*t1), vreinterpretq_u64_s32(*t3));
  x = vtrnq_s32(vreinterpretq_s32_u64(a.val[0]),
		vreinterpretq_s32_u64(a.val[1]));
  *t0 = x.val[0];
  *t1 = x.val[1];
  x = vtrnq_s32(vreinterpretq_s32_u64(b.val[0]),
		vreinterpretq_s32_u64(b.val[1]));
  *t2 = x.val[0];
  *t3 = x.val[1];
}

OD_SIMD_INLINE void od_load4(const od_coeff *x, int xstride,
 int32x4_t *t0, int32x4_t *t1, int32x4_t *t2, int32x4_t *t3) {
  *t0 = vld1q_s32((const int *)(x + 0*xstride));
  *t1 = vld1q_s32((const int *)(x + 1*xstride));
  *t2 = vld1q_s32((const int *)(x + 2*xstride));
  *t3 = vld1q_s32((const int *)(x + 3*xstride));
}

OD_SIMD_INLINE void od_store4(od_coeff *x, int xstride,
 int32x4_t t0, int32x4_t t1, int32x4_t t2, int32x4_t t3) {
  vst1q_s32((int *)(x + 0*xstride), t0);
  vst1q_s32((int *)(x + 1*xstride), t1);
  vst1q_s32((int *)(x + 2*xstride), t2);
  vst1q_s32((int *)(x + 3*xstride), t3);
}

OD_SIMD_INLINE void od_fdct4_kernel(int32x4_t *x0, int32x4_t *x1,
 int32x4_t *x2, int32x4_t *x3) {
  /*9 adds, 2 shifts, 3 "muls".*/
  int32x4_t t0 = *x0;
  int32x4_t t2 = *x1;
  int32x4_t t1 = *x2;
  int32x4_t t3 = *x3;
  int32x4_t t2h;
  /*+1/-1 butterflies:*/
  t3 = vsubq_s32(t0, t3);
  t2 = vaddq_s32(t2, t1);
  t2h = OD_UNBIASED_RSHIFT(t2, 1);
  t1 = vsubq_s32(t2h, t1);
  t0 = vsubq_s32(t0, OD_UNBIASED_RSHIFT(t3, 1));
  /*+ Embedded 2-point type-II DCT.*/
  t0 = vaddq_s32(t0, t2h);
  t2 = vsubq_s32(t0, t2);
  /*+ Embedded 2-point type-IV DST.*/
  /*23013/32768 ~= 4*sin(\frac{\pi}{8}) - 2*tan(\frac{\pi}{8}) ~=
     0.70230660471416898931046248770220*/
  od_overflow_check_epi32(t1, 23013, 16384, 0);
  t3 = vsubq_s32(t3, OD_DCT_MUL(t1, 23013, 16384, 15));
  /*21407/32768~=\sqrt{1/2}*cos(\frac{\pi}{8}))
     ~=0.65328148243818826392832158671359*/
  od_overflow_check_epi32(t3, 21407, 16384, 1);
  t1 = vaddq_s32(t1, OD_DCT_MUL(t3, 21407, 16384, 15));
  /*18293/16384 ~= 4*sin(\frac{\pi}{8}) - tan(\frac{\pi}{8}) ~=
     1.1165201670872640381121512119119*/
  od_overflow_check_epi32(t3, 18293, 8192, 2);
  t3 = vsubq_s32(t3, OD_DCT_MUL(t1, 18293, 8192, 14));
  od_transpose4(&t0, &t1, &t2, &t3);
  *x0 = t0;
  *x1 = t1;
  *x2 = t2;
  *x3 = t3;
}

void od_bin_fdct4x4_neon(od_coeff *y, int ystride,
 const od_coeff *x, int xstride) {
  int32x4_t t0;
  int32x4_t t1;
  int32x4_t t2;
  int32x4_t t3;
#if defined(OD_CHECKASM)
  od_coeff ref[4*4];
  od_bin_fdct4x4(ref, 4, x, xstride);
#endif
  od_load4(x, xstride, &t0, &t1, &t2, &t3);
  od_fdct4_kernel(&t0, &t1, &t2, &t3);
  od_fdct4_kernel(&t0, &t1, &t2, &t3);
  od_store4(y, ystride, t0, t1, t2, t3);
#if defined(OD_CHECKASM)
  od_dct_check(0, ref, y, ystride);
#endif
}

OD_SIMD_INLINE void od_idct4_kernel(int32x4_t *y0, int32x4_t *y1,
 int32x4_t *y2, int32x4_t *y3) {
  int32x4_t t0 = *y0;
  int32x4_t t1 = *y1;
  int32x4_t t2 = *y2;
  int32x4_t t3 = *y3;
  int32x4_t t2h;
  od_transpose4(&t0, &t1, &t2, &t3);
  t3 = vaddq_s32(t3, OD_DCT_MUL(t1, 18293, 8192, 14));
  t1 = vsubq_s32(t1, OD_DCT_MUL(t3, 21407, 16384, 15));
  t3 = vaddq_s32(t3, OD_DCT_MUL(t1, 23013, 16384, 15));
  t2 = vsubq_s32(t0, t2);
  t2h = OD_UNBIASED_RSHIFT(t2, 1);
  t0 = vsubq_s32(t0, vsubq_s32(t2h, OD_UNBIASED_RSHIFT(t3, 1)));
  t1 = vsubq_s32(t2h, t1);
  *y0 = t0;
  *y1 = vsubq_s32(t2, t1);
  *y2 = t1;
  *y3 = vsubq_s32(t0, t3);
}

void od_bin_idct4x4_neon(od_coeff *x, int xstride,
 const od_coeff *y, int ystride) {
  int32x4_t t0;
  int32x4_t t1;
  int32x4_t t2;
  int32x4_t t3;
#if defined(OD_CHECKASM)
  od_coeff ref[4*4];
  od_bin_idct4x4(ref, 4, y, ystride);
#endif
  od_load4(y, ystride, &t0, &t1, &t2, &t3);
  od_idct4_kernel(&t0, &t1, &t2, &t3);
  od_idct4_kernel(&t0, &t1, &t2, &t3);
  od_store4(x, xstride, t0, t1, t2, t3);
#if defined(OD_CHECKASM)
  od_dct_check(0, ref, x, xstride);
#endif
}
