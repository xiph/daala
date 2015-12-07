/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

#if !defined(_x86_x86int_H)
# define _x86_x86int_H (1)
# include "../state.h"

# if OD_GNUC_PREREQ(3, 0, 0)
#  define OD_SIMD_INLINE static __inline __attribute__((always_inline))
# else
#  define OD_SIMD_INLINE static
# endif

void od_state_opt_vtbl_init_x86(od_state *_state);
extern const od_filter_dering_direction_func
 OD_DERING_DIRECTION_SSE2[OD_DERINGSIZES];
extern const od_filter_dering_orthogonal_func
 OD_DERING_ORTHOGONAL_SSE2[OD_DERINGSIZES];
void od_mc_predict1fmv8_sse2(od_state *state, unsigned char *_dst,
 const unsigned char *_src, int _systride, int32_t _mvx, int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz);
void od_mc_predict1fmv16_sse2(od_state *state, unsigned char *_dst,
 const unsigned char *_src, int _systride, int32_t _mvx, int32_t _mvy,
 int _log_xblk_sz, int _log_yblk_sz);
void od_mc_blend_full8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz);
void od_mc_blend_full_split8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz);
void od_bin_fdct4x4_sse2(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_fdct4x4_sse41(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct4x4_sse2(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct4x4_sse41(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_fdct8x8_sse2(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_fdct8x8_sse41(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct8x8_sse2(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct8x8_sse41(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_fdct8x8_avx2(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct8x8_avx2(od_coeff *x, int xstride,
 const od_coeff *y, int ystride);
void od_copy_16x16_8_sse2(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride);
void od_copy_32x32_8_sse2(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride);
void od_copy_64x64_8_sse2(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride);
void od_filter_dering_direction_4x4_sse2(int16_t *y, int ystride,
 int16_t *in, int threshold, int dir);
void od_filter_dering_direction_8x8_sse2(int16_t *y, int ystride,
 int16_t *in, int threshold, int dir);
void od_filter_dering_orthogonal_4x4_sse2(int16_t *y, int ystride,
 int16_t *in, int16_t *x, int xstride, int threshold, int dir);
void od_filter_dering_orthogonal_8x8_sse2(int16_t *y, int ystride,
 int16_t *in, int16_t *x, int xstride, int threshold, int dir);
#endif
