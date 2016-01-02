/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_x86_x86enc_H)
# define _x86_x86enc_H (1)
# include "../encint.h"

void od_enc_opt_vtbl_init_x86(od_enc_ctx *enc);

int32_t od_mc_compute_sad8_4x4_sse(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride);
int32_t od_mc_compute_sad8_8x8_sse(const unsigned char *_src,
 int systride, const unsigned char *ref, int dystride);
int32_t od_mc_compute_sad8_16x16_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride);
int32_t od_mc_compute_sad8_32x32_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride);
int32_t od_mc_compute_sad8_64x64_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride);

int32_t od_mc_compute_satd8_4x4_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd8_8x8_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd8_16x16_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd8_32x32_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd8_64x64_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);

int32_t od_mc_compute_satd16_4x4_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd16_8x8_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd16_16x16_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd16_32x32_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);
int32_t od_mc_compute_satd16_64x64_sse2(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride);

#endif
