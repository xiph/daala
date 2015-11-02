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

#if defined(HAVE_CONFIG_H)
# include "config.h"
#endif

#include "x86enc.h"
#include "cpu.h"

#if defined(OD_X86ASM)

#include <stdio.h>

void od_enc_opt_vtbl_init_x86(od_enc_ctx *enc) {
  od_enc_opt_vtbl_init_c(enc);
#if defined(OD_GCC_INLINE_ASSEMBLY)
  if (enc->state.cpu_flags & OD_CPU_X86_SSE) {
    enc->opt_vtbl.mc_compute_sad_4x4 =
     od_mc_compute_sad8_4x4_sse;
    enc->opt_vtbl.mc_compute_sad_8x8 =
     od_mc_compute_sad8_8x8_sse;
  }
  if (enc->state.cpu_flags & OD_CPU_X86_SSE2) {
    enc->opt_vtbl.mc_compute_sad_16x16 =
     od_mc_compute_sad8_16x16_sse2;
    enc->opt_vtbl.mc_compute_sad_32x32 =
     od_mc_compute_sad8_32x32_sse2;
    enc->opt_vtbl.mc_compute_satd_4x4 =
     od_mc_compute_satd8_4x4_sse2;
    enc->opt_vtbl.mc_compute_satd_8x8 =
     od_mc_compute_satd8_8x8_sse2;
    enc->opt_vtbl.mc_compute_satd_16x16 =
     od_mc_compute_satd8_16x16_sse2;
    enc->opt_vtbl.mc_compute_satd_32x32 =
     od_mc_compute_satd8_32x32_sse2;
    enc->opt_vtbl.mc_compute_satd_64x64 =
     od_mc_compute_satd8_64x64_sse2;
  }
#endif
}

#endif
