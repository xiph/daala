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

#if defined(OD_X86ASM)

# if defined(OD_CHECKASM)

#  include <stdio.h>

void od_mc_compute_sad_check(const unsigned char *src, int systride,
 const unsigned char *ref, int dystride, int dxstride, int w, int h, int sad) {
  int c_sad;
  c_sad = od_mc_compute_sad_c(src, systride, ref, dystride, dxstride, w, h);
  if (sad != c_sad) {
    fprintf(stderr, "od_mc_compute_sad %ix%i check failed: %i!=%i\n",
     w, h, sad, c_sad);
  }
  OD_ASSERT(sad == c_sad);
}
# endif

/*Handle one 4x4 block with dxstride == 1.*/
int od_mc_compute_sad_4x4_xstride_1_sse(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride){
  ptrdiff_t srow;
  ptrdiff_t drow;
  int ret;
  __asm__ __volatile__(
    "movd (%[src]), %%mm0\n"
    "movd (%[src], %[systride]), %%mm2\n"
    "punpckldq %%mm2, %%mm0\n"
    "movd (%[ref]), %%mm1\n"
    "movd (%[ref], %[dystride]), %%mm3\n"
    "punpckldq %%mm3, %%mm1\n"
    "psadbw %%mm1, %%mm0\n"
    "lea (%[src], %[systride], 2), %[srow]\n"
    "lea (%[ref], %[dystride], 2), %[drow]\n"
    : [srow]"=&r"(srow), [drow]"=r"(drow)
    : [src]"r"(src), [systride]"r"((ptrdiff_t)systride), [ref]"r"(ref),
     [dystride]"r"((ptrdiff_t)dystride)
  );
  __asm__ __volatile__(
    "movd (%[srow]), %%mm4\n"
    "movd (%[srow], %[systride]), %%mm5\n"
    "punpckldq %%mm5, %%mm4\n"
    "movd (%[drow]), %%mm6\n"
    "movd (%[drow], %[dystride]), %%mm7\n"
    "punpckldq %%mm7, %%mm6\n"
    "psadbw %%mm6, %%mm4\n"
    "paddd %%mm4, %%mm0\n"
    "movd %%mm0, %[ret]\n"
    : [ret]"=r"(ret)
    : [systride]"r"((ptrdiff_t)systride),
     [dystride]"r"((ptrdiff_t)dystride), [srow]"r"(srow), [drow]"r"(drow)
  );
# if defined(OD_CHECKASM)
  od_mc_compute_sad_check(src, systride, ref, dystride, 1, 4, 4, ret);
# endif
  return ret;
}

/*Handle one 8x8 block with dxstride == 1.*/
int od_mc_compute_sad_8x8_xstride_1_sse(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride){
  const unsigned char *srow;
  const unsigned char *drow;
  int ret;
  int i;
  srow = src;
  drow = ref;
  __asm__ __volatile__(
    "pxor %mm2, %mm2\n"
  );
  for (i = 0; i < 8; i++) {
    __asm__ __volatile__(
      "movq (%[srow]), %%mm0\n"
      "movq (%[drow]), %%mm1\n"
      "psadbw %%mm1, %%mm0\n"
      "paddd %%mm0, %%mm2\n"
      :
      : [srow]"r"(srow), [drow]"r"(drow)
    );
    srow += systride;
    drow += dystride;
  }
  __asm__ __volatile__(
    "movd %%mm2, %[ret]\n"
    :[ret]"=r"(ret)
  );
# if defined(OD_CHECKASM)
  od_mc_compute_sad_check(src, systride, ref, dystride, 1, 8, 8, ret);
# endif
  return ret;
}

/*Handle one 16x16 block with dxstride == 1.*/
int od_mc_compute_sad_16x16_xstride_1_sse2(const unsigned char *src,
 int systride, const unsigned char *ref, int dystride){
  const unsigned char *srow;
  const unsigned char *drow;
  int ret;
  int i;
  srow = src;
  drow = ref;
  __asm__ __volatile__(
    "pxor %xmm2, %xmm2\n\t"
  );
  for (i = 0; i < 16; i++) {
    __asm__ __volatile__(
      "movdqu (%[drow]), %%xmm1\n"
      "movdqu (%[srow]), %%xmm0\n"
      "psadbw %%xmm1, %%xmm0\n"
      "paddq %%xmm0, %%xmm2\n"
      :
      : [srow]"r"(srow), [drow]"r"(drow)
    );
    srow += systride;
    drow += dystride;
  }
  __asm__ __volatile__(
    "movdqa %%xmm2, %%xmm0\n"
    "punpckhqdq %%xmm2, %%xmm0\n"
    "paddq %%xmm2, %%xmm0\n"
    "movd %%xmm0, %[ret]\n"
    : [ret]"=r"(ret)
  );
# if defined(OD_CHECKASM)
  od_mc_compute_sad_check(src, systride, ref, dystride, 1, 16, 16, ret);
# endif
  return ret;
}

#endif
