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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "x86int.h"
#include "cpu.h"

#if defined(OD_X86ASM)

/*So here are the constraints we have:
  We want (as much as possible) the same code to compile on x86-32 and x86-64.
  We want the same code to be used with and without -fomit-frame-pointer and
   with and without -fPIC.
  These are, in fact, the reason we're using inline gcc assembly instead of
   separate asm files.
  Separate files using the Intel syntax might have a chance of compiling on
   win32 without requiring mingw32, but would introduce a dependency on a
   separate assembler like nasm everywhere else, would require us to maintain
   multiple versions of the code for, e.g., 32- and 64-bit processors, and
   would require horrible macro garbage to correct for things like name
   mangling, which are not consistent between platforms.
  gcc can take care of all that crap for us, _and_ I can cross-compile for
   win32 from Linux.
  That certainly made the decision easy for me.

  These are the rules:
  We cannot use either %ebp (not available without -fomit-frame-pointer) or
   %ebx (not available with -fPIC) as general purpose registers, nor can we
   rely on gcc being able to use them when it constructs inputs/outputs for us.
  This gives the register-starved x86 architecture just 5 general purpose
   registers.
  Compile without -fomit-frame-pointer and with -fPIC during testing to make
   sure you don't break register allocation.
  To make matters worse, we can't combine static variables in "m" inputs with
   an index register or scale, because with -fPIC they might be using %ebx as a
   base register, and there's no syntactically valid way to combine the two.
  Therefore we have to load them into a register first, and then index them.
  Finally, to use scratch registers that contain pointers or pointer offsets,
   we must declare a local C variable of type ptrdiff_t and let gcc allocate it
   to a register for us.
  This will automatically select a 32- or 64-bit register as appropriate for
   the target platform.
  Type long is not good enough thanks to win64, where long is still 32 bits.*/

/*Constants used for the bilinear blending weights.*/
static const unsigned short __attribute__((aligned(16),used)) OD_BIL4H[8]={
  0x0,0x1,0x2,0x3,0x0,0x1,0x2,0x3
};
static const unsigned short __attribute__((aligned(16),used)) OD_BILH[16]={
  0x0,0x1,0x2,0x3,0x4,0x5,0x6,0x7,
  0x8,0x9,0xA,0xB,0xC,0xD,0xE,0xF
};
static const unsigned short __attribute__((aligned(16),used)) OD_BIL4V[64]={
  0x0,0x0,0x0,0x0,0x1,0x1,0x1,0x1,
  0x2,0x2,0x2,0x2,0x3,0x3,0x3,0x3,
  0x4,0x4,0x4,0x4,0x5,0x5,0x5,0x5,
  0x6,0x6,0x6,0x6,0x7,0x7,0x7,0x7,
  0x8,0x8,0x8,0x8,0x9,0x9,0x9,0x9,
  0xA,0xA,0xA,0xA,0xB,0xB,0xB,0xB,
  0xC,0xC,0xC,0xC,0xD,0xD,0xD,0xD,
  0xE,0xE,0xE,0xE,0xF,0xF,0xF,0xF
};
static const unsigned short __attribute__((aligned(16),used)) OD_BILV[128]={
  0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,
  0x1,0x1,0x1,0x1,0x1,0x1,0x1,0x1,
  0x2,0x2,0x2,0x2,0x2,0x2,0x2,0x2,
  0x3,0x3,0x3,0x3,0x3,0x3,0x3,0x3,
  0x4,0x4,0x4,0x4,0x4,0x4,0x4,0x4,
  0x5,0x5,0x5,0x5,0x5,0x5,0x5,0x5,
  0x6,0x6,0x6,0x6,0x6,0x6,0x6,0x6,
  0x7,0x7,0x7,0x7,0x7,0x7,0x7,0x7,
  0x8,0x8,0x8,0x8,0x8,0x8,0x8,0x8,
  0x9,0x9,0x9,0x9,0x9,0x9,0x9,0x9,
  0xA,0xA,0xA,0xA,0xA,0xA,0xA,0xA,
  0xB,0xB,0xB,0xB,0xB,0xB,0xB,0xB,
  0xC,0xC,0xC,0xC,0xC,0xC,0xC,0xC,
  0xD,0xD,0xD,0xD,0xD,0xD,0xD,0xD,
  0xE,0xE,0xE,0xE,0xE,0xE,0xE,0xE,
  0xF,0xF,0xF,0xF,0xF,0xF,0xF,0xF
};

#if defined(OD_CHECKASM)
void od_mc_predict1fmv8_check(unsigned char *_dst,const unsigned char *_src,
 int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz){
  unsigned char dst[16*16];
  int           xblk_sz;
  int           yblk_sz;
  int           failed;
  int           i;
  int           j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  failed=0;
  od_mc_predict1fmv8_c(dst,_src,_systride,_mvx,_mvy,
   _log_xblk_sz,_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      if(_dst[i+(j<<_log_xblk_sz)]!=dst[i+(j<<_log_xblk_sz)]){
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         _dst[i+(j<<_log_xblk_sz)],dst[i+(j<<_log_xblk_sz)],i,j);
        failed=1;
      }
    }
  }
  if(failed){
    fprintf(stderr,"od_mc_predict1fmv8 %ix%i check failed.\n",
     (1<<_log_xblk_sz),(1<<_log_yblk_sz));
  }
  OD_ASSERT(!failed);
}
#endif

/*Initializes %%xmm7 to contain the mask 0x00FF00FF00FF00FF00FF00FF00FF00FF,
   %xmm5 to contain 8 replicated copies of %[hscale], and %xmm6 to contain 8
   replicated copies of %[vscale].
  We can't keep both %xmm5 and %xmm6 free throughout the HV blend, so we spill
   %xmm6 to %[vscale8] as it is used less often, and at the end, where there is
   less register pressure.*/
#define OD_MC_PREDICT1FMV8HV_PROLOG \
  "#OD_MC_PREDICT1FMV8HV_PROLOG\n\t" \
  "movd %[hscale],%%xmm5\n\t" \
  "pcmpeqw %%xmm7,%%xmm7\n\t" \
  "movd %[vscale],%%xmm6\n\t" \
  "pshuflw $0x00,%%xmm5,%%xmm5\n\t" \
  "pshuflw $0x00,%%xmm6,%%xmm6\n\t" \
  "pshufd $0x00,%%xmm5,%%xmm5\n\t" \
  "pshufd $0x00,%%xmm6,%%xmm6\n\t" \
  "psrlw $8,%%xmm7\n\t" \
  "movdqa %%xmm6,%[vscale8]\n\t" \

/*Initializes %%xmm7 to contain the mask 0x00FF00FF00FF00FF00FF00FF00FF00FF and
   %xmm5 to contain 8 replicated copies of %[hscale].*/
#define OD_MC_PREDICT1FMV8H_PROLOG \
  "#OD_MC_PREDICT1FMV8H_PROLOG\n\t" \
  "movd %[hscale],%%xmm5\n\t" \
  "pcmpeqw %%xmm7,%%xmm7\n\t" \
  "pshuflw $0x00,%%xmm5,%%xmm5\n\t" \
  "psrlw $8,%%xmm7\n\t" \
  "pshufd $0x00,%%xmm5,%%xmm5\n\t" \

/*Initializes %%xmm7 to contain the mask 0x00FF00FF00FF00FF00FF00FF00FF00FF and
   %xmm6 to contain 8 replicated copies of %[vscale].*/
#define OD_MC_PREDICT1FMV8V_PROLOG \
  "#OD_MC_PREDICT1FMV8V_PROLOG\n\t" \
  "movd %[vscale],%%xmm6\n\t" \
  "pcmpeqw %%xmm7,%%xmm7\n\t" \
  "pshuflw $0x00,%%xmm6,%%xmm6\n\t" \
  "psrlw $8,%%xmm7\n\t" \
  "pshufd $0x00,%%xmm6,%%xmm6\n\t" \

/*Initializes %%xmm7 to contain the mask 0x00FF00FF00FF00FF00FF00FF00FF00F.*/
#define OD_MC_PREDICT1FMV8_PROLOG \
  "#OD_MC_PREDICT1FMV8_PROLOG\n\t" \
  "pcmpeqw %%xmm7,%%xmm7\n\t" \
  "psrlw $8,%%xmm7\n\t" \

#define OD_MC_PREDICT1FMV8HV_4x4 \
  "#OD_MC_PREDICT1FMV8HV_4x4\n\t" \
  "movq (%[src],%[systride],2),%%xmm1\n\t" \
  "lea (%[src],%[systride],2),%[a]\n\t" \
  "movq (%[a],%[systride],4),%%xmm3\n\t" \
  "pslldq $8,%%xmm1\n\t" \
  "movq (%[src]),%%xmm0\n\t" \
  "pslldq $8,%%xmm3\n\t" \
  "movq (%[src],%[systride],4),%%xmm2\n\t" \
  "por %%xmm1,%%xmm0\n\t" \
  "por %%xmm3,%%xmm2\n\t" \
  "movdqa %%xmm0,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "movdqa %%xmm2,%%xmm3\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "psrlw $8,%%xmm1\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  /*The blending here is subtle. \
    Like OD_IMV_SAMPLE16, we want to compute (a*((1<<16)-w)+b*w)>>16. \
    However, for fixed motion vectors, the weights only ever use the top 5 \
     bits of the word: 2 for 1/8th pel resolution+2 for chroma decimation+1 \
     for split edges). \
    But we're truncating to 8 bits anyway, so we can ignore the extra copy of \
     this weight that may be added to the result. \
    We already have a suitable mask in xmm7 to remove them, which we used to \
     split apart the even and odd pixel values. \
    This is simply used to clear these bits before additional computations, \
     and before packing (since SSE does not have a non-saturating pack \
     instruction), saving considerable computation over the general case and \
     using one less register.*/ \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm1\n\t" \
  "lea (%[src],%[systride]),%[a]\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "movq (%[a]),%%xmm6\n\t" \
  "paddw %%xmm3,%%xmm2\n\t" \
  "movq (%[a],%[systride],2),%%xmm1\n\t" \
  "pslldq $8,%%xmm1\n\t" \
  "lea (%[a],%[systride],4),%[a]\n\t" \
  "por %%xmm6,%%xmm1\n\t" \
  "movq (%[a],%[systride],2),%%xmm4\n\t" \
  "movdqa %%xmm1,%%xmm6\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "movq (%[a]),%%xmm3\n\t" \
  "pslldq $8,%%xmm4\n\t" \
  "psrlw $8,%%xmm6\n\t" \
  "por %%xmm3,%%xmm4\n\t" \
  "psubw %%xmm1,%%xmm6\n\t" \
  "movdqa %%xmm4,%%xmm3\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "pmulhuw %%xmm5,%%xmm6\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "paddw %%xmm6,%%xmm1\n\t" \
  "psubw %%xmm4,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "movdqa %[vscale8],%%xmm6\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "paddw %%xmm3,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm4\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "paddw %%xmm4,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8H_4x4 \
  "#OD_MC_PREDICT1FMV8H_4x4\n\t" \
  "movq (%[src],%[systride],2),%%xmm1\n\t" \
  "lea (%[src],%[systride],2),%[a]\n\t" \
  "pslldq $8,%%xmm1\n\t" \
  "movq (%[a],%[systride],4),%%xmm3\n\t" \
  "pslldq $8,%%xmm3\n\t" \
  "movq (%[src]),%%xmm0\n\t" \
  "por %%xmm1,%%xmm0\n\t" \
  "movq (%[src],%[systride],4),%%xmm2\n\t" \
  "por %%xmm3,%%xmm2\n\t" \
  "movdqa %%xmm0,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "movdqa %%xmm2,%%xmm3\n\t" \
  "psrlw $8,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "paddw %%xmm3,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8V_4x4 \
  "#OD_MC_PREDICT1FMV8V_4x4\n\t" \
  "movq (%[src],%[systride],2),%%xmm1\n\t" \
  "lea (%[src],%[systride],2),%[a]\n\t" \
  "pslldq $8,%%xmm1\n\t" \
  "movq (%[a],%[systride],4),%%xmm3\n\t" \
  "pslldq $8,%%xmm3\n\t" \
  "movq (%[src]),%%xmm0\n\t" \
  "por %%xmm1,%%xmm0\n\t" \
  "movq (%[src],%[systride],4),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "lea (%[src],%[systride]),%[a]\n\t" \
  "por %%xmm3,%%xmm2\n\t" \
  "movq (%[a],%[systride],2),%%xmm3\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "movq (%[a],%[systride],4),%%xmm4\n\t" \
  "pslldq $8,%%xmm3\n\t" \
  "movq (%[a]),%%xmm1\n\t" \
  "lea (%[a],%[systride],2),%[a]\n\t" \
  "movq (%[a],%[systride],4),%%xmm5\n\t" \
  "por %%xmm3,%%xmm1\n\t" \
  "pslldq $8,%%xmm5\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "por %%xmm5,%%xmm4\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm4\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "paddw %%xmm4,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8_4x4 \
  "#OD_MC_PREDICT1FMV8V_4x4\n\t" \
  "movq (%[src],%[systride],2),%%xmm1\n\t" \
  "lea (%[src],%[systride],2),%[a]\n\t" \
  "pslldq $8,%%xmm1\n\t" \
  "movq (%[a],%[systride],4),%%xmm3\n\t" \
  "pslldq $8,%%xmm3\n\t" \
  "movq (%[src]),%%xmm0\n\t" \
  "por %%xmm1,%%xmm0\n\t" \
  "movq (%[src],%[systride],4),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "por %%xmm3,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

/*TODO:
  SSE3's lddqu instruction would actually help out here, but it is not
   supported on P4s.
  I'm not sure if it's worth adding a whole extra version of all this code to
   the binary for the small advantage, but it could be done by passing in
   "lddqu" instead of "movdqu" for the *_8x2 and *_16x1 macros below.*/

#define OD_MC_PREDICT1FMV8HV_8x2(_lddqu) \
  "#OD_MC_PREDICT1FMV8HV_8x2\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  "lea (%[src],%[systride]),%[a]\n\t" \
  "movdqa %%xmm0,%%xmm1\n\t" \
  _lddqu " (%[src],%[systride],2),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "movdqa %%xmm2,%%xmm3\n\t" \
  "psrlw $8,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "paddw %%xmm3,%%xmm2\n\t" \
  _lddqu " (%[a]),%%xmm6\n\t" \
  _lddqu " (%[a],%[systride],2),%%xmm3\n\t" \
  "movdqa %%xmm6,%%xmm1\n\t" \
  "psrlw $8,%%xmm6\n\t" \
  "movdqa %%xmm3,%%xmm4\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "psubw %%xmm1,%%xmm6\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "pmulhuw %%xmm5,%%xmm6\n\t" \
  "psubw %%xmm4,%%xmm3\n\t" \
  "paddw %%xmm6,%%xmm1\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm3,%%xmm4\n\t" \
  "movdqa %[vscale8],%%xmm6\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm4\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm6,%%xmm4\n\t" \
  "paddw %%xmm4,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8H_8x2(_lddqu) \
  "#OD_MC_PREDICT1FMV8H_8x2\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  "movdqa %%xmm0,%%xmm1\n\t" \
  _lddqu " (%[src],%[systride],2),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "movdqa %%xmm2,%%xmm3\n\t" \
  "psrlw $8,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "paddw %%xmm3,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8V_8x2(_lddqu) \
  "#OD_MC_PREDICT1FMV8V_8x2\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  "lea (%[src],%[systride]),%[a]\n\t" \
  _lddqu " (%[src],%[systride],2),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  _lddqu " (%[a]),%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  _lddqu " (%[a],%[systride],2),%%xmm4\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm6,%%xmm4\n\t" \
  "paddw %%xmm4,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8_8x2(_lddqu) \
  "#OD_MC_PREDICT1FMV8_8x2\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  _lddqu " (%[src],%[systride],2),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8HV_16x1(_lddqu) \
  "#OD_MC_PREDICT1FMV8HV_16x1\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  "movdqa %%xmm0,%%xmm1\n\t" \
  _lddqu " 0x10(%[src]),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "movdqa %%xmm2,%%xmm3\n\t" \
  "psrlw $8,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "paddw %%xmm3,%%xmm2\n\t" \
  _lddqu " (%[src],%[systride]),%%xmm6\n\t" \
  "movdqa %%xmm6,%%xmm1\n\t" \
  _lddqu " 0x10(%[src],%[systride]),%%xmm3\n\t" \
  "psrlw $8,%%xmm6\n\t" \
  "movdqa %%xmm3,%%xmm4\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "psubw %%xmm1,%%xmm6\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "pmulhuw %%xmm5,%%xmm6\n\t" \
  "psubw %%xmm4,%%xmm3\n\t" \
  "paddw %%xmm6,%%xmm1\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm3,%%xmm4\n\t" \
  "movdqa %[vscale8],%%xmm6\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm4\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm6,%%xmm4\n\t" \
  "paddw %%xmm4,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8H_16x1(_lddqu) \
  "#OD_MC_PREDICT1FMV8H_16x1\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  "movdqa %%xmm0,%%xmm1\n\t" \
  _lddqu " 0x10(%[src]),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "movdqa %%xmm2,%%xmm3\n\t" \
  "psrlw $8,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "psrlw $8,%%xmm3\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm3\n\t" \
  "pmulhuw %%xmm5,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm5,%%xmm3\n\t" \
  "paddw %%xmm3,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8V_16x1(_lddqu) \
  "#OD_MC_PREDICT1FMV8V_16x1\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  "lea (%[src],%[systride]),%[a]\n\t" \
  _lddqu " 0x10(%[src]),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  _lddqu " (%[a]),%%xmm1\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  _lddqu " 0x10(%[a]),%%xmm4\n\t" \
  "pand %%xmm7,%%xmm1\n\t" \
  "pand %%xmm7,%%xmm4\n\t" \
  "psubw %%xmm0,%%xmm1\n\t" \
  "psubw %%xmm2,%%xmm4\n\t" \
  "pmulhuw %%xmm6,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm0\n\t" \
  "pmulhuw %%xmm6,%%xmm4\n\t" \
  "paddw %%xmm4,%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

#define OD_MC_PREDICT1FMV8_16x1(_lddqu) \
  "#OD_MC_PREDICT1FMV8_16x1\n\t" \
  _lddqu " (%[src]),%%xmm0\n\t" \
  _lddqu " 0x10(%[src]),%%xmm2\n\t" \
  "pand %%xmm7,%%xmm0\n\t" \
  "pand %%xmm7,%%xmm2\n\t" \
  "packuswb %%xmm2,%%xmm0\n\t" \

/*Defines a pure-C implementation with hard-coded loop limits for block sizes
   we don't want to implement manually (e.g., that have fewer than 16 bytes,
   require byte-by-byte unaligned loads, etc.).
  This should let the compiler aggressively unroll loops, etc.
  It can't vectorize it itself because of the difference in operand sizes.*/
#if 0
#define OD_MC_PREDICT1FMV8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_predict1fmv8_##_n##x##_m(unsigned char *_dst, \
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){ \
  int i; \
  int j; \
  if(_mvxf!=0){ \
    if(_mvyf!=0){ \
      for(j=0;j<(_m);j++){ \
        for(i=0;i<(_n);i++){ \
          ogg_uint32_t a; \
          ogg_uint32_t b; \
          a=_src[i<<1]+((_src[i<<1|1]-_src[i<<1])*_mvxf>>16); \
          b=(_src+_systride)[i<<1]+ \
           (((_src+_systride)[i<<1|1]-(_src+_systride)[i<<1])*_mvxf>>16); \
          _dst[i]=(unsigned char)(a+((b-a)*_mvyf>>16)); \
        } \
        _src+=_systride<<1; \
        _dst+=(_n); \
      } \
    } \
    else{ \
      for(j=0;j<(_m);j++){ \
        for(i=0;i<(_n);i++){ \
          _dst[i]=(unsigned char) \
           (_src[i<<1]+((_src[i<<1|1]-_src[i<<1])*_mvxf>>16)); \
        } \
        _src+=_systride<<1; \
        _dst+=(_n); \
      } \
    } \
  } \
  else{ \
    if(_mvyf!=0){ \
      for(j=0;j<(_m);j++){ \
        for(i=0;i<(_n);i++){ \
          _dst[i]=(unsigned char) \
           (_src[i<<1]+(((_src+_systride)[(i<<1)]-_src[i<<1])*_mvyf>>16)); \
        } \
        _src+=_systride<<1; \
        _dst+=(_n); \
      } \
    } \
    else{ \
      for(j=0;j<(_m);j++){ \
        for(i=0;i<(_n);i++)_dst[i]=_src[i<<1]; \
        _src+=_systride<<1; \
        _dst+=(_n); \
      } \
    } \
  } \
} \

#else
/*The above is great and all, but not really worth the 24K of code it generates
   considering how seldom most of it is used.
  Even gcc won't inline these itself.*/
#define OD_MC_PREDICT1FMV8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_predict1fmv8_##_n##x##_m(unsigned char *_dst, \
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){ \
  od_mc_predict1fmv8_c(_dst,_src,_systride, \
   (ogg_int32_t)_mvxf,(ogg_int32_t)_mvyf,_log_xblk_sz,_log_yblk_sz); \
} \

#endif

OD_MC_PREDICT1FMV8_C(1,1,0,0)
OD_MC_PREDICT1FMV8_C(1,2,0,1)
OD_MC_PREDICT1FMV8_C(1,4,0,2)
OD_MC_PREDICT1FMV8_C(1,8,0,3)
OD_MC_PREDICT1FMV8_C(1,16,0,4)

OD_MC_PREDICT1FMV8_C(2,1,1,0)
OD_MC_PREDICT1FMV8_C(2,2,1,1)
OD_MC_PREDICT1FMV8_C(2,4,1,2)
OD_MC_PREDICT1FMV8_C(2,8,1,3)
OD_MC_PREDICT1FMV8_C(2,16,1,4)

OD_MC_PREDICT1FMV8_C(4,1,2,0)
OD_MC_PREDICT1FMV8_C(4,2,2,1)

static void od_mc_predict1fmv8_4x4(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        OD_MC_PREDICT1FMV8HV_4x4
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        OD_MC_PREDICT1FMV8H_4x4
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [hscale]"m"(_mvxf)
      );
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        OD_MC_PREDICT1FMV8V_4x4
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [vscale]"m"(_mvyf)
      );
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        OD_MC_PREDICT1FMV8_4x4
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride)
      );
    }
  }
}

static void od_mc_predict1fmv8_4x8(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_4x16(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_4x4
          "lea (%[src],%[systride],8),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

OD_MC_PREDICT1FMV8_C(8,1,3,0)

static void od_mc_predict1fmv8_8x2(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        OD_MC_PREDICT1FMV8HV_8x2("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        OD_MC_PREDICT1FMV8H_8x2("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [hscale]"m"(_mvxf)
      );
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        OD_MC_PREDICT1FMV8V_8x2("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [vscale]"m"(_mvyf)
      );
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        OD_MC_PREDICT1FMV8_8x2("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride)
      );
    }
  }
}

static void od_mc_predict1fmv8_8x4(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_8x8(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_8x16(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_8x2("movdqu")
          "lea (%[src],%[systride],4),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_16x1(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        OD_MC_PREDICT1FMV8HV_16x1("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        OD_MC_PREDICT1FMV8H_16x1("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [hscale]"m"(_mvxf)
      );
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        OD_MC_PREDICT1FMV8V_16x1("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride),
         [vscale]"m"(_mvyf)
      );
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        OD_MC_PREDICT1FMV8_16x1("movdqu")
        "movdqa %%xmm0,(%[dst])\n\t"
        :[a]"=&r"(a)
        :[dst]"r"(_dst),[src]"r"(_src),[systride]"r"((ptrdiff_t)_systride)
      );
    }
  }
}

static void od_mc_predict1fmv8_16x2(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x20;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_16x4(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x40;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_16x8(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x80;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_16x1("movdqu")
          "lea (%[src],%[systride]),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

static void od_mc_predict1fmv8_16x16(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf){
  unsigned short __attribute__((aligned(16))) mvyf[8];
  ptrdiff_t                                   a;
  ptrdiff_t                                   row;
  if(_mvxf!=0){
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8HV_PROLOG
        :
        :[hscale]"m"(_mvxf),[vscale]"m"(_mvyf),[vscale8]"m"(*mvyf)
      );
      for(row=0;row<0x100;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8HV_16x1("movdqu")
          "lea (%[src],%[systride],2),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride),
           [vscale8]"m"(*mvyf)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8H_PROLOG
        :
        :[hscale]"m"(_mvxf)
      );
      for(row=0;row<0x100;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8H_16x1("movdqu")
          "lea (%[src],%[systride],2),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
  else{
    if(_mvyf!=0){
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8V_PROLOG
        :
        :[vscale]"m"(_mvyf)
      );
      for(row=0;row<0x100;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8V_16x1("movdqu")
          "lea (%[src],%[systride],2),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
    else{
      __asm__ __volatile__(
        OD_MC_PREDICT1FMV8_PROLOG
        ::
      );
      for(row=0;row<0x100;row+=0x10){
        __asm__ __volatile__(
          OD_MC_PREDICT1FMV8_16x1("movdqu")
          "lea (%[src],%[systride],2),%[src]\n\t"
          "movdqa %%xmm0,(%[dst],%[row])\n\t"
          :[src]"+r"(_src),[a]"=&r"(a),[row]"+r"(row)
          :[dst]"r"(_dst),[systride]"r"((ptrdiff_t)_systride)
        );
      }
    }
  }
}

typedef void (*od_mc_predict1fmv8_fixed_func)(unsigned char *_dst,
 const unsigned char *_src,int _systride,unsigned _mvxf,unsigned _mvyf);

void od_mc_predict1fmv8_sse2(unsigned char *_dst,const unsigned char *_src,
 int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz){
  static const od_mc_predict1fmv8_fixed_func VTBL[5][5]={
    {
      od_mc_predict1fmv8_1x1,
      od_mc_predict1fmv8_1x2,
      od_mc_predict1fmv8_1x4,
      od_mc_predict1fmv8_1x8,
      od_mc_predict1fmv8_1x16
    },
    {
      od_mc_predict1fmv8_2x1,
      od_mc_predict1fmv8_2x2,
      od_mc_predict1fmv8_2x4,
      od_mc_predict1fmv8_2x8,
      od_mc_predict1fmv8_2x16
    },
    {
      od_mc_predict1fmv8_4x1,
      od_mc_predict1fmv8_4x2,
      od_mc_predict1fmv8_4x4,
      od_mc_predict1fmv8_4x8,
      od_mc_predict1fmv8_4x16
    },
    {
      od_mc_predict1fmv8_8x1,
      od_mc_predict1fmv8_8x2,
      od_mc_predict1fmv8_8x4,
      od_mc_predict1fmv8_8x8,
      od_mc_predict1fmv8_8x16
    },
    {
      od_mc_predict1fmv8_16x1,
      od_mc_predict1fmv8_16x2,
      od_mc_predict1fmv8_16x4,
      od_mc_predict1fmv8_16x8,
      od_mc_predict1fmv8_16x16
    }
  };
  (*VTBL[_log_xblk_sz][_log_yblk_sz])(_dst,
   _src+(_mvx>>16)+_systride*(_mvy>>16),_systride,_mvx&0xFFFFU,_mvy&0xFFFFU);
#if defined(OD_CHECKASM)
  od_mc_predict1fmv8_check(_dst,_src,_systride,_mvx,_mvy,
   _log_xblk_sz,_log_yblk_sz);
  /*fprintf(stderr,"od_mc_predict1fmv8 %ix%i check finished.\n",
   1<<_log_xblk_sz,1<<_log_yblk_sz);*/
#endif
}



#if defined(OD_CHECKASM)
static void od_mc_blend_full8_check(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  unsigned char  dst[16*16];
  int            xblk_sz;
  int            yblk_sz;
  int            failed;
  int            i;
  int            j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  failed=0;
  od_mc_blend_full8_c(dst,xblk_sz,_src,_log_xblk_sz,_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      if(dst[i+(j<<_log_xblk_sz)]!=(_dst+j*_dystride)[i]){
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         dst[i+(j<<_log_xblk_sz)],(_dst+j*_dystride)[i],i,j);
        failed=1;
      }
    }
  }
  if(failed){
    fprintf(stderr,"od_mc_predict1fmv8 %ix%i check failed.\n",
     (1<<_log_xblk_sz),(1<<_log_yblk_sz));
  }
  OD_ASSERT(!failed);
}
#endif

/*Loads a block of 16 bytes from each of the 4 images into xmm0...xmm3.
  We swap images 2 and 3 here, so that the order more closely follows the
   natural rectilinear indexing, instead of the circular indexing the rest of
   the code uses.*/
#define OD_IM_LOAD16 \
  "#OD_IM_LOAD16\n\t" \
  "mov (%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm0\n\t" \
  /*The "c" prefix here means to leave off the leading $ normally used for \
     immediates (since we want to move from an address, not an immediate). \
    Currently, this is not documented in the gcc manual, and my editor doesn't \
     even syntax highlight it correctly; who knows what versions of gcc \
     support it. \
    It _is_ documented in the internals manual: \
     http://gcc.gnu.org/onlinedocs/gccint/Output-Template.html#Output-Template \
    I found it on this thread: \
     http://gcc.gnu.org/ml/gcc-help/2006-09/msg00301.html*/ \
  "mov %c[pstride](%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm1\n\t" \
  "mov %c[pstride]*3(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm2\n\t" \
  "mov %c[pstride]*2(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm3\n\t" \

/*Unpacks an 8-bit register into two 16-bit registers.
  _rega: The input register, which will contain the low-order output.
  _regb: Will contain the high-order output.
  _zero: A register containing the value zero.*/
#define OD_IM_UNPACK(_rega,_regb,_zero) \
  "#OD_IM_UNPACK\n\t" \
  "movdqa " _rega "," _regb "\n\t" \
  "punpcklbw " _zero "," _rega "\n\t" \
  "punpckhbw " _zero "," _regb "\n\t" \

/*Bilinearly blends 2 pairs of 16-bit registers.
  _reg0a:  The low-order register of the first pair, which will contain the
            output.
  _reg0b:  The high-order register of the first pair, which will contain the
            output.
  _reg1a:  The low-order register of the second pair.
  _reg1b:  The high-order register of the second pair.
  _shift:  The number of bits of precision the blending adds (e.g., the
            precision of the weights).
  _scalea: The weights to apply to the low-order register in the second pair.
           The weights applied to the low-order register in the first pair are
            (1<<_shift)-_scalea.
  _scaleb: The weights to apply to the high-order register in the second pair.
           The weights applied to the high-order register in the first pair are
            (1<<_shift)-_scaleb.*/
#define OD_IM_BLEND(_reg0a,_reg0b,_reg1a,_reg1b,_shift,_scalea,_scaleb) \
  "#OD_IM_BLEND\n\t" \
  "psubw " _reg0a "," _reg1a "\n\t" \
  "psubw " _reg0b "," _reg1b "\n\t" \
  "psllw " _shift "," _reg0a "\n\t" \
  "pmullw " _scalea "," _reg1a "\n\t" \
  "psllw " _shift "," _reg0b "\n\t" \
  "pmullw " _scaleb "," _reg1b "\n\t" \
  "paddw " _reg1a "," _reg0a "\n\t" \
  "paddw " _reg1b "," _reg0b "\n\t" \

/*Rounds, shifts, and packs two 16-bit registers into one 8 bit register.
  _rega:  The low-order register, which will contain the output.
  _regb:  The high-order register.
  _round: The register containing the rounding offset.
  _shift: The immediate that is the amount to shift.*/
#define OD_IM_PACK(_rega,_regb,_round,_shift) \
  "#OD_IM_PACK\n\t" \
  "paddw " _round "," _rega "\n\t" \
  "paddw " _round "," _regb "\n\t" \
  "psrlw " _shift "," _rega "\n\t" \
  "psrlw " _shift "," _regb "\n\t" \
  "packuswb " _regb "," _rega "\n\t" \

/*Note: We lea the address of our blending weights into a register before
   invoking this macro because of PIC.
  We are CAPABLE of referencing the memory directly in all cases, but
    a) There's no way to tell gcc to combine the base register %ebx it uses for
     PIC with our constant offset, index register, and scale values (because
     gcc just concatenates strings; it does no assembly parsing itself),
    b) gcc is the only one who knows if we're using PIC or not,
    c) gcc is the only one who knows the real name of the (static) symbol.
  The resulting performance difference is unmeasurable.*/

/*Blends 4 rows of a 4xN block (N up to 64).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL8_4x4(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BIL4H],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa (%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$2","%%xmm6","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$2","%%xmm6","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "lea %[OD_BIL4V],%[a]\n\t" \
  "psllw $" #_log_yblk_sz "+1,%%xmm7\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row])","0x10(%[a],%[row])") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+2") \
  /*Get it back out to memory. \
    We have to do this 4 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movdqa %%xmm0,%%xmm1\n\t" \
  "lea (%[dst],%[dystride]),%[a]\n\t" \
  "psrldq $4,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst])\n\t" \
  "movd %%xmm1,(%[a])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "psrldq $8,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst],%[dystride],2)\n\t" \
  "movd %%xmm1,(%[a],%[dystride],2)\n\t" \

/*Blends 2 rows of an 8xN block (N up to 32).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL8_8x2(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa (%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$3","%%xmm6","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$3","%%xmm6","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "lea %[OD_BILV],%[a]\n\t" \
  "psllw $" #_log_yblk_sz "+2,%%xmm7\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row],2)","0x10(%[a],%[row],2)") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+3") \
  /*Get it back out to memory. \
    We have to do this 8 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movq %%xmm0,(%[dst])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "movq %%xmm0,(%[dst],%[dystride])\n\t" \

/*Blends 1 row of a 16xN block (N up to 16).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL8_16x1(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the 4 images to blend.*/ \
  OD_IM_LOAD16 \
  /*Unpack and blend the 0 and 1 images.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  "movdqa 0x10(%[a]),%%xmm6\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$4","(%[a])","%%xmm6") \
  /*Unpack and blend the 2 and 3 images.*/ \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm1","%%xmm7") \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm1","$4","(%[a])","%%xmm6") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "lea %[OD_BILV],%[a]\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "movdqa (%[a],%[row]),%%xmm6\n\t" \
  "psllw $" #_log_yblk_sz "+3,%%xmm7\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "%%xmm6","%%xmm6") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+4") \
  /*Get it back out to memory.*/ \
  "movdqa %%xmm0,(%[dst])\n\t" \

#if 0
/*Defines a pure-C implementation with hard-coded loop limits for block sizes
   we don't want to implement manually (e.g., that have fewer than 16 bytes,
   require byte-by-byte unaligned loads, etc.).
  This should let the compiler aggressively unroll loops, etc.
  It can't vectorize it itself because of the difference in operand sizes.*/
#define OD_MC_BLEND_FULL8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full8_##_n##x##_m(unsigned char *_dst,int _dystride, \
 const unsigned char *_src[4]){ \
  int      o; \
  unsigned a; \
  unsigned b; \
  int      i; \
  int      j; \
  o=0; \
  for(j=0;j<(_m);j++){ \
    for(i=0;i<(_n);i++){ \
      a=(_src[0][o+i]<<(_log_xblk_sz))+(_src[1][o+i]-_src[0][o+i])*i; \
      b=(_src[3][o+i]<<(_log_xblk_sz))+(_src[2][o+i]-_src[3][o+i])*i; \
      _dst[i]=(unsigned char)((a<<(_log_yblk_sz))+(b-a)*j+ \
       (1<<(_log_xblk_sz)+(_log_yblk_sz))/2>>(_log_xblk_sz)+(_log_yblk_sz)); \
    } \
    o+=(_m); \
    _dst+=_dystride; \
  } \
} \

#else
/*With -O3 and inter-module optimization, gcc inlines these anyway.
  I'd rather leave the choice to the compiler.*/
#define OD_MC_BLEND_FULL8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full8_##_n##x##_m(unsigned char *_dst,int _dystride, \
 const unsigned char *_src[4]){ \
  od_mc_blend_full8_c(_dst,_dystride,_src,_log_xblk_sz,_log_yblk_sz); \
} \

#endif

OD_MC_BLEND_FULL8_C(1,1,0,0)
OD_MC_BLEND_FULL8_C(1,2,0,1)
OD_MC_BLEND_FULL8_C(1,4,0,2)
OD_MC_BLEND_FULL8_C(1,8,0,3)
OD_MC_BLEND_FULL8_C(1,16,0,4)

OD_MC_BLEND_FULL8_C(2,1,1,0)
OD_MC_BLEND_FULL8_C(2,2,1,1)
OD_MC_BLEND_FULL8_C(2,4,1,2)
OD_MC_BLEND_FULL8_C(2,8,1,3)
OD_MC_BLEND_FULL8_C(2,16,1,4)

OD_MC_BLEND_FULL8_C(4,1,2,0)
OD_MC_BLEND_FULL8_C(4,2,2,1)

static void od_mc_blend_full8_4x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL8_4x4(2)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    /*Note that we pass the constant 0 for [row] here.
      We'll still use it in indexing expression in the asm, but the overhead is
       negligible, and it's easier than writing a special case of
       OD_MC_BLEND_FULL8_4x4 for it.*/
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full8_4x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  /*We use loops like these so gcc can decide to unroll them if it wants (and
     so that it can do the jumps, etc., for us, instead of trying to figure out
     how to put that in a macro).
    row can't count by 1, because (in the 8x2 versions below) we will need to
     scale it by either 16 or 32, and an index register can only be scaled by
     16.
    Therefore we pre-scale it by 16, and do so everywhere (even for the 4x4 and
     8x2 versions) so that things are consistent.*/
  for(row=0;row<0x20;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_4x4(3)
      "lea (%[dst],%[dystride],4),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_4x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x40;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_4x4(4)
      "lea (%[dst],%[dystride],4),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

OD_MC_BLEND_FULL8_C(8,1,3,0)

static void od_mc_blend_full8_8x2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL8_8x2(1)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full8_8x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x20;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_8x2(2)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_8x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x40;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_8x2(3)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_8x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x80;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_8x2(4)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x1(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL8_16x1(0)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full8_16x2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x20;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(1)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x40;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(2)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x80;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(3)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full8_16x16(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x100;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL8_16x1(4)
      "lea (%[dst],%[dystride]),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

typedef void (*od_mc_blend_full8_fixed_func)(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4]);


/*Perform normal bilinear blending.*/
void od_mc_blend_full8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  static const od_mc_blend_full8_fixed_func VTBL[5][5]={
    {
      od_mc_blend_full8_1x1,od_mc_blend_full8_1x2,
      od_mc_blend_full8_1x4,od_mc_blend_full8_1x8,
      od_mc_blend_full8_1x16
    },
    {
      od_mc_blend_full8_2x1,od_mc_blend_full8_2x2,
      od_mc_blend_full8_2x4,od_mc_blend_full8_2x8,
      od_mc_blend_full8_2x16
    },
    {
      od_mc_blend_full8_4x1,od_mc_blend_full8_4x2,
      od_mc_blend_full8_4x4,od_mc_blend_full8_4x8,
      od_mc_blend_full8_4x16
    },
    {
      od_mc_blend_full8_8x1,od_mc_blend_full8_8x2,
      od_mc_blend_full8_8x4,od_mc_blend_full8_8x8,
      od_mc_blend_full8_8x16
    },
    {
      od_mc_blend_full8_16x1,od_mc_blend_full8_16x2,
      od_mc_blend_full8_16x4,od_mc_blend_full8_16x8,
      od_mc_blend_full8_16x16
    }
  };
  (*VTBL[_log_xblk_sz][_log_yblk_sz])(_dst,_dystride,_src);
#if defined(OD_CHECKASM)
  od_mc_blend_full8_check(_dst,_dystride,_src,_log_xblk_sz,_log_yblk_sz);
  /*fprintf(stderr,"od_mc_blend_full8 %ix%i check finished.\n",
   1<<_log_xblk_sz,1<<_log_yblk_sz);*/
#endif
}



#if defined(OD_CHECKASM)
void od_mc_blend_full_split8_check(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz){
  unsigned char dst[16][16];
  int           xblk_sz;
  int           yblk_sz;
  int           failed;
  int           i;
  int           j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  failed=0;
  od_mc_blend_full_split8_c(dst[0],sizeof(dst[0]),_src,_c,_s,
   _log_xblk_sz,_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      if((_dst+j*_dystride)[i]!=dst[j][i]){
        fprintf(stderr,"ASM mismatch: 0x%02X!=0x%02X @ (%2i,%2i)\n",
         (_dst+j*_dystride)[i],dst[j][i],i,j);
        failed=1;
      }
    }
  }
  if(failed){
    fprintf(stderr,"od_mc_blend_full_split8 %ix%i check failed.\n",
     (1<<_log_xblk_sz),(1<<_log_yblk_sz));
  }
  OD_ASSERT(!failed);
}
#endif


/*Loads a block of 16 bytes from each the first 2 images into xmm0...xmm3.
  xmm2 and xmm3 contain duplicate copies of xmm0 and xmm1, or not, depending on
   whether the block edges are split or not.*/
#define OD_IM_LOAD16A \
  "#OD_IM_LOAD16A\n\t" \
  "mov (%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm0\n\t" \
  "mov %c[pstride](%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm1\n\t" \
  "mov %c[pstride]*4(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm2\n\t" \
  "mov %c[pstride]*5(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm3\n\t" \

/*Loads a block of 16 bytes from the third image into xmm2 and xmm1.
  xmm1 contains a duplicate copy of xmm2, or not, depending on whether the
   block edge is split or not.*/
#define OD_IM_LOAD16B \
  "#OD_IM_LOAD16B\n\t" \
  "mov %c[pstride]*3(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm2\n\t" \
  "mov %c[pstride]*7(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm1\n\t" \

/*Loads a block of 16 bytes from the fourth image into xmm3 and xmm1.
  xmm1 contains a duplicate copy of xmm3, or not, depending on whether the
   block edge is split or not.*/
#define OD_IM_LOAD16C \
  "#OD_IM_LOAD16C\n\t" \
  "mov %c[pstride]*2(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm3\n\t" \
  "mov %c[pstride]*6(%[src]),%[a]\n\t" \
  "movdqa (%[a],%[row]),%%xmm1\n\t" \

/*Blends 4 rows of a 4xN block with split edges (N up to 32).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL_SPLIT8_4x4(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the first two images to blend.*/ \
  OD_IM_LOAD16A \
  /*Unpack and merge the 0 image.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  OD_IM_UNPACK("%%xmm2","%%xmm6","%%xmm7") \
  "paddw %%xmm2,%%xmm0\n\t" \
  "paddw %%xmm6,%%xmm4\n\t" \
  /*Unpack and merge the 1 image.*/ \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  "lea %[OD_BIL4H],%[a]\n\t" \
  "paddw %%xmm3,%%xmm1\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  "movdqa (%[a]),%%xmm6\n\t" \
  /*Blend the 0 and 1 images.*/ \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$2","%%xmm6","%%xmm6") \
  /*Load, unpack, and merge the 2 image.*/ \
  OD_IM_LOAD16B \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm1","%%xmm6","%%xmm7") \
  "paddw %%xmm1,%%xmm2\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  /*Load, unpack, and merge the 3 image.*/ \
  OD_IM_LOAD16C \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  /*We run out of registers here, and have to overwrite our zero. \
    Fortunately, the first row of OD_BILV is also all zeros. \
  "lea %[OD_BILV],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm7","(%[a])") \
  "paddw %%xmm1,%%xmm3\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "paddw %%xmm7,%%xmm6\n\t" \
  "pxor %%xmm7,%%xmm7\n\t"*/ \
  /*Alternate version: Saves 1 memory reference, but has a longer dependency \
     chain.*/ \
  "movdqa %%xmm1,%%xmm7\n\t" \
  "punpcklbw %[OD_BILV],%%xmm7\n\t" \
  "paddw %%xmm7,%%xmm3\n\t" \
  "pxor %%xmm7,%%xmm7\n\t" \
  "punpckhbw %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm6\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  /*End alternate version.*/ \
  "lea %[OD_BIL4H],%[a]\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "movdqa (%[a]),%%xmm1\n\t" \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm6","$2","%%xmm1","%%xmm1") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "psllw $" #_log_yblk_sz "+2,%%xmm7\n\t" \
  "lea %[OD_BIL4V],%[a]\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row])","0x10(%[a],%[row])") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+3") \
  /*Get it back out to memory. \
    We have to do this 4 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movdqa %%xmm0,%%xmm1\n\t" \
  "lea (%[dst],%[dystride]),%[a]\n\t" \
  "psrldq $4,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst])\n\t" \
  "movd %%xmm1,(%[a])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "psrldq $8,%%xmm1\n\t" \
  "movd %%xmm0,(%[dst],%[dystride],2)\n\t" \
  "movd %%xmm1,(%[a],%[dystride],2)\n\t" \

/*Blends 2 rows of an 8xN block with split edges (N up to 16).
  %[dst] must be manually advanced to the proper row beforehand because of its
   stride.*/
#define OD_MC_BLEND_FULL_SPLIT8_8x2(_log_yblk_sz) \
  "pxor %%xmm7,%%xmm7\n\t" \
  /*Load the first two images to blend.*/ \
  OD_IM_LOAD16A \
  /*Unpack and merge the 0 image.*/ \
  OD_IM_UNPACK("%%xmm0","%%xmm4","%%xmm7") \
  OD_IM_UNPACK("%%xmm2","%%xmm6","%%xmm7") \
  "paddw %%xmm2,%%xmm0\n\t" \
  "paddw %%xmm6,%%xmm4\n\t" \
  /*Unpack and merge the 1 image.*/ \
  OD_IM_UNPACK("%%xmm1","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  "lea %[OD_BILH],%[a]\n\t" \
  "paddw %%xmm3,%%xmm1\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  "movdqa (%[a]),%%xmm6\n\t" \
  /*Blend the 0 and 1 images.*/ \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm1","%%xmm5","$3","%%xmm6","%%xmm6") \
  /*Load, unpack, and merge the 2 image.*/ \
  OD_IM_LOAD16B \
  OD_IM_UNPACK("%%xmm2","%%xmm5","%%xmm7") \
  OD_IM_UNPACK("%%xmm1","%%xmm6","%%xmm7") \
  "paddw %%xmm1,%%xmm2\n\t" \
  "paddw %%xmm6,%%xmm5\n\t" \
  /*Load, unpack, and merge the 3 image.*/ \
  OD_IM_LOAD16C \
  OD_IM_UNPACK("%%xmm3","%%xmm6","%%xmm7") \
  /*"lea %[OD_BILV],%[a]\n\t" \
  OD_IM_UNPACK("%%xmm1","%%xmm7","(%[a])") \
  "paddw %%xmm1,%%xmm3\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  "paddw %%xmm7,%%xmm6\n\t" \
  "pxor %%xmm7,%%xmm7\n\t"*/ \
  /*Alternate version: Saves 1 memory reference, but has a longer dependency \
     chain.*/ \
  "movdqa %%xmm1,%%xmm7\n\t" \
  "punpcklbw %[OD_BILV],%%xmm7\n\t" \
  "paddw %%xmm7,%%xmm3\n\t" \
  "pxor %%xmm7,%%xmm7\n\t" \
  "punpckhbw %%xmm7,%%xmm1\n\t" \
  "paddw %%xmm1,%%xmm6\n\t" \
  "pcmpeqw %%xmm1,%%xmm1\n\t" \
  /*End alternate version.*/ \
  "lea %[OD_BILH],%[a]\n\t" \
  "psubw %%xmm1,%%xmm7\n\t" \
  "movdqa (%[a]),%%xmm1\n\t" \
  OD_IM_BLEND("%%xmm2","%%xmm5","%%xmm3","%%xmm6","$3","%%xmm1","%%xmm1") \
  /*Blend, shift, and re-pack images 0+1 and 2+3.*/ \
  "psllw $" #_log_yblk_sz "+3,%%xmm7\n\t" \
  "lea %[OD_BILV],%[a]\n\t" \
  OD_IM_BLEND("%%xmm0","%%xmm4","%%xmm2","%%xmm5","$" #_log_yblk_sz, \
   "(%[a],%[row],2)","0x10(%[a],%[row],2)") \
  OD_IM_PACK("%%xmm0","%%xmm4","%%xmm7","$" #_log_yblk_sz "+4") \
  /*Get it back out to memory. \
    We have to do this 8 bytes at a time because the destination will not in \
     general be packed, nor aligned.*/ \
  "movq %%xmm0,(%[dst])\n\t" \
  "psrldq $8,%%xmm0\n\t" \
  "movq %%xmm0,(%[dst],%[dystride])\n\t" \

#if 0
/*Defines a pure-C implementation with hard-coded loop limits for block sizes
   we don't want to implement manually (e.g., that have fewer than 16 bytes,
   require byte-by-byte unaligned loads, etc.).
  This should let the compiler aggressively unroll loops, etc.
  It can't vectorize it itself because of the difference in operand sizes.*/
/*TODO: This approach (using pointer aliasing to allow us to use normal
   bilinear weights) might actually be faster than the pure-C routine we're
   currently using, which adjusts the weights.
  This should be investigated.*/
#define OD_MC_BLEND_FULL_SPLIT8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full_split8_##_n##x##_m(unsigned char *_dst, \
 int _dystride,const unsigned char *_src[8]){ \
  int      o; \
  unsigned a; \
  unsigned b; \
  int      i; \
  int      j; \
  o=0; \
  for(j=0;j<(_m);j++){ \
    for(i=0;i<(_n);i++){ \
      a=(_src[0][o+i]+_src[4+0][o+i]<<(_log_xblk_sz))+ \
       (_src[1][o+i]-_src[0][o+i]+_src[4+1][o+i]-_src[4+0][o+i])*i; \
      b=(_src[3][o+i]+_src[4+3][o+i]<<(_log_xblk_sz))+ \
       (_src[2][o+i]-_src[3][o+i]+_src[4+2][o+i]-_src[4+3][o+i])*i; \
      _dst[i]=(unsigned char)((a<<(_log_yblk_sz))+(b-a)*j+ \
       (1<<(_log_xblk_sz)+(_log_yblk_sz))>>(_log_xblk_sz)+(_log_yblk_sz)+1); \
    } \
    o+=(_m); \
    _dst+=_dystride; \
  } \
} \

#else
/*TODO: This approach (using pointer aliasing to allow us to use normal
   bilinear weights) might actually be faster than the pure-C routine we're
   currently using, which adjusts the weights.
  This should be investigated.*/
static void od_mc_blend_full_split8_bil_c(unsigned char *_dst,
 int _dystride,const unsigned char *_src[8],int _log_xblk_sz,int _log_yblk_sz){
  int      xblk_sz;
  int      yblk_sz;
  int      round;
  int      o;
  unsigned a;
  unsigned b;
  int      i;
  int      j;
  o=0;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  round=1<<(_log_xblk_sz+_log_yblk_sz);
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      a=((_src[0][o+i]+_src[4+0][o+i])<<_log_xblk_sz)+
       (_src[1][o+i]-_src[0][o+i]+_src[4+1][o+i]-_src[4+0][o+i])*i;
      b=((_src[3][o+i]+_src[4+3][o+i])<<_log_xblk_sz)+
       (_src[2][o+i]-_src[3][o+i]+_src[4+2][o+i]-_src[4+3][o+i])*i;
      _dst[i]=(unsigned char)(((a<<_log_yblk_sz)+(b-a)*j+
       round)>>(_log_xblk_sz+_log_yblk_sz+1));
    }
    o+=xblk_sz;
    _dst+=_dystride;
  }
}

/*With -O3 and inter-module optimization, gcc inlines these anyway.
  I'd rather leave the choice to the compiler.*/
#define OD_MC_BLEND_FULL_SPLIT8_C(_n,_m,_log_xblk_sz,_log_yblk_sz) \
static void od_mc_blend_full_split8_##_n##x##_m(unsigned char *_dst, \
 int _dystride,const unsigned char *_src[8]){ \
  od_mc_blend_full_split8_bil_c(_dst,_dystride,_src, \
  _log_xblk_sz,_log_yblk_sz); \
} \

#endif

OD_MC_BLEND_FULL_SPLIT8_C(1,1,0,0)
OD_MC_BLEND_FULL_SPLIT8_C(1,2,0,1)
OD_MC_BLEND_FULL_SPLIT8_C(1,4,0,2)
OD_MC_BLEND_FULL_SPLIT8_C(1,8,0,3)

OD_MC_BLEND_FULL_SPLIT8_C(2,1,1,0)
OD_MC_BLEND_FULL_SPLIT8_C(2,2,1,1)
OD_MC_BLEND_FULL_SPLIT8_C(2,4,1,2)
OD_MC_BLEND_FULL_SPLIT8_C(2,8,1,3)

OD_MC_BLEND_FULL_SPLIT8_C(4,1,2,0)
OD_MC_BLEND_FULL_SPLIT8_C(4,2,2,1)

static void od_mc_blend_full_split8_4x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL_SPLIT8_4x4(2)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full_split8_4x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x20;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_4x4(3)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BIL4H]"m"(*OD_BIL4H),[OD_BIL4V]"m"(*OD_BIL4V),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

OD_MC_BLEND_FULL_SPLIT8_C(8,1,3,0)

static void od_mc_blend_full_split8_8x2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  __asm__ __volatile__(
    OD_MC_BLEND_FULL_SPLIT8_8x2(1)
    :[dst]"+r"(_dst),[a]"=&r"(a)
    :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),[row]"r"((ptrdiff_t)0),
     [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
     [pstride]"i"(sizeof(*_src))
  );
}

static void od_mc_blend_full_split8_8x4(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x20;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_8x2(2)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}

static void od_mc_blend_full_split8_8x8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[8]){
  ptrdiff_t a;
  ptrdiff_t row;
  for(row=0;row<0x40;row+=0x10){
    __asm__ __volatile__(
      OD_MC_BLEND_FULL_SPLIT8_8x2(3)
      "lea (%[dst],%[dystride],2),%[dst]\t\n"
      :[dst]"+r"(_dst),[row]"+r"(row),[a]"=&r"(a)
      :[src]"r"(_src),[dystride]"r"((ptrdiff_t)_dystride),
       [OD_BILH]"m"(*OD_BILH),[OD_BILV]"m"(*OD_BILV),
       [pstride]"i"(sizeof(*_src))
    );
  }
}


typedef void (*od_mc_blend_full_split8_fixed_func)(unsigned char *_dst,
 int _dystride,const unsigned char *_src[8]);



/*Sets up a second set of image pointers based on the given split state to
   properly shift weight from one image to another.*/
static void od_mc_setup_split_ptrs(const unsigned char *_drc[4],
 const unsigned char *_src[4],int _c,int _s){
  int j;
  int k;
  _drc[_c]=_src[_c];
  j=(_c+(_s&1))&3;
  k=(_c+1)&3;
  _drc[k]=_src[j];
  j=(_c+(_s&2)+((_s&2)>>1))&3;
  k=(_c+3)&3;
  _drc[k]=_src[j];
  k=_c^2;
  _drc[k]=_src[k];
}

/*Perform normal bilinear blending.*/
void od_mc_blend_full_split8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz){
  static const od_mc_blend_full_split8_fixed_func VTBL[4][4]={
    {
      od_mc_blend_full_split8_1x1,od_mc_blend_full_split8_1x2,
      od_mc_blend_full_split8_1x4,od_mc_blend_full_split8_1x8,
    },
    {
      od_mc_blend_full_split8_2x1,od_mc_blend_full_split8_2x2,
      od_mc_blend_full_split8_2x4,od_mc_blend_full_split8_2x8,
    },
    {
      od_mc_blend_full_split8_4x1,od_mc_blend_full_split8_4x2,
      od_mc_blend_full_split8_4x4,od_mc_blend_full_split8_4x8,
    },
    {
      od_mc_blend_full_split8_8x1,od_mc_blend_full_split8_8x2,
      od_mc_blend_full_split8_8x4,od_mc_blend_full_split8_8x8,
    }
  };
  /*We pack all the image pointers in one array to save a register.*/
  const unsigned char *drc[8];
  memcpy(drc,_src,sizeof(*drc)*4);
  od_mc_setup_split_ptrs(drc+4,drc,_c,_s);
  (*VTBL[_log_xblk_sz][_log_yblk_sz])(_dst,_dystride,drc);
#if defined(OD_CHECKASM)
  od_mc_blend_full_split8_check(_dst,_dystride,_src,_c,_s,
   _log_xblk_sz,_log_yblk_sz);
  /*fprintf(stderr,"od_mc_blend_full_split8 %ix%i check finished.\n",
   1<<_log_xblk_sz,1<<_log_yblk_sz);*/
#endif
}

#endif
