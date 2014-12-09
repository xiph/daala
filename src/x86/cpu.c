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

#include "cpu.h"
#include "x86int.h"
#if defined(OD_X86ASM)

# if defined(__amd64__)||defined(__x86_64__)
/*On x86-64, gcc seems to be able to figure out how to save %rbx for us when
 *    compiling with -fPIC.*/
#  define cpuid(_op,_subop,_eax,_ebx,_ecx,_edx) \
    __asm__ __volatile__( \
           "cpuid\n\t" \
           :[eax]"=a"(_eax),[ebx]"=b"(_ebx),[ecx]"=c"(_ecx),[edx]"=d"(_edx) \
           :"a"(_op),"c"(_subop) \
          )
# else
/*On x86-32, not so much.*/
#  define cpuid(_op,_subop,_eax,_ebx,_ecx,_edx) \
    __asm__ __volatile__( \
           "xchgl %%ebx,%[ebx]\n\t" \
           "cpuid\n\t" \
           "xchgl %%ebx,%[ebx]\n\t" \
           :[eax]"=a"(_eax),[ebx]"=r"(_ebx),[ecx]"=c"(_ecx),[edx]"=d"(_edx) \
           :"a"(_op),"c"(_subop) \
          )
# endif

ogg_uint32_t od_cpu_flags_get(void){
  ogg_uint32_t eax;
  ogg_uint32_t ebx;
  ogg_uint32_t ecx;
  ogg_uint32_t edx;
  ogg_uint32_t flags;
#if !defined(__amd64__)&&!defined(__x86_64__)
  /*x86-32: Check to see if we have the cpuid instruction.
    This is done by attempting to flip the ID bit in the rFLAGS register.
    If the bit is not writable, this processor doesn't support cpuid.*/
  __asm__(
    "pushfl\n\t"
    "pushfl\n\t"
    "popl %[eax]\n\t"
    "movl %[eax],%[ebx]\n\t"
    "xorl $0x200000,%[eax]\n\t"
    "pushl %[eax]\n\t"
    "popfl\n\t"
    "pushfl\n\t"
    "popl %[eax]\n\t"
    "popfl\n\t"
    :[eax]"=r"(eax),[ebx]"=r"(ebx)
    :
    :"cc"
  );
  if(eax==ebx)return 0;
  /*x86-64: All CPUs support cpuid, so there's no need to check.*/
#endif
  cpuid(0,0,eax,ebx,ecx,edx);
  /*         l e t n          I e n i          u n e G*/
  if(ecx==0x6C65746E&&edx==0x49656E69&&ebx==0x756E6547){
    /*Intel:*/
    cpuid(1,0,eax,ebx,ecx,edx);
    /*If there isn't even MMX, give up.*/
    if(!(edx&0x00800000))return 0;
    flags=OD_CPU_X86_MMX;
    if(edx&0x02000000)flags|=OD_CPU_X86_MMXEXT|OD_CPU_X86_SSE;
    if(edx&0x04000000)flags|=OD_CPU_X86_SSE2;
    if(ecx&0x00000001)flags|=OD_CPU_X86_PNI;
    if(ecx&0x00080000)flags|=OD_CPU_X86_SSE4_1;
    cpuid(7,0,eax,ebx,ecx,edx);
    if(ebx&0x00000020)flags|=OD_CPU_X86_AVX2;
  }
  /*Also          R E T T          E B S I            D M A
     is found in some engineering samples c. 1994.
    We don't bother to check for that here.*/
  /*              D M A c          i t n e          h t u A*/
  else if((ecx==0x444D4163&&edx==0x69746E65&&ebx==0x68747541)||
   /*      C S N            y b   e          d o e G*/
   (ecx==0x43534E20&&edx==0x79622065&&ebx==0x646F6547)){
    /*AMD:*/
    cpuid(0x80000000,0,eax,ebx,ecx,edx);
    if(eax<=0x80000000){
      /*No extended functions supported.
        Use normal cpuid flags.*/
      cpuid(1,0,eax,ebx,ecx,edx);
      /*If there isn't even MMX, give up.*/
      if(!(edx&0x00800000))return 0;
      flags=OD_CPU_X86_MMX;
      if(edx&0x02000000)flags|=OD_CPU_X86_MMXEXT|OD_CPU_X86_SSE;
    }
    else{
      cpuid(0x80000001,0,eax,ebx,ecx,edx);
      /*If there isn't even MMX, give up.*/
      if(!(edx&0x00800000))return 0;
      flags=OD_CPU_X86_MMX;
      if(edx&0x80000000)flags|=OD_CPU_X86_3DNOW;
      if(edx&0x40000000)flags|=OD_CPU_X86_3DNOW2;
      if(edx&0x00400000)flags|=OD_CPU_X86_MMXEXT;
      /*Also check for SSE.*/
      cpuid(1,0,eax,ebx,ecx,edx);
      if(edx&0x02000000)flags|=OD_CPU_X86_SSE;
    }
    if(edx&0x04000000)flags|=OD_CPU_X86_SSE2;
    if(ecx&0x00000001)flags|=OD_CPU_X86_PNI;
  }
  /*Unhandled processor manufacturers.*/
  /*Transmeta:
                  6 8 x M          T e n i          u n e G*/
  /*VIA:
                  d a e t          s n I x          i r y C*/
  /*NexGen Inc. (acquired by AMD):
                  n e v i          r D n e          G x e N*/
  /*Rise:
                  e s i R          e s i R          e s i R*/
  /*IDT/Centaur (now VIA):
                  s l u a          H r u a          t n e C
    These can actually be configured to return any string, and require a
     special detection method.*/
  /*United Microelectronics Corporation:
                    C M U            C M U            C M U*/
  else flags=0;
  return flags;
}

#endif
