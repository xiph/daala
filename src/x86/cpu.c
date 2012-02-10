/*
    Daala video codec
    Copyright (C) 2006-2010 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


#include "cpu.h"
#include "x86int.h"
#if defined(OD_X86ASM)

# if defined(__amd64__)||defined(__x86_64__)
/*On x86-64, gcc seems to be able to figure out how to save %rbx for us when
 *    compiling with -fPIC.*/
#  define cpuid(_op,_eax,_ebx,_ecx,_edx) \
    __asm__ __volatile__( \
           "cpuid\n\t" \
           :[eax]"=a"(_eax),[ebx]"=b"(_ebx),[ecx]"=c"(_ecx),[edx]"=d"(_edx) \
           :"a"(_op) \
           :"cc" \
          )
# else
/*On x86-32, not so much.*/
#  define cpuid(_op,_eax,_ebx,_ecx,_edx) \
    __asm__ __volatile__( \
           "xchgl %%ebx,%[ebx]\n\t" \
           "cpuid\n\t" \
           "xchgl %%ebx,%[ebx]\n\t" \
           :[eax]"=a"(_eax),[ebx]"=r"(_ebx),[ecx]"=c"(_ecx),[edx]"=d"(_edx) \
           :"a"(_op) \
           :"cc" \
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
  cpuid(0,eax,ebx,ecx,edx);
  /*         l e t n          I e n i          u n e G*/
  if(ecx==0x6C65746E&&edx==0x49656E69&&ebx==0x756E6547){
    /*Intel:*/
    cpuid(1,eax,ebx,ecx,edx);
    /*If there isn't even MMX, give up.*/
    if(!(edx&0x00800000))return 0;
    flags=OD_CPU_X86_MMX;
    if(edx&0x02000000)flags|=OD_CPU_X86_MMXEXT|OD_CPU_X86_SSE;
    if(edx&0x04000000)flags|=OD_CPU_X86_SSE2;
    if(ecx&0x00000001)flags|=OD_CPU_X86_PNI;
  }
  /*Also          R E T T          E B S I            D M A
     is found in some engineering samples c. 1994.
    We don't bother to check for that here.*/
  /*              D M A c          i t n e          h t u A*/
  else if(ecx==0x444D4163&&edx==0x69746E65&&ebx==0x68747541||
   /*      C S N            y b   e          d o e G*/
   ecx==0x43534E20&&edx==0x79622065&&ebx==0x646F6547){
    /*AMD:*/
    cpuid(0x80000000,eax,ebx,ecx,edx);
    if(eax<=0x80000000){
      /*No extended functions supported.
        Use normal cpuid flags.*/
      cpuid(1,eax,ebx,ecx,edx);
      /*If there isn't even MMX, give up.*/
      if(!(edx&0x00800000))return 0;
      flags=OD_CPU_X86_MMX;
      if(edx&0x02000000)flags|=OD_CPU_X86_MMXEXT|OD_CPU_X86_SSE;
    }
    else{
      cpuid(0x80000001,eax,ebx,ecx,edx);
      /*If there isn't even MMX, give up.*/
      if(!(edx&0x00800000))return 0;
      flags=OD_CPU_X86_MMX;
      if(edx&0x80000000)flags|=OD_CPU_X86_3DNOW;
      if(edx&0x40000000)flags|=OD_CPU_X86_3DNOW2;
      if(edx&0x00400000)flags|=OD_CPU_X86_MMXEXT;
      /*Also check for SSE.*/
      cpuid(1,eax,ebx,ecx,edx);
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
