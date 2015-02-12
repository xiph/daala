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
#include "armint.h"
#if defined(OD_ARMASM)

#if defined(_MSC_VER)
/*For GetExceptionCode() and EXCEPTION_ILLEGAL_INSTRUCTION.*/
# define WIN32_LEAN_AND_MEAN
# define WIN32_EXTRA_LEAN
# include <windows.h>

static OD_INLINE ogg_uint32_t od_cpu_flags_get(void){
  ogg_uint32_t flags;
  flags=0;
  /* MSVC has no OD_INLINE __asm support for ARM, but it does let you __emit
   * instructions via their assembled hex code.
   * All of these instructions should be essentially nops. */
# if defined(OD_ARM_MAY_HAVE_EDSP)
  __try{
    /*PLD [r13]*/
    __emit(0xF5DDF000);
    flags|=OD_CPU_ARM_EDSP;
  }
  __except(GetExceptionCode()==EXCEPTION_ILLEGAL_INSTRUCTION){
    /*Ignore exception.*/
  }
#  if defined(OD_ARM_MAY_HAVE_MEDIA)
  __try{
    /*SHADD8 r3,r3,r3*/
    __emit(0xE6333F93);
    flags|=OD_CPU_ARM_MEDIA;
  }
  __except(GetExceptionCode()==EXCEPTION_ILLEGAL_INSTRUCTION){
    /*Ignore exception.*/
  }
#   if defined(OD_ARM_MAY_HAVE_NEON)
  __try{
    /*VORR q0,q0,q0*/
    __emit(0xF2200150);
    flags|=OD_CPU_ARM_NEON;
  }
  __except(GetExceptionCode()==EXCEPTION_ILLEGAL_INSTRUCTION){
    /*Ignore exception.*/
  }
#   endif
#  endif
# endif
  return flags;
}

#elif defined(__linux__)
/* Linux based */
ogg_uint32_t od_cpu_flags_get(void)
{
  ogg_uint32_t flags = 0;
  FILE *cpuinfo;

  /* Reading /proc/self/auxv would be easier, but that doesn't work reliably on
   * Android */
  cpuinfo = fopen("/proc/cpuinfo", "r");

  if(cpuinfo != NULL)
  {
    /* 512 should be enough for anybody (it's even enough for all the flags that
     * x86 has accumulated... so far). */
    char buf[512];

    while(fgets(buf, 512, cpuinfo) != NULL)
    {
# if defined(OD_ARM_MAY_HAVE_EDSP) || defined(OD_ARM_MAY_HAVE_NEON)
      /* Search for edsp and neon flag */
      if(memcmp(buf, "Features", 8) == 0)
      {
        char *p;
#  if defined(OD_ARM_MAY_HAVE_EDSP)
        p = strstr(buf, " edsp");
        if(p != NULL && (p[5] == ' ' || p[5] == '\n'))
          flags |= OD_CPU_ARM_EDSP;
#  endif

#  if defined(OD_ARM_MAY_HAVE_NEON)
        p = strstr(buf, " neon");
        if(p != NULL && (p[5] == ' ' || p[5] == '\n'))
          flags |= OD_CPU_ARM_NEON;
#  endif
      }
# endif

# if defined(OD_ARM_MAY_HAVE_MEDIA)
      /* Search for media capabilities (>= ARMv6) */
      if(memcmp(buf, "CPU architecture:", 17) == 0)
      {
        int version;
        version = atoi(buf+17);

        if(version >= 6)
          flags |= OD_CPU_ARM_MEDIA;
      }
# endif
    }

    fclose(cpuinfo);
  }
  return flags;
}
#else
/* The feature registers which can tell us what the processor supports are
 * accessible in priveleged modes only, so we can't have a general user-space
 * detection method like on x86.*/
# error "Configured to use ARM asm but no CPU detection method available for " \
   "your platform.  Reconfigure with --disable-asm (or send patches)."
#endif


#endif
