#if !defined(_x86_cpu_H)
# define _x86_cpu_H (1)
#include "../internal.h"

#define OD_CPU_X86_MMX    (1<<0)
#define OD_CPU_X86_MMXEXT (1<<1)
#define OD_CPU_X86_3DNOW  (1<<2)
#define OD_CPU_X86_3DNOW2 (1<<3)
#define OD_CPU_X86_SSE    (1<<4)
#define OD_CPU_X86_SSE2   (1<<5)
/*Prescott New Instructions, also known as SSE3.*/
#define OD_CPU_X86_PNI    (1<<6)

ogg_uint32_t od_cpu_flags_get(void);

#endif
