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
