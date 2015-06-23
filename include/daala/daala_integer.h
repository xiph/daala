/*
 *  Original copyright (c) 2010 The WebM project authors. All Rights Reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  * Neither the name of Google, nor the WebM Project, nor the names
 *    of its contributors may be used to endorse or promote products
 *    derived from this software without specific prior written
 *    permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


#ifndef DAALA_DAALA_INTEGER_H_
#define DAALA_DAALA_INTEGER_H_

/* get ptrdiff_t, size_t, wchar_t, NULL */
#include <stddef.h>

#if defined(_MSC_VER)
#define DAALA_FORCE_INLINE __forceinline
#define DAALA_INLINE __inline
#else
#define DAALA_FORCE_INLINE __inline__ __attribute__(always_inline)
/* TODO(jbb): Allow a way to force inline off for older compilers. */
#define DAALA_INLINE inline
#endif

#if (defined(_MSC_VER) && (_MSC_VER < 1600)) || defined(DAALA_EMULATE_INTTYPES)
typedef signed char  int8_t;
typedef signed short int16_t;
typedef signed int   int32_t;

typedef unsigned char  uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int   uint32_t;

#if (defined(_MSC_VER) && (_MSC_VER < 1600))
typedef signed __int64   int64_t;
typedef unsigned __int64 uint64_t;
#define INT64_MAX _I64_MAX
#define INT32_MAX _I32_MAX
#define INT32_MIN _I32_MIN
#define INT16_MAX _I16_MAX
#define INT16_MIN _I16_MIN
#endif

#ifndef _UINTPTR_T_DEFINED
typedef size_t uintptr_t;
#endif

#else

/* Most platforms have the C99 standard integer types. */

#if defined(__cplusplus)
# if !defined(__STDC_FORMAT_MACROS)
#  define __STDC_FORMAT_MACROS
# endif
# if !defined(__STDC_LIMIT_MACROS)
#  define __STDC_LIMIT_MACROS
# endif
/* __cplusplus */
#endif

#include <stdint.h>

#endif

/* VS2010 defines stdint.h, but not inttypes.h */
#if defined(_MSC_VER) && _MSC_VER < 1800
#define PRId64 "I64d"
#else
#include <inttypes.h>
#endif

/* DAALA_DAALA_INTEGER_H_ */
#endif
