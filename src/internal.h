/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#if !defined(_internal_H)
# define _internal_H (1)
# include <limits.h>
# include "../include/daala/codec.h"
# include "odintrin.h"

# include <stddef.h>

# if defined(_MSC_VER)
#  define _USE_MATH_DEFINES
# elif OD_GNUC_PREREQ(4, 2, 0)
#  pragma GCC diagnostic ignored "-Wlong-long"
#  pragma GCC diagnostic ignored "-Woverlength-strings"
# endif

# define OD_VERSION_MAJOR (0)
# define OD_VERSION_MINOR (0)
# define OD_VERSION_SUB   (0)

/* PACKAGE_STRING needs to be defined for deterministic builds */
# if defined(PACKAGE_STRING)
#  define OD_VENDOR_STRING "Xiph's experimental encoder library " \
    PACKAGE_STRING
# else
#  define OD_VENDOR_STRING "Xiph's experimental encoder library " __DATE__
# endif

/*Smallest blocks are 4x4*/
# define OD_LOG_BSIZE0 (2)
/*There are 4 block sizes total (4x4, 8x8, 16x16, 32x32).*/
# define OD_NBSIZES    (4)
/*The log of the maximum length of the side of a block.*/
# define OD_LOG_BSIZE_MAX (OD_LOG_BSIZE0 + OD_NBSIZES - 1)
/*The maximum length of the side of a block.*/
# define OD_BSIZE_MAX     (1 << OD_LOG_BSIZE_MAX)

/*Largest motion compensation partition sizes are 64x64.*/
# define OD_LOG_MVBSIZE_MAX (6)
# define OD_MVBSIZE_MAX (1 << OD_LOG_MVBSIZE_MAX)
/*Smallest motion compensation partition sizes are 8x8.*/
# define OD_LOG_MVBSIZE_MIN (3)
# define OD_MVBSIZE_MIN (1 << OD_LOG_MVBSIZE_MIN)
/*log(2) of the spacing between level-0 MV grid points.*/
# define OD_LOG_MVB_DELTA0 (OD_LOG_MVBSIZE_MAX - OD_LOG_MVBSIZE_MIN)
/*The number of different MV block sizes.*/
# define OD_NMVBSIZES (OD_LOG_MVB_DELTA0 + 1)
/*The spacing between level-0 MV grid points.*/
# define OD_MVB_DELTA0 (1 << OD_LOG_MVB_DELTA0)
/*A mask used to get an MV grid point's offset relative to a level-0 grid
   point.*/
# define OD_MVB_MASK (OD_MVB_DELTA0 - 1)
/*The largest subdivision level.*/
# define OD_MC_LEVEL_MAX (2*OD_LOG_MVB_DELTA0)
/*The total number of subdivision levels.*/
# define OD_MC_NLEVELS (OD_MC_LEVEL_MAX + 1)

# define OD_LIMIT_BSIZE_MIN (OD_BLOCK_4X4)
# define OD_LIMIT_BSIZE_MAX (OD_BLOCK_32X32)
# if OD_LIMIT_BSIZE_MIN > OD_BLOCK_32X32 || OD_LIMIT_BSIZE_MAX > OD_BLOCK_32X32
#  error "block sizes above 32x32 not supported"
# endif
# define OD_DISABLE_FILTER (0)
# define OD_DEBLOCKING (0)
# define OD_DISABLE_CFL (0)

# define OD_ROBUST_STREAM (1)

# define OD_USE_HAAR_WAVELET (0)

# define OD_SIGNAL_Q_SCALING (1)

# define OD_COEFF_SHIFT (4)
/*OD_QUALITY_SHIFT specifies the number of fractional bits in a
   passed in 'quality' parameter.
  For example, an OD_QUALITY_SHIFT of (4) specifies the quality parameter is
   in Q4 format.*/
# define OD_QUALITY_SHIFT (4)

# if defined(OD_ENABLE_ASSERTIONS)
#  if OD_GNUC_PREREQ(2, 5, 0)
__attribute__((noreturn))
#  endif
void od_fatal_impl(const char *_str, const char *_file, int _line);

#  define OD_FATAL(_str) (od_fatal_impl(_str, __FILE__, __LINE__))

#  define OD_ASSERT(_cond) \
  do { \
    if (!(_cond)) { \
      OD_FATAL("assertion failed: " # _cond); \
    } \
  } \
  while (0)

#  define OD_ASSERT2(_cond, _message) \
  do { \
    if (!(_cond)) { \
      OD_FATAL("assertion failed: " # _cond "\n" _message); \
    } \
  } \
  while (0)

#  define OD_ALWAYS_TRUE(_cond) OD_ASSERT(_cond)

# else
#  define OD_ASSERT(_cond)
#  define OD_ASSERT2(_cond, _message)
#  define OD_ALWAYS_TRUE(_cond) ((void)(_cond))
# endif

/*Like OD_ASSERT, but return an error code instead of terminating.*/
# define OD_RETURN_CHECK(_cond, _err) \
  do { \
    if (!(_cond)) { \
      return (_err); \
    } \
  } while(0)

# define OD_MEM_SIZE_MAX (~(size_t)0 >> 1)
# define OD_MEM_DIFF_MAX ((ptrdiff_t)OD_MEM_SIZE_MAX)

# if OD_GNUC_PREREQ(3, 0, 0)
/*Another alternative is
    (__builtin_constant_p(_x)?!!(_x):__builtin_expect(!!(_x),1))
   but that evaluates _x multiple times, which may be bad.*/
#  define OD_LIKELY(_x) (__builtin_expect(!!(_x),1))
#  define OD_UNLIKELY(_x) (__builtin_expect(!!(_x),0))
# else
#  define OD_LIKELY(_x)   (!!(_x))
#  define OD_UNLIKELY(_x) (!!(_x))
# endif

/*Currently this structure is only in Tremor, and is read-only.*/
typedef struct oggbyte_buffer oggbyte_buffer;

/*Simple libogg1-style buffer.*/
struct oggbyte_buffer{
  unsigned char *buf;
  unsigned char *ptr;
  ptrdiff_t      storage;
};

/*Encoding functions.*/
void oggbyte_writeinit(oggbyte_buffer *_b);
void oggbyte_writetrunc(oggbyte_buffer *_b, ptrdiff_t _bytes);
void oggbyte_write1(oggbyte_buffer *_b, unsigned _value);
void oggbyte_write4(oggbyte_buffer *_b, uint32_t _value);
void oggbyte_writecopy(oggbyte_buffer *_b, const void *_source,
 ptrdiff_t _bytes);
void oggbyte_writeclear(oggbyte_buffer *_b);
/*Decoding functions.*/
void oggbyte_readinit(oggbyte_buffer *_b, unsigned char *_buf,
 ptrdiff_t _bytes);
int oggbyte_look1(oggbyte_buffer *_b);
int oggbyte_look4(oggbyte_buffer *_b, uint32_t *_val);
void oggbyte_adv1(oggbyte_buffer *_b);
void oggbyte_adv4(oggbyte_buffer *_b);
int oggbyte_read1(oggbyte_buffer *_b);
int oggbyte_read4(oggbyte_buffer *_b, uint32_t *_val);
int oggbyte_readcopy(oggbyte_buffer *_b, void *_dest, uint32_t _bytes);

/*Shared functions.*/
void oggbyte_reset(oggbyte_buffer *_b);
ptrdiff_t oggbyte_bytes(oggbyte_buffer *_b);
unsigned char *oggbyte_get_buffer(oggbyte_buffer *_b);
ptrdiff_t oggbyte_bytes_left(oggbyte_buffer *_b);

int od_ilog(uint32_t _v);
void **od_malloc_2d(size_t _height, size_t _width, size_t _sz);
void **od_calloc_2d(size_t _height, size_t _width, size_t _sz);
void od_free_2d(void *_ptr);

# define OD_DIVU_DMAX (64)

extern uint32_t OD_DIVU_SMALL_CONSTS[OD_DIVU_DMAX][2];

/*Translate unsigned division by small divisors into multiplications.*/
# define OD_DIVU_SMALL(_x, _d) \
  ((uint32_t)((OD_DIVU_SMALL_CONSTS[(_d)-1][0]* \
  (unsigned long long)(_x)+OD_DIVU_SMALL_CONSTS[(_d)-1][1])>>32)>> \
  (OD_ILOG(_d)-1))

#endif
