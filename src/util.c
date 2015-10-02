/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdlib.h>

#include "internal.h"
#include "util.h"
#include "x86/x86int.h"

#define OD_COPY_C(_blk_sz) \
void od_copy_##_blk_sz##x##_blk_sz##_c(unsigned char *_dst, int _dstride, \
 const unsigned char *_src, int _sstride) { \
  int y; \
  for (y = 0; y < _blk_sz; y++) { \
    memcpy(_dst, _src, _blk_sz); \
    _dst += _dstride; \
    _src += _sstride; \
  } \
} \

OD_COPY_C(2)
OD_COPY_C(4)
OD_COPY_C(8)
OD_COPY_C(16)
OD_COPY_C(32)
OD_COPY_C(64)

/*Block copy functions. Copying of overlapping regions has undefined
   behavior. Only 16x16, 32x32, 64x64 show a SIMD improvement.*/

void od_copy_16x16(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride) {
#if defined(OD_SSE2_INTRINSICS)
  od_copy_16x16_sse2(_dst, _dstride, _src, _sstride);
#elif
  od_copy_16x16_c(_dst, _dstride, _src, _sstride);
#endif
}

void od_copy_32x32(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride) {
#if defined(OD_SSE2_INTRINSICS)
  od_copy_32x32_sse2(_dst, _dstride, _src, _sstride);
#elif
  od_copy_32x32_c(_dst, _dstride, _src, _sstride);
#endif
}

void od_copy_64x64(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride) {
#if defined(OD_SSE2_INTRINSICS)
  od_copy_64x64_sse2(_dst, _dstride, _src, _sstride);
#elif
  od_copy_64x64_c(_dst, _dstride, _src, _sstride);
#endif
}

typedef void (*od_copy_fixed_func)(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride);

void od_copy_log_nxn(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride, int _log_n) {
  static const od_copy_fixed_func
   VTBL[7] = {
    NULL,
    od_copy_2x2_c,
    od_copy_4x4_c,
    od_copy_8x8_c,
    od_copy_16x16,
    od_copy_32x32,
    od_copy_64x64
  };
  OD_ASSERT(_log_n > 0 || _log_n <= 6);
  (*VTBL[_log_n])(_dst, _dstride, _src, _sstride);
}

void od_copy_log_nxm(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride, int _log_n, int _log_m) {
  int j;
  if (_log_n == _log_m) {
    od_copy_log_nxn(_dst, _dstride, _src, _sstride, _log_n);
    return;
  }
  for (j = 0; j < 1 << _log_m; j++) {
    OD_COPY(_dst, _src, 1 << _log_n);
    _dst += _dstride;
    _src += _sstride;
  }
}
