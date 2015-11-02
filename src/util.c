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

#define OD_COPY8_C(_blk_sz) \
void od_copy_##_blk_sz##x##_blk_sz##_8_c(unsigned char *_dst, int _dstride, \
 const unsigned char *_src, int _sstride) { \
  int y; \
  for (y = 0; y < _blk_sz; y++) { \
    memcpy(_dst, _src, _blk_sz); \
    _dst += _dstride; \
    _src += _sstride; \
  } \
} \

#define OD_COPY16_C(_blk_sz) \
void od_copy_##_blk_sz##x##_blk_sz##_16_c(unsigned char *_dst, int _dstride, \
 const unsigned char *_src, int _sstride) { \
  int y; \
  for (y = 0; y < _blk_sz; y++) { \
    memcpy(_dst, _src, _blk_sz << 1); \
    _dst += _dstride; \
    _src += _sstride; \
  } \
} \

OD_COPY8_C(2)
OD_COPY8_C(4)
OD_COPY8_C(8)
OD_COPY8_C(16)
OD_COPY8_C(32)
OD_COPY8_C(64)

OD_COPY16_C(2)
OD_COPY16_C(4)
OD_COPY16_C(8)
OD_COPY16_C(16)
OD_COPY16_C(32)
OD_COPY16_C(64)

const od_copy_nxn_func OD_COPY_NXN_8_C[OD_LOG_COPYBSIZE_MAX + 1] = {
  NULL,
  od_copy_2x2_8_c,
  od_copy_4x4_8_c,
  od_copy_8x8_8_c,
  od_copy_16x16_8_c,
  od_copy_32x32_8_c,
  od_copy_64x64_8_c
};

const od_copy_nxn_func OD_COPY_NXN_16_C[OD_LOG_COPYBSIZE_MAX + 1] = {
  NULL,
  od_copy_2x2_16_c,
  od_copy_4x4_16_c,
  od_copy_8x8_16_c,
  od_copy_16x16_16_c,
  od_copy_32x32_16_c,
  od_copy_64x64_16_c
};

/*This will work for 8 or 16 bit, but 16 bit copies must add one to _log_n.*/
void od_copy_nxm(unsigned char *_dst, int _dstride,
 const unsigned char *_src, int _sstride, int _log_n, int _log_m) {
  int j;
  OD_ASSERT(_log_n != _log_m);
  for (j = 0; j < 1 << _log_m; j++) {
    OD_COPY(_dst, _src, 1 << _log_n);
    _dst += _dstride;
    _src += _sstride;
  }
}
