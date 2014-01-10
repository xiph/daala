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

#if !defined(_block_size_h)
# define _block_size_h

# include "odintrin.h"
# include <ogg/ogg.h>
# include "entenc.h"

extern const ogg_uint16_t od_switch_size32_cdf[][3];
extern const ogg_uint16_t od_switch_size16_cdf[][8];
extern const ogg_uint16_t od_switch_size8_cdf[][16];

# define OD_BLOCK_SIZE4x4(bsize, bstride, bx, by)\
   ((bsize)[((by)>>1)*(bstride) + ((bx)>>1)])
# define OD_BLOCK_SIZE8x8(bsize, bstride, bx, by)\
   ((bsize)[(by)*(bstride) + (bx)])

int od_block_size_prob32(const unsigned char *bsize, int stride);

int od_block_size_prob16(const unsigned char *bsize, int stride);

int od_block_size_cdf8_id(const unsigned char *bsize, int stride);

#endif
