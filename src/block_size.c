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

#include "block_size.h"
#include <string.h>
#include "internal.h"
#include "entenc.h"


int od_block_size_prob32(const unsigned char *bsize, int stride) {
  int i;
  int sum32;
  sum32 = 0;
  for (i = -1; i < 4; i++) sum32 += bsize[-stride + i];
  for (i = 0; i < 4; i++) sum32 += bsize[stride*i - 1];
  return sum32;
}

int od_block_size_prob16(const unsigned char *bsize, int stride) {
  return 16*bsize[-stride + 2] + 4*bsize[2*stride - 1] + bsize[stride + 1];
}

int od_block_size_cdf8_id(const unsigned char *bsize, int stride) {
  int upleft;
  int up;
  int left;
  int id;
  upleft = bsize[-stride - 1];
  /* For up and left, count only the combinations that are possible */
  up = bsize[-stride] < 2 ? bsize[-stride]*2 + bsize[-stride + 1] :
   bsize[-stride] + 2;
  left = bsize[-1] < 2 ? bsize[-1]*2 + bsize[stride - 1] : bsize[-1] + 2;
  id = up*24 + left*4 + upleft;
  return id;
}
