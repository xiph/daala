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
#include "block_size_dec.h"

void od_block_size_decode(od_ec_dec *dec, unsigned char *bsize, int stride) {
  int i, j;
  int inefficient;
  inefficient = od_ec_dec_uint(dec, 17);
  if (inefficient == 16) {
    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++) {
        bsize[i*stride + j] = 3;
      }
    }
  }
  else {
    /*bsize[0] = (inefficient&1) ? 2 : 0;
    bsize[2] = (inefficient&2) ? 2 : 0;
    bsize[2*stride] = (inefficient&4) ? 2 : 0;
    bsize[2*stride + 2] = (inefficient&8) ? 2 : 0;*/
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        unsigned char *bsize16 = &bsize[2*i*stride + 2*j];
        if ((inefficient & (1 << 2*i << j)) == 0) {
          const ogg_uint16_t *cdf;
          int split;
          int cdf_id;
          cdf_id = od_block_size_cdf16_id(&bsize16[0], stride);
          cdf = od_switch_size8_cdf[cdf_id];
          split = od_ec_decode_cdf_q15(dec, cdf, 16);
          bsize16[0] = (split&1) ? 1 : 0;
          bsize16[1] = (split&2) ? 1 : 0;
          bsize16[stride] = (split&4) ? 1 : 0;
          bsize16[stride + 1] = (split&8) ? 1 : 0;
        }
        else {
          bsize16[0] = bsize16[1] = bsize16[stride] = bsize16[stride + 1] = 2;
        }
      }
    }
  }
}
