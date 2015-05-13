/*Daala video codec
Copyright (c) 2014 Daala project contributors.  All rights reserved.

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

#include <stdlib.h>

#include "entenc.h"
#include "generic_code.h"
#include "laplace_code.h"
#include "intra_paint.h"
#include "state.h"

void od_paint_dering(od_adapt_ctx *adapt, od_ec_enc *enc, unsigned char *paint, const unsigned char *img,
 int w, int h, int stride, const unsigned char *dec8, int bstride,
 unsigned char *mode, int mstride, int *edge_sum, int *edge_sum2, int *edge_count, int q) {
  int i, j;
  for(i = 0; i < 32*w; i++) paint[-stride + i] = paint[i];
  for(i = 0; i < 32*h; i++) paint[stride*i - 1] = paint[stride*i];
  paint[-stride - 1] = paint[0];
  /* First pass on the image: collect the stats. */
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      /* Computes the direction for each block and accumulates sums for each edge. */
      od_intra_paint_analysis(paint, stride, dec8, bstride, mode,
       mstride, edge_sum, edge_sum2, edge_count, j, i, 3);
    }
  }
#if 0
  for(i = 0; i < h*8; i++) {
    for(j = 0; j < w*8; j++) {
      printf("%d ", mode[i*mstride+j]);
      /*printf("%d ", mode[i*mstride+j]<<3>>dec8[(i>>1)*bstride + (j>>1)]);*/
    }
    printf("\n");
  }
#endif
  /* Second pass on the image: do the painting and compute the gains. */
  for(i = 0; i < h; i++) {
    for(j = 0; j < w; j++) {
      int idx;
      int gi;
      int k;
      int m;
      int dist;
      int best_dist;
      int best_gain;
      /* Computes both the painted image and the Wiener filter gains for each
         block, by first computing edges, then painting. */
      od_paint_compute_edge_mask(adapt, enc, paint_out, paint, paint_mask, stride, dec8, bstride, mode,
        mstride, edge_sum, edge_sum2, edge_count, q, j, i, 3);
      best_gain = 0;
      dist = 0;
      /* Now we find the optimal deringing strength. We start by computing the
         distortion of the coded image without deringing. */
      for (k = 0; k < 32; k++) {
        for (m = 0; m < 32; m++) {
          idx = (32*i+k)*stride + 32*j + m;
          dist += (img[idx] - paint[idx])*(int)(img[idx] - paint[idx]);
        }
      }
      best_dist = dist;
      /* We compute the distortion for 4 different deringing strengths. */
      for (gi = 1; gi <= 4; gi++) {
        dist = 0;
        for (k = 0; k < 32; k++) {
          for (m = 0; m < 32; m++) {
            int x;
            int y;
            idx = (32*i+k)*stride + 32*j + m;
            x = img[idx];
            y = OD_CLAMPI(0, paint[idx] + ((OD_MINI((int)paint_mask[idx]<<gi>>1, 255)*(paint_out[idx] - paint[idx]) + 128) >> 8), 255);
            dist += (x - y)*(int)(x - y);
          }
        }
        if (dist < best_dist) {
          best_dist = dist;
          best_gain = gi;
          for (k = 0; k < 32; k++) {
            for (m = 0; m < 32; m++) {
              int y;
              idx = (32*i+k)*stride + 32*j + m;
              y = OD_CLAMPI(0, paint[idx] + ((OD_MINI((int)paint_mask[idx]<<gi>>1, 255)*(paint_out[idx] - paint[idx]) + 128) >> 8), 255);
              paint[idx] = y;
            }
          }
        }
      }
      od_encode_cdf_adapt(enc, best_gain, gain_cdf, 9, 128);
      /*printf("%d ", best_gain);*/
    }
    /*printf("\n");*/
  }
}
