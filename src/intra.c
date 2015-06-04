/*Daala video codec
Copyright (c) 2012-2013 Daala project contributors.  All rights reserved.

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
#include <limits.h>
#include "block_size.h"
#include "filter.h"
#include "intra.h"
#include "tf.h"
#include "state.h"

void od_hv_intra_pred(od_coeff *pred, od_coeff *d, int w, int bx, int by,
 unsigned char *bsize, int bstride, int bs) {
  int i;
  od_coeff *t;
  double g1;
  double g2;
  int top;
  int left;
  int n;
  n = 1 << (bs + OD_LOG_BSIZE0);
  top = by > 0 && OD_BLOCK_SIZE4x4(bsize, bstride, bx, by - 1) == bs;
  left = bx > 0 && OD_BLOCK_SIZE4x4(bsize, bstride, bx - 1, by) == bs;
  t = &d[((by << OD_LOG_BSIZE0))*w + (bx << OD_LOG_BSIZE0)];
  g1 = g2 = 0;
  if (top) for (i = 1; i < 4; i++) g1 += t[-n*w + i]*(double)t[-n*w + i];
  if (left) for (i = 1; i < 4; i++) g2 += t[-n + i*w]*(double)t[-n + i*w];
  if (top) for (i = 4; i < n; i++) pred[i] = t[-n*w + i];
  if (left) for (i = 4; i < n; i++) pred[i*n] = t[-n + i*w];
  if (g1 > g2) {
    if (top) for (i = 1; i < 4; i++) pred[i] = t[-n*w + i];
  }
  else {
    if (left) for (i = 1; i < 4; i++) pred[i*n] = t[-n + i*w];
  }
}

/*Trained using a linear regression on subset3.
  See the dump_cfl_scaling4 branch.*/
static const ogg_int16_t OD_CFL_SCALING4[4][4] = {
  { 128, 128, 100, 36 },
  { 128, 80, 71, 35 },
  { 100, 71, 35, 31 },
  { 36, 35, 31, 18 },
};

void od_resample_luma_coeffs(od_coeff *chroma_pred, int cpstride,
 const od_coeff *decoded_luma, int dlstride, int xdec, int ydec, int bs,
 int chroma_bs) {
  int n;
  n = 4 << bs;
  if (chroma_bs == 0 && (xdec || ydec)) {
    if (xdec) {
      if (ydec) {
        int i;
        od_tf_up_hv_lp(chroma_pred, cpstride, decoded_luma, dlstride, n, n, n);
        for (i = 0; i < 4; i++) {
          int j;
          for (j = 0; j < 4; j++) {
            chroma_pred[i*cpstride + j] =
             (OD_CFL_SCALING4[j][i] * chroma_pred[i*cpstride + j] + 64) >> 7;
          }
        }
      }
      else od_tf_up_h_lp(chroma_pred, cpstride, decoded_luma, dlstride, n, n);
    }
    else {
      OD_ASSERT(ydec);
      od_tf_up_v_lp(chroma_pred, cpstride, decoded_luma, dlstride, n, n);
    }
  }
  else {
    int x;
    int y;
    OD_ASSERT(xdec == ydec);
    /*When the transform we code chroma with is smaller than the luma one,
       downsampling just requires copying the upper-left quarter coeffs.*/
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        chroma_pred[y*cpstride + x] = decoded_luma[y*dlstride + x];
      }
    }
  }
}
