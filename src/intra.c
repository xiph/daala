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
 unsigned char *bsize, int bstride, int ln) {
  int i;
  od_coeff *t;
  double g1;
  double g2;
  int top;
  int left;
  int n;
  n = 1 << (ln + OD_LOG_BSIZE0);
  top = by > 0 && OD_BLOCK_SIZE4x4(bsize, bstride, bx, by - 1) == ln;
  left = bx > 0 && OD_BLOCK_SIZE4x4(bsize, bstride, bx - 1, by) == ln;
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

/* Trained using a linear regression on subset3. See dump_cfl_scaling4.*/
static ogg_int16_t od_cfl_scaling4[4][4] = {
  { 128, 128, 100, 36},
  { 128, 80, 71, 35},
  { 100, 71, 35, 31},
  { 36, 35, 31, 18},
};

void od_resample_luma_coeffs(od_coeff *l, int lstride,
 const od_coeff *c, int cstride, int xdec, int ydec, int ln, int cln) {
  int n;
  n = 4 << ln;
  if (cln == 0 && (xdec || ydec)) {
    if (xdec) {
      if (ydec) {
        int i;
        od_tf_up_hv_lp(l, lstride, c, cstride, n, n, n);
        for (i = 0; i < 4; i++) {
          int j;
          for (j = 0; j < 4; j++) {
            l[i*lstride + j] = (od_cfl_scaling4[j][i] * l[i*lstride + j]
             + 64) >> 7;
          }
        }
      }
      else od_tf_up_h_lp(l, lstride, c, cstride, n, n);
    }
    else {
      OD_ASSERT(ydec);
      od_tf_up_v_lp(l, lstride, c, cstride, n, n);
    }
  }
  else {
    /*When the transform we code chroma with is smaller than the luma one,
       downsampling just requires copying the upper lower quarter coeffs.*/
    int x;
    int y;
    OD_ASSERT(xdec == ydec);
    for (y = 0; y < n; y++) {
      for (x = 0; x < n; x++) {
        l[y*lstride + x] = c[y*cstride + x];
      }
    }
  }
}
