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

#include "dct.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "odintrin.h"
#include "state.h"
#include "intra_paint.h"

/* All modes are numbered clockwise starting from mode 0 oriented 45 degrees
   up-right. For an NxN block, mode N is horizontal, mode 2*N is 45-degrees
   down-right, mode 3*N is vertical, and mode 4*N-1 is just next to mode 0.
   Mode 4*N means DC/gradient.  */

/* This function computes the position of the four points used to interpolate
   pixel (i,j) within an block of size 2^ln. The weights w[] are in Q7. We
   always use two edges when interpolating so we can make a smooth transition
   across edges. */
void pixel_interp(int pi[4], int pj[4], int w[4], int m, int i, int j,
 int ln) {
  int k;
  int n;
  int rev;
  int r;
  int y0;
  int y1;
  int x0;
  int x1;
  int f0;
  int f1;
  int d0;
  int d1;
  int dir;
  n = 1 << ln;
  /* The the top and left edges, we reuse the points as-is. */
  if (i == n - 1 || j == n - 1) {
    pi[0] = i;
    pj[0] = j;
    w[0] = 128;
    for (k = 1; k < 4; k++) pi[k] = pj[k] = w[k] = 0;
    if (m == 3*n && i < n - 1) {
      pi[1] = n - 1;
      pj[1] = j;
      pi[2] = -1;
      pj[2] = j;
      w[0] = 32;
      w[1] = 96*(i + 1)/n;
      w[2] = 96 - w[1];
    }
    if (m == n && j < n - 1) {
      pi[1] = i;
      pj[1] = n - 1;
      pi[2] = i;
      pj[2] = -1;
      w[0] = 32;
      w[1] = 96*(j + 1)/n;
      w[2] = 96 - w[1];
    }
    return;
  }
  i++;
  j++;
  /* DC/Gradient mode, weights are proportional to 1/distance. The 255 alias
     is a temporary fix because mode is sometimes unsigned. */
  if (m == 4*n) {
    pi[0] = -1;
    pj[0] = j - 1;
    pi[1] = n - 1;
    pj[1] = j - 1;
    pi[2] = i - 1;
    pj[2] = -1;
    pi[3] = i - 1;
    pj[3] = n - 1;
    if (0) {
    w[0] = (n - i) << 6 >> ln;
    w[1] = i << 6 >> ln;
    w[2] = (n - j) << 6 >> ln;
    w[3] = j << 6 >> ln;
    }
    else {
      int sum;
      w[0] = 1024/i;
      w[1] = 1024/(n - i);
      w[2] = 1024/j;
      w[3] = 1024/(n - j);
      sum = w[0]+w[1]+w[2]+w[3];
      w[0] = w[0]*128/sum;
      w[1] = w[1]*128/sum;
      w[2] = w[2]*128/sum;
      w[3] = 128-w[0]-w[1]-w[2];
    }
    return;
  }
  if (m > 2*n) {
    int tmp;
    tmp = i;
    i = j;
    j = tmp;
    m = 4*n - m;
    rev = 1;
  }
  else {
    rev = 0;
  }
  dir = n - m;
  r = dir << 7 >> ln;
  y0 = (i << 7) + j*r;
  if (y0 >= 0 && y0 < (n << 7)) {
    pi[0] = y0 >> 7;
    f0 = y0 & 0x7f;
    pi[1] = OD_MINI(n, pi[0] + 1);
    pj[0] = pj[1] = 0;
    d0 = j*sqrt(128*128 + r*r);
  }
  else {
    int r_1;
    r_1 = (1 << 7 << ln) / dir;
    if (dir > 0) {
      x0 = (j << 7) - (n - i)*r_1;
      pi[0] = pi[1] = n;
      d0 = n - i;
    }
    else {
      x0 = (j << 7) + i*r_1;
      pi[0] = pi[1] = 0;
      d0 = i;
    }
    d0 = d0 * sqrt(128*128 + r_1*r_1);
    pj[0] = x0 >> 7;
    f0 = x0 & 0x7f;
    pj[1] = OD_MINI(n, pj[0] + 1);
  }

  y1 = (i << 7) - (n - j)*r;
  if (y1 >= 0 && y1 < (n << 7)) {
    pi[2] = y1 >> 7;
    f1 = y1 & 0x7f;
    pi[3] = OD_MINI(n, pi[2] + 1);
    pj[2] = pj[3] = n;
    d1 = (n - j)*sqrt(128*128 + r*r);
  }
  else {
    int r_1;
    r_1 = (1 << 7 << ln) / dir;
    if (dir > 0) {
      x1 = (j << 7) + i*r_1;
      pi[2] = pi[3] = 0;
      d1 = i;
    }
    else {
      x1 = (j << 7) - (n - i)*r_1;
      pi[2] = pi[3] = n;
      d1 = n - i;
    }
    d1 = d1 * sqrt(128*128 + r_1*r_1);
    pj[2] = x1 >> 7;
    f1 = x1 & 0x7f;
    pj[3] = OD_MINI(n, pj[2] + 1);
  }
  if (1) {
    w[0] = (128-f0)*d1/(d0+d1);
    w[1] = (f0)*d1/(d0+d1);
    w[2] = (128-f1)*d0/(d0+d1);
    w[3] = 128-w[0]-w[1]-w[2];
  }
  else {
    /* Pseudo-hanning blending -- doesn't seem to help. */
    double h = (double)d1/(double)(d0+d1);
    h = .5 - .5*cos(2*M_PI*h);
    w[0] = (128-f0)*h;
    w[1] = (f0)*h;
    w[2] = (128-f1)*(1-h);
    w[3] = 128-w[0]-w[1]-w[2];
  }
  if (rev) {
    for (k = 0; k < 4; k++) {
      int tmp = pi[k];
      pi[k] = pj[k];
      pj[k] = tmp;
    }
  }
  for (k = 0; k < 4; k++) pi[k]--;
  for (k = 0; k < 4; k++) pj[k]--;
}

/* Compute the actual interpolation (reconstruction) for a block. */
void interp_block(unsigned char *img, const unsigned char *edge_accum, int n,
 int stride, int m) {
  int i;
  int j;
  int pi[4];
  int pj[4];
  int w[4];
  int ln;
  ln = 0;
  while (1 << ln < n) ln++;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      int k;
      int sum;
      pixel_interp(pi, pj, w, m, i, j, ln);
      sum = 0;
      for (k = 0; k < 4; k++) sum += edge_accum[pi[k]*stride + pj[k]]*w[k];
      img[i*stride + j] = OD_CLAMPI(0, (sum + 64) >> 7, 255);
    }
  }
}

