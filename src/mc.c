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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stddef.h>
#include "logging.h"
#include "mc.h"

/*Motion compensation routines shared between the encoder and decoder.*/

/*A table of indices used to set up the rotated versions of each vector
   interpolation formula.*/
static const int MIDXS[][4] = {
  { 0, 1, 2, 3 },
  { 0, 0, 2, 3 },
  { 1, 1, 2, 3 },
  { 0, 0, 3, 3 },
  { 1, 1, 2, 2 },
  { 0, 1, 0, 0 },
  { 2, 1, 2, 2 },
  { 0, 1, 1, 1 },
};

/*Set up the finite differences needed to interpolate a motion vector
   component.
  dmv: Returns the motion vector deltas.
  dmv[0]: The initial value.
  dmv[1]: The initial amount to increment per unit change in i.
  dmv[2]: The amount to increment per unit change in j.
  dmv[3]: The amount to increment _dxdi by per unit change in j.
  mvx: The component value of the 4 motion vectors.
  m: The index of the motion vector to use for each corner in the
      base orientation.
  r: The amount to rotate (clockwise) the formulas by (0...3).
  log_xblk_sz: The log base 2 of the horizontal block dimension.
  log_yblk_sz: The log base 2 of the vertical block dimension.*/
void od_mc_setup_mvc(ogg_int32_t dmv[4], const ogg_int32_t mvs[4],
 const int m[4], int r, int log_xblk_sz, int log_yblk_sz) {
  int c0;
  int c1;
  int c2;
  int c3;
  c0 = (m[(0 - r) & 3] + r) & 3;
  c1 = (m[(1 - r) & 3] + r) & 3;
  c2 = (m[(2 - r) & 3] + r) & 3;
  c3 = (m[(3 - r) & 3] + r) & 3;
  dmv[0] = mvs[c0];
  dmv[1] = (mvs[c1] - dmv[0]) >> log_xblk_sz;
  dmv[2] = (mvs[c3] - dmv[0]) >> log_yblk_sz;
  dmv[3] = (mvs[c0] + mvs[c2] - mvs[c1] - mvs[c3]) >>
   (log_xblk_sz + log_yblk_sz);
  /*Advance the vector to the (0.5, 0.5) position.*/
  dmv[0] += dmv[2] >> 1;
  dmv[1] += dmv[3] >> 1;
  dmv[0] += dmv[1] >> 1;
  dmv[2] += dmv[3] >> 1;
}

/*Form the prediction given by one fixed motion vector.
  dst: The destination buffer (xstride must be 1).
  src: The source buffer (xstride must be 1).
  systride: The byte offset between source pixel rows.
  mvx: The X component of the motion vector.
  mvy: The Y component of the motion vector.
  log_xblk_sz: The log base 2 of the horizontal block dimension.
  log_yblk_sz: The log base 2 of the vertical block dimension.*/
void od_mc_predict1fmv8_c(unsigned char *dst, const unsigned char *src,
 int systride, ogg_int32_t mvx, ogg_int32_t mvy,
 int log_xblk_sz, int log_yblk_sz) {
  ogg_uint32_t mvxf;
  ogg_uint32_t mvyf;
  int xblk_sz;
  int yblk_sz;
  int p00;
  int p01;
  int p10;
  int i;
  int j;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  src += (mvx >> 16) + (mvy >> 16)*systride;
  mvxf = (ogg_uint32_t)(mvx & 0xFFFF);
  mvyf = (ogg_uint32_t)(mvy & 0xFFFF);
  if (mvxf != 0) {
    if (mvyf != 0) {
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i++) {
          ogg_uint32_t a;
          ogg_uint32_t b;
          int p11;
          /*printf("<%16.12f, %16.12f>%s", mvx/(double)0x40000,
           mvy/(double)0x40000, i + 1 < xblk_sz ? "::" : "\n");*/
          p00 = src[i<<1];
          p01 = src[i<<1 | 1];
          p10 = (src + systride)[i<<1];
          p11 = (src + systride)[i<<1 | 1];
          a = (((ogg_uint32_t)p00 << 16) + (p01 - p00)*mvxf) >> 16;
          b = (((ogg_uint32_t)p10 << 16) + (p11 - p10)*mvxf) >> 16;
          dst[j*xblk_sz + i] = (unsigned char)(((a<<16) + (b - a)*mvyf) >> 16);
        }
        src += systride << 1;
      }
    }
    else {
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i++) {
          /*printf("<%16.12f, %16.12f>%s", mvx/(double)0x40000,
           mvy/(double)0x40000, i + 1 < xblk_sz ? "::" : "\n");*/
          p00 = src[i<<1];
          p01 = src[i<<1 | 1];
          dst[j*xblk_sz + i] = (unsigned char)(
           (((ogg_uint32_t)p00 << 16) + (p01 - p00)*mvxf) >> 16);
        }
        src += systride << 1;
      }
    }
  }
  else {
    if (mvyf != 0) {
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i++) {
          /*printf("<%16.12f, %16.12f>%s", mvx/(double)0x40000,
           mvy/(double)0x40000, i + 1 < xblk_sz ? "::" : "\n");*/
          p00 = src[i<<1];
          p10 = (src + systride)[i<<1];
          dst[j*xblk_sz + i] = (unsigned char)(
           (((ogg_uint32_t)p00 << 16) + (p10 - p00)*mvyf) >> 16);
        }
        src += systride << 1;
      }
    }
    else {
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i++) {
          /*printf("<%16.12f, %16.12f>%s", mvx/(double)0x40000,
           mvy/(double)0x40000, i + 1 < xblk_sz ? "::" : "\n");*/
          dst[j*xblk_sz + i] = src[i<<1];
        }
        src += systride << 1;
      }
    }
  }
  /*dst -= xblk_sz*yblk_sz;
  for (j = 0; j < yblk_sz; j++) {
    for (i = 0; i < xblk_sz; i++) printf("%2X ", *(dst + i + j*xblk_sz));
    printf("\n");
  }*/
}

static void od_mc_predict1fmv8(od_state *state, unsigned char *dst,
 const unsigned char *src, int systride, ogg_int32_t mvx, ogg_int32_t mvy,
 int log_xblk_sz, int log_yblk_sz) {
  (*state->opt_vtbl.mc_predict1fmv8)(dst, src, systride, mvx, mvy,
   log_xblk_sz, log_yblk_sz);
}

/*Perform normal bilinear blending.*/
void od_mc_blend_full8_c(unsigned char *dst, int dystride,
 const unsigned char *src[4], int log_xblk_sz, int log_yblk_sz) {
  int log_blk_sz2;
  int xblk_sz;
  int yblk_sz;
  int round;
  int i;
  int j;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  log_blk_sz2 = log_xblk_sz + log_yblk_sz;
  round = 1 << (log_blk_sz2 - 1);
  for (j = 0; j < yblk_sz; j++) {
    for (i = 0; i < xblk_sz; i++) {
      unsigned a;
      unsigned b;
      a = src[0][j*xblk_sz + i];
      b = src[3][j*xblk_sz + i];
      a = (a << log_xblk_sz) + (src[1][j*xblk_sz + i] - a)*i;
      b = (b << log_xblk_sz) + (src[2][j*xblk_sz + i] - b)*i;
      dst[i] = (unsigned char)(((a << log_yblk_sz) + (b - a)*j + round) >>
       log_blk_sz2);
    }
    dst += dystride;
  }
}

/*Perform normal bilinear blending.*/
static void od_mc_blend_full8(od_state *state, unsigned char *dst,
 int dystride, const unsigned char *src[4], int log_xblk_sz, int log_yblk_sz) {
  (*state->opt_vtbl.mc_blend_full8)(dst, dystride, src,
   log_xblk_sz, log_yblk_sz);
}

/* Pulled aut of mcenc so it can be used in decoder as well. */
/* maybe call od_mv */
void od_state_mvs_clear(od_state *state) {
  int vx;
  int vy;
  int nhmvbs;
  int nvmvbs;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  for (vy = 0; vy < nvmvbs; vy++) {
    od_mv_grid_pt *grid;
    grid = state->mv_grid[vy];
    for (vx = 0; vx < nhmvbs; vx++) {
      grid[vx].valid = 0;
      grid[vx].right = 0;
      grid[vx].down = 0;
      grid[vx].mv[0] = 0;
      grid[vx].mv[1] = 0;
    }
  }
}
#if 0
/*Perform multiresolution bilinear blending.*/
static void od_mc_blend_multi8(unsigned char *dst, int dystride,
 const unsigned char *src[4], int log_xblk_sz, int log_yblk_sz) {
  const unsigned char *p;
  unsigned char *dst0;
  ptrdiff_t o;
  ptrdiff_t o0;
  int ll[4];
  int lh;
  int hl;
  int hh;
  int a;
  int b;
  int c;
  int d;
  int log_blk_sz2;
  int xblk_sz;
  int yblk_sz;
  int xblk_sz_2;
  int yblk_sz_2;
  int i;
  int j;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  log_blk_sz2 = log_xblk_sz + log_yblk_sz;
  o0 = 0;
  dst0 = dst;
  /*Perform multiresolution blending.*/
  xblk_sz_2 = xblk_sz >> 1;
  yblk_sz_2 = yblk_sz >> 1;
  for (j = 1; j < yblk_sz_2; j += 2) {
    o = o0;
    dst = dst0;
    /*Upper-left quadrant.*/
    for (i = 1; i < xblk_sz_2; i += 2) {
      p = src[0] + o;
      /*Forward Haar wavelet.*/
      ll[0] = p[0] + p[1];
      lh = p[0] - p[1];
      hl = (p + xblk_sz)[0] + (p + xblk_sz)[1];
      hh = (p + xblk_sz)[0] - (p + xblk_sz)[1];
      c = ll[0] - hl;
      ll[0] += hl;
      hl = c;
      /*No need to finish the transform; we'd just invert it later.*/
      p = src[1] + o;
      ll[1] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[2] + o;
      ll[2] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[3] + o;
      ll[3] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      /*LL blending.*/
      a = (ll[0] << log_xblk_sz) + (ll[1] - ll[0])*i;
      b = (ll[3] << log_xblk_sz) + (ll[2] - ll[3])*i;
      a = (int)((((ogg_int32_t)a << log_yblk_sz) +
       (ogg_int32_t)(b - a)*j) >> log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c = (a - hl + 1) >> 1;
      a = (a + hl + 1) >> 1;
      d = (c - hh + 1) >> 1;
      c = (c + hh + 1) >> 1;
      b = (a - lh + 1) >> 1;
      a = (a + lh + 1) >> 1;
      dst[0] = OD_CLAMP255(a);
      dst[1] = OD_CLAMP255(b);
      (dst + dystride)[0] = OD_CLAMP255(c);
      (dst + dystride)[1] = OD_CLAMP255(d);
      o += 2;
      dst += 2;
    }
    /*Upper-right quadrant.*/
    for (; i < xblk_sz; i += 2) {
      p = src[0] + o;
      ll[0] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[1] + o;
      /*Forward Haar wavelet.*/
      ll[1] = p[0] + p[1];
      lh = p[0] - p[1];
      hl = (p + xblk_sz)[0] + (p + xblk_sz)[1];
      hh = (p + xblk_sz)[0] - (p + xblk_sz)[1];
      c = ll[1] - hl;
      ll[1] += hl;
      hl = c;
      /*No need to finish the transform; we'd just invert it later.*/
      p = src[2] + o;
      ll[2] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[3] + o;
      ll[3] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      /*LL blending.*/
      a = (ll[0] << log_xblk_sz) + (ll[1] - ll[0])*i;
      b = (ll[3] << log_xblk_sz) + (ll[2] - ll[3])*i;
      a = (int)((((ogg_int32_t)a << log_yblk_sz)
       + (ogg_int32_t)(b - a)*j) >> log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c = (a - hl + 1) >> 1;
      a = (a + hl + 1) >> 1;
      d = (c - hh + 1) >> 1;
      c = (c + hh + 1) >> 1;
      b = (a - lh + 1) >> 1;
      a = (a + lh + 1) >> 1;
      dst[0] = OD_CLAMP255(a);
      dst[1] = OD_CLAMP255(b);
      (dst + dystride)[0] = OD_CLAMP255(c);
      (dst + dystride)[1] = OD_CLAMP255(d);
      o += 2;
      dst += 2;
    }
    o0 += xblk_sz << 1;
    dst0 += dystride << 1;
  }
  for (; j < yblk_sz; j += 2) {
    o = o0;
    dst = dst0;
    /*Lower-left quadrant.*/
    for (i = 1; i < xblk_sz_2; i += 2) {
      p = src[0] + o;
      ll[0] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[1] + o;
      ll[1] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[2] + o;
      ll[2] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[3] + o;
      /*Forward Haar wavelet.*/
      ll[3] = p[0] + p[1];
      lh = p[0] - p[1];
      hl = (p + xblk_sz)[0] + (p + xblk_sz)[1];
      hh = (p + xblk_sz)[0] - (p + xblk_sz)[1];
      c = ll[3] - hl;
      ll[3] += hl;
      hl = c;
      /*No need to finish the transform; we'd just invert it later.*/
      /*LL blending.*/
      a = (ll[0] << log_xblk_sz) + (ll[1] - ll[0])*i;
      b = (ll[3] << log_xblk_sz) + (ll[2] - ll[3])*i;
      a = (int)(((ogg_int32_t)a << log_yblk_sz)
       + (ogg_int32_t)(b - a)*j >> log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c = (a - hl + 1) >> 1;
      a = (a + hl + 1) >> 1;
      d = (c - hh + 1) >> 1;
      c = (c + hh + 1) >> 1;
      b = (a - lh + 1) >> 1;
      a = (a + lh + 1) >> 1;
      dst[0] = OD_CLAMP255(a);
      dst[1] = OD_CLAMP255(b);
      (dst + dystride)[0] = OD_CLAMP255(c);
      (dst + dystride)[1] = OD_CLAMP255(d);
      o += 2;
      dst += 2;
    }
    /*Lower-right quadrant.*/
    for (; i < xblk_sz; i += 2) {
      p = src[0] + o;
      ll[0] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[1] + o;
      ll[1] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[2] + o;
      /*Forward Haar wavelet.*/
      ll[2] = p[0] + p[1];
      lh = p[0] - p[1];
      hl = (p + xblk_sz)[0] + (p + xblk_sz)[1];
      hh = (p + xblk_sz)[0] - (p + xblk_sz)[1];
      c = ll[2] - hl;
      ll[2] += hl;
      hl = c;
      /*No need to finish the transform; we'd just invert it later.*/
      p = src[3] + o;
      ll[3] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      /*LL blending.*/
      a = (ll[0] << log_xblk_sz) + (ll[1] - ll[0])*i;
      b = (ll[3] << log_xblk_sz) + (ll[2] - ll[3])*i;
      a = (int)((((ogg_int32_t)a << log_yblk_sz)
       + (ogg_int32_t)(b - a)*j) >> log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c = (a - hl + 1) >> 1;
      a = (a + hl + 1) >> 1;
      d = (c - hh + 1) >> 1;
      c = (c + hh + 1) >> 1;
      b = (a - lh + 1) >> 1;
      a = (a + lh + 1) >> 1;
      dst[0] = OD_CLAMP255(a);
      dst[1] = OD_CLAMP255(b);
      (dst + dystride)[0] = OD_CLAMP255(c);
      (dst + dystride)[1] = OD_CLAMP255(d);
      o += 2;
      dst += 2;
    }
    o0 += xblk_sz << 1;
    dst0 += dystride << 1;
  }
}
#else

/*Perform multiresolution bilinear blending.*/
static void od_mc_blend_multi8(unsigned char *dst, int dystride,
 const unsigned char *src[4], int log_xblk_sz, int log_yblk_sz) {
  unsigned char src_ll[4][8][8];
  int dst_ll[8][8];
  const unsigned char *p;
  ptrdiff_t o;
  int a;
  int b;
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int log_blk_sz2;
  int xblk_sz;
  int yblk_sz;
  int xblk_sz_2;
  int yblk_sz_2;
  int xblk_sz_4;
  int yblk_sz_4;
  int round;
  int i;
  int i2;
  int j;
  int j2;
  int k;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  /*Perform multiresolution blending.*/
  xblk_sz_2 = xblk_sz >> 1;
  yblk_sz_2 = yblk_sz >> 1;
  log_blk_sz2 = log_xblk_sz + log_yblk_sz;
  round = 1 << (log_blk_sz2 - 1);
  /*Compute the low-pass band for each src block.*/
  for (k = 0; k < 4; k++) {
    unsigned lh[4][8];
    p = src[k];
    src_ll[k][0][0] = p[0];
    for (i = 1; i < xblk_sz_2; i++) {
      i2 = i << 1;
      src_ll[k][0][i] =
       (unsigned char)((p[i2 - 1] + 2*p[i2] + p[i2+1] + 2) >> 2);
    }
    p += xblk_sz;
    lh[1][0] = p[0] << 2;
    for (i = 1; i < xblk_sz_2; i++) {
      i2 = i << 1;
      lh[1][i] = p[i2-1] + 2*p[i2] + p[i2+1];
    }
    p += xblk_sz;
    for (j = 1; j < yblk_sz_2; j++) {
      j2 = j << 1 & 3;
      lh[j2][0] = p[0] << 2;
      for (i = 1; i < xblk_sz_2; i++) {
        i2 = i << 1;
        lh[j2][i] = p[i2 - 1] + 2*p[i2] + p[i2 + 1];
      }
      p += xblk_sz;
      lh[j2 + 1][0] = p[0] << 2;
      for (i = 1; i < xblk_sz_2; i++) {
        i2 = i << 1;
        lh[j2 + 1][i] = p[i2 - 1] + 2*p[i2] + p[i2 + 1];
      }
      p += xblk_sz;
      for (i = 0; i < xblk_sz_2; i++) {
        src_ll[k][j][i] = (unsigned char)(
         (lh[(j2 - 1) & 3][i] + 2*lh[j2][i] + lh[j2 + 1][i] + 8) >> 4);
      }
    }
  }
  /*Blend the low-pass bands.*/
  for (j = 0; j < xblk_sz_2; j++) {
    for (i = 0; i < xblk_sz_2; i++) {
      a = (src_ll[0][j][i] << (log_xblk_sz - 1))
       + (src_ll[1][j][i] - src_ll[0][j][i])*i;
      b = (src_ll[3][j][i] << (log_xblk_sz - 1))
       + (src_ll[2][j][i] - src_ll[3][j][i])*i;
      dst_ll[j][i] = (a << (log_yblk_sz - 1)) + (b - a)*j;
    }
  }
  /*Perform the high-pass filtering for each quadrant.*/
  xblk_sz_4 = xblk_sz >> 2;
  yblk_sz_4 = yblk_sz >> 2;
  o = 0;
  for (j = 0; j < yblk_sz_4; j++) {
    /*Upper-left quadrant.*/
    for (i = 0; i < xblk_sz_4; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[0][j][i] << log_blk_sz2;
      f = (src_ll[0][j][i] + src_ll[0][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[0][j][i] + src_ll[0][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[0][j][i] + src_ll[0][j][i + 1]
       + src_ll[0][j + 1][i] + src_ll[0][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       (((src[0] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
      dst[i2 + 1] = OD_CLAMP255(
       (((src[0] + o)[i2 + 1] <<log_blk_sz2) + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(
       (((src[0] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
       log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(
       (((src[0] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
       log_blk_sz2);
    }
    /*Upper-right quadrant.*/
    for (; i < xblk_sz_2 - 1; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[1][j][i] << log_blk_sz2;
      f = (src_ll[1][j][i] + src_ll[1][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[1][j][i] + src_ll[1][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[1][j][i] + src_ll[1][j][i + 1]
       + src_ll[1][j + 1][i] + src_ll[1][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       (((src[1] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
      dst[i2 + 1] = OD_CLAMP255(
       (((src[1] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(
       (((src[1] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
       log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(
       (((src[1] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
       log_blk_sz2);
    }
    /*Upper-right quadrant, last column.*/
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (3*dst_ll[j][i] - dst_ll[j][i - 1]) << 1;
    c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j + 1][i])
     - (dst_ll[j][i - 1] + dst_ll[j + 1][i - 1]);
    e = src_ll[1][j][i] << log_blk_sz2;
    f = (3*src_ll[1][j][i] - src_ll[1][j][i - 1]) << (log_blk_sz2 - 1);
    g = (src_ll[1][j][i] + src_ll[1][j+1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[1][j][i] + src_ll[1][j + 1][i])
     - (src_ll[1][j][i - 1] + src_ll[1][j + 1][i - 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     (((src[1] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
    dst[i2 + 1] = OD_CLAMP255(
     (((src[1] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(
     (((src[1] + o + xblk_sz)[i2]<< log_blk_sz2) + c - g + round) >>
     log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(
     (((src[1] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
     log_blk_sz2);
    o += xblk_sz << 1;
    dst += dystride << 1;
  }
  for (; j < yblk_sz_2 - 1; j++) {
    /*Lower-left quadrant.*/
    for (i = 0; i < xblk_sz_4; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = (dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1]);
      e = src_ll[3][j][i] << log_blk_sz2;
      f = (src_ll[3][j][i] + src_ll[3][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[3][j][i] + src_ll[3][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[3][j][i] + src_ll[3][j][i + 1]
       + src_ll[3][j + 1][i] + src_ll[3][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       (((src[3] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
      dst[i2 + 1] = OD_CLAMP255(
       (((src[3] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(
       (((src[3] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
       log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(
       (((src[3] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
       log_blk_sz2);
    }
    /*Lower-right quadrant.*/
    for (; i < xblk_sz_2 - 1; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[2][j][i] << log_blk_sz2;
      f = (src_ll[2][j][i] + src_ll[2][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[2][j][i] + src_ll[2][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[2][j][i] + src_ll[2][j][i + 1]
       + src_ll[2][j + 1][i] + src_ll[2][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       (((src[2] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
      dst[i2+1] = OD_CLAMP255(
       (((src[2] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(
       (((src[2] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
       log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(
       (((src[2] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
       log_blk_sz2);
    }
    /*Lower-right quadrant, last column.*/
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (3*dst_ll[j][i] - dst_ll[j][i - 1]) << 1;
    c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j + 1][i])
     - (dst_ll[j][i - 1] + dst_ll[j + 1][i - 1]);
    e = src_ll[2][j][i] << log_blk_sz2;
    f = (3*src_ll[2][j][i] - src_ll[2][j][i - 1]) << (log_blk_sz2 - 1);
    g = (src_ll[2][j][i] + src_ll[2][j + 1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[2][j][i] + src_ll[2][j + 1][i])
     - (src_ll[2][j][i - 1] + src_ll[2][j + 1][i - 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     (((src[2] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
    dst[i2 + 1] = OD_CLAMP255(
     (((src[2] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(
     (((src[2] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
     log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(
     (((src[2] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
     log_blk_sz2);
    o += xblk_sz << 1;
    dst += dystride << 1;
  }
  /*Lower-left quadrant, last row.*/
  for (i = 0; i < xblk_sz_4; i++) {
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (dst_ll[j][i] + dst_ll[j][i+1]) << 1;
    c = (3*dst_ll[j][i] - dst_ll[j - 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j][i + 1])
     - (dst_ll[j - 1][i] + dst_ll[j - 1][i + 1]);
    e = src_ll[3][j][i] << log_blk_sz2;
    f = (src_ll[3][j][i] + src_ll[3][j][i + 1]) << (log_blk_sz2 - 1);
    g = (3*src_ll[3][j][i] - src_ll[3][j - 1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[3][j][i] + src_ll[3][j][i + 1])
     - (src_ll[3][j - 1][i] + src_ll[3][j - 1][i + 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     (((src[3] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
    dst[i2 + 1] = OD_CLAMP255(
     (((src[3] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(
     (((src[3] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
     log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(
     (((src[3] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
     log_blk_sz2);
  }
  /*Lower-right quadrant, last row.*/
  for (; i < xblk_sz_2 - 1; i++) {
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
    c = (3*dst_ll[j][i] - dst_ll[j - 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j][i + 1])
     - (dst_ll[j - 1][i] + dst_ll[j - 1][i + 1]);
    e = src_ll[2][j][i] << log_blk_sz2;
    f = (src_ll[2][j][i] + src_ll[2][j][i + 1]) << (log_blk_sz2 - 1);
    g = (3*src_ll[2][j][i] - src_ll[2][j - 1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[2][j][i] + src_ll[2][j][i + 1])
     - (src_ll[2][j - 1][i] + src_ll[2][j - 1][i + 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     (((src[2] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
    dst[i2+1] = OD_CLAMP255(
     (((src[2] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(
     (((src[2] + o + xblk_sz)[i2] << log_blk_sz2) + c - g + round) >>
     log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(
     (((src[2] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
     log_blk_sz2);
  }
  /*Lower-right quadrant, last row and column.*/
  i2 = i << 1;
  a = dst_ll[j][i] << 2;
  b = (3*dst_ll[j][i] - dst_ll[j][i-1]) << 1;
  c = (3*dst_ll[j][i] - dst_ll[j-1][i]) << 1;
  d = 9*dst_ll[j][i] - 3*(dst_ll[j - 1][i] + dst_ll[j][i - 1])
   + dst_ll[j - 1][i - 1];
  e = src_ll[2][j][i] << log_blk_sz2;
  f = (3*src_ll[2][j][i] - src_ll[2][j][i - 1]) << (log_blk_sz2 - 1);
  g = (3*src_ll[2][j][i] - src_ll[2][j - 1][i]) << (log_blk_sz2 - 1);
  h = (9*src_ll[2][j][i] - 3*(src_ll[2][j][i - 1] + src_ll[2][j - 1][i])
   + src_ll[2][j - 1][i - 1]) << (log_blk_sz2 - 2);
  dst[i2] = OD_CLAMP255(
   (((src[2] + o)[i2] << log_blk_sz2) + a - e + round) >> log_blk_sz2);
  dst[i2+1] = OD_CLAMP255(
   (((src[2] + o)[i2 + 1] << log_blk_sz2) + b - f + round) >> log_blk_sz2);
  (dst + dystride)[i2] = OD_CLAMP255(
   (((src[2] + o + xblk_sz)[i2]<<log_blk_sz2) + c - g + round) >> log_blk_sz2);
  (dst + dystride)[i2 + 1] = OD_CLAMP255(
   (((src[2] + o + xblk_sz)[i2 + 1] << log_blk_sz2) + d - h + round) >>
   log_blk_sz2);
}
#endif

static void od_mc_setup_s_split(int s0[4], int dsdi[4], int dsdj[4],
 int ddsdidj[4], int oc, int s, int log_xblk_sz, int log_yblk_sz) {
  int log_blk_sz2;
  int k;
  log_blk_sz2 = log_xblk_sz + log_yblk_sz;
  s0[0] = 2 << log_blk_sz2;
  s0[1] = s0[2] = s0[3] = 0;
  dsdi[0] = -2 << log_xblk_sz;
  dsdi[1] = 2 << log_xblk_sz;
  dsdi[2] = dsdi[3] = 0;
  dsdj[0] = -2 << log_yblk_sz;
  dsdj[1] = dsdj[2] = 0;
  dsdj[3] = 2 << log_yblk_sz;
  ddsdidj[0] = ddsdidj[2] = 2;
  ddsdidj[1] = ddsdidj[3] = -2;
  if (!(s & 1)) {
    k = (oc + 1) & 3;
    s0[k] >>= 1;
    s0[oc] += s0[k];
    dsdi[k] >>= 1;
    dsdi[oc] += dsdi[k];
    dsdj[k] >>= 1;
    dsdj[oc] += dsdj[k];
    ddsdidj[k] >>= 1;
    ddsdidj[oc] += ddsdidj[k];
  }
  if (!(s & 2)) {
    k = (oc + 3) & 3;
    s0[k] >>= 1;
    s0[oc] += s0[k];
    dsdi[k] >>= 1;
    dsdi[oc] += dsdi[k];
    dsdj[k] >>= 1;
    dsdj[oc] += dsdj[k];
    ddsdidj[k] >>= 1;
    ddsdidj[oc] += ddsdidj[k];
  }
  /*Advance the weights to the (0.5, 0.5) position.
    LOOP VECTORIZES.
  for (k = 0; k < 4; k++) {
    s0[k] += dsdj[k] >> 1;
    dsdi[k] += ddsdidj[k] >> 1;
    s0[k] += dsdi[k] >> 1;
    dsdj[k] += ddsdidj[k] >> 1;
  }*/
}

/*Perform normal blending with bilinear weights modified for unsplit edges.*/
void od_mc_blend_full_split8_c(unsigned char *dst, int dystride,
 const unsigned char *src[4], int oc, int s,
 int log_xblk_sz, int log_yblk_sz) {
  int sw[4];
  int s0[4];
  int dsdi[4];
  int dsdj[4];
  int ddsdidj[4];
  int xblk_sz;
  int yblk_sz;
  int log_blk_sz2p1;
  int round;
  int i;
  int j;
  int k;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  /*The block is too small; perform normal blending.*/
  log_blk_sz2p1 = log_xblk_sz + log_yblk_sz + 1;
  round = 1 << (log_blk_sz2p1 - 1);
  od_mc_setup_s_split(s0, dsdi, dsdj, ddsdidj,
   oc, s, log_xblk_sz, log_yblk_sz);
  /*LOOP VECTORIZES.*/
  for (k = 0; k < 4; k++) sw[k] = s0[k];
  for (j = 0; j < yblk_sz; j++) {
    for (i = 0; i < xblk_sz; i++) {
      int a;
      int b;
      int c;
      int d;
      a = src[0][j*xblk_sz + i];
      b = (src[1][j*xblk_sz + i] - a)*sw[1];
      c = (src[2][j*xblk_sz + i] - a)*sw[2];
      d = (src[3][j*xblk_sz + i] - a)*sw[3];
      dst[i] = (unsigned char)(((a << log_blk_sz2p1)
       + b + c + d + round) >> log_blk_sz2p1);
      /*LOOP VECTORIZES.*/
      for (k = 0; k < 4; k++) sw[k] += dsdi[k];
    }
    dst += dystride;
    /*LOOP VECTORIZES.*/
    for (k = 0; k < 4; k++) {
      s0[k] += dsdj[k];
      sw[k] = s0[k];
      dsdi[k] += ddsdidj[k];
    }
  }
}

/*Perform normal blending with bilinear weights modified for unsplit edges.*/
void od_mc_blend_full_split8(od_state *state, unsigned char *dst,
 int dystride, const unsigned char *src[4], int oc, int s,
 int log_xblk_sz, int log_yblk_sz) {
  (*state->opt_vtbl.mc_blend_full_split8)(dst, dystride, src, oc, s,
   log_xblk_sz, log_yblk_sz);
}

#if 0
/*There are other ways to implement multiresolution blending for the modified
   bilinear weights, but using a table lookup to select which predictor to
   draw the high-frequency coefficients from moves all the complexity into the
   tables, and leaves the code dead simple.*/

/*The MV from which to use the high-frequency coefficients for a 2x2 LL band.*/
static const unsigned char OD_MC_SIDXS_22[3][4][4] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0,
      0, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1,
      3, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 2,
      2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 1,
      3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 1,
      0, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1,
      3, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 2,
      3, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1,
      3, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0,
      3, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 1,
      3, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 1,
      2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 1,
      3, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 2x4 LL band.*/
static const unsigned char OD_MC_SIDXS_24[3][4][8] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0, 0, 0,
      0, 0, 2, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1, 1, 1,
      3, 3, 1, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 0, 2, 2,
      2, 2, 2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 3, 1, 1,
      3, 3, 3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 0, 1, 1,
      0, 0, 2, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1, 1, 1,
      3, 3, 2, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 0, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1, 1, 1,
      3, 3, 2, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0, 0, 0,
      3, 3, 2, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 0, 1, 1,
      3, 3, 1, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 0, 1, 1,
      2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 3, 1, 1,
      3, 3, 2, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 2x8 LL band.*/
static const unsigned char OD_MC_SIDXS_28[3][4][16] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 0, 0, 0, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 3, 3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 0, 0, 0, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 2, 2, 2, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0, 0, 0, 0, 0, 0, 0,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 0, 0, 0, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 0, 0, 0, 1, 1, 1, 1,
      2, 2, 2, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 2, 2, 2, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 4x2 LL band.*/
static const unsigned char OD_MC_SIDXS_42[3][4][8] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0,
      0, 0,
      0, 2,
      0, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1,
      1, 1,
      3, 1,
      3, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 2,
      0, 2,
      2, 2,
      2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 1,
      3, 1,
      3, 3,
      3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 1,
      0, 1,
      0, 2,
      0, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1,
      1, 1,
      3, 2,
      3, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 2,
      0, 2,
      3, 2,
      3, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1,
      1, 1,
      3, 2,
      3, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0,
      0, 0,
      3, 2,
      3, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 1,
      0, 1,
      3, 1,
      3, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 1,
      0, 1,
      2, 2,
      2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 1,
      3, 1,
      3, 2,
      3, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 4x4 LL band.*/
static const unsigned char OD_MC_SIDXS_44[3][4][16] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 2, 2,
      0, 0, 2, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1, 1, 1,
      1, 1, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 0, 2, 2,
      0, 0, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 3, 3,
      3, 3, 3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 2, 2,
      0, 0, 2, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1, 1, 1,
      1, 1, 1, 1,
      3, 3, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 0, 2, 2,
      0, 0, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1, 1, 1,
      1, 1, 1, 1,
      3, 3, 2, 2,
      3, 3, 2, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0, 0, 0,
      0, 0, 0, 0,
      3, 3, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 0, 1, 1,
      0, 0, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 0, 1, 1,
      0, 0, 1, 1,
      2, 2, 2, 2,
      2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 2, 2,
      3, 3, 2, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 4x8 LL band.*/
static const unsigned char OD_MC_SIDXS_48[3][4][32] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2,
      0, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1,
      3, 3, 3, 3, 3, 3, 3, 1,
      3, 3, 3, 3, 3, 3, 3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      3, 3, 3, 2, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for an 8x2 LL
  band.*/
static const unsigned char OD_MC_SIDXS_82[3][4][16] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0,
      0, 0,
      0, 0,
      0, 0,
      0, 2,
      0, 2,
      0, 2,
      0, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1,
      1, 1,
      1, 1,
      1, 1,
      3, 1,
      3, 1,
      3, 1,
      3, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 2,
      0, 2,
      0, 2,
      0, 2,
      2, 2,
      2, 2,
      2, 2,
      2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 1,
      3, 1,
      3, 1,
      3, 1,
      3, 3,
      3, 3,
      3, 3,
      3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 1,
      0, 1,
      0, 1,
      0, 1,
      0, 2,
      0, 2,
      0, 2,
      0, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1,
      1, 1,
      1, 1,
      1, 1,
      3, 2,
      3, 2,
      3, 2,
      3, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 2,
      0, 2,
      0, 2,
      0, 2,
      3, 2,
      3, 2,
      3, 2,
      3, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1,
      1, 1,
      1, 1,
      1, 1,
      3, 2,
      3, 2,
      3, 2,
      3, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0,
      0, 0,
      0, 0,
      0, 0,
      3, 2,
      3, 2,
      3, 2,
      3, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 1,
      0, 1,
      0, 1,
      0, 1,
      3, 1,
      3, 1,
      3, 1,
      3, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 1,
      0, 1,
      0, 1,
      0, 1,
      0, 1,
      2, 2,
      2, 2,
      2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 1,
      3, 1,
      3, 1,
      3, 1,
      3, 2,
      3, 2,
      3, 2,
      3, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for an 8x4 LL
  band.*/
static const unsigned char OD_MC_SIDXS_84[3][4][32] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 2,
      0, 0, 0, 2,
      0, 0, 2, 2,
      0, 0, 2, 2,
      0, 2, 2, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      3, 1, 1, 1,
      3, 1, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 3, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 0, 0, 2,
      0, 0, 2, 2,
      0, 0, 2, 2,
      0, 2, 2, 2,
      0, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 1, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 3, 1,
      3, 3, 3, 1,
      3, 3, 3, 3,
      3, 3, 3, 3,
      3, 3, 3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 2, 2,
      0, 0, 2, 2,
      0, 0, 2, 2,
      0, 2, 2, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      3, 1, 1, 1,
      3, 3, 1, 2,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 0, 0, 2,
      0, 0, 2, 2,
      0, 0, 2, 2,
      0, 0, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1, 1, 1,
      1, 1, 1, 1,
      1, 1, 1, 1,
      3, 1, 1, 1,
      3, 3, 1, 2,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 2,
      3, 0, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 3, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 1, 1,
      0, 0, 2, 1,
      0, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2,
      2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 1, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 1, 1,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2,
      3, 3, 2, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for an 8x8 LL
  band.*/
static const unsigned char OD_MC_SIDXS_88[3][4][64] = {
  {
    {
      /*Corner: 0; split: none*/
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 2,
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: none*/
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1
    },
    {
      /*Corner: 2; split: none*/
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2,
      0, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: none*/
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1,
      3, 3, 3, 3, 3, 3, 3, 1,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3,
      3, 3, 3, 3, 3, 3, 3, 3
    }
  },
  {
    {
      /*Corner: 0; split: 1*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 0, 1, 1, 1,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2,
      0, 0, 0, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: 2*/
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 2; split: 3*/
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 0, 0, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      0, 0, 0, 0, 2, 2, 2, 2,
      3, 3, 3, 2, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 0*/
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      3, 3, 1, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    }
  },
  {
    {
      /*Corner: 0; split: 3*/
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 2, 2,
      3, 3, 3, 0, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    },
    {
      /*Corner: 1; split: 0*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1,
      3, 3, 3, 3, 3, 1, 1, 1
    },
    {
      /*Corner: 2; split: 1*/
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 1, 1, 1, 1,
      0, 0, 0, 0, 2, 1, 1, 1,
      0, 0, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 2, 2, 2, 2
    },
    {
      /*Corner: 3; split: 2*/
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 1, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 1, 1, 1, 1,
      3, 3, 3, 3, 3, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2,
      3, 3, 3, 3, 2, 2, 2, 2
    }
  }
};

/*The MV from which to use the high-frequency coefficients, indexed by:
   [log_yblk_sz - 2][log_xblk_sz - 2]
   [!!s3 << 1 | !!s1][oc][j << (log_xblk_sz - 1) | i].*/
static const unsigned char *OD_MC_SIDXS[3][3][3][4] = {
  {
    {
      {
        OD_MC_SIDXS_22[0][0], OD_MC_SIDXS_22[0][1],
        OD_MC_SIDXS_22[0][2], OD_MC_SIDXS_22[0][3]
      },
      {
        OD_MC_SIDXS_22[1][0], OD_MC_SIDXS_22[1][1],
        OD_MC_SIDXS_22[1][2], OD_MC_SIDXS_22[1][3]
      },
      {
        OD_MC_SIDXS_22[2][0], OD_MC_SIDXS_22[2][1],
        OD_MC_SIDXS_22[2][2], OD_MC_SIDXS_22[2][3]
      }
    },
    {
      {
        OD_MC_SIDXS_24[0][0], OD_MC_SIDXS_24[0][1],
        OD_MC_SIDXS_24[0][2], OD_MC_SIDXS_24[0][3]
      },
      {
        OD_MC_SIDXS_24[1][0], OD_MC_SIDXS_24[1][1],
        OD_MC_SIDXS_24[1][2], OD_MC_SIDXS_24[1][3]
      },
      {
        OD_MC_SIDXS_24[2][0], OD_MC_SIDXS_24[2][1],
        OD_MC_SIDXS_24[2][2], OD_MC_SIDXS_24[2][3]
      }
    },
    {
      {
        OD_MC_SIDXS_28[0][0], OD_MC_SIDXS_28[0][1],
        OD_MC_SIDXS_28[0][2], OD_MC_SIDXS_28[0][3]
      },
      {
        OD_MC_SIDXS_28[1][0], OD_MC_SIDXS_28[1][1],
        OD_MC_SIDXS_28[1][2], OD_MC_SIDXS_28[1][3]
      },
      {
        OD_MC_SIDXS_28[2][0], OD_MC_SIDXS_28[2][1],
        OD_MC_SIDXS_28[2][2], OD_MC_SIDXS_28[2][3]
      }
    }
  },
  {
    {
      {
        OD_MC_SIDXS_42[0][0], OD_MC_SIDXS_42[0][1],
        OD_MC_SIDXS_42[0][2], OD_MC_SIDXS_42[0][3]
      },
      {
        OD_MC_SIDXS_42[1][0], OD_MC_SIDXS_42[1][1],
        OD_MC_SIDXS_42[1][2], OD_MC_SIDXS_42[1][3]
      },
      {
        OD_MC_SIDXS_42[2][0], OD_MC_SIDXS_42[2][1],
        OD_MC_SIDXS_42[2][2], OD_MC_SIDXS_42[2][3]
      }
    },
    {
      {
        OD_MC_SIDXS_44[0][0], OD_MC_SIDXS_44[0][1],
        OD_MC_SIDXS_44[0][2], OD_MC_SIDXS_44[0][3]
      },
      {
        OD_MC_SIDXS_44[1][0], OD_MC_SIDXS_44[1][1],
        OD_MC_SIDXS_44[1][2], OD_MC_SIDXS_44[1][3]
      },
      {
        OD_MC_SIDXS_44[2][0], OD_MC_SIDXS_44[2][1],
        OD_MC_SIDXS_44[2][2], OD_MC_SIDXS_44[2][3]
      }
    },
    {
      {
        OD_MC_SIDXS_48[0][0], OD_MC_SIDXS_48[0][1],
        OD_MC_SIDXS_48[0][2], OD_MC_SIDXS_48[0][3]
      },
      {
        OD_MC_SIDXS_48[1][0], OD_MC_SIDXS_48[1][1],
        OD_MC_SIDXS_48[1][2], OD_MC_SIDXS_48[1][3]
      },
      {
        OD_MC_SIDXS_48[2][0], OD_MC_SIDXS_48[2][1],
        OD_MC_SIDXS_48[2][2], OD_MC_SIDXS_48[2][3]
      }
    }
  },
  {
    {
      {
        OD_MC_SIDXS_82[0][0], OD_MC_SIDXS_82[0][1],
        OD_MC_SIDXS_82[0][2], OD_MC_SIDXS_82[0][3]
      },
      {
        OD_MC_SIDXS_82[1][0], OD_MC_SIDXS_82[1][1],
        OD_MC_SIDXS_82[1][2], OD_MC_SIDXS_82[1][3]
      },
      {
        OD_MC_SIDXS_82[2][0], OD_MC_SIDXS_82[2][1],
        OD_MC_SIDXS_82[2][2], OD_MC_SIDXS_82[2][3]
      }
    },
    {
      {
        OD_MC_SIDXS_84[0][0], OD_MC_SIDXS_84[0][1],
        OD_MC_SIDXS_84[0][2], OD_MC_SIDXS_84[0][3]
      },
      {
        OD_MC_SIDXS_84[1][0], OD_MC_SIDXS_84[1][1],
        OD_MC_SIDXS_84[1][2], OD_MC_SIDXS_84[1][3]
      },
      {
        OD_MC_SIDXS_84[2][0], OD_MC_SIDXS_84[2][1],
        OD_MC_SIDXS_84[2][2], OD_MC_SIDXS_84[2][3]
      }
    },
    {
      {
        OD_MC_SIDXS_88[0][0], OD_MC_SIDXS_88[0][1],
        OD_MC_SIDXS_88[0][2], OD_MC_SIDXS_88[0][3]
      },
      {
        OD_MC_SIDXS_88[1][0], OD_MC_SIDXS_88[1][1],
        OD_MC_SIDXS_88[1][2], OD_MC_SIDXS_88[1][3]
      },
      {
        OD_MC_SIDXS_88[2][0], OD_MC_SIDXS_88[2][1],
        OD_MC_SIDXS_88[2][2], OD_MC_SIDXS_88[2][3]
      }
    }
  }
};

/*Perform multiresolution blending with bilinear weights modified for unsplit
   edges.*/
static void od_mc_blend_multi_split8(unsigned char *dst, int dystride,
 const unsigned char *src[4], int oc, int s,
 int log_xblk_sz, int log_yblk_sz) {
  const unsigned char *p;
  const unsigned char *sidx0;
  const unsigned char *sidx;
  unsigned char *dst0;
  ptrdiff_t o;
  ptrdiff_t o0;
  int ll[4];
  int lh;
  int hl;
  int hh;
  int sw[4];
  int s0[4];
  int dsdi[4];
  int dsdj[4];
  int ddsdidj[4];
  int xblk_sz;
  int yblk_sz;
  int log_blk_sz2p1;
  int i;
  int j;
  int k;
  /*Perform multiresolution blending.*/
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  o0 = 0;
  dst0 = dst;
  log_blk_sz2p1 = log_xblk_sz + log_yblk_sz + 1;
  od_mc_setup_s_split(s0, dsdi, dsdj, ddsdidj, oc, s,
   log_xblk_sz - 1, log_yblk_sz - 1);
  sidx0 = OD_MC_SIDXS[log_yblk_sz - 2][log_xblk_sz - 2][s][oc];
  for (k = 0; k < 4; k++) sw[k] = s0[k];
  for (j = 1; j < yblk_sz; j += 2) {
    o = o0;
    dst = dst0;
    sidx = sidx0;
    /*Upper-left quadrant.*/
    for (i = 1; i < xblk_sz; i += 2) {
      int a;
      int b;
      int c;
      int d;
      k = *sidx++;
      p = src[k] + o;
      /*Forward Haar wavelet.*/
      ll[k] = p[0] + p[1];
      lh = p[0] - p[1];
      hl = (p + xblk_sz)[0] + (p + xblk_sz)[1];
      hh = (p + xblk_sz)[0] - (p + xblk_sz)[1];
      c = ll[k] - hl;
      ll[k] += hl;
      hl = c;
      /*No need to finish the transform; we'd just invert it later.*/
      p = src[(k + 1) & 3] + o;
      ll[(k + 1) & 3] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[(k + 2)&3] + o;
      ll[(k + 2) & 3] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      p = src[(k + 3) & 3] + o;
      ll[(k + 3) & 3] = p[0] + p[1] + (p + xblk_sz)[0] + (p + xblk_sz)[1];
      /*LL blending.*/
      a = (int)(((ogg_int32_t)ll[0]*sw[0]
       + (ogg_int32_t)ll[1]*sw[1] + (ogg_int32_t)ll[2]*sw[2]
       + (ogg_int32_t)ll[3]*sw[3]) >> log_blk_sz2p1);
      /*Inverse Haar wavelet.*/
      c = (a - hl + 1) >> 1;
      a = (a + hl + 1) >> 1;
      d = (c - hh + 1) >> 1;
      c = (c + hh + 1) >> 1;
      b = (a - lh + 1) >> 1;
      a = (a + lh + 1) >> 1;
      dst[0] = OD_CLAMP255(a);
      dst[1] = OD_CLAMP255(b);
      (dst + dystride)[0] = OD_CLAMP255(c);
      (dst + dystride)[1] = OD_CLAMP255(d);
      o += 2;
      dst += 2;
      for (k = 0; k < 4; k++) sw[k] += dsdi[k];
    }
    o0 += xblk_sz << 1;
    dst0 += dystride << 1;
    sidx0 += xblk_sz >> 1;
    for (k = 0; k < 4; k++) {
      s0[k] += dsdj[k];
      sw[k] = s0[k];
      dsdi[k] += ddsdidj[k];
    }
  }
}
#else
/*Sets up a second set of image pointers based on the given split state to
   properly shift weight from one image to another.*/
static void od_mc_setup_split_ptrs(const unsigned char *drc[4],
 const unsigned char *src[4], int oc, int s) {
  int j;
  int k;
  drc[oc] = src[oc];
  j = (oc + (s & 1)) & 3;
  k = (oc + 1) & 3;
  drc[k] = src[j];
  j = (oc + (s & 2) + ((s & 2) >> 1)) & 3;
  k = (oc + 3) & 3;
  drc[k] = src[j];
  k = oc ^ 2;
  drc[k] = src[k];
}

/*Perform multiresolution bilinear blending.*/
static void od_mc_blend_multi_split8(unsigned char *dst, int dystride,
 const unsigned char *src[4], int oc, int s,
 int log_xblk_sz, int log_yblk_sz) {
  unsigned char src_ll[4][8][8];
  int dst_ll[8][8];
  const unsigned char *drc[4];
  const unsigned char *p;
  const unsigned char *q;
  ptrdiff_t o;
  int a;
  int b;
  int c;
  int d;
  int e;
  int f;
  int g;
  int h;
  int log_blk_sz2;
  int xblk_sz;
  int yblk_sz;
  int xblk_sz_2;
  int yblk_sz_2;
  int xblk_sz_4;
  int yblk_sz_4;
  int round;
  int i;
  int i2;
  int j;
  int j2;
  int k;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  od_mc_setup_split_ptrs(drc, src, oc, s);
  /*Perform multiresolution blending.*/
  xblk_sz_2 = xblk_sz >> 1;
  yblk_sz_2 = yblk_sz >> 1;
  log_blk_sz2 = log_xblk_sz + log_yblk_sz;
  round = 1 << (log_blk_sz2 - 1);
  /*Compute the low-pass band for each src block.*/
  for (k = 0; k < 4; k++) {
    unsigned lh[4][8];
    p = src[k];
    q = drc[k];
    src_ll[k][0][0] = (p[0] + q[0] + 1) >> 1;
    for (i = 1; i < xblk_sz_2; i++) {
      i2 = i << 1;
      src_ll[k][0][i] = (unsigned char)((p[i2 - 1] + q[i2 - 1]
       + 2*(p[i2] + q[i2]) + p[i2 + 1] + q[i2 + 1] + 4) >> 3);
    }
    p += xblk_sz;
    q += xblk_sz;
    lh[1][0] = (p[0] + q[0]) << 2;
    for (i = 1; i < xblk_sz_2; i++) {
      i2 = i << 1;
      lh[1][i] = p[i2 - 1] + q[i2 - 1]
       + 2*(p[i2] + q[i2]) + p[i2 + 1] + q[i2 + 1];
    }
    p += xblk_sz;
    q += xblk_sz;
    for (j = 1; j < yblk_sz_2; j++) {
      j2 = j << 1 & 3;
      lh[j2][0] = (p[0] + q[0]) << 2;
      for (i = 1; i < xblk_sz_2; i++) {
        i2 = i << 1;
        lh[j2][i] = p[i2-1] + q[i2-1] +
         2*(p[i2] + q[i2]) + p[i2 + 1] + q[i2 + 1];
      }
      p += xblk_sz;
      q += xblk_sz;
      lh[j2 + 1][0] = p[0] << 2;
      for (i = 1; i < xblk_sz_2; i++) {
        i2 = i << 1;
        lh[j2 + 1][i] = p[i2 - 1] + q[i2 - 1]
         + 2*(p[i2] + q[i2]) + p[i2 + 1] + q[i2 + 1];
      }
      p += xblk_sz;
      q += xblk_sz;
      for (i = 0; i < xblk_sz_2; i++) {
        src_ll[k][j][i] = (unsigned char)(
         (lh[(j2 - 1) & 3][i] + 2*lh[j2][i] + lh[j2 + 1][i] + 16) >> 5);
      }
    }
  }
  /*Blend the low-pass bands.*/
  for (j = 0; j < xblk_sz_2; j++) {
    for (i = 0; i < xblk_sz_2; i++) {
      a = (src_ll[0][j][i] << (log_xblk_sz - 1))
       + (src_ll[1][j][i] - src_ll[0][j][i])*i;
      b = (src_ll[3][j][i] << (log_xblk_sz - 1))
       + (src_ll[2][j][i] - src_ll[3][j][i])*i;
      dst_ll[j][i] = (a << (log_yblk_sz - 1)) + (b - a)*j;
    }
  }
  /*Perform the high-pass filtering for each quadrant.*/
  xblk_sz_4 = xblk_sz >> 2;
  yblk_sz_4 = yblk_sz >> 2;
  o = 0;
  for (j = 0; j < yblk_sz_4; j++) {
    /*Upper-left quadrant.*/
    for (i = 0; i < xblk_sz_4; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] +dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] +dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[0][j][i] << log_blk_sz2;
      f = (src_ll[0][j][i] + src_ll[0][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[0][j][i] + src_ll[0][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[0][j][i] + src_ll[0][j][i + 1]
       + src_ll[0][j + 1][i] + src_ll[0][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       ((((src[0] + o)[i2] + (drc[0] + o)[i2]) << (log_blk_sz2 - 1))
       + a - e + round) >> log_blk_sz2);
      dst[i2+1] = OD_CLAMP255(
       ((((src[0] + o)[i2 + 1] + (drc[0] + o)[i2 + 1]) << (log_blk_sz2 - 1))
       + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(((((src[0] + o + xblk_sz)[i2]
       + (drc[0] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
       + c - g + round) >> log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[0] + o + xblk_sz)[i2 + 1]
       + (drc[0] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
       + d - h + round) >> log_blk_sz2);
    }
    /*Upper-right quadrant.*/
    for (; i < xblk_sz_2 - 1; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[1][j][i] << log_blk_sz2;
      f = (src_ll[1][j][i] + src_ll[1][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[1][j][i] + src_ll[1][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[1][j][i] + src_ll[1][j][i + 1]
       + src_ll[1][j + 1][i] + src_ll[1][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       ((((src[1] + o)[i2] + (drc[1] + o)[i2]) << (log_blk_sz2 - 1))
       + a - e + round) >> log_blk_sz2);
      dst[i2 + 1] = OD_CLAMP255(
       ((((src[1] + o)[i2 + 1] + (drc[1] + o)[i2 + 1]) << (log_blk_sz2 - 1))
       + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(((((src[1] + o + xblk_sz)[i2]
       + (drc[1] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
       + c - g + round) >> log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[1] + o + xblk_sz)[i2 + 1]
       + (drc[1] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
       + d - h + round) >> log_blk_sz2);
    }
    /*Upper-right quadrant, last column.*/
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (3*dst_ll[j][i] - dst_ll[j][i - 1]) << 1;
    c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j + 1][i])
     - (dst_ll[j][i - 1] + dst_ll[j + 1][i - 1]);
    e = src_ll[1][j][i] << log_blk_sz2;
    f = (3*src_ll[1][j][i] - src_ll[1][j][i - 1]) << (log_blk_sz2 - 1);
    g = (src_ll[1][j][i] + src_ll[1][j + 1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[1][j][i] + src_ll[1][j + 1][i])
     - (src_ll[1][j][i - 1] + src_ll[1][j + 1][i - 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     ((((src[1] + o)[i2] + (drc[1] + o)[i2]) << (log_blk_sz2 - 1))
     + a - e + round) >> log_blk_sz2);
    dst[i2 + 1] = OD_CLAMP255(
     ((((src[1] + o)[i2 + 1] + (drc[1] + o)[i2 + 1]) << (log_blk_sz2 - 1))
     + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(((((src[1] + o + xblk_sz)[i2]
     + (drc[1] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
     + c - g + round) >> log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[1] + o + xblk_sz)[i2 + 1]
     + (drc[1] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
     + d - h + round) >> log_blk_sz2);
    o += xblk_sz << 1;
    dst += dystride << 1;
  }
  for (; j < yblk_sz_2 - 1; j++) {
    /*Lower-left quadrant.*/
    for (i = 0; i < xblk_sz_4; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[3][j][i] << log_blk_sz2;
      f = (src_ll[3][j][i] + src_ll[3][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[3][j][i] + src_ll[3][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[3][j][i] + src_ll[3][j][i + 1]
       + src_ll[3][j + 1][i] + src_ll[3][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       ((((src[3] + o)[i2] + (drc[3] + o)[i2]) << (log_blk_sz2 - 1))
       + a - e + round) >> log_blk_sz2);
      dst[i2+1] = OD_CLAMP255(
       ((((src[3] + o)[i2 + 1] + (drc[3] + o)[i2 + 1]) << (log_blk_sz2 - 1))
       + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(((((src[3] + o + xblk_sz)[i2]
       + (drc[3] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
       + c - g + round) >> log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[3] + o + xblk_sz)[i2 + 1]
       + (drc[3] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
       + d - h + round) >> log_blk_sz2);
    }
    /*Lower-right quadrant.*/
    for (; i < xblk_sz_2 - 1; i++) {
      i2 = i << 1;
      a = dst_ll[j][i] << 2;
      b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
      c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
      d = dst_ll[j][i] + dst_ll[j][i + 1]
       + dst_ll[j + 1][i] + dst_ll[j + 1][i + 1];
      e = src_ll[2][j][i] << log_blk_sz2;
      f = (src_ll[2][j][i] + src_ll[2][j][i + 1]) << (log_blk_sz2 - 1);
      g = (src_ll[2][j][i] + src_ll[2][j + 1][i]) << (log_blk_sz2 - 1);
      h = (src_ll[2][j][i] + src_ll[2][j][i + 1]
       + src_ll[2][j + 1][i] + src_ll[2][j + 1][i + 1]) << (log_blk_sz2 - 2);
      dst[i2] = OD_CLAMP255(
       ((((src[2] + o)[i2] + (drc[2] + o)[i2]) << (log_blk_sz2 - 1))
       + a - e + round) >> log_blk_sz2);
      dst[i2 + 1] = OD_CLAMP255(
       ((((src[2] + o)[i2 + 1] + (drc[2] + o)[i2 + 1]) << (log_blk_sz2 - 1))
       + b - f + round) >> log_blk_sz2);
      (dst + dystride)[i2] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2]
       + (drc[2] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
       + c - g + round) >> log_blk_sz2);
      (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2 + 1]
       + (drc[2] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
       + d - h + round) >> log_blk_sz2);
    }
    /*Lower-right quadrant, last column.*/
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (3*dst_ll[j][i] - dst_ll[j][i - 1]) << 1;
    c = (dst_ll[j][i] + dst_ll[j + 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j + 1][i])
     - (dst_ll[j][i - 1] + dst_ll[j + 1][i - 1]);
    e = src_ll[2][j][i] << log_blk_sz2;
    f = (3*src_ll[2][j][i] - src_ll[2][j][i - 1]) << (log_blk_sz2 - 1);
    g = (src_ll[2][j][i] + src_ll[2][j+1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[2][j][i] + src_ll[2][j + 1][i])
     - (src_ll[2][j][i - 1] + src_ll[2][j + 1][i - 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     ((((src[2] + o)[i2] + (drc[2] + o)[i2 + 1]) << (log_blk_sz2 - 1))
     + a - e + round) >> log_blk_sz2);
    dst[i2 + 1] = OD_CLAMP255(
     ((((src[2] + o)[i2 + 1] + (drc[2] + o)[i2]) << (log_blk_sz2 - 1))
     + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2] +
     (drc[2] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
     + c - g + round) >> log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2 + 1]
     + (drc[2] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
     + d - h + round) >> log_blk_sz2);
    o += xblk_sz << 1;
    dst += dystride << 1;
  }
  /*Lower-left quadrant, last row.*/
  for (i = 0; i < xblk_sz_4; i++) {
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
    c = (3*dst_ll[j][i] - dst_ll[j - 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j][i + 1])
     - (dst_ll[j - 1][i] + dst_ll[j - 1][i + 1]);
    e = src_ll[3][j][i] << log_blk_sz2;
    f = (src_ll[3][j][i] + src_ll[3][j][i + 1]) << (log_blk_sz2 - 1);
    g = (3*src_ll[3][j][i] - src_ll[3][j - 1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[3][j][i] + src_ll[3][j][i + 1])
     - (src_ll[3][j - 1][i] + src_ll[3][j - 1][i + 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     ((((src[3] + o)[i2] + (drc[3] + o)[i2]) << (log_blk_sz2 - 1))
     + a - e + round) >> log_blk_sz2);
    dst[i2+1] = OD_CLAMP255(
     ((((src[3] + o)[i2 + 1] + (drc[3] + o)[i2 + 1]) << (log_blk_sz2 - 1))
     + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(((((src[3] + o + xblk_sz)[i2]
     + (drc[3] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
     + c - g + round) >> log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[3] + o + xblk_sz)[i2 + 1]
     + (drc[3] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
     + d - h + round) >> log_blk_sz2);
  }
  /*Lower-right quadrant, last row.*/
  for (; i < xblk_sz_2 - 1; i++) {
    i2 = i << 1;
    a = dst_ll[j][i] << 2;
    b = (dst_ll[j][i] + dst_ll[j][i + 1]) << 1;
    c = (3*dst_ll[j][i] - dst_ll[j - 1][i]) << 1;
    d = 3*(dst_ll[j][i] + dst_ll[j][i + 1])
     - (dst_ll[j - 1][i] + dst_ll[j - 1][i + 1]);
    e = src_ll[2][j][i] << log_blk_sz2;
    f = (src_ll[2][j][i] + src_ll[2][j][i + 1]) << (log_blk_sz2 - 1);
    g = (3*src_ll[2][j][i] - src_ll[2][j - 1][i]) << (log_blk_sz2 - 1);
    h = (3*(src_ll[2][j][i] + src_ll[2][j][i + 1])
     - (src_ll[2][j - 1][i] + src_ll[2][j - 1][i + 1])) << (log_blk_sz2 - 2);
    dst[i2] = OD_CLAMP255(
     ((((src[2] + o)[i2] + (drc[2] + o)[i2]) << (log_blk_sz2 - 1))
     + a - e + round) >> log_blk_sz2);
    dst[i2 + 1] = OD_CLAMP255(
     ((((src[2] + o)[i2 + 1] + (drc[2] + o)[i2 + 1]) << (log_blk_sz2 - 1))
     + b - f + round) >> log_blk_sz2);
    (dst + dystride)[i2] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2]
     + (drc[2] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
     + c - g + round) >> log_blk_sz2);
    (dst + dystride)[i2 + 1] = OD_CLAMP255(
     ((((src[2] + o + xblk_sz)[i2 + 1] +
     (drc[2] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
     + d - h + round) >>log_blk_sz2);
  }
  /*Lower-right quadrant, last row and column.*/
  i2 = i << 1;
  a = dst_ll[j][i] << 2;
  b = (3*dst_ll[j][i] - dst_ll[j][i - 1]) << 1;
  c = (3*dst_ll[j][i] - dst_ll[j - 1][i]) << 1;
  d = 9*dst_ll[j][i] - 3*(dst_ll[j - 1][i]
   + dst_ll[j][i - 1]) + dst_ll[j - 1][i - 1];
  e = src_ll[2][j][i] << log_blk_sz2;
  f = (3*src_ll[2][j][i] - src_ll[2][j][i - 1]) << (log_blk_sz2 - 1);
  g = (3*src_ll[2][j][i] - src_ll[2][j - 1][i]) << (log_blk_sz2 - 1);
  h = (9*src_ll[2][j][i] - 3*(src_ll[2][j][i - 1] + src_ll[2][j - 1][i])
   + src_ll[2][j - 1][i - 1]) << (log_blk_sz2 - 2);
  dst[i2] = OD_CLAMP255(
   ((((src[2] + o)[i2] + (drc[2] + o)[i2]) << (log_blk_sz2 - 1))
   + a - e + round) >> log_blk_sz2);
  dst[i2 + 1] = OD_CLAMP255(
   ((((src[2] + o)[i2 + 1] + (drc[2] + o)[i2 + 1]) << (log_blk_sz2 - 1))
   + b - f + round) >> log_blk_sz2);
  (dst + dystride)[i2] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2]
   + (drc[2] + o + xblk_sz)[i2]) << (log_blk_sz2 - 1))
   + c - g + round) >> log_blk_sz2);
  (dst + dystride)[i2 + 1] = OD_CLAMP255(((((src[2] + o + xblk_sz)[i2 + 1]
   + (drc[2] + o + xblk_sz)[i2 + 1]) << (log_blk_sz2 - 1))
   + d - h + round) >> log_blk_sz2);
}
#endif

static void od_mc_blend8(od_state *state, unsigned char *dst, int dystride,
 const unsigned char *src[4], int oc, int s,
 int log_xblk_sz, int log_yblk_sz) {
  if (0 && log_xblk_sz > 1 && log_yblk_sz > 1) {
    /*Perform multiresolution blending.*/
    if (s == 3) {
      od_mc_blend_multi8(dst, dystride, src, log_xblk_sz, log_yblk_sz);
    }
    else {
      od_mc_blend_multi_split8(dst, dystride, src,
       oc, s, log_xblk_sz, log_yblk_sz);
    }
  }
  else {
    /*The block is too small; perform normal blending.*/
    if (s == 3) {
      od_mc_blend_full8(state, dst, dystride, src,
       log_xblk_sz, log_yblk_sz);
    }
    else {
      od_mc_blend_full_split8(state, dst, dystride, src,
       oc, s, log_xblk_sz, log_yblk_sz);
    }
  }
}

void od_mc_predict8(od_state *state, unsigned char *dst, int dystride,
 const unsigned char *src, int systride,
                    const ogg_int32_t mvx[4], /*this is x coord for the four
                                               motion vectors of the four
                                               corners (in rotation not
                                               raster order) */
 const ogg_int32_t mvy[4],
                    int interp_type, /* type of interpolation (bbbb for now) */
                    int oc,  /* index of outside corner  */
                    int s, /* two split flags that indicate if the corners are split*/
                    int log_xblk_sz,   /* log 2 of block size */
                    int log_yblk_sz
) {
  const unsigned char *pred[4];
  unsigned char __attribute__((aligned(16))) buf[4][16*16];
  OD_ASSERT(interp_type == OD_MC_INTERP_BBBB);
  (void)interp_type;
  od_mc_predict1fmv8(state, buf[0], src, systride,
   mvx[0], mvy[0], log_xblk_sz, log_yblk_sz);
  pred[0] = buf[0];
  if (mvx[1] == mvx[0] && mvy[1] == mvy[1]) pred[1] = pred[0];
  else {
    od_mc_predict1fmv8(state, buf[1], src, systride,
     mvx[1], mvy[1], log_xblk_sz, log_yblk_sz);
    pred[1] = buf[1];
  }
  if (mvx[2] == mvx[0] && mvy[2] == mvy[0]) pred[2] = pred[0];
  else if (mvx[2] == mvx[1] && mvy[2] == mvy[1]) pred[2] = pred[1];
  else {
    od_mc_predict1fmv8(state, buf[2], src, systride,
     mvx[2], mvy[2], log_xblk_sz, log_yblk_sz);
    pred[2] = buf[2];
  }
  if (mvx[3] == mvx[0] && mvy[3] == mvy[0]) pred[3] = pred[0];
  else if (mvx[3] == mvx[1] && mvy[3] == mvy[1]) pred[3] = pred[1];
  else if (mvx[3] == mvx[2] && mvy[3] == mvy[2]) pred[3] = pred[2];
  else {
    od_mc_predict1fmv8(state, buf[3], src, systride,
     mvx[3], mvy[3], log_xblk_sz, log_yblk_sz);
    pred[3] = buf[3];
  }
  od_mc_blend8(state, dst, dystride, pred,
   oc, s, log_xblk_sz, log_yblk_sz);
}

/*Gets the predictor for a given MV node at the given MV resolution.*/
void od_state_get_predictor(od_state *state,
 int pred[2], int vx, int vy, int level, int mv_res) {
  int nhmvbs;
  int nvmvbs;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  if (vx < 2 || vy < 2 || vx > nhmvbs - 2 || vy > nvmvbs - 2) {
    pred[0] = pred[1] = 0;
  }
  else {
    int mvb_sz;
    od_mv_grid_pt *cneighbors[4];
    int a[4][2];
    int ncns;
    int ci;
    ncns = 4;
    mvb_sz = 1 << ((4 - level) >> 1);
    if (level == 0) {
      cneighbors[0] = state->mv_grid[vy - 4] + vx - 4;
      cneighbors[1] = state->mv_grid[vy - 4] + vx;
      cneighbors[2] = state->mv_grid[vy - 4] + vx + 4;
      cneighbors[3] = state->mv_grid[vy] + vx - 4;
    }
    else {
      if (level & 1) {
        cneighbors[0] = state->mv_grid[vy - mvb_sz] + vx - mvb_sz;
        cneighbors[1] = state->mv_grid[vy - mvb_sz] + vx + mvb_sz;
        cneighbors[2] = state->mv_grid[vy + mvb_sz] + vx - mvb_sz;
        cneighbors[3] = state->mv_grid[vy + mvb_sz] + vx + mvb_sz;
      }
      else {
        cneighbors[0] = state->mv_grid[vy - mvb_sz] + vx;
        cneighbors[1] = state->mv_grid[vy] + vx - mvb_sz;
        /*NOTE: Only one of these candidates can be excluded at a time, so
           there will always be at least 3.*/
        if (vx + mvb_sz > ((vx + 3) & ~3)) ncns--;
        else cneighbors[2] = state->mv_grid[vy] + vx + mvb_sz;
        if (vy + mvb_sz > ((vy + 3) & ~3)) ncns--;
        else cneighbors[ncns - 1] = state->mv_grid[vy + mvb_sz] + vx;
      }
    }
    for (ci = 0; ci < ncns; ci++) {
      a[ci][0] = cneighbors[ci]->mv[0];
      a[ci][1] = cneighbors[ci]->mv[1];
    }
    /*Median-of-4.*/
    if (ncns > 3) {
      OD_LOG((OD_LOG_MOTION_COMPENSATION, OD_LOG_DEBUG, "Median of 4:"));
      OD_LOG((OD_LOG_MOTION_COMPENSATION, OD_LOG_DEBUG,
       "(%i, %i) (%i, %i) (%i, %i) (%i, %i)", a[0][0], a[0][1],
       a[1][0], a[1][1], a[2][0], a[2][1], a[3][0], a[3][1]));
/*Sorting network for 4 elements:
0000 0001 0010 0011 0100 0101 0110 0111 1000 1001 1010 1011 1100 1101 1110 1111
0001 0010 0011 0100 0101 0110 0111 1001 1010 1011 1101
0:1
0010 0010 0011 0100 0110 0110 0111 1010 1010 1011 1110
0010 0011 0100 0110 0111 1010 1011
2:3
0010 0011 1000 1010 1011 1100 1110
0010 0011 1010 1011
0:2
0010 0110 1010 1110
0010 0110 1010
1:3
1000 1100 1010
1010
This last compare is unneeded for a median:
1:2
1100*/
      OD_SORT2I(a[0][0], a[1][0]);
      OD_SORT2I(a[0][1], a[1][1]);
      OD_SORT2I(a[2][0], a[3][0]);
      OD_SORT2I(a[2][1], a[3][1]);
      OD_SORT2I(a[0][0], a[2][0]);
      OD_SORT2I(a[0][1], a[2][1]);
      OD_SORT2I(a[1][0], a[3][0]);
      OD_SORT2I(a[1][1], a[3][1]);
      OD_LOG((OD_LOG_MOTION_COMPENSATION, OD_LOG_DEBUG,
       "(%i, %i) (%i, %i) (%i, %i) (%i, %i)", a[0][0], a[0][1],
       a[1][0], a[1][1], a[2][0], a[2][1], a[3][0], a[3][1]));
      pred[0] = OD_DIV_POW2_RE(a[1][0] + a[2][0], mv_res + 1);
      pred[1] = OD_DIV_POW2_RE(a[1][1] + a[2][1], mv_res + 1);
    }
    /*Median-of-3.*/
    else {
      OD_LOG((OD_LOG_MOTION_COMPENSATION, OD_LOG_DEBUG, "Median of 3:"));
      OD_LOG((OD_LOG_MOTION_COMPENSATION, OD_LOG_DEBUG,
       "(%i, %i) (%i, %i) (%i, %i)",
       a[0][0], a[0][1], a[1][0], a[1][1], a[2][0], a[2][1]));
      OD_SORT2I(a[0][0], a[1][0]);
      OD_SORT2I(a[0][1], a[1][1]);
      OD_SORT2I(a[1][0], a[2][0]);
      OD_SORT2I(a[1][1], a[2][1]);
      OD_SORT2I(a[0][0], a[1][0]);
      OD_SORT2I(a[0][1], a[1][1]);
      OD_LOG((OD_LOG_MOTION_COMPENSATION, OD_LOG_DEBUG,
       "(%i, %i) (%i, %i) (%i, %i)",
       a[0][0], a[0][1], a[1][0], a[1][1], a[2][0], a[2][1]));
      pred[0] = OD_DIV_POW2_RE(a[1][0], mv_res);
      pred[1] = OD_DIV_POW2_RE(a[1][1], mv_res);
    }
  }
}

/*Probability that a given level 1 MV is coded, given the nearest
   left and up l1 neighbors and the consistency of the motion field
   down and to the righ.*/
int od_mv_level1_prob(od_mv_grid_pt **grid, int vx, int vy) {
  const int probs[3][3] =
   {{28323, 30610, 32128}, {18468, 21082, 24253}, {10799, 12839, 14144}};
  od_mv_grid_pt *vur;
  od_mv_grid_pt *vdr;
  od_mv_grid_pt *vdl;
  int lf;
  int uf;
  int rf;
  int bf;
  vur = &(grid[vy - 2][vx + 2]);
  vdr = &(grid[vy + 2][vx + 2]);
  vdl = &(grid[vy + 2][vx - 2]);
  lf = vx > 3 ? grid[vy][vx - 4].valid : 0;
  uf = vy > 3 ? grid[vy - 4][vx].valid : 0;
  rf = (vur->mv[0] == vdr->mv[0]) && (vur->mv[1] == vdr->mv[1]);
  bf = (vdr->mv[0] == vdl->mv[0]) && (vdr->mv[1] == vdl->mv[1]);
  return probs[lf + uf][rf + bf];
}

#if 0
# include <stdio.h>
# include <string.h>

static unsigned char mask[4][4] = {
  {  0,  8,  2, 10 },
  { 12,  4, 14,  6 },
  {  3, 11,  1,  9 },
  { 15,  7, 13,  5 }
};

static unsigned char img[16*7][16*7];

static ogg_int32_t mvs[4][4][2];
static ogg_int32_t mvs2[4][4][2];

static int edge_types[4][4] = {
  { 0x0, 0x1, 0x2, 0x8 },
  { 0x4, 0x6, 0xA, 0xC },
  { 0x5, 0x7, 0xE, 0xD },
  { 0x3, 0xF, 0xB, 0x9 }
};

static int edge_types2[8][8] = {
  { 0x0, 0x4, 0x1, 0x5, 0x0, 0x2, 0xA, 0x8 },
  { 0x2, 0x9, 0x0, 0x1, 0x2, 0xA, 0x8, 0x0 },
  { 0x6, 0x8, 0x2, 0xA, 0xC, 0x6, 0xC, 0x4 },
  { 0x5, 0x4, 0x6, 0xE, 0x9, 0x3, 0xF, 0xD },
  { 0x1, 0x5, 0x3, 0xF, 0xA, 0xE, 0xB, 0x9 },
  { 0x6, 0xD, 0x4, 0x7, 0xE, 0xF, 0xE, 0xC },
  { 0x5, 0x7, 0xB, 0xF, 0xD, 0x7, 0xB, 0xD },
  { 0x3, 0xB, 0xC, 0x7, 0x9, 0x3, 0x8, 0x1 }
};

static void fill_mvs(int log_blk_sz) {
  int i;
  int j;
  for (j = 0; j < 4; j++) {
    for (i = 0; i < 4; i++) {
      mvs[j][i][0] =
       (mask[(j + 0) & 3][(i + 0) & 3] - 8) << (log_blk_sz + 15) ^
       mask[(j + 2) & 3][(i + 2) & 3] << (log_blk_sz + 12);
      mvs[j][i][1] =
       (mask[(j + 1) & 3][(i + 1) & 3] - 8) << (log_blk_sz + 15) ^
       mask[(j + 3) & 3][(i + 3) & 3] << (log_blk_sz + 12);
      mvs2[j][i][0] =
       (mask[(j + 3) & 3][(i + 3) & 3] - 8) << (log_blk_sz + 15) ^
       mask[(j + 1) & 3][(i + 1) & 3] << (log_blk_sz + 12);
      mvs2[j][i][1] =
       (mask[(j + 2) & 3][(i + 2) & 3] - 8) << (log_blk_sz + 15) ^
       mask[(j + 0) & 3][(i + 0) & 3] << (log_blk_sz + 12);
    }
  }
}

static void fill_img(int log_blk_sz) {
  int i;
  int j;
  for (j = 0; j < 7; j++) {
    for (i = 0; i < 7; i++) {
      int c;
      int x;
      int y;
      c = mask[j & 3][i & 3];
      c = c << 4 | 15;
      for (y = j << log_blk_sz; y < (j + 1) << log_blk_sz; y++) {
        for (x = i << log_blk_sz; x < (i + 1) << log_blk_sz; x++) {
          img[y][x] = c;
        }
      }
    }
  }
  for (j = 0; j < 6 << log_blk_sz; j++) {
    for (i = 0; i < 6 << log_blk_sz; i++) {
      printf("%2X%c", img[j][i], i + 1 < 6 << log_blk_sz ? ' ' : '\n');
    }
  }
}

int main(void) {
  int log_blk_sz;
  for (log_blk_sz = 2; log_blk_sz <= 4; log_blk_sz++) {
    int blk_sz;
    int i;
    int j;
    blk_sz = 1 << log_blk_sz;
    fill_img(log_blk_sz);
    fill_mvs(log_blk_sz);
    for (j = 0; j < 4; j++) {
      for (i = 0; i < 4; i++) {
        ogg_int32_t mvx[4];
        ogg_int32_t mvy[4];
        unsigned char dst[17][17];
        unsigned char dst2[4][9][9];
        unsigned mismatch[4][9][9];
        int etype;
        int x;
        int y;
        int c;
        mvx[0] = mvs[j][i][0];
        mvy[0] = mvs[j][i][1];
        mvx[1] = mvs[j][(i + 1) & 3][0];
        mvy[1] = mvs[j][(i + 1) & 3][1];
        mvx[2] = mvs[(j + 1) & 3][(i + 1) & 3][0];
        mvy[2] = mvs[(j + 1) & 3][(i + 1) & 3][1];
        mvx[3] = mvs[(j + 1) & 3][i][0];
        mvy[3] = mvs[(j + 1) & 3][i][1];
        etype = edge_types[j][i];
        printf("Block (%i, %i): size %i, "
         "interpolation type: %c%c%c%c (0x%X)\n", i, j, 1 << log_blk_sz,
         etype & 1 ? 'V' : 'B', etype & 2 ? 'V' : 'B',
         etype & 4 ? 'V' : 'B', etype & 8 ? 'V' : 'B', etype);
        printf("<%8.4f, %8.4f> <%8.4f, %8.4f>\n",
         mvx[0]/(double)0x40000, mvy[0]/(double)0x40000,
         mvx[1]/(double)0x40000, mvy[1]/(double)0x40000);
        printf("<%8.4f, %8.4f> <%8.4f, %8.4f>\n",
         mvx[3]/(double)0x40000, mvy[3]/(double)0x40000,
         mvx[2]/(double)0x40000, mvy[2]/(double)0x40000);
        od_mc_predict8(dst[0], sizeof(dst[0]),
         img[(j + 1) << log_blk_sz] + ((i + 1) << log_blk_sz),
         sizeof(img[0]), mvx, mvy, etype, 0, 0, 0, log_blk_sz, log_blk_sz);
        for (y = 0; y < blk_sz; y++) {
          for (x = 0; x < blk_sz; x++) {
            printf("%2X%c", dst[y][x], x + 1 < blk_sz ? ' ' : '\n');
          }
        }
        printf("\n");
        for (c = 0; c < 4; c++) {
          int s1;
          int s3;
          mvx[0] = mvs[j][i][0];
          mvy[0] = mvs[j][i][1];
          mvx[1] = mvs[j][(i + 1) & 3][0];
          mvy[1] = mvs[j][(i + 1) & 3][1];
          mvx[2] = mvs[(j+1) & 3][(i+1) & 3][0];
          mvy[2] = mvs[(j+1) & 3][(i+1) & 3][1];
          mvx[3] = mvs[(j+1) & 3][i][0];
          mvy[3] = mvs[(j+1) & 3][i][1];
          mvx[(c + 2) & 3] = mvs2[j][i][0];
          mvy[(c + 2) & 3] = mvs2[j][i][1];
          etype = edge_types2[j << 1 | (c >> 1)][i << 1 | (((c + 1)&3) >> 1)];
          if (1 || !(c & 1)) {
            etype &= ~(1 << ((c + 1) & 3));
            etype |= (etype | etype << 4) >> 3 & 1 << ((c + 1) & 3);
            s1 = 0;
          }
          else s1 = 1;
          if (s1 || (etype >> c & 1)) {
            mvx[(c + 1) & 3] = (mvx[c] + mvx[(c + 1) & 3]) >> 1;
            mvy[(c + 1) & 3] = (mvy[c] + mvy[(c + 1) & 3]) >> 1;
          }
          if (1 || (c & 1)) {
            etype &= ~(1 << ((c + 2) &3));
            etype |= (etype | etype << 4) >> 1 & 1 << ((c + 2) & 3);
            s3 = 0;
          }
          else s3 = 1;
          if (s3 || (etype << (-c & 3) & 8)) {
            mvx[(c + 3) & 3] = (mvx[c] + mvx[(c + 3) & 3]) >> 1;
            mvy[(c + 3) & 3] = (mvy[c] + mvy[(c + 3) & 3]) >> 1;
          }
          printf("Block (%i.%i, %i.%i): size %i, "
           "interpolation type: %c%c%c%c (0x%X)\n",
           i, (((c + 1) & 3) >> 1)*5, j, (c >> 1)*5, 1 << (log_blk_sz - 1),
           etype & 1 ? 'V' : 'B', etype & 2 ? 'V' : 'B',
           etype & 4 ? 'V' : 'B', etype & 8 ? 'V' : 'B', etype);
          printf("<%9.5f, %9.5f> <%9.5f, %9.5f>\n",
           mvx[0]/(double)0x40000, mvy[0]/(double)0x40000,
           mvx[1]/(double)0x40000, mvy[1]/(double)0x40000);
          printf("<%9.5f, %9.5f> <%9.5f, %9.5f>\n",
           mvx[3]/(double)0x40000, mvy[3]/(double)0x40000,
           mvx[2]/(double)0x40000, mvy[2]/(double)0x40000);
          od_mc_predict8(dst2[c][0], sizeof(dst2[c][0]),
           img[((j + 1) << 1 | (c >> 1)) << (log_blk_sz - 1)]
           + (((i + 1) << 1 | (((c + 1) & 3) >> 1)) << (log_blk_sz - 1)),
           sizeof(img[0]), mvx, mvy, etype, c, s1, s3,
           log_blk_sz - 1, log_blk_sz - 1);
          memset(mismatch[c][0], 0, sizeof(mismatch[c]));
          switch(c) {
            case 0: {
              for (x = 0; x < blk_sz >> 1; x++) {
                if (dst2[c][0][x] != dst[0][x]) mismatch[c][0][x]++;
              }
              for (y = 1; y < blk_sz >> 1; y++) {
                if (dst2[c][y][0] != dst[y][0]) mismatch[c][y][0]++;
              }
              break;
            }
            case 1: {
              for (x = 0; x < blk_sz >> 1; x++) {
                if (dst2[c][0][x] != dst[0][x + (blk_sz >> 1)]) {
                  mismatch[c][0][x]++;
                }
              }
              for (y = 1; y < blk_sz >> 1; y++) {
                if (dst2[c][y][blk_sz >> 1] != dst[y][blk_sz]) {
                  mismatch[c][y][blk_sz >> 1]++;
                }
              }
              for (y = 0; y < blk_sz >> 1; y++) {
                if (dst2[c][y][0] != dst2[0][y][blk_sz >> 1]) {
                  mismatch[c][y][0]++;
                }
              }
              break;
            }
            case 2: {
              for (x = 0; x < blk_sz >> 1; x++) {
                if (dst2[c][0][x] != dst2[1][blk_sz >> 1][x]) {
                  mismatch[c][0][x]++;
                }
              }
              for (y = 0; y < blk_sz >> 1; y++) {
                if (dst2[c][y][blk_sz >> 1] !=
                 dst[y + (blk_sz >> 1)][blk_sz]) {
                  mismatch[c][y][blk_sz >> 1]++;
                }
              }
              for (x = 0; x < blk_sz >> 1; x++) {
                if (dst2[c][blk_sz >> 1][x] !=
                 dst[blk_sz][x + (blk_sz >> 1)]) {
                  mismatch[c][blk_sz >> 1][x]++;
                }
              }
              break;
            }
            case 3: {
              for (x = 0; x < blk_sz >> 1; x++) {
                if (dst2[c][0][x] != dst2[0][ blk_sz >> 1][x]) {
                  mismatch[c][0][x]++;
                }
              }
              for (y = 1; y < blk_sz >> 1; y++) {
                if (dst2[c][y][blk_sz >> 1] != dst2[2][y][0]) {
                  mismatch[c][y][blk_sz >> 1]++;
                }
              }
              for (x = 0; x < blk_sz >> 1; x++) {
                if (dst2[c][blk_sz >> 1][x] != dst[blk_sz][x]) {
                  mismatch[c][blk_sz >> 1][x]++;
                }
              }
              for (y = 0; y < blk_sz >> 1; y++) {
                if (dst2[c][y][0] != dst[y + (blk_sz >> 1)][0]) {
                  mismatch[c][y][0]++;
                }
              }
              break;
            }
          }
          for (y = 0; y < blk_sz >> 1; y++) {
            for (x = 0; x < blk_sz >> 1; x++) {
              printf("%c%2X", mismatch[c][y][x] ? '!' : ' ', dst2[c][y][x]);
            }
            printf("\n");
          }
          printf("\n");
        }
      }
    }
  }
  return 0;
}

#endif
