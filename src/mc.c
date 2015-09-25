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
#include "state.h"

/*Motion compensation routines shared between the encoder and decoder.*/

/*ME/MC Subpel interpolation filter set.*/
/*Based on fractional part of MV, which is 000, 001, ..., 111
   for 1/8 pel precision, apply corresponding filter.
  For this, we need 8 kinds of 1D filter, each corresponds to
   one of eight fractional positions.*/
# if 0
  /*6-tap filter #1 :
    Windowed-sinc 6-tap for 1/2 pel, and bilinear for 1/4 and 1/8 pel.
     (Same design as adopted in frame level upsampling function,
     which is currently used in master branch.)
    Filter coefficient is scaled up by 7 bit.*/
  const OD_ALIGN16(int16_t) OD_SUBPEL_FILTER_SET[8][8] = {
  /*Padded with extra 0's to make each row have 8 elements.*/
  /* -2   -1  [ 0    1]   2  3  : pixel position in support region.*/
    { 0,   0, 128,   0,   0, 0, 0, 0 },
    { 1,  -5, 116,  20,  -5, 1, 0, 0 },
    { 2, -10, 104,  40, -10, 2, 0, 0 },
    { 3, -15,  92,  60, -15, 3, 0, 0 },
    { 4, -20,  80,  80, -20, 4, 0, 0 },
    { 3, -15,  60,  92, -15, 3, 0, 0 },
    { 2, -10,  40, 104, -10, 2, 0, 0 },
    { 1,  -5,  20, 116,  -5, 1, 0, 0 }
  };
# else
  /*6-tap filter #2 : sinc(x) * sinc(x/3)
    Windowed-sinc 6-tap: sinc(x) * sinc(x/3) for 1/4 pel and 1/8 pel.
    We keep the 1/2 pel filter from set #1 above.
    Filter coefficient is scaled up by 7 bit.*/
  const OD_ALIGN16(int16_t) OD_SUBPEL_FILTER_SET[8][8] = {
  /*Extra 0's are padded to make each row have 8 elements.*/
  /* -2   -1  [ 0    1]   2  3  : pixel position in support region.*/
    { 0, 0, 128, 0, 0, 0, 0, 0},
    { 1, -9, 122, 18, -5, 1, 0, 0},
    { 3, -15, 112, 37, -11, 2, 0, 0},
    { 3, -18, 97, 58, -15, 3, 0, 0},
    { 4, -20, 80, 80, -20, 4, 0, 0},
    { 3, -15, 58, 97, -18, 3, 0, 0},
    { 2, -11, 37, 112, -15, 3, 0, 0},
    { 1, -5, 18, 122, -9, 1, 0, 0},
  };
# endif

/*Form the prediction given by one fixed motion vector.
  dst: The destination buffer (xstride must be 1).
  src: The source buffer (xstride must be 1).
  systride: The byte offset between source pixel rows.
  mvx: The X component of the motion vector.
  mvy: The Y component of the motion vector.
  log_xblk_sz: The log base 2 of the horizontal block dimension.
  log_yblk_sz: The log base 2 of the vertical block dimension.*/
/*TODO: this code does not currently special-case the boundaries of the buffer.
  It simply assumes there's sufficient extra padding to not cause a bounds
   overrun.
  For now, we add the extra padding, but we should eventually add
   special-casing to the code to avoid the additional ~useless buffer
   overhead. */
void od_mc_predict1fmv8_c(unsigned char *dst, const unsigned char *src,
 int systride, int32_t mvx, int32_t mvy,
 int log_xblk_sz, int log_yblk_sz) {
  int mvxf;
  int mvyf;
  int xblk_sz;
  int yblk_sz;
  int i;
  int j;
  /*Pointer to the start of an image block in local buffer (defined
     below, buff[]), where the buffer contans the top and bottom apron
     area of the image block.
    Used as output for 1st stage horizontal filtering then as input for
     2nd stage vertical filtering.*/
  int16_t *buff_p;
  /*A pointer to input row for both 1st and 2nd stage filtering*/
  const unsigned char *src_p;
  unsigned char *dst_p;
  /*1D filter chosen for the current fractional position of x mv.*/
  const int16_t *fx;
  const int16_t *fy;
  /*2D buffer to store the result of 1st stage (i.e. horizontal) 1D filtering
     of a block. The 1st stage filtering requires to output results for
     top and bottom aprons of input image block, because the 2nd stage
     filtering (i.e vertical) requires support region on those apron pixels.
    The size of the buffer is :
     wxh = OD_MVBSIZE_MAX x (OD_MVBSIZE_MAX + OD_SUBPEL_BUFF_APRON_SZ).*/
  int16_t buff[(OD_MVBSIZE_MAX + OD_SUBPEL_BUFF_APRON_SZ)
   *OD_MVBSIZE_MAX];
  int32_t sum;
  int k;
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  src_p = src + (mvx >> 3) + (mvy >> 3)*systride;
  dst_p = dst;
  /*Fetch LSB 3 bits, i.e. fractional MV.*/
  mvxf = mvx & 0x07;
  mvyf = mvy & 0x07;
  /*Check whether mvxf and mvyf are in the range [0...7],
     i.e. downto 1/8 precision.*/
  OD_ASSERT(mvxf <= 7);
  fx = OD_SUBPEL_FILTER_SET[mvxf];
  OD_ASSERT(mvyf <= 7);
  fy = OD_SUBPEL_FILTER_SET[mvyf];
  /*MC with subpel MV?*/
  if (mvxf || mvyf) {
    buff_p = buff;
    src_p -= systride*OD_SUBPEL_TOP_APRON_SZ;
    /*1st stage 1D filtering, Horizontal.*/
    if (mvxf) {
      for (j = -OD_SUBPEL_TOP_APRON_SZ;
       j < yblk_sz + OD_SUBPEL_BOTTOM_APRON_SZ; j++) {
        for (i = 0; i < xblk_sz; i++) {
          sum = 0;
          for (k = 0; k < OD_SUBPEL_FILTER_TAP_SIZE; k++) {
            sum += src_p[i + k - OD_SUBPEL_TOP_APRON_SZ]*fx[k];
          }
          buff_p[i] = sum - (128 << OD_SUBPEL_COEFF_SCALE);
        }
        src_p += systride;
        buff_p += xblk_sz;
      }
    }
    /*The mvx is in integer position.*/
    else {
      for (j = -OD_SUBPEL_TOP_APRON_SZ;
       j < yblk_sz + OD_SUBPEL_BOTTOM_APRON_SZ; j++) {
        for (i = 0; i < xblk_sz; i++) {
          buff_p[i] = (src_p[i] << OD_SUBPEL_COEFF_SCALE)
           - OD_SUBPEL_COEFF_NORMALIZE;
        }
        src_p += systride;
        buff_p += xblk_sz;
      }
    }
    /*2nd stage 1D filtering, Vertical.*/
    buff_p = buff + xblk_sz*OD_SUBPEL_TOP_APRON_SZ;
    if (mvyf) {
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i++) {
          sum = 0;
          for (k = 0; k < OD_SUBPEL_FILTER_TAP_SIZE; k++) {
            sum += buff_p[i + (k - OD_SUBPEL_TOP_APRON_SZ)*xblk_sz] * fy[k];
          }
          dst_p[i] = OD_CLAMP255((sum + OD_SUBPEL_RND_OFFSET3)
           >> OD_SUBPEL_COEFF_SCALE2);
        }
        buff_p += xblk_sz;
        dst_p += xblk_sz;
      }
    }
    /*The mvy is in integer position.*/
    else {
      for (j = 0; j < yblk_sz; j++) {
        for (i = 0; i < xblk_sz; i++) {
          dst_p[i] = OD_CLAMP255((buff_p[i] + OD_SUBPEL_RND_OFFSET4)
           >> OD_SUBPEL_COEFF_SCALE);
        }
        buff_p += xblk_sz;
        dst_p += xblk_sz;
      }
    }
  }
  /*MC with full-pel MV, i.e. integer position.*/
  else {
    for (j = 0; j < yblk_sz; j++) {
      OD_COPY(dst_p, src_p, xblk_sz);
      src_p += systride;
      dst_p += xblk_sz;
    }
  }
}

static void od_mc_predict1fmv8(od_state *state, unsigned char *dst,
 const unsigned char *src, int systride, int32_t mvx, int32_t mvy,
 int log_xblk_sz, int log_yblk_sz) {
  OD_ASSERT(OD_SUBPEL_TOP_APRON_SZ <= OD_RESAMPLE_PADDING
   && OD_SUBPEL_BOTTOM_APRON_SZ <= OD_RESAMPLE_PADDING);
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
      int32_t a;
      int32_t b;
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
  nhmvbs = state->nhmvbs;
  nvmvbs = state->nvmvbs;
  for (vy = 0; vy <= nvmvbs; vy++) {
    od_mv_grid_pt *grid;
    grid = state->mv_grid[vy];
    for (vx = 0; vx <= nhmvbs; vx++) {
      grid[vx].valid = 0;
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
      a = (int)((((int32_t)a << log_yblk_sz) +
       (int32_t)(b - a)*j) >> log_blk_sz2);
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
      a = (int)((((int32_t)a << log_yblk_sz)
       + (int32_t)(b - a)*j) >> log_blk_sz2);
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
      a = (int)(((int32_t)a << log_yblk_sz)
       + (int32_t)(b - a)*j >> log_blk_sz2);
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
      a = (int)((((int32_t)a << log_yblk_sz)
       + (int32_t)(b - a)*j) >> log_blk_sz2);
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
      int32_t a;
      int32_t b;
      int32_t c;
      int32_t d;
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
  int log_blk_sz2m1;
  int round;
  int i;
  int j;
  int k;
  /*Perform multiresolution blending.*/
  xblk_sz = 1 << log_xblk_sz;
  yblk_sz = 1 << log_yblk_sz;
  o0 = 0;
  dst0 = dst;
  log_blk_sz2m1 = log_xblk_sz + log_yblk_sz - 1;
  round = 1 << (log_blk_sz2m1 - 1);
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
      a = ll[0];
      b = (ll[1] - a)*sw[1];
      c = (ll[2] - a)*sw[2];
      d = (ll[3] - a)*sw[3];
      a = (int)(((a << log_blk_sz2m1) + b + c + d + round) >> log_blk_sz2m1);
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


void od_mc_predict8_singleref(od_state *state, unsigned char *dst,
 int dystride, const unsigned char *src, int systride,
 const int32_t mvx[4], /* This is x coord for the four
                            motion vectors of the four corners
                            (in rotation not raster order). */
 const int32_t mvy[4],
 int oc, /* Index of outside corner. */
 int s, /* Two split flags that indicate if the corners are split. */
 int log_xblk_sz,   /* Log 2 of block size. */
 int log_yblk_sz
) {
  const unsigned char *pred[4];
  od_mc_predict1fmv8(state, state->mc_buf[0], src, systride,
   mvx[0], mvy[0], log_xblk_sz, log_yblk_sz);
  pred[0] = state->mc_buf[0];
  if (mvx[1] == mvx[0] && mvy[1] == mvy[0]) pred[1] = pred[0];
  else {
    od_mc_predict1fmv8(state, state->mc_buf[1], src, systride,
     mvx[1], mvy[1], log_xblk_sz, log_yblk_sz);
    pred[1] = state->mc_buf[1];
  }
  if (mvx[2] == mvx[0] && mvy[2] == mvy[0]) pred[2] = pred[0];
  else if (mvx[2] == mvx[1] && mvy[2] == mvy[1]) pred[2] = pred[1];
  else {
    od_mc_predict1fmv8(state, state->mc_buf[2], src, systride,
     mvx[2], mvy[2], log_xblk_sz, log_yblk_sz);
    pred[2] = state->mc_buf[2];
  }
  if (mvx[3] == mvx[0] && mvy[3] == mvy[0]) pred[3] = pred[0];
  else if (mvx[3] == mvx[1] && mvy[3] == mvy[1]) pred[3] = pred[1];
  else if (mvx[3] == mvx[2] && mvy[3] == mvy[2]) pred[3] = pred[2];
  else {
    od_mc_predict1fmv8(state, state->mc_buf[3], src, systride,
     mvx[3], mvy[3], log_xblk_sz, log_yblk_sz);
    pred[3] = state->mc_buf[3];
  }
  od_mc_blend8(state, dst, dystride, pred,
   oc, s, log_xblk_sz, log_yblk_sz);
}

/* TODO: Identical MV optimizations like od_mc_predict8_singleref. */
void od_mc_predict8(od_state *state, unsigned char *dst,
 int dystride, const unsigned char *src[4], int systride,
 const int32_t mvx[4], const int32_t mvy[4],
 int oc,  /* index of outside corner  */
 int s, /* two split flags that indicate if the corners are split*/
 int log_xblk_sz,   /* log 2 of block size */
 int log_yblk_sz) {
  const unsigned char *pred[4];
  if ((src[0] == src[1]) && (src[0] == src[2]) && (src[0] == src[3])) {
    od_mc_predict8_singleref(state, dst, dystride, src[0], systride, mvx, mvy,
     oc, s, log_xblk_sz, log_yblk_sz);
    return;
  }
  od_mc_predict1fmv8(state, state->mc_buf[0], src[0], systride,
   mvx[0], mvy[0], log_xblk_sz, log_yblk_sz);
  pred[0] = state->mc_buf[0];
  od_mc_predict1fmv8(state, state->mc_buf[1], src[1], systride,
   mvx[1], mvy[1], log_xblk_sz, log_yblk_sz);
  pred[1] = state->mc_buf[1];
  od_mc_predict1fmv8(state, state->mc_buf[2], src[2], systride,
   mvx[2], mvy[2], log_xblk_sz, log_yblk_sz);
  pred[2] = state->mc_buf[2];
  od_mc_predict1fmv8(state, state->mc_buf[3], src[3], systride,
   mvx[3], mvy[3], log_xblk_sz, log_yblk_sz);
  pred[3] = state->mc_buf[3];

  od_mc_blend8(state, dst, dystride, pred,
   oc, s, log_xblk_sz, log_yblk_sz);
}

int od_mc_get_ref_predictor(od_state *state, int vx, int vy, int level) {
  static const od_mv_grid_pt ZERO_GRID_PT = { {0, 0}, 1, OD_FRAME_PREV};
  int mvb_sz;
  const od_mv_grid_pt *cneighbors[4];
  int ncns;
  int ci;
  int hist[2] = {0, 0};
  int max_count = 0;
  int max_ref = OD_FRAME_PREV;
  ncns = 4;
  mvb_sz = 1 << ((4 - level) >> 1);
  mvb_sz = 1 << ((OD_MC_LEVEL_MAX - level) >> 1);
  if (level == 0) {
    if (vy >= mvb_sz) {
      cneighbors[0] = vx >= mvb_sz ?
       state->mv_grid[vy - mvb_sz] + vx - mvb_sz : &ZERO_GRID_PT;
      cneighbors[1] = state->mv_grid[vy - mvb_sz] + vx;
      cneighbors[2] = vx + mvb_sz <= state->nhmvbs ?
       state->mv_grid[vy - mvb_sz] + vx + mvb_sz : &ZERO_GRID_PT;
    }
    else cneighbors[2] = cneighbors[1] = cneighbors[0] = &ZERO_GRID_PT;
    cneighbors[3] = vx >= mvb_sz ?
     state->mv_grid[vy] + vx - mvb_sz : &ZERO_GRID_PT;
  }
  else {
    if (level & 1) {
      cneighbors[0] = state->mv_grid[vy - mvb_sz] + vx - mvb_sz;
      cneighbors[1] = state->mv_grid[vy - mvb_sz] + vx + mvb_sz;
      cneighbors[2] = state->mv_grid[vy + mvb_sz] + vx - mvb_sz;
      cneighbors[3] = state->mv_grid[vy + mvb_sz] + vx + mvb_sz;
    }
    else {
      cneighbors[0] = vy >= mvb_sz ?
       state->mv_grid[vy - mvb_sz] + vx : &ZERO_GRID_PT;
      cneighbors[1] = vx >= mvb_sz ?
       state->mv_grid[vy] + vx - mvb_sz : &ZERO_GRID_PT;
      /*NOTE: Only one of these candidates can be excluded at a time, so
         there will always be at least 3.*/
      if (vx > 0 && vx + mvb_sz > ((vx + OD_MVB_MASK) & ~OD_MVB_MASK)) ncns--;
      else cneighbors[2] = state->mv_grid[vy] + vx + mvb_sz;
      if (vy > 0 && vy + mvb_sz > ((vy + OD_MVB_MASK) & ~OD_MVB_MASK)) ncns--;
      else cneighbors[ncns - 1] = state->mv_grid[vy + mvb_sz] + vx;
    }
  }
  for (ci = 0; ci < ncns; ci++) {
    int ref = cneighbors[ci]->ref;
    OD_ASSERT(ref < 2);
    hist[ref]++;
    if (hist[ref] > max_count) {
      max_ref = ref;
      max_count = hist[ref];
    }
  }
  return max_ref;
}

/*Gets the predictor for a given MV node at the given MV resolution.*/
int od_state_get_predictor(od_state *state,
 int pred[2], int vx, int vy, int level, int mv_res, int ref) {
  static const od_mv_grid_pt ZERO_GRID_PT = { {0, 0}, 1, OD_FRAME_PREV};
  const od_mv_grid_pt *cneighbors[4];
  int a[4][2];
  int equal_mvs;
  int mvb_sz;
  int ncns;
  int ci;
  ncns = 4;
  mvb_sz = 1 << ((OD_MC_LEVEL_MAX - level) >> 1);
  if (level == 0) {
    if (vy >= mvb_sz) {
      cneighbors[0] = vx >= mvb_sz ?
       state->mv_grid[vy - mvb_sz] + vx - mvb_sz : &ZERO_GRID_PT;
      cneighbors[1] = state->mv_grid[vy - mvb_sz] + vx;
      cneighbors[2] = vx + mvb_sz <= state->nhmvbs ?
       state->mv_grid[vy - mvb_sz] + vx + mvb_sz : &ZERO_GRID_PT;
    }
    else cneighbors[2] = cneighbors[1] = cneighbors[0] = &ZERO_GRID_PT;
    cneighbors[3] = vx >= mvb_sz ?
     state->mv_grid[vy] + vx - mvb_sz : &ZERO_GRID_PT;
  }
  else {
    if (level & 1) {
      cneighbors[0] = state->mv_grid[vy - mvb_sz] + vx - mvb_sz;
      cneighbors[1] = state->mv_grid[vy - mvb_sz] + vx + mvb_sz;
      cneighbors[2] = state->mv_grid[vy + mvb_sz] + vx - mvb_sz;
      cneighbors[3] = state->mv_grid[vy + mvb_sz] + vx + mvb_sz;
    }
    else {
      cneighbors[0] = vy >= mvb_sz ?
       state->mv_grid[vy - mvb_sz] + vx : &ZERO_GRID_PT;
      cneighbors[1] = vx >= mvb_sz ?
       state->mv_grid[vy] + vx - mvb_sz : &ZERO_GRID_PT;
      /*NOTE: Only one of these candidates can be excluded at a time, so
         there will always be at least 3.*/
      if (vx > 0 && vx + mvb_sz > ((vx + OD_MVB_MASK) & ~OD_MVB_MASK)) ncns--;
      else cneighbors[2] = state->mv_grid[vy] + vx + mvb_sz;
      if (vy > 0 && vy + mvb_sz > ((vy + OD_MVB_MASK) & ~OD_MVB_MASK)) ncns--;
      else cneighbors[ncns - 1] = state->mv_grid[vy + mvb_sz] + vx;
    }
  }
  for (ci = 0; ci < ncns; ci++) {
    a[ci][0] = cneighbors[ci]->mv[0];
    a[ci][1] = cneighbors[ci]->mv[1];
#if defined(OD_ENABLE_LOGGING)
    if (!cneighbors[ci]->valid) {
      OD_LOG((OD_LOG_MOTION_ESTIMATION, OD_LOG_ERR,
       "Failure in MV pred: predictor (%i, %i) for (%i, %i) is not valid",
	     a[ci][0], a[ci][1], vx, vy));
    }
#endif
    OD_ASSERT(cneighbors[ci]->valid);
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
    OD_ASSERT(ncns == 3);
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
  equal_mvs = 0;
  for (ci = 0; ci < ncns; ci++) {
    if (pred[0] == OD_DIV_POW2_RE(cneighbors[ci]->mv[0], mv_res) &&
     pred[1] == OD_DIV_POW2_RE(cneighbors[ci]->mv[1], mv_res)) {
      equal_mvs++;
    }
  }
  return equal_mvs;
}

int od_mv_split_flag_ctx(od_mv_grid_pt **grid, int vx, int vy,int level) {
  od_mv_grid_pt *v1;
  od_mv_grid_pt *v2;
  od_mv_grid_pt *v3;
  int split1;
  int split2;
  int same1;
  int same2;
  int mvb_sz;
  mvb_sz = 1 << ((OD_MC_LEVEL_MAX - level) >> 1);
  if (level & 1) {
    v1 = grid[vy - mvb_sz] + vx + mvb_sz;
    v2 = grid[vy + mvb_sz] + vx + mvb_sz;
    v3 = grid[vy + mvb_sz] + vx - mvb_sz;
  }
  else {
    v1 = vy >= mvb_sz ? grid[vy - mvb_sz] + vx : NULL;
    v2 = vx >= mvb_sz ? grid[vy] + vx - mvb_sz : NULL;
    v3 = vx & mvb_sz ? grid[vy] + vx + mvb_sz : grid[vy + mvb_sz] + vx;
  }
  split1 = vx >= 2*mvb_sz ? grid[vy][vx - 2*mvb_sz].valid : 0;
  split2 = vy >= 2*mvb_sz ? grid[vy - 2*mvb_sz][vx].valid : 0;
  same1 = v1 != NULL && v2 != NULL
   && (v1->mv[0] == v2->mv[0]) && (v1->mv[1] == v2->mv[1]);
  same2 = v2 != NULL
   && (v2->mv[0] == v3->mv[0]) && (v2->mv[1] == v3->mv[1]);
  return 3*(split1 + split2) + same1 + same2;
}

uint16_t *od_mv_split_flag_cdf(od_state *state,
 int vx, int vy, int level) {
  int ctx;
  ctx = od_mv_split_flag_ctx(state->mv_grid, vx, vy, level);
  OD_ASSERT(0 < level && level <= OD_MC_LEVEL_MAX);
  return state->adapt.split_flag_cdf[level - 1][ctx];
}
