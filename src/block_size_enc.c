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
#include "block_size_enc.h"
#include <string.h>
#include <stdlib.h>

/* Actual 2D coding gains of lapped transforms (the 32x32 one is made-up).
   We divide by 6 to get bits from dB. */
#define CG4 (15.943/6)
#define CG8 (16.7836/6)
#define CG16 (16.9986/6)
#define CG32 (17.1/6)

/* Weighting of the 8x8 masking compared to 4x4 */
#define PSY8_FUDGE (.5f)

/* This advanced macro computes the product of x by itself, otherwise known as
   raising to the power of two, squaring, or inverse square-root. It can be
   used to square integers, but not circles. */
#define SQUARE(x) ((int)(x)*(int)(x))

static void od_compute_stats(const signed char *img, int stride,
 BlockStats *stats) {
  const signed char *x;
  int i;
  int j;
  int off8;
  ogg_int32_t (*Sx2)[SIZE2_SUMS];
  ogg_int32_t (*Sxx2)[SIZE2_SUMS];
  ogg_int32_t (*Sx4)[SIZE4_SUMS];
  ogg_int32_t (*Sxx4)[SIZE4_SUMS];
  ogg_int32_t (*Sx8)[SIZE8_SUMS];
  ogg_int32_t (*Sxx8)[SIZE8_SUMS];
  ogg_int32_t (*Var4)[SIZE4_SUMS];
  ogg_int32_t (*invVar4)[SIZE4_SUMS];
  ogg_int32_t (*Var8)[SIZE8_SUMS];
  ogg_int32_t (*invVar8)[SIZE8_SUMS];

  Sx2 = stats->Sx2;
  Sxx2 = stats->Sxx2;
  Sx4 = stats->Sx4;
  Sxx4 = stats->Sxx4;
  Sx8 = stats->Sx8;
  Sxx8 = stats->Sxx8;
  Var4 = stats->Var4;
  invVar4 = stats->invVar4;
  Var8 = stats->Var8;
  invVar8 = stats->invVar8;

  x = img - BLOCK_OFFSET(stride);
  for (i = 0; i < SIZE2_SUMS; i++) {
    for (j = 0; j < SIZE2_SUMS; j++) {
      Sx2[i][j] = x[2*j] + x[2*j + 1] + x[2*j + stride] + x[2*j+stride + 1];
      Sxx2[i][j] = SQUARE(x[2*j]) + SQUARE(x[2*j + 1])
       + SQUARE(x[2*j + stride]) + SQUARE(x[2*j + stride + 1]);
    }
    x += 2*stride;
  }
  for (i = 0; i < SIZE4_SUMS; i++) {
    for (j = 0; j < SIZE4_SUMS; j++) {
      Sx4[i][j] = Sx2[i][j] + Sx2[i][j + 1]
       + Sx2[i + 1][j] + Sx2[i + 1][j + 1];
      Sxx4[i][j] = Sxx2[i][j] + Sxx2[i][j + 1]
       + Sxx2[i + 1][j] + Sxx2[i + 1][j + 1];
    }
  }
  off8 = OFF32 - 2*OFF8_32;
  OD_ASSERT(off8 >= 0);
  for (i = 0; i < SIZE8_SUMS; i++) {
    for (j = 0; j < SIZE8_SUMS; j++) {
      Sx8[i][j] = Sx4[2*i + off8][2*j+off8] + Sx4[2*i + off8][2*j + 2 + off8]
       + Sx4[2*i + 2 + off8][2*j + off8] + Sx4[2*i + 2 + off8][2*j + 2 + off8];
      Sxx8[i][j] = Sxx4[2*i + off8][2*j + off8] + Sxx4[2*i + off8][2*j + 2 + off8]
       + Sxx4[2*i + 2 + off8][2*j + off8] + Sxx4[2*i + 2 + off8][2*j + 2 + off8];
    }
  }
  for (i = 0; i < SIZE4_SUMS; i++) {
    for (j = 0; j < SIZE4_SUMS; j++) {
      ogg_int32_t var_floor;
      Var4[i][j] = (Sxx4[i][j] - (SQUARE(Sx4[i][j]) >> 4)) >> 5;
      var_floor = 4 + ((Sx4[i][j]+(128<<4)) >> 8);
      if (Var4[i][j] < var_floor) Var4[i][j] = var_floor;
      invVar4[i][j] = 16384/Var4[i][j];
    }
  }
  for (i = 0; i < SIZE8_SUMS; i++) {
    for (j = 0; j < SIZE8_SUMS; j++) {
      ogg_int32_t var_floor;
      Var8[i][j] = (Sxx8[i][j] - (SQUARE(Sx8[i][j]) >> 6)) >> 5;
      var_floor = 4 + ((Sx8[i][j]+(128<<6)) >> 8);
      if (Var8[i][j] < var_floor) Var8[i][j] = var_floor;
      invVar8[i][j] = 16384/Var8[i][j];
    }
  }
}

/* This function decides how to partition a 32x32 superblock based on a simple
 * activity masking model. The masking at any given point is assumed to be
 * proportional to the local variance. The decision is made using a simple
 * dynamic programming algorithm, working from 8x8 decisions up to 32x32.
 * @param [scratch] bs Sratch space for computation
 * @param [in]      psy_img Image on which to compute the psy model (should not be a residual)
 * @param [in]      img     Image on which to compute the noise model (i.e. what is to be coded)
 * @param [in]      stride  Image stride
 * @param [out]     dec     Decision for each 8x8 block in the image. 0=4x4, 1=8x8, 2=16x16, 3=32x32
 */
void process_block_size32(BlockSizeComp *bs, const unsigned char *psy_img,
 int stride, const unsigned char *pred, int pred_stride, int bsize[4][4],
 int q) {
  int i;
  int j;
  /* Tuning parameter for block decision (higher values results in smaller
      blocks) */
  double psy_lambda;
  const unsigned char *x0;
  x0 = psy_img - BLOCK_OFFSET(stride);
  /* The passed in q value is now a quantizer with the same scaling as
     the coefficients. */
  psy_lambda = q ? 6*sqrt((double)(1<<OD_COEFF_SHIFT)/q) : 6;
  for (i = 0; i < 2*SIZE2_SUMS; i++) {
    for (j = 0; j < 2*SIZE2_SUMS; j++) {
      bs->res[i][j] = (int)x0[i*stride + j] - 128;
    }
  }
  od_compute_stats(&bs->res[2*OFF32][2*OFF32], 2*SIZE2_SUMS, &bs->psy_stats);
  if (psy_img == pred || pred == NULL) {
    OD_COPY(&bs->img_stats, &bs->psy_stats, 1);
  }
  else {
    const unsigned char *p0;
    p0 = pred - BLOCK_OFFSET(pred_stride);
    for (i = 0; i < 2*SIZE2_SUMS; i++) {
      for (j = 0; j < 2*SIZE2_SUMS; j++) {
        bs->res[i][j] = OD_CLAMPI(-128, (int)x0[i*stride + j]
         - (int)p0[i*pred_stride + j], 127);
      }
    }
    od_compute_stats(&bs->res[2*OFF32][2*OFF32], 2*SIZE2_SUMS,
     &bs->img_stats);
  }
  /* Compute 4x4 masking */
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      int k;
      int m;
      ogg_int32_t sum_var = 0;
      float psy = 0;
      /* Masking based on 4x4 variances */
      for (k = 0; k < 3; k++) {
        for (m = 0; m < 3; m++) {
          sum_var +=
           bs->img_stats.Var4[2*i + k + OFF32 - 1][2*j + m + OFF32 - 1];
        }
      }
      bs->noise4_4[i][j] = sum_var/(3*3);
      for (k = 0; k < 3; k++) {
        for (m = 0; m < 3; m++) {
          psy += OD_LOG2(1 + bs->noise4_4[i][j]*bs->psy_stats.invVar4
           [2*i + k + OFF32 - 1][2*j + m + OFF32 - 1]/16384.);
        }
      }
      bs->psy4[i][j] = OD_MAXF(psy/(3*3) - 1., 0);
    }
  }
  /* Compute 8x8 masking and make 4x4 vs 8x8 decision */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      int k;
      int m;
      float gain4, gain8;
      float psy4_avg;
      ogg_int32_t sum_var = 0;
      float psy = 0;
      /* Masking based on 4x4 variances */
      for (k = 0; k < COUNT8; k++) {
        for (m = 0; m < COUNT8; m++) {
          sum_var +=
           bs->img_stats.Var4[4*i + k + OFF32 - OFF8][4*j + m + OFF32 - OFF8];
        }
      }
      bs->noise4_8[i][j] = sum_var/(COUNT8*COUNT8);
      for (k = 0; k < COUNT8; k++) {
        for (m = 0; m < COUNT8; m++) {
          psy += OD_LOG2(1 + bs->noise4_8[i][j]*bs->psy_stats.invVar4
           [4*i + k + OFF32 - OFF8]
           [4*j + m + OFF32 - OFF8]/16384.);
        }
      }
      bs->psy8[i][j] = OD_MAXF(psy/(COUNT8*COUNT8) - 1., 0);
      psy4_avg = .25*(bs->psy4[2*i][2*j] + bs->psy4[2*i][2*j + 1]
       + bs->psy4[2*i + 1][2*j] + bs->psy4[2*i + 1][2*j + 1]);
      gain4 = CG4 - psy_lambda*(psy4_avg);
      gain8 = CG8 - psy_lambda*(bs->psy8[i][j]);
      if (gain8 >= gain4) {
        bsize[i][j] = 1;
        bs->dec_gain8[i][j] = gain8;
      }
      else {
        bsize[i][j] = 0;
        bs->dec_gain8[i][j] = gain4;
      }
    }
  }
#if (OD_LIMIT_LOG_BSIZE_MAX >= 4) && (OD_LIMIT_LOG_BSIZE_MIN <= 4)
  /* Compute 16x16 masking and make 4x4/8x8 vs 16x16 decision */
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      int k;
      int m;
      float gain16;
      float gain8_avg;
      ogg_int32_t sum_var = 0;
      float psy = 0;
      float psy8 = 0;
      /* Masking based on 4x4 variances */
      for (k = 0; k < COUNT16; k++) {
        for (m = 0; m < COUNT16; m++) {
          sum_var += bs->img_stats.Var4[8*i + k + OFF32 - OFF16]
           [8*j + m + OFF32 - OFF16];
        }
      }
      bs->noise4_16[i][j] = sum_var/(COUNT16*COUNT16);
      for (k = 0; k < COUNT16; k++) {
        for (m = 0; m < COUNT16; m++) {
          psy += OD_LOG2(1 + bs->noise4_16[i][j]*bs->psy_stats.invVar4
           [8*i + k + OFF32 - OFF16][8*j + m + OFF32 - OFF16]/16384.);
        }
      }
      bs->psy16[i][j] = OD_MAXF(psy/(COUNT16*COUNT16) - 1., 0);
      sum_var = 0;
      /* Use 8x8 variances */
      for (k = 0; k < COUNT8_16; k++) {
        for (m = 0; m < COUNT8_16; m++) {
          sum_var += bs->img_stats.Var8[4*i + k + OFF8_32 - OFF8_16]
           [4*j + m + OFF8_32 - OFF8_16];
        }
      }
      bs->noise8_16[i][j] = sum_var/(COUNT8_16*COUNT8_16);
      for (k = 0; k < COUNT8_16; k++) {
        for (m = 0; m < COUNT8_16; m++) {
          psy8 += OD_LOG2(1 + bs->noise8_16[i][j]*bs->psy_stats.invVar8
           [4*i + k + OFF8_32 - OFF8_16][4*j + m + OFF8_32 - OFF8_16]/16384.);
        }
      }
      bs->psy16[i][j] = OD_MAXF(bs->psy16[i][j], PSY8_FUDGE*
       OD_MAXF(psy8/(COUNT8_16*COUNT8_16) - 1., 0));
      gain8_avg = .25*(bs->dec_gain8[2*i][2*j] + bs->dec_gain8[2*i][2*j + 1]
       + bs->dec_gain8[2*i + 1][2*j] + bs->dec_gain8[2*i + 1][2*j + 1]);
      gain16 = CG16 - psy_lambda*(bs->psy16[i][j]);
      if (gain16 >= gain8_avg) {
        bsize[2*i][2*j] = bsize[2*i][2*j + 1]
         = bsize[2*i + 1][2*j] = bsize[2*i + 1][2*j + 1] = 2;
        bs->dec_gain16[i][j] = gain16;
      }
      else {
        bs->dec_gain16[i][j] = gain8_avg;
      }
    }
  }
#endif
#if OD_LIMIT_LOG_BSIZE_MAX>=5
  /* Compute 32x32 masking and make final 4x4/8x8/16x16 vs 32x32 decision */
  {
    int k;
    int m;
    float gain32;
    float gain16_avg;
    ogg_int32_t sum_var = 0;
    float psy = 0;
    float psy8 = 0;
    /* Masking based on 4x4 variances */
    for (k = 0; k < COUNT32; k++) {
      for (m = 0; m < COUNT32; m++) {
        sum_var += bs->img_stats.Var4[k][m];
      }
    }
    bs->noise4_32 = sum_var/(COUNT32*COUNT32);
    for (k = 0; k < COUNT32; k++) {
      for (m = 0; m < COUNT32; m++) {
        psy += OD_LOG2(1 + bs->noise4_32*bs->psy_stats.invVar4[k][m]/16384.);
      }
    }
    bs->psy32 = OD_MAXF(psy/(COUNT32*COUNT32) - 1., 0);
    sum_var = 0;
    /* Use 8x8 variances */
    for (k = 0; k < COUNT8_32; k++) {
      for (m = 0; m < COUNT8_32; m++) {
        sum_var += bs->img_stats.Var8[k][m];
      }
    }
    bs->noise8_32 = sum_var/(COUNT8_32*COUNT8_32);
    for (k = 0; k < COUNT8_32; k++) {
      for (m = 0; m < COUNT8_32; m++) {
        psy8 += OD_LOG2(1 + bs->noise8_32*bs->psy_stats.invVar8[k][m]/16384.);
      }
    }
    bs->psy32 = OD_MAXF(bs->psy32,
     PSY8_FUDGE*OD_MAXF(psy8/(COUNT8_32*COUNT8_32) - 1., 0));
    gain16_avg = .25*(bs->dec_gain16[0][0] + bs->dec_gain16[0][1]
     + bs->dec_gain16[1][0] + bs->dec_gain16[1][1]);
    gain32 = CG32 - psy_lambda*(bs->psy32);
    if (gain32 >= gain16_avg) {
      for (k = 0; k < 4; k++)
        for (m = 0; m < 4; m++) bsize[k][m] = 3;
    }
  }
#endif
}

void od_block_size_encode(od_ec_enc *enc, od_adapt_ctx *adapt,
 const unsigned char *bsize, int stride) {
  int i, j;
  int min_size;
  int max_size;
  int range_id;
  int split16;
  int min_ctx_size;
  min_ctx_size = bsize[-1 - stride];
  for (i = 0; i < 4; i++) {
    min_ctx_size = OD_MINI(min_ctx_size, bsize[i*stride - 1]);
    min_ctx_size = OD_MINI(min_ctx_size, bsize[-stride + i]);
  }
  min_size = max_size = bsize[0];
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      min_size = OD_MINI(min_size, bsize[i*stride + j]);
      max_size = OD_MAXI(max_size, bsize[i*stride + j]);
    }
  }
  range_id = -1;
  for (i = 0; i < 7; i++) {
    if (od_range_ids[i][0] == min_size && od_range_ids[i][1] == max_size) {
      range_id = i;
      break;
    }
  }
  OD_ASSERT(range_id != -1);
  /* The first symbol we encode is the min/max range of the block sizes. We
     use the min neighbor block size as context. */
  od_encode_cdf_adapt(enc, range_id, adapt->bsize_range_cdf[min_ctx_size], 7,
   adapt->bsize_range_increment);
  split16 = 8*(bsize[0] == 2) + 4*(bsize[2] == 2) + 2*(bsize[2*stride] == 2)
   + (bsize[2 + 2*stride] == 2);
  /* Encode whether the 16x16 blocks is split only when the range allows for
     both 8x8 and 16x16.
     TODO: Add some context. */
  if (max_size >= 2 && min_size < 2) {
    od_encode_cdf_adapt(enc, split16, adapt->bsize16_cdf[min_size], 16,
     adapt->bsize16_increment);
  }
  /* Encode the 8x8 block splits only when the range allows both 4x4 and 8x8.*/
  if (min_size == 0 && max_size > 0) {
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        if (bsize[2*i*stride + 2*j] < 2) {
          int split8;
          /* TODO: Consider the fact that we can't have all 8x8 if the
             min is 4x4. */
          split8 = bsize[2*i*stride + 2*j] + 2*bsize[2*i*stride + 2*j + 1]
           + 4*bsize[(2*i + 1)*stride + 2*j] + 8*bsize[(2*i + 1)*stride + 2*j + 1];
          od_encode_cdf_adapt(enc, split8, adapt->bsize8_cdf, 16,
            adapt->bsize8_increment);
        }
      }
    }
  }
}
