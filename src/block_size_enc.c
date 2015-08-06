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
#define OD_CG4 (15.943f/6)
#define OD_CG8 (16.7836f/6)
#define OD_CG16 (16.9986f/6)
#define OD_CG32 (17.1f/6)

/* Weighting of the 8x8 masking compared to 4x4 */
#define PSY8_FUDGE (.5f)

/* This advanced macro computes the product of x by itself, otherwise known as
   raising to the power of two, squaring, or inverse square-root. It can be
   used to square integers, but not circles. */
#define OD_SQUARE(x) ((int)(x)*(int)(x))

/* Compute statistics used to determine how to split a superblock.
 * @param [in]      img    Image on which to compute the statistics, a margin
 *                          of OD_MAX_OVERLAP pixels is required around the
 *                          superblock
 * @param [in]      stride Image stride
 * @param [out]     stats  Computed statistics
 */
static void od_compute_stats(const signed char *img, int stride,
 od_superblock_stats *stats) {
  const signed char *x;
  int i;
  int j;
  int off8;
  int32_t (*Sx2)[OD_SIZE2_SUMS];
  int32_t (*Sxx2)[OD_SIZE2_SUMS];
  int32_t (*Sx4)[OD_SIZE4_SUMS];
  int32_t (*Sxx4)[OD_SIZE4_SUMS];
  int32_t (*Sx8)[OD_SIZE8_SUMS];
  int32_t (*Sxx8)[OD_SIZE8_SUMS];
  int32_t (*Var4)[OD_SIZE4_SUMS];
  int32_t (*invVar4)[OD_SIZE4_SUMS];
  int32_t (*Var8)[OD_SIZE8_SUMS];
  int32_t (*invVar8)[OD_SIZE8_SUMS];

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

  x = img - OD_BLOCK_OFFSET(stride);
  for (i = 0; i < OD_SIZE2_SUMS; i++) {
    for (j = 0; j < OD_SIZE2_SUMS; j++) {
      Sx2[i][j] = x[2*j] + x[2*j + 1] + x[2*j + stride] + x[2*j+stride + 1];
      Sxx2[i][j] = OD_SQUARE(x[2*j]) + OD_SQUARE(x[2*j + 1])
       + OD_SQUARE(x[2*j + stride]) + OD_SQUARE(x[2*j + stride + 1]);
    }
    x += 2*stride;
  }
  for (i = 0; i < OD_SIZE4_SUMS; i++) {
    for (j = 0; j < OD_SIZE4_SUMS; j++) {
      Sx4[i][j] = Sx2[i][j] + Sx2[i][j + 1]
       + Sx2[i + 1][j] + Sx2[i + 1][j + 1];
      Sxx4[i][j] = Sxx2[i][j] + Sxx2[i][j + 1]
       + Sxx2[i + 1][j] + Sxx2[i + 1][j + 1];
    }
  }
  off8 = OD_MAX_OVERLAP - 2*OD_MAX_OVERLAP_8;
  OD_ASSERT(off8 >= 0);
  for (i = 0; i < OD_SIZE8_SUMS; i++) {
    for (j = 0; j < OD_SIZE8_SUMS; j++) {
      Sx8[i][j] = Sx4[2*i + off8][2*j+off8] + Sx4[2*i + off8][2*j + 2 + off8]
       + Sx4[2*i + 2 + off8][2*j + off8] + Sx4[2*i + 2 + off8][2*j + 2 + off8];
      Sxx8[i][j] = Sxx4[2*i + off8][2*j + off8]
       + Sxx4[2*i + off8][2*j + 2 + off8]
       + Sxx4[2*i + 2 + off8][2*j + off8]
       + Sxx4[2*i + 2 + off8][2*j + 2 + off8];
    }
  }
  for (i = 0; i < OD_SIZE4_SUMS; i++) {
    for (j = 0; j < OD_SIZE4_SUMS; j++) {
      int32_t var_floor;
      Var4[i][j] = (Sxx4[i][j] - (OD_SQUARE(Sx4[i][j]) >> 4)) >> 5;
      var_floor = 4 + ((Sx4[i][j]+(128<<4)) >> 8);
      if (Var4[i][j] < var_floor) Var4[i][j] = var_floor;
      invVar4[i][j] = 16384/Var4[i][j];
    }
  }
  for (i = 0; i < OD_SIZE8_SUMS; i++) {
    for (j = 0; j < OD_SIZE8_SUMS; j++) {
      int32_t var_floor;
      Var8[i][j] = (Sxx8[i][j] - (OD_SQUARE(Sx8[i][j]) >> 6)) >> 5;
      var_floor = 4 + ((Sx8[i][j]+(128<<6)) >> 8);
      if (Var8[i][j] < var_floor) Var8[i][j] = var_floor;
      invVar8[i][j] = 16384/Var8[i][j];
    }
  }
}

/* Number of overlapping 4x4 blocks along one direction of the block. */
static int od_count_overlapping4x4(int bsize) {
  int non_overlapped_count_4x4;
  non_overlapped_count_4x4 = 1 << bsize;
  /*Each overlapping 4x4 block starts at a 2x2 offset, so along one direction
     there is one overlapped 4x4 block between each pair of non-overlapped 4x4
     blocks.*/
  return 2*non_overlapped_count_4x4 - 1;
}

/* Margin in units of 2x2 blocks used for the masking computations in
 *  `od_noise_var4x4` and `od_psy_var4x4`.
 * The margin must be less than or equal to OD_MAX_OVERLAP.
 * Because our transform is lapped, we need to take into account some pixels
 *  outside the current block for the masking computations, the amount of
 *  lapping used in the transform depends on the size of the block and of the
 *  neighboring blocks, since we can't know the size of all the neighboring
 *  blocks, we use an approximation.
 */
static unsigned int od_overlap_var4x4[OD_BLOCK_SIZES] = { 1, 1, 2, 3 };

/* Calculate the noise of a block based on the variances of overlapping 4x4
 *  blocks.
 * Let X represent a 2x2 block, suppose we want to find the noise of a 4x4
 *  block:
 *   +--+
 *   |XX|
 *   |XX|
 *   +--+
 *
 * We first use `od_count_overlapping4x4` to find the amount of overlapping 4x4
 *  blocks which are part of this block, for a 4x4 it's 1.
 * We then use `od_overlap_var4x4` to find the amount of 2x2 blocks around the
 *  block that we will take into account (because our transform is lapped), for
 *  a 4x4 it's 1:
 *   +----+
 *   |XXXX|
 *   |XXXX|
 *   |XXXX|
 *   |XXXX|
 *   +----+
 *
 * Then the noise is defined as the average variance of the 4x4 blocks at every
 *  2x2 offset:
 *   +----+  +----+  +----+  +----+  +----+  +----+  +----+  +----+  +----+
 *   |XX  |  | XX |  |  XX|  |    |  |    |  |    |  |    |  |    |  |    |
 *   |XX  |  | XX |  |  XX|  |XX  |  | XX |  |  XX|  |    |  |    |  |    |
 *   |    |  |    |  |    |  |XX  |  | XX |  |  XX|  |XX  |  | XX |  |  XX|
 *   |    |  |    |  |    |  |    |  |    |  |    |  |XX  |  | XX |  |  XX|
 *   +----+  +----+  +----+  +----+  +----+  +----+  +----+  +----+  +----+
 *
 * @param img_stats Residual if we're in inter or the image itself if we're in
 *                   intra
 * @param bsize     Size of the block (see `od_block_size`)
 * @param y         y offset of the block inside the 32x32 superblock
 * @param x         x offset of the block inside the 32x32 superblock
 * @retval block noise
 */
static int od_noise_var4x4(od_superblock_stats *img_stats,
 int bsize, int y, int x) {
  int i;
  int j;
  int length;
  int overlap;
  int sum_var;
  int count;
  OD_ASSERT(y % 4 == 0);
  OD_ASSERT(x % 4 == 0);
  length = od_count_overlapping4x4(bsize);
  overlap = od_overlap_var4x4[bsize];
  count = length + 2*overlap;
  sum_var = 0;
  for (i = -overlap; i < length + overlap; i++) {
    for (j = -overlap; j < length + overlap; j++) {
      sum_var += img_stats->Var4
       [OD_MAX_OVERLAP + y/2 + i][OD_MAX_OVERLAP + x/2 + j];
    }
  }
  return sum_var/(count*count);
}

/* Masking based on 4x4 variances assuming that noise is spatially flat and
 * proportional to the amount of residual we have to code (or the image itself
 * if we're in intra).
 *
 * @param psy_stats Image on which to compute the psy model (should not be a
 *                   residual)
 * @param bsize     Size of the block (see `od_block_size`)
 * @param y         y offset of the block inside the 32x32 superblock
 * @param x         x offset of the block inside the 32x32 superblock
 * @param noise     Noise of the block
 */
static float od_psy_var4x4(od_superblock_stats *psy_stats,
 int bsize, int y, int x, int noise) {
  int i;
  int j;
  int length;
  int overlap;
  float psy;
  int count;
  OD_ASSERT(y % 4 == 0);
  OD_ASSERT(x % 4 == 0);
  length = od_count_overlapping4x4(bsize);
  overlap = od_overlap_var4x4[bsize];
  count = length + 2*overlap;
  psy = 0;
  for (i = -overlap; i < length + overlap; i++) {
    for (j = -overlap; j < length + overlap; j++) {
      psy += OD_LOG2(1 + noise*psy_stats->invVar4
       [OD_MAX_OVERLAP + y/2 + i][OD_MAX_OVERLAP + x/2 + j]/16384.f);
    }
  }
  return OD_MAXF(psy/(count*count) - 1.f, 0);
}

/* Number of overlapping 8x8 blocks along one direction of the block. */
static int od_count_overlapping8x8(int bsize) {
  int non_overlapped_count_8x8;
  OD_ASSERT(bsize >= OD_BLOCK_16X16);
  non_overlapped_count_8x8 = 1 << (bsize - 1);
  /*Each overlapping 8x8 block starts at a 4x4 offset, so along one direction
     there is one overlapped 8x8 block between each pair of non-overlapped 8x8
     blocks.*/
  return 2*non_overlapped_count_8x8 - 1;
}

/* Margin in units of 4x4 blocks used for the masking computations in
 * `od_noise_var8x8` and `od_psy_var8x8`.
 * The margin must be less than or equal to OD_MAX_OVERLAP_8.
 */
static unsigned int od_overlap_var8x8[OD_BLOCK_SIZES] = { 0, 0, 1, 1 };

/* Same as `od_noise_var4x4` but using overlapping 8x8 blocks. */
static int od_noise_var8x8(od_superblock_stats *img_stats,
 int bsize, int y, int x) {
  int i;
  int j;
  int length;
  int overlap;
  int sum_var;
  int count;
  OD_ASSERT(y % 8 == 0);
  OD_ASSERT(x % 8 == 0);
  OD_ASSERT(bsize >= OD_BLOCK_16X16);
  length = od_count_overlapping8x8(bsize);
  overlap = od_overlap_var8x8[bsize];
  count = length + 2*overlap;
  sum_var = 0;
  for (i = -overlap; i < length + overlap; i++) {
    for (j = -overlap; j < length + overlap; j++) {
      sum_var += img_stats->Var8
       [OD_MAX_OVERLAP_8 + y/4 + i][OD_MAX_OVERLAP_8 + x/4 + j];
    }
  }
  return sum_var/(count*count);
}

/* Same as `od_psy_var4x4` but using overlapping 8x8 blocks. */
static float od_psy_var8x8(od_superblock_stats *psy_stats,
 int bsize, int y, int x, int noise) {
  int i;
  int j;
  int length;
  int overlap;
  float psy;
  int count;
  OD_ASSERT(y % 8 == 0);
  OD_ASSERT(x % 8 == 0);
  OD_ASSERT(bsize >= OD_BLOCK_16X16);
  length = od_count_overlapping8x8(bsize);
  overlap = od_overlap_var8x8[bsize];
  count = length + 2*overlap;
  psy = 0;
  for (i = -overlap; i < length + overlap; i++) {
    for (j = -overlap; j < length + overlap; j++) {
      psy += OD_LOG2(1 + noise*psy_stats->invVar8
       [OD_MAX_OVERLAP_8 + y/4 + i][OD_MAX_OVERLAP_8 + x/4 + j]/16384.f);
    }
  }
  return OD_MAXF(psy/(count*count) - 1.f, 0);
}

/* This function decides how to split a 32x32 superblock based on a simple
 * activity masking model. The masking at any given point is assumed to be
 * proportional to the local variance. The decision is made using a simple
 * dynamic programming algorithm, working from 8x8 decisions up to 32x32.
 * @param [scratch] bs          Scratch space for computation
 * @param [in]      psy_img     Image on which to compute the psy model
 *                               (should not be a residual)
 * @param [in]      stride      Image stride
 * @param [in]      pred        Prediction input (NULL means no prediction
 *                               available)
 * @param [in]      pred_stride Prediction input stride
 * @param [out]     bsize       Decision for each 8x8 block in the image
 *                               (see OD_BLOCK_* macros in block_size.h for
 *                               possible values)
 * @param [in]      q           Quality tuning parameter
 */
void od_split_superblock(od_block_size_comp *bs,
 const unsigned char *psy_img, int stride,
 const unsigned char *pred, int pred_stride,
  int bsize[OD_BSIZE_GRID][OD_BSIZE_GRID], int q) {
  int i;
  int j;
  /* Tuning parameter for block decision (higher values results in smaller
      blocks) */
  double psy_lambda;
  const unsigned char *x0;
  double cg4;
  double cg8;
  x0 = psy_img - OD_BLOCK_OFFSET(stride);
  /* The passed in q value is now a quantizer with the same scaling as
     the coefficients. */
  psy_lambda = q ? 6*sqrt((double)(1<<OD_COEFF_SHIFT)/q) : 6;
  for (i = 0; i < 2*OD_SIZE2_SUMS; i++) {
    for (j = 0; j < 2*OD_SIZE2_SUMS; j++) {
      bs->res[i][j] = (int)x0[i*stride + j] - 128;
    }
  }
  cg4 = OD_CG4;
  cg8 = OD_CG8;
  od_compute_stats(&bs->res[2*OD_MAX_OVERLAP][2*OD_MAX_OVERLAP],
   2*OD_SIZE2_SUMS, &bs->psy_stats);
  if (psy_img == pred || pred == NULL) {
    OD_COPY(&bs->img_stats, &bs->psy_stats, 1);
  }
  else {
    const unsigned char *p0;
    cg4 -= .01*OD_MAXI((q >> OD_COEFF_SHIFT) - 40, 0);
    cg8 -= .005*OD_MAXI(((q >> OD_COEFF_SHIFT) - 40), 0);
    p0 = pred - OD_BLOCK_OFFSET(pred_stride);
    for (i = 0; i < 2*OD_SIZE2_SUMS; i++) {
      for (j = 0; j < 2*OD_SIZE2_SUMS; j++) {
        bs->res[i][j] = OD_CLAMPI(-128, (int)x0[i*stride + j]
         - (int)p0[i*pred_stride + j], 127);
      }
    }
    od_compute_stats(&bs->res[2*OD_MAX_OVERLAP][2*OD_MAX_OVERLAP],
     2*OD_SIZE2_SUMS, &bs->img_stats);
  }
  /* Compute 4x4 masking */
  for (i = 0; i < 8; i++) {
    for (j = 0; j < 8; j++) {
      bs->noise4_4[i][j] = od_noise_var4x4(&bs->img_stats,
       OD_BLOCK_4X4, 4*i, 4*j);
      bs->psy4[i][j] = od_psy_var4x4(&bs->psy_stats,
       OD_BLOCK_4X4, 4*i, 4*j, bs->noise4_4[i][j]);
    }
  }
  /* Compute 8x8 masking and make 4x4 vs 8x8 decision */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      float gain4, gain8;
      float psy4_avg;
      bs->noise4_8[i][j] = od_noise_var4x4(&bs->img_stats,
        OD_BLOCK_8X8, 8*i, 8*j);
      bs->psy8[i][j] = od_psy_var4x4(&bs->psy_stats,
       OD_BLOCK_8X8, 8*i, 8*j, bs->noise4_8[i][j]);
      psy4_avg = .25f*(bs->psy4[2*i][2*j] + bs->psy4[2*i][2*j + 1]
       + bs->psy4[2*i + 1][2*j] + bs->psy4[2*i + 1][2*j + 1]);
      gain4 = cg4 - psy_lambda*(psy4_avg);
      gain8 = cg8 - psy_lambda*(bs->psy8[i][j]);
      if (gain8 >= gain4) {
        bsize[i][j] = OD_BLOCK_8X8;
        bs->dec_gain8[i][j] = gain8;
      }
      else {
        bsize[i][j] = OD_BLOCK_4X4;
        bs->dec_gain8[i][j] = gain4;
      }
    }
  }
#if (OD_LIMIT_BSIZE_MAX >= OD_BLOCK_16X16) && (OD_LIMIT_BSIZE_MIN <= OD_BLOCK_16X16)
  /* Compute 16x16 masking and make 4x4/8x8 vs 16x16 decision */
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      float gain16;
      float gain8_avg;
      bs->noise4_16[i][j] = od_noise_var4x4(&bs->img_stats,
       OD_BLOCK_16X16, 16*i, 16*j);
      bs->psy16[i][j] = od_psy_var4x4(&bs->psy_stats,
       OD_BLOCK_16X16, 16*i, 16*j, bs->noise4_16[i][j]);
      bs->noise8_16[i][j] = od_noise_var8x8(&bs->img_stats,
       OD_BLOCK_16X16, 16*i, 16*j);
      bs->psy16[i][j] = OD_MAXF(bs->psy16[i][j], PSY8_FUDGE*
       od_psy_var8x8(&bs->psy_stats, OD_BLOCK_16X16, 16*i, 16*j,
       bs->noise8_16[i][j]));
      gain8_avg = .25*(bs->dec_gain8[2*i][2*j] + bs->dec_gain8[2*i][2*j + 1]
       + bs->dec_gain8[2*i + 1][2*j] + bs->dec_gain8[2*i + 1][2*j + 1]);
      gain16 = OD_CG16 - psy_lambda*(bs->psy16[i][j]);
      if (gain16 >= gain8_avg) {
        bsize[2*i][2*j] = bsize[2*i][2*j + 1]
         = bsize[2*i + 1][2*j] = bsize[2*i + 1][2*j + 1] = OD_BLOCK_16X16;
        bs->dec_gain16[i][j] = gain16;
      }
      else {
        bs->dec_gain16[i][j] = gain8_avg;
      }
    }
  }
#endif
#if OD_LIMIT_BSIZE_MAX >= OD_BLOCK_32X32
  /* Compute 32x32 masking and make final 4x4/8x8/16x16 vs 32x32 decision */
  {
    int k;
    int m;
    float gain32;
    float gain16_avg;
    bs->noise4_32 = od_noise_var4x4(&bs->img_stats, OD_BLOCK_32X32, 0, 0);
    bs->psy32 = od_psy_var4x4(&bs->psy_stats,
     OD_BLOCK_32X32, 0, 0, bs->noise4_32);
    bs->noise8_32 = od_noise_var8x8(&bs->img_stats, OD_BLOCK_32X32, 0, 0);
    bs->psy32 = OD_MAXF(bs->psy32, PSY8_FUDGE*
     od_psy_var8x8(&bs->psy_stats, OD_BLOCK_32X32, 0, 0, bs->noise8_32));
    gain16_avg = .25f*(bs->dec_gain16[0][0] + bs->dec_gain16[0][1]
     + bs->dec_gain16[1][0] + bs->dec_gain16[1][1]);
    gain32 = OD_CG32 - psy_lambda*(bs->psy32);
    if (gain32 >= gain16_avg) {
      for (k = 0; k < 4; k++)
        for (m = 0; m < 4; m++) bsize[k][m] = OD_BLOCK_32X32;
    }
  }
#endif
}
