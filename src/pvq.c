/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

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

#include "odintrin.h"
#include "pvq.h"
#include "partition.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "filter.h"

/*These tables were generated using compute_basis.c. If OD_FILT_SIZE is
   changed, they have to be regenerated.*/
static const double MAG4[] = {
  0.870774, 0.872037, 0.949493, 0.947936
};
static const double MAG8[] = {
  0.936496, 0.892830, 0.938452, 0.970087,
  0.974272, 0.967954, 0.974035, 0.990480
};
static const double MAG16[] = {
  0.968807, 0.940969, 0.947977, 0.957741,
  0.969762, 0.978644, 0.984885, 0.988009,
  0.987424, 0.985569, 0.984215, 0.984462,
  0.987205, 0.991415, 0.994985, 0.998237
};
static const double MAG32[] = {
  0.985068, 0.970006, 0.969893, 0.973192,
  0.973444, 0.975881, 0.979601, 0.981070,
  0.984989, 0.987520, 0.988830, 0.990983,
  0.992376, 0.992884, 0.993447, 0.993381,
  0.993712, 0.994060, 0.993294, 0.992392,
  0.991338, 0.992410, 0.992051, 0.993874,
  0.993488, 0.994162, 0.995318, 0.995925,
  0.997475, 0.999027, 0.998303, 1.001413,
};
static const double MAG64[] = {
  0.992453, 0.984930, 0.985137, 0.985029,
  0.985514, 0.985784, 0.986269, 0.986854,
  0.989932, 0.987780, 0.988269, 0.989175,
  0.989951, 0.990466, 0.991145, 0.991839,
  0.990773, 0.993191, 0.993618, 0.994221,
  0.994662, 0.995259, 0.995826, 0.995996,
  0.999070, 0.996624, 0.996835, 0.996948,
  0.997022, 0.996973, 0.996993, 0.996996,
  0.996871, 0.996828, 0.996598, 0.996688,
  0.996845, 0.996407, 0.996327, 0.996435,
  0.999173, 0.996216, 0.995981, 0.996173,
  0.996595, 0.996334, 0.996512, 0.996627,
  0.994976, 0.997113, 0.997248, 0.997548,
  0.997943, 0.998121, 0.998291, 0.998687,
  1.001696, 0.999133, 0.999315, 0.999621,
  0.999745, 0.999905, 0.999936, 1.000075
};

static const double MAG4_CHROMA_420[] = {
  0.870774, 0.872037, 0.949493, 0.947936
};
static const double MAG8_CHROMA_420[] = {
  0.936496, 0.892830, 0.938452, 0.970087,
  0.974272, 0.967954, 0.974035, 0.990480
};
static const double MAG16_CHROMA_420[] = {
  0.968807, 0.940969, 0.947977, 0.957741,
  0.969762, 0.978644, 0.984885, 0.988009,
  0.987424, 0.985569, 0.984215, 0.984462,
  0.987205, 0.991415, 0.994985, 0.998237
};
static const double MAG32_CHROMA_420[] = {
  0.985068, 0.970006, 0.969893, 0.973192,
  0.973444, 0.975881, 0.979601, 0.981070,
  0.984989, 0.987520, 0.988830, 0.990983,
  0.992376, 0.992884, 0.993447, 0.993381,
  0.993712, 0.994060, 0.993294, 0.992392,
  0.991338, 0.992410, 0.992051, 0.993874,
  0.993488, 0.994162, 0.995318, 0.995925,
  0.997475, 0.999027, 0.998303, 1.001413
};
static const double MAG64_CHROMA_420[] = {
  0.992453, 0.984930, 0.985137, 0.985029,
  0.985514, 0.985784, 0.986269, 0.986854,
  0.989932, 0.987780, 0.988269, 0.989175,
  0.989951, 0.990466, 0.991145, 0.991839,
  0.990773, 0.993191, 0.993618, 0.994221,
  0.994662, 0.995259, 0.995826, 0.995996,
  0.999070, 0.996624, 0.996835, 0.996948,
  0.997022, 0.996973, 0.996993, 0.996996,
  0.996871, 0.996828, 0.996598, 0.996688,
  0.996845, 0.996407, 0.996327, 0.996435,
  0.999173, 0.996216, 0.995981, 0.996173,
  0.996595, 0.996334, 0.996512, 0.996627,
  0.994976, 0.997113, 0.997248, 0.997548,
  0.997943, 0.998121, 0.998291, 0.998687,
  1.001696, 0.999133, 0.999315, 0.999621,
  0.999745, 0.999905, 0.999936, 1.000075
};

const double *OD_BASIS_MAG[2][OD_NBSIZES + 1] = {
  {
    MAG4, MAG8, MAG16, MAG32, MAG64
  },
  {
    MAG4_CHROMA_420, MAG8_CHROMA_420, MAG16_CHROMA_420, MAG32_CHROMA_420,
    MAG64_CHROMA_420
  }
};

/* Quantization matrices for 8x8. For other block sizes, we currently just do
   resampling. */
/* Flat quantization, i.e. optimize for PSNR. */
const int OD_QM8_Q4_FLAT[] = {
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16
};
# if 0
/* M1: MPEG2 matrix for inter (which has a dead zone). */
const int OD_QM8_Q4[] = {
  16, 17, 18, 19, 20, 21, 22, 23,
  17, 18, 19, 20, 21, 22, 23, 24,
  18, 19, 20, 21, 22, 23, 24, 25,
  19, 20, 21, 22, 23, 24, 26, 27,
  20, 21, 22, 23, 25, 26, 27, 28,
  21, 22, 23, 24, 26, 27, 28, 30,
  22, 23, 24, 26, 27, 28, 30, 31,
  23, 24, 25, 27, 28, 30, 31, 33};
# endif
# if 0
/* M2: MPEG2 matrix for intra (no dead zone). */
const int OD_QM8_Q4[] = {
  16, 16, 19, 22, 22, 26, 26, 27,
  16, 16, 22, 22, 26, 27, 27, 29,
  19, 22, 26, 26, 27, 29, 29, 35,
  22, 24, 27, 27, 29, 32, 34, 38,
  26, 27, 29, 29, 32, 35, 38, 46,
  27, 29, 34, 34, 35, 40, 46, 56,
  29, 34, 34, 37, 40, 48, 56, 69,
  34, 37, 38, 40, 48, 58, 69, 83
};
# endif
# if 0
/* M3: Taken from dump_psnrhvs. */
const int OD_QM8_Q4[] = {
  16, 16, 17, 20, 24, 29, 36, 42,
  16, 17, 17, 19, 22, 26, 31, 37,
  17, 17, 21, 23, 26, 30, 34, 40,
  20, 19, 23, 28, 31, 35, 39, 45,
  24, 22, 26, 31, 36, 41, 46, 51,
  29, 26, 30, 35, 41, 47, 52, 58,
  36, 31, 34, 39, 46, 52, 59, 66,
  42, 37, 40, 45, 51, 58, 66, 73
};
# endif
# if 1
/* M4: a compromise equal to .5*(M3 + .5*(M2+transpose(M2))) */
const int OD_QM8_Q4_HVS[] = {
  16, 16, 18, 21, 24, 28, 32, 36,
  16, 17, 20, 21, 24, 27, 31, 35,
  18, 20, 24, 25, 27, 31, 33, 38,
  21, 21, 25, 28, 30, 34, 37, 42,
  24, 24, 27, 30, 34, 38, 43, 49,
  28, 27, 31, 34, 38, 44, 50, 58,
  32, 31, 33, 37, 43, 50, 58, 68,
  36, 35, 38, 42, 49, 58, 68, 78
};
#endif

/* Constants for the beta parameter, which controls how activity masking is
   used.
   beta = 1 / (1 - alpha), so when beta is 1, alpha is 0 and activity
   masking is disabled. When beta is 1.5, activity masking is used. Note that
   activity masking is neither used for 4x4 blocks nor for chroma. */
#define OD_BETA(b) OD_QCONST32(b, OD_BETA_SHIFT)
static const od_val16 OD_PVQ_BETA4_LUMA[1] = {OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA8_LUMA[4] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA16_LUMA[7] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA32_LUMA[10] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA64_LUMA[13] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.)};

static const od_val16 OD_PVQ_BETA4_LUMA_MASKING[1] = {OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA8_LUMA_MASKING[4] = {OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5)};
static const od_val16 OD_PVQ_BETA16_LUMA_MASKING[7] = {OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5),
 OD_BETA(1.5)};
static const od_val16 OD_PVQ_BETA32_LUMA_MASKING[10] = {OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5)};
static const od_val16 OD_PVQ_BETA64_LUMA_MASKING[13] = {OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5), OD_BETA(1.5),
 OD_BETA(1.5), OD_BETA(1.5)};

static const od_val16 OD_PVQ_BETA4_CHROMA[1] = {OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA8_CHROMA[4] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA16_CHROMA[7] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA32_CHROMA[10] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.)};
static const od_val16 OD_PVQ_BETA64_CHROMA[13] = {OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.),
 OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.), OD_BETA(1.)};

const od_val16 *const OD_PVQ_BETA[2][OD_NPLANES_MAX][OD_NBSIZES + 1] = {
 {{OD_PVQ_BETA4_LUMA, OD_PVQ_BETA8_LUMA,
   OD_PVQ_BETA16_LUMA, OD_PVQ_BETA32_LUMA,
   OD_PVQ_BETA64_LUMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA}},
 {{OD_PVQ_BETA4_LUMA_MASKING, OD_PVQ_BETA8_LUMA_MASKING,
   OD_PVQ_BETA16_LUMA_MASKING, OD_PVQ_BETA32_LUMA_MASKING,
   OD_PVQ_BETA64_LUMA_MASKING},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA},
  {OD_PVQ_BETA4_CHROMA, OD_PVQ_BETA8_CHROMA,
   OD_PVQ_BETA16_CHROMA, OD_PVQ_BETA32_CHROMA,
   OD_PVQ_BETA64_CHROMA}}
};

void od_adapt_pvq_ctx_reset(od_pvq_adapt_ctx *state, int is_keyframe) {
  od_pvq_codeword_ctx *ctx;
  int i;
  int pli;
  int bs;
  ctx = &state->pvq_codeword_ctx;
  generic_model_init(&state->pvq_param_model[0]);
  generic_model_init(&state->pvq_param_model[1]);
  generic_model_init(&state->pvq_param_model[2]);
  for (i = 0; i < 2*OD_NBSIZES; i++) {
    ctx->pvq_adapt[4*i + OD_ADAPT_K_Q8] = 384;
    ctx->pvq_adapt[4*i + OD_ADAPT_SUM_EX_Q8] = 256;
    ctx->pvq_adapt[4*i + OD_ADAPT_COUNT_Q8] = 104;
    ctx->pvq_adapt[4*i + OD_ADAPT_COUNT_EX_Q8] = 128;
  }
  ctx->pvq_k1_increment = 128;
  OD_CDFS_INIT(ctx->pvq_k1_cdf, ctx->pvq_k1_increment);
  for (pli = 0; pli < OD_NPLANES_MAX; pli++) {
    for (bs = 0; bs < OD_NBSIZES; bs++)
    for (i = 0; i < PVQ_MAX_PARTITIONS; i++) {
      state->pvq_exg[pli][bs][i] = 2 << 16;
    }
  }
  for (i = 0; i < OD_NBSIZES*PVQ_MAX_PARTITIONS; i++) {
    state->pvq_ext[i] = is_keyframe ? 24576 : 2 << 16;
  }
  state->pvq_gaintheta_increment = 128;
  OD_CDFS_INIT(state->pvq_gaintheta_cdf, state->pvq_gaintheta_increment >> 2);
  state->pvq_skip_dir_increment = 128;
  OD_CDFS_INIT(state->pvq_skip_dir_cdf, state->pvq_skip_dir_increment >> 2);
  ctx->pvq_split_increment = 128;
  OD_CDFS_INIT(ctx->pvq_split_cdf, ctx->pvq_split_increment >> 1);
}

/* QMs are arranged from smallest to largest blocksizes, first for
   blocks with decimation=0, followed by blocks with decimation=1.*/
int od_qm_offset(int bs, int xydec)
{
    return xydec*OD_QM_STRIDE + OD_QM_OFFSET(bs);
}

#if defined(OD_FLOAT_PVQ)
#define OD_DEFAULT_MAG 1.0
#else
#define OD_DEFAULT_MAG OD_QM_SCALE
#endif

/* Initialize the quantization matrix with the magnitude compensation applied.
   We need to compensate for the magnitude because lapping causes some basis
   functions to be smaller, so they would end up being quantized too finely
   (the same error in the quantized domain would result in a smaller pixel
   domain error). */
void od_init_qm(int16_t *x, int16_t *x_inv, const int *qm) {
  int i;
  int j;
  int16_t y[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int16_t y_inv[OD_BSIZE_MAX*OD_BSIZE_MAX];
  int16_t *x1;
  int16_t *x1_inv;
  int off;
  int bs;
  int xydec;
  for (bs = 0; bs < OD_NBSIZES; bs++) {
    for (xydec = 0; xydec < 2; xydec++) {
      off = od_qm_offset(bs, xydec);
      x1 = x + off;
      x1_inv = x_inv + off;
      for (i = 0; i < 4 << bs; i++) {
        for (j = 0; j < 4 << bs; j++) {
          /*This will ultimately be clamped to fit in 16 bits.*/
          od_val32 mag;
          int16_t ytmp;
#if OD_DEBLOCKING || OD_DISABLE_FILTER
          mag = OD_DEFAULT_MAG;
#else
          /*FIXME: Do this rounding when generating OD_BASIS_MAG.*/
          mag = OD_ROUND32(OD_DEFAULT_MAG*OD_BASIS_MAG[xydec][bs][i]*
           OD_BASIS_MAG[xydec][bs][j]);
#endif
          if (i == 0 && j == 0) {
            mag = OD_DEFAULT_MAG;
          }
          else {
#if defined(OD_FLOAT_PVQ)
            mag /= 0.0625*qm[(i << 1 >> bs)*8 + (j << 1 >> bs)];
#else
            int qmv;
            qmv = qm[(i << 1 >> bs)*8 + (j << 1 >> bs)];
            mag *= 16;
            mag = (mag + (qmv >> 1))/qmv;
#endif
            OD_ASSERT(mag > 0.0);
          }
          /*Convert to fit in 16 bits.*/
#if defined(OD_FLOAT_PVQ)
          y[i*(4 << bs) + j] = (int16_t)OD_MINI(OD_QM_SCALE_MAX,
           (int32_t)floor(.5 + mag*OD_QM_SCALE));
          y_inv[i*(4 << bs) + j] = (int16_t)floor(.5
           + OD_QM_SCALE*OD_QM_INV_SCALE/(double)y[i*(4 << bs) + j]);
#else
          y[i*(4 << bs) + j] = (int16_t)OD_MINI(OD_QM_SCALE_MAX, mag);
          ytmp = y[i*(4 << bs) + j];
          y_inv[i*(4 << bs) + j] = (int16_t)((OD_QM_SCALE*OD_QM_INV_SCALE
           + (ytmp >> 1))/ytmp);
#endif
        }
      }
      od_raster_to_coding_order_16(x1, 4 << bs, y, 4 << bs);
      od_raster_to_coding_order_16(x1_inv, 4 << bs, y_inv, 4 << bs);
    }
  }
}

/* Maps each possible size (n) in the split k-tokenizer to a different value.
   Possible values of n are:
   2, 3, 4, 7, 8, 14, 15, 16, 31, 32, 63, 64, 127, 128
   Since we don't care about the order (even in the bit-stream) the simplest
   ordering (implemented here) is:
   14, 2, 3, 4, 7, 8, 15, 16, 31, 32, 63, 64, 127, 128 */
int od_pvq_size_ctx(int n) {
  int logn;
  int odd;
  logn = OD_ILOG(n - 1);
  odd = n & 1;
  return 2*logn - 1 - odd - 7*(n == 14);
}

/* Maps a length n to a context for the (k=1, n<=16) coder, with a special
   case when n is the original length (orig_length=1) of the vector (i.e. we
   haven't split it yet). For orig_length=0, we use the same mapping as
   od_pvq_size_ctx() up to n=16. When orig_length=1, we map lengths
   7, 8, 14, 15 to contexts 8 to 11. */
int od_pvq_k1_ctx(int n, int orig_length) {
  if (orig_length) return 8 + 2*(n > 8) + (n & 1);
  else return od_pvq_size_ctx(n);
}

/* Indexing for the packed quantization matrices. */
int od_qm_get_index(int bs, int band) {
  /* The -band/3 term is due to the fact that we force corresponding horizontal
     and vertical bands to have the same quantization. */
  OD_ASSERT(bs >= 0 && bs < OD_NBSIZES);
  return bs*(bs + 1) + band - band/3;
}

#if !defined(OD_FLOAT_PVQ)
/*See celt/mathops.c in Opus and tools/cos_search.c.*/
static int16_t od_pvq_cos_pi_2(int16_t x)
{
  int16_t x2;
  x2 = OD_MULT16_16_Q15(x, x);
  return OD_MINI(32767, (1073758164 - x*x + x2*(-7654 + OD_MULT16_16_Q16(x2,
   16573 + OD_MULT16_16_Q16(-2529, x2)))) >> 15);
}
#endif

/*Approximates cos(x) for -pi < x < pi.
  Input is in OD_THETA_SCALE.*/
od_val16 od_pvq_cos(od_val32 x) {
#if defined(OD_FLOAT_PVQ)
  return cos(x);
#else
  /*Wrap x around by masking, since cos is periodic.*/
  x = x & 0x0001ffff;
  if (x > (1 << 16)) {
    x = (1 << 17) - x;
  }
  if (x & 0x00007fff) {
    if (x < (1 << 15)) {
       return od_pvq_cos_pi_2((int16_t)x);
    }
    else {
      return -od_pvq_cos_pi_2((int16_t)(65536 - x));
    }
  }
  else {
    if (x & 0x0000ffff) {
      return 0;
    }
    else if (x & 0x0001ffff) {
      return -32767;
    }
    else {
      return 32767;
    }
  }
#endif
}

/*Approximates sin(x) for 0 <= x < pi.
  Input is in OD_THETA_SCALE.*/
od_val16 od_pvq_sin(od_val32 x) {
#if defined(OD_FLOAT_PVQ)
  return sin(x);
#else
  return od_pvq_cos(32768 - x);
#endif
}

#if !defined(OD_FLOAT_PVQ)
/* Computes an upper-bound on the number of bits required to store the L2 norm
   of a vector (excluding sign). */
int od_vector_log_mag(const od_coeff *x, int n) {
  int i;
  int32_t sum;
  sum = 0;
  for (i = 0; i < n; i++) {
    int16_t tmp;
    tmp = x[i] >> 8;
    sum += tmp*(int32_t)tmp;
  }
  /* We add one full bit (instead of rounding OD_ILOG() up) for safety because
     the >> 8 above causes the sum to be slightly underestimated. */
  return 8 + 1 + OD_ILOG(n + sum)/2;
}
#endif

/** Computes Householder reflection that aligns the reference r to the
 *  dimension in r with the greatest absolute value. The reflection
 *  vector is returned in r.
 *
 * @param [in,out]  r      reference vector to be reflected, reflection
 *                         also returned in r
 * @param [in]      n      number of dimensions in r
 * @param [in]      gr     gain of reference vector
 * @param [out]     sign   sign of reflection
 * @return                 dimension number to which reflection aligns
 **/
int od_compute_householder(od_val16 *r, int n, od_val32 gr, int *sign,
 int shift) {
  int m;
  int i;
  int s;
  od_val16 maxr;
  OD_UNUSED(shift);
  /* Pick component with largest magnitude. Not strictly
   * necessary, but it helps numerical stability */
  m = 0;
  maxr = 0;
  for (i = 0; i < n; i++) {
    if (OD_ABS(r[i]) > maxr) {
      maxr = OD_ABS(r[i]);
      m = i;
    }
  }
  s = r[m] > 0 ? 1 : -1;
  /* This turns r into a Householder reflection vector that would reflect
   * the original r[] to e_m */
  r[m] += OD_SHR_ROUND(gr*s, shift);
  *sign = s;
  return m;
}

#if !defined(OD_FLOAT_PVQ)
#define OD_RCP_INSHIFT 15
#define OD_RCP_OUTSHIFT 14
static od_val16 od_rcp(od_val16 x)
{
  int i;
  od_val16 n;
  od_val16 r;
  i = OD_ILOG(x) - 1;
  /*n is Q15 with range [0,1).*/
  n = OD_VSHR_ROUND(x, i - OD_RCP_INSHIFT) - (1 << OD_RCP_INSHIFT);
  /*Start with a linear approximation:
    r = 1.8823529411764706-0.9411764705882353*n.
    The coefficients and the result are Q14 in the range [15420,30840].*/
  r = 30840 + OD_MULT16_16_Q15(-15420, n);
  /*Perform two Newton iterations:
    r -= r*((r*n)-1.Q15)
       = r*((r*n)+(r-1.Q15)).*/
  r = r - OD_MULT16_16_Q15(r, (OD_MULT16_16_Q15(r, n) + r - 32768));
  /*We subtract an extra 1 in the second iteration to avoid overflow; it also
     neatly compensates for truncation error in the rest of the process.*/
  r = r - (1 + OD_MULT16_16_Q15(r, OD_MULT16_16_Q15(r, n) + r - 32768));
  /*r is now the Q15 solution to 2/(n+1), with a maximum relative error
     of 7.05346E-5, a (relative) RMSE of 2.14418E-5, and a peak absolute
     error of 1.24665/32768.*/
  return OD_VSHR_ROUND(r, i - OD_RCP_OUTSHIFT);
}
#endif

/** Applies Householder reflection from compute_householder(). The
 * reflection is its own inverse.
 *
 * @param [out]     out    reflected vector
 * @param [in]      x      vector to be reflected
 * @param [in]      r      reflection
 * @param [in]      n      number of dimensions in x,r
 */
void od_apply_householder(od_val16 *out, const od_val16 *x, const od_val16 *r,
 int n) {
  int i;
  od_val32 proj;
  od_val16 proj_1;
  od_val32 l2r;
#if !defined(OD_FLOAT_PVQ)
  od_val16 proj_norm;
  od_val16 l2r_norm;
  od_val16 rcp;
  int proj_shift;
  int l2r_shift;
  int outshift;
#endif
  /*FIXME: Can we get l2r and/or l2r_shift from an earlier computation?*/
  l2r = 0;
  for (i = 0; i < n; i++) {
    l2r += OD_MULT16_16(r[i], r[i]);
  }
  /* Apply Householder reflection */
  proj = 0;
  for (i = 0; i < n; i++) {
    proj += OD_MULT16_16(r[i], x[i]);
  }
#if defined(OD_FLOAT_PVQ)
  proj_1 = proj*2./(1e-100 + l2r);
  for (i = 0; i < n; i++) {
    out[i] = x[i] - r[i]*proj_1;
  }
#else
  /*l2r_norm is [0.5, 1.0[ in Q15.*/
  l2r_shift = (OD_ILOG(l2r) - 1) - 14;
  l2r_norm = OD_VSHR_ROUND(l2r, l2r_shift);
  rcp = od_rcp(l2r_norm);
  proj_shift = (OD_ILOG(abs(proj)) - 1) - 14;
  /*proj_norm is [0.5, 1.0[ in Q15.*/
  proj_norm = OD_VSHR_ROUND(proj, proj_shift);
  proj_1 = OD_MULT16_16_Q15(proj_norm, rcp);
  /*The proj*2. in the float code becomes -1 in the final outshift.
    The sign of l2r_shift is positive since we're taking the reciprocal of
     l2r_norm and this is a right shift.*/
  outshift = OD_MINI(30, OD_RCP_OUTSHIFT - proj_shift - 1 + l2r_shift);
  if (outshift >= 0) {
    for (i = 0; i < n; i++) {
      int32_t tmp;
      tmp = OD_MULT16_16(r[i], proj_1);
      tmp = OD_SHR_ROUND(tmp, outshift);
      out[i] = x[i] - tmp;
    }
  }
  else {
    /*FIXME: Can we make this case impossible?
      Right now, if r[] is all zeros except for 1, 2, or 3 ones, and
       if x[] is all zeros except for large values at the same position as the
       ones in r[], then we can end up with a shift of -1.*/
    for (i = 0; i < n; i++) {
      int32_t tmp;
      tmp = OD_MULT16_16(r[i], proj_1);
      tmp = OD_SHL(tmp, -outshift);
      out[i] = x[i] - tmp;
    }
  }
#endif
}

#if !defined(OD_FLOAT_PVQ)
static od_val16 od_beta_rcp(od_val16 beta){
  if (beta == OD_BETA(1.))
    return OD_BETA(1.);
  else if (beta == OD_BETA(1.5))
    return OD_BETA(1./1.5);
  else {
    od_val16 rcp_beta;
    /*Shift by 1 less, transposing beta to range [.5, .75] and thus < 32768.*/
    rcp_beta = od_rcp(beta << (OD_RCP_INSHIFT - 1 - OD_BETA_SHIFT));
    return OD_SHR_ROUND(rcp_beta, OD_RCP_OUTSHIFT + 1 - OD_BETA_SHIFT);
  }
}

#define OD_EXP2_INSHIFT 15
#define OD_EXP2_FRACSHIFT 15
#define OD_EXP2_OUTSHIFT 15
static const int32_t OD_EXP2_C[5] = {32768, 22709, 7913, 1704, 443};
/*Output is [1.0, 2.0) in Q(OD_EXP2_FRACSHIFT).
  It does not include the integer offset, which is added in od_exp2 after the
   final shift).*/
static int32_t od_exp2_frac(int32_t x)
{
  return OD_MULT16_16_Q15(x, (OD_EXP2_C[1] + OD_MULT16_16_Q15(x,
   (OD_EXP2_C[2] + OD_MULT16_16_Q15(x, (OD_EXP2_C[3]
   + OD_MULT16_16_Q15(x, OD_EXP2_C[4])))))));
}

/** Base-2 exponential approximation (2^x) with Q15 input and output.*/
static int32_t od_exp2(int32_t x)
{
  int integer;
  int32_t frac;
  integer = x >> OD_EXP2_INSHIFT;
  if (integer > 14)
    return 0x7f000000;
  else if (integer < -15)
    return 0;
  frac = od_exp2_frac(x - OD_SHL(integer, OD_EXP2_INSHIFT));
  return OD_VSHR_ROUND(OD_EXP2_C[0] + frac, -integer) + 1;
}

#define OD_LOG2_INSHIFT 15
#define OD_LOG2_OUTSHIFT 15
#define OD_LOG2_INSCALE_1 (1./(1 << OD_LOG2_INSHIFT))
#define OD_LOG2_OUTSCALE (1 << OD_LOG2_OUTSHIFT)
static int16_t od_log2(int16_t x)
{
  return x + OD_MULT16_16_Q15(x, (14482 + OD_MULT16_16_Q15(x, (-23234
   + OD_MULT16_16_Q15(x, (13643 + OD_MULT16_16_Q15(x, (-6403
   + OD_MULT16_16_Q15(x, 1515)))))))));
}

static int32_t od_pow(int32_t x, od_val16 beta)
{
  int16_t t;
  int xshift;
  int log2_x;
  od_val32 logr;
  /*FIXME: this conditional is to avoid doing log2(0).*/
  if (x == 0)
    return 0;
  log2_x = (OD_ILOG(x) - 1);
  xshift = log2_x - OD_LOG2_INSHIFT;
  /*t should be in range [0.0, 1.0[ in Q(OD_LOG2_INSHIFT).*/
  t = OD_VSHR(x, xshift) - (1 << OD_LOG2_INSHIFT);
  /*log2(g/OD_COMPAND_SCALE) = log2(x) - OD_COMPAND_SHIFT in
     Q(OD_LOG2_OUTSHIFT).*/
  logr = od_log2(t) + (log2_x - OD_COMPAND_SHIFT)*OD_LOG2_OUTSCALE;
  logr = OD_MULT16_32_QBETA(beta, logr);
  return od_exp2(logr);
}
#endif

/** Gain companding: raises gain to the power 1/beta for activity masking.
 *
 * @param [in]  g     real (uncompanded) gain
 * @param [in]  q0    uncompanded quality parameter
 * @param [in]  beta  activity masking beta param (exponent)
 * @return            g^(1/beta)
 */
static od_val32 od_gain_compand(od_val32 g, int q0, od_val16 beta) {
#if defined(OD_FLOAT_PVQ)
  if (beta == 1) return OD_CGAIN_SCALE*g/(double)q0;
  else {
    return OD_CGAIN_SCALE*OD_COMPAND_SCALE*pow(g*OD_COMPAND_SCALE_1,
     1./beta)/(double)q0;
  }
#else
  if (beta == OD_BETA(1)) return (OD_CGAIN_SCALE*g + (q0 >> 1))/q0;
  else {
    int32_t expr;
    expr = od_pow(g, od_beta_rcp(beta));
    expr <<= OD_CGAIN_SHIFT + OD_COMPAND_SHIFT - OD_EXP2_OUTSHIFT;
    return (expr + (q0 >> 1))/q0;
  }
#endif
}

#if !defined(OD_FLOAT_PVQ)
#define OD_SQRT_INSHIFT 16
#define OD_SQRT_OUTSHIFT 15
static int16_t od_rsqrt_norm(int16_t x);

static int16_t od_sqrt_norm(int32_t x)
{
  OD_ASSERT(x < 65536);
  return OD_MINI(OD_SHR_ROUND(x*od_rsqrt_norm(x), OD_SQRT_OUTSHIFT), 32767);
}

static int16_t od_sqrt(int32_t x, int *sqrt_shift)
{
  int k;
  int s;
  int32_t t;
  if (x == 0) {
    *sqrt_shift = 0;
     return 0;
  }
  OD_ASSERT(x < (1 << 30));
  k = ((OD_ILOG(x) - 1) >> 1);
  /*t is x in the range [0.25, 1) in QINSHIFT, or x*2^(-s).
    Shift by log2(x) - log2(0.25*(1 << INSHIFT)) to ensure 0.25 lower bound.*/
  s = 2*k - (OD_SQRT_INSHIFT - 2);
  t = OD_VSHR(x, s);
  /*We want to express od_sqrt() in terms of od_sqrt_norm(), which is
     defined as (2^OUTSHIFT)*sqrt(t*(2^-INSHIFT)) with t=x*(2^-s).
    This simplifies to 2^(OUTSHIFT-(INSHIFT/2)-(s/2))*sqrt(x), so the caller
     needs to shift right by OUTSHIFT - INSHIFT/2 - s/2.*/
  *sqrt_shift = OD_SQRT_OUTSHIFT - ((s + OD_SQRT_INSHIFT) >> 1);
  return od_sqrt_norm(t);
}
#endif

/** Gain expanding: raises gain to the power beta for activity masking.
 *
 * @param [in]  cg    companded gain
 * @param [in]  q0    uncompanded quality parameter
 * @param [in]  beta  activity masking beta param (exponent)
 * @return            g^beta
 */
od_val32 od_gain_expand(od_val32 cg0, int q0, od_val16 beta) {
  if (beta == OD_BETA(1)) {
    /*The multiply fits into 28 bits because the expanded gain has a range from
       0 to 2^20.*/
    return OD_SHR_ROUND(cg0*q0, OD_CGAIN_SHIFT);
  }
  else if (beta == OD_BETA(1.5)) {
#if defined(OD_FLOAT_PVQ)
    double cg;
    cg = cg0*OD_CGAIN_SCALE_1;
    cg *= q0*OD_COMPAND_SCALE_1;
    return OD_COMPAND_SCALE*cg*sqrt(cg);
#else
    int32_t irt;
    int64_t tmp;
    int sqrt_inshift;
    int sqrt_outshift;
    /*cg0 is in Q(OD_CGAIN_SHIFT) and we need to divide it by
       2^OD_COMPAND_SHIFT.*/
    irt = od_sqrt(cg0*q0, &sqrt_outshift);
    sqrt_inshift = (OD_CGAIN_SHIFT + OD_COMPAND_SHIFT) >> 1;
    /*tmp is in Q(OD_CGAIN_SHIFT + OD_COMPAND_SHIFT).*/
    tmp = cg0*q0*(int64_t)irt;
    /*Expanded gain must be in Q(OD_COMPAND_SHIFT), thus OD_COMPAND_SHIFT is
       not included here.*/
    return OD_VSHR_ROUND(tmp, OD_CGAIN_SHIFT + sqrt_outshift + sqrt_inshift);
#endif
  }
  else {
#if defined(OD_FLOAT_PVQ)
    /*Expanded gain must be in Q(OD_COMPAND_SHIFT), hence the multiply by
       OD_COMPAND_SCALE.*/
    double cg;
    cg = cg0*OD_CGAIN_SCALE_1;
    return OD_COMPAND_SCALE*pow(cg*q0*OD_COMPAND_SCALE_1, beta);
#else
    int32_t expr;
    int32_t cg;
    cg = OD_SHR_ROUND(cg0*q0, OD_CGAIN_SHIFT);
    expr = od_pow(cg, beta);
    /*Expanded gain must be in Q(OD_COMPAND_SHIFT), hence the subtraction by
       OD_COMPAND_SHIFT.*/
    return OD_SHR_ROUND(expr, OD_EXP2_OUTSHIFT - OD_COMPAND_SHIFT);
#endif
  }
}

/** Computes the raw and quantized/companded gain of a given input
 * vector
 *
 * @param [in]      x      vector of input data
 * @param [in]      n      number of elements in vector x
 * @param [in]      q0     quantizer
 * @param [out]     g      raw gain
 * @param [in]      beta   activity masking beta param
 * @param [in]      bshift shift to be applied to raw gain
 * @return                 quantized/companded gain
 */
od_val32 od_pvq_compute_gain(const od_val16 *x, int n, int q0, od_val32 *g,
 od_val16 beta, int bshift) {
  int i;
  od_val32 acc;
#if !defined(OD_FLOAT_PVQ)
  od_val32 irt;
  int sqrt_shift;
#else
  OD_UNUSED(bshift);
#endif
  acc = 0;
  for (i = 0; i < n; i++) {
    acc += x[i]*(od_val32)x[i];
  }
#if defined(OD_FLOAT_PVQ)
  *g = sqrt(acc);
#else
  irt = od_sqrt(acc, &sqrt_shift);
  *g = OD_VSHR_ROUND(irt, sqrt_shift - bshift);
#endif
  /* Normalize gain by quantization step size and apply companding
     (if ACTIVITY != 1). */
  return od_gain_compand(*g, q0, beta);
}

/** Compute theta quantization range from quantized/companded gain
 *
 * @param [in]      qcg    quantized companded gain value
 * @param [in]      beta   activity masking beta param
 * @return                 max theta value
 */
int od_pvq_compute_max_theta(od_val32 qcg, od_val16 beta){
  /* Set angular resolution (in ra) to match the encoded gain */
#if defined(OD_FLOAT_PVQ)
  int ts = (int)floor(.5 + qcg*OD_CGAIN_SCALE_1*M_PI/(2*beta));
#else
  int ts = OD_SHR_ROUND(qcg*OD_MULT16_16_QBETA(OD_QCONST32(M_PI/2,
   OD_CGAIN_SHIFT), od_beta_rcp(beta)), OD_CGAIN_SHIFT*2);
#endif
  /* Special case for low gains -- will need to be tuned anyway */
  if (qcg < OD_QCONST32(1.4, OD_CGAIN_SHIFT)) ts = 1;
  return ts;
}

/** Decode quantized theta value from coded value
 *
 * @param [in]      t          quantized companded gain value
 * @param [in]      max_theta  maximum theta value
 * @return                     decoded theta value
 */
od_val32 od_pvq_compute_theta(int t, int max_theta) {
  if (max_theta != 0) {
#if defined(OD_FLOAT_PVQ)
    return OD_MINI(t, max_theta - 1)*.5*M_PI/max_theta;
#else
    return (OD_MAX_THETA_SCALE*OD_MINI(t, max_theta - 1)
     + (max_theta >> 1))/max_theta;
#endif
  }
  else return 0;
}

#define OD_SQRT_TBL_SHIFT (10)

#define OD_ITHETA_SHIFT 15
/** Compute the number of pulses used for PVQ encoding a vector from
 * available metrics (encode and decode side)
 *
 * @param [in]      qcg        quantized companded gain value
 * @param [in]      itheta     quantized PVQ error angle theta
 * @param [in]      theta      PVQ error angle theta
 * @param [in]      noref      indicates present or lack of reference
 *                             (prediction)
 * @param [in]      n          number of elements to be coded
 * @param [in]      beta       activity masking beta param
 * @param [in]      nodesync   do not use info that depends on the reference
 * @return                     number of pulses to use for coding
 */
int od_pvq_compute_k(od_val32 qcg, int itheta, od_val32 theta, int noref, int n,
 od_val16 beta, int nodesync) {
#if !defined(OD_FLOAT_PVQ)
  /*Lookup table for sqrt(n+3/2) and sqrt(n+2/2) in Q10.
    Real max values are 32792 and 32784, but clamped to stay within 16 bits.
    Update with tools/gen_sqrt_tbl if needed.*/
  static const od_val16 od_sqrt_table[2][13] = {
   {0, 0, 0, 0, 2290, 2985, 4222, 0, 8256, 0, 16416, 0, 32767},
   {0, 0, 0, 0, 2401, 3072, 4284, 0, 8287, 0, 16432, 0, 32767}};
#endif
  if (noref) {
    if (qcg == 0) return 0;
    if (n == 15 && qcg == OD_CGAIN_SCALE) {
      return 1;
    }
    else {
#if defined(OD_FLOAT_PVQ)
      return OD_MAXI(1, (int)floor(.5 + (qcg*OD_CGAIN_SCALE_1 - .2)*
       sqrt((n + 3)/2)/beta));
#else
      od_val16 rt;
      OD_ASSERT(OD_ILOG(n + 1) < 13);
      rt = od_sqrt_table[1][OD_ILOG(n + 1)];
      /*FIXME: get rid of 64-bit mul.*/
      return OD_MAXI(1, OD_SHR_ROUND((int64_t)((qcg
       - (int64_t)OD_QCONST32(.2, OD_CGAIN_SHIFT))*
       OD_MULT16_16_QBETA(od_beta_rcp(beta), rt)), OD_CGAIN_SHIFT
       + OD_SQRT_TBL_SHIFT));
#endif
    }
  }
  else {
    if (itheta == 0) return 0;
    /* Sets K according to gain and theta, based on the high-rate
       PVQ distortion curves (see PVQ document). Low-rate will have to be
       perceptually tuned anyway. We subtract 0.2 from the radius as an
       approximation for the fact that the coefficients aren't identically
       distributed within a band so at low gain the number of dimensions that
       are likely to have a pulse is less than n. */
    if (n == 15 && itheta == 1 && beta == OD_BETA(1)) return 1;
    if (nodesync) {
#if defined(OD_FLOAT_PVQ)
      return OD_MAXI(1, (int)floor(.5 + (itheta - .2)*sqrt((n + 2)/2)));
#else
      od_val16 rt;
      OD_ASSERT(OD_ILOG(n + 1) < 13);
      rt = od_sqrt_table[0][OD_ILOG(n + 1)];
      /*FIXME: get rid of 64-bit mul.*/
      return OD_MAXI(1, OD_VSHR_ROUND(((OD_SHL(itheta, OD_ITHETA_SHIFT)
       - OD_QCONST32(.2, OD_ITHETA_SHIFT)))*(int64_t)rt,
       OD_SQRT_TBL_SHIFT + OD_ITHETA_SHIFT));
#endif
    }
    else {
      return OD_MAXI(1, (int)floor(.5 + (qcg*OD_CGAIN_SCALE_1*
       od_pvq_sin(theta)*OD_TRIG_SCALE_1 - .2)*sqrt((n
       + 2)/2)/(beta*OD_BETA_SCALE_1)));
    }
  }
}

#if !defined(OD_FLOAT_PVQ)
#define OD_RSQRT_INSHIFT 16
#define OD_RSQRT_OUTSHIFT 14
/** Reciprocal sqrt approximation where the input is in the range [0.25,1) in
     Q16 and the output is in the range (1.0, 2.0] in Q14).
    Error is always within +/1 of round(1/sqrt(t))*/
static int16_t od_rsqrt_norm(int16_t t)
{
  int16_t n;
  int32_t r;
  int32_t r2;
  int32_t ry;
  int32_t y;
  int32_t ret;
  /* Range of n is [-16384,32767] ([-0.5,1) in Q15).*/
  n = t - 32768;
  OD_ASSERT(n >= -16384);
  /*Get a rough initial guess for the root.
    The optimal minimax quadratic approximation (using relative error) is
     r = 1.437799046117536+n*(-0.823394375837328+n*0.4096419668459485).
    Coefficients here, and the final result r, are Q14.*/
  r = (23565 + OD_MULT16_16_Q15(n, (-13481 + OD_MULT16_16_Q15(n, 6711))));
  /*We want y = t*r*r-1 in Q15, but t is 32-bit Q16 and r is Q14.
    We can compute the result from n and r using Q15 multiplies with some
     adjustment, carefully done to avoid overflow.*/
  r2 = r*r;
  y = (((r2 >> 15)*n + r2) >> 12) - 131077;
  ry = r*y;
  /*Apply a 2nd-order Householder iteration: r += r*y*(y*0.375-0.5).
    This yields the Q14 reciprocal square root of the Q16 t, with a maximum
     relative error of 1.04956E-4, a (relative) RMSE of 2.80979E-5, and a peak
     absolute error of 2.26591/16384.*/
  ret = r + ((((ry >> 16)*(3*y) >> 3) - ry) >> 18);
  OD_ASSERT(ret >= 16384 && ret < 32768);
  return (int16_t)ret;
}

static int16_t od_rsqrt(int32_t x, int *rsqrt_shift)
{
   int k;
   int s;
   int16_t t;
   k = (OD_ILOG(x) - 1) >> 1;
  /*t is x in the range [0.25, 1) in QINSHIFT, or x*2^(-s).
    Shift by log2(x) - log2(0.25*(1 << INSHIFT)) to ensure 0.25 lower bound.*/
   s = 2*k - (OD_RSQRT_INSHIFT - 2);
   t = OD_VSHR(x, s);
   /*We want to express od_rsqrt() in terms of od_rsqrt_norm(), which is
      defined as (2^OUTSHIFT)/sqrt(t*(2^-INSHIFT)) with t=x*(2^-s).
     This simplifies to 2^(OUTSHIFT+(INSHIFT/2)+(s/2))/sqrt(x), so the caller
      needs to shift right by OUTSHIFT + INSHIFT/2 + s/2.*/
   *rsqrt_shift = OD_RSQRT_OUTSHIFT + ((s + OD_RSQRT_INSHIFT) >> 1);
   return od_rsqrt_norm(t);
}
#endif

/** Synthesizes one parition of coefficient values from a PVQ-encoded
 * vector.  This 'partial' version is called by the encode loop where
 * the Householder reflection has already been computed and there's no
 * need to recompute it.
 *
 * @param [out]     xcoeff  output coefficient partition (x in math doc)
 * @param [in]      ypulse  PVQ-encoded values (y in the math doc); in
 *                          the noref case, this vector has n entries,
 *                          in the reference case it contains n-1 entries
 *                          (the m-th entry is not included)
 * @param [in]      r       reference vector (prediction)
 * @param [in]      n       number of elements in this partition
 * @param [in]      noref   indicates presence or lack of prediction
 * @param [in]      g       decoded quantized vector gain
 * @param [in]      theta   decoded theta (prediction error)
 * @param [in]      m       alignment dimension of Householder reflection
 * @param [in]      s       sign of Householder reflection
 * @param [in]      qm_inv  inverse of the QM with magnitude compensation
 */
void od_pvq_synthesis_partial(od_coeff *xcoeff, const od_coeff *ypulse,
 const od_val16 *r16, int n, int noref, od_val32 g, od_val32 theta, int m, int s,
 const int16_t *qm_inv) {
  int i;
  int yy;
  od_val32 scale;
  int nn;
#if !defined(OD_FLOAT_PVQ)
  int gshift;
  int qshift;
#endif
  OD_ASSERT(g != 0);
  nn = n-(!noref); /* when noref==0, vector in is sized n-1 */
  yy = 0;
  for (i = 0; i < nn; i++)
    yy += ypulse[i]*(int32_t)ypulse[i];
#if !defined(OD_FLOAT_PVQ)
  /* Shift required for the magnitude of the pre-qm synthesis to be guaranteed
     to fit in 16 bits. In practice, the range will be 8192-16384 after scaling
     most of the time. */
  gshift = OD_MAXI(0, OD_ILOG(g) - 14);
#endif
  /*scale is g/sqrt(yy) in Q(16-gshift) so that x[]*scale has a norm that fits
     in 16 bits.*/
  if (yy == 0) scale = 0;
#if defined(OD_FLOAT_PVQ)
  else {
    scale = g/sqrt(yy);
  }
#else
  else {
    int rsqrt_shift;
    int16_t rsqrt;
    /*FIXME: should be < int64_t*/
    int64_t tmp;
    rsqrt = od_rsqrt(yy, &rsqrt_shift);
    tmp = rsqrt*(int64_t)g;
    scale = OD_VSHR_ROUND(tmp, rsqrt_shift + gshift - 16);
  }
  /* Shift to apply after multiplying by the inverse QM, taking into account
     gshift. */
  qshift = OD_QM_INV_SHIFT - gshift;
#endif
  if (noref) {
    for (i = 0; i < n; i++) {
      od_val32 x;
      /* This multiply doesn't round, so it introduces some bias.
         It would be nice (but not critical) to fix this. */
      x = OD_MULT16_32_Q16(ypulse[i], scale);
#if defined(OD_FLOAT_PVQ)
      xcoeff[i] = (od_coeff)floor(.5
       + x*(qm_inv[i]*OD_QM_INV_SCALE_1));
#else
      xcoeff[i] = OD_SHR_ROUND(x*qm_inv[i], qshift);
#endif
    }
  }
  else{
    od_val16 x[MAXN];
    scale = OD_ROUND32(scale*OD_TRIG_SCALE_1*od_pvq_sin(theta));
    /* The following multiply doesn't round, but it's probably OK since
       the Householder reflection is likely to undo most of the resulting
       bias. */
    for (i = 0; i < m; i++)
      x[i] = OD_MULT16_32_Q16(ypulse[i], scale);
    x[m] = OD_ROUND16(-s*(OD_SHR_ROUND(g, gshift))*OD_TRIG_SCALE_1*
     od_pvq_cos(theta));
    for (i = m; i < nn; i++)
      x[i+1] = OD_MULT16_32_Q16(ypulse[i], scale);
    od_apply_householder(x, x, r16, n);
    for (i = 0; i < n; i++) {
#if defined(OD_FLOAT_PVQ)
      xcoeff[i] = (od_coeff)floor(.5 + (x[i]*(qm_inv[i]*OD_QM_INV_SCALE_1)));
#else
      xcoeff[i] = OD_SHR_ROUND(x[i]*qm_inv[i], qshift);
#endif
    }
  }
}

