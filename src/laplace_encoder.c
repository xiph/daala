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

#include <stdio.h>

#include "laplace_code.h"
#include "entenc.h"
#include "entdec.h"
#include "logging.h"
#include "odintrin.h"

/** Encodes the tail of a Laplace-distributed variable, i.e. it doesn't
 * do anything special for the zero case.
 *
 * @param [in,out] enc     range encoder
 * @param [in]     x       variable to encode (has to be positive)
 * @param [in]     decay   decay factor of the distribution in Q8 format,
 * i.e. pdf ~= decay^x
 * @param [in]     max     maximum possible value of x (used to truncate
 * the pdf)
 */
void laplace_encode_special(od_ec_enc *enc, int x, unsigned decay, int max) {
  int shift;
  int xs;
  int ms;
  int sym;
  const uint16_t *cdf;
  shift = 0;
  if (max == 0) return;
  /* We don't want a large decay value because that would require too many
     symbols. However, it's OK if the max is below 15. */
  while (((max >> shift) >= 15 || max == -1) && decay > 235) {
    decay = (decay*decay + 128) >> 8;
    shift++;
  }
  OD_ASSERT(x <= max || max == -1);
  decay = OD_MINI(decay, 254);
  decay = OD_MAXI(decay, 2);
  xs = x >> shift;
  ms = max >> shift;
  cdf = EXP_CDF_TABLE[(decay + 1) >> 1];
  OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "decay = %d", decay));
  do {
    sym = OD_MINI(xs, 15);
    {
      int i;
      OD_LOG((OD_LOG_PVQ, OD_LOG_DEBUG, "%d %d %d %d %d\n", x, xs, shift,
       sym, max));
      for (i = 0; i < 16; i++) {
        OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "%d ", cdf[i]));
      }
      OD_LOG_PARTIAL((OD_LOG_PVQ, OD_LOG_DEBUG, "\n"));
    }
    if (ms > 0 && ms < 15) {
      /* Simple way of truncating the pdf when we have a bound */
      od_ec_encode_cdf_unscaled(enc, sym, cdf, ms + 1);
    }
    else {
      od_ec_encode_cdf_q15(enc, sym, cdf, 16);
    }
    xs -= 15;
    ms -= 15;
  }
  while (sym >= 15 && ms != 0);
  if (shift) od_ec_enc_bits(enc, x & ((1 << shift) - 1), shift);
}

/** Encodes a Laplace-distributed variable for use in PVQ
 *
 * @param [in,out] enc  range encoder
 * @param [in]     x    variable to encode (including sign)
 * @param [in]     ExQ8 expectation of the absolute value of x in Q8
 * @param [in]     K    maximum value of |x|
 */
void laplace_encode(od_ec_enc *enc, int x, int ex_q8, int k) {
  int j;
  int shift;
  int xs;
  uint16_t cdf[16];
  int sym;
  int decay;
  int offset;
  /* shift down x if expectation is too high */
  shift = od_ilog(ex_q8) - 11;
  if (shift < 0) shift = 0;
  /* Apply the shift with rounding to Ex, K and xs */
  ex_q8 = (ex_q8 + (1 << shift >> 1)) >> shift;
  k = (k + (1 << shift >> 1)) >> shift;
  xs = (x + (1 << shift >> 1)) >> shift;
  decay = OD_MINI(254, 256*ex_q8/(ex_q8 + 256));
  offset = LAPLACE_OFFSET[(decay + 1) >> 1];
  for (j = 0; j < 16; j++) {
    cdf[j] = EXP_CDF_TABLE[(decay + 1) >> 1][j] - offset;
  }
  sym = xs;
  if (sym > 15) sym = 15;
  /* Simple way of truncating the pdf when we have a bound */
  if (k != 0) od_ec_encode_cdf_unscaled(enc, sym, cdf, OD_MINI(k + 1, 16));
  if (shift) {
    int special;
    /* Because of the rounding, there's only half the number of possibilities
       for xs=0 */
    special = xs == 0;
    if (shift - special > 0) {
      od_ec_enc_bits(enc, x - (xs << shift) + (!special << (shift - 1)),
       shift - special);
    }
  }
  /* Handle the exponentially-decaying tail of the distribution */
  OD_ASSERT(xs - 15 <= k - 15);
  if (xs >= 15) laplace_encode_special(enc, xs - 15, decay, k - 15);
}

static void laplace_encode_vector_delta(od_ec_enc *enc, const od_coeff *y, int n, int k,
                                        int32_t *curr, const int32_t *means) {
  int i;
  int prev;
  int sum_ex;
  int sum_c;
  int first;
  int k_left;
  int coef;
  prev = 0;
  sum_ex = 0;
  sum_c = 0;
  first = 1;
  k_left = k;
  coef = 256*means[OD_ADAPT_COUNT_Q8]/
   (1 + means[OD_ADAPT_COUNT_EX_Q8]);
  coef = OD_MAXI(coef, 1);
  for (i = 0; i < n; i++) {
    if (y[i] != 0) {
      int j;
      int count;
      int mag;
      mag = abs(y[i]);
      count = i - prev;
      if (first) {
        int decay;
        int ex = coef*(n - prev)/k_left;
        if (ex > 65280) decay = 255;
        else {
          decay = OD_MINI(255,
           (int)((256*ex/(ex + 256) + (ex>>5)*ex/((n + 1)*(n - 1)*(n - 1)))));
        }
        /*Update mean position.*/
        OD_ASSERT(count <= n - 1);
        laplace_encode_special(enc, count, decay, n - 1);
        first = 0;
      }
      else laplace_encode(enc, count, coef*(n - prev)/k_left, n - prev - 1);
      sum_ex += 256*(n - prev);
      sum_c += count*k_left;
      od_ec_enc_bits(enc, y[i] < 0, 1);
      for (j = 0; j < mag - 1; j++) {
        laplace_encode(enc, 0, coef*(n - i)/(k_left - 1 - j), n - i - 1);
        sum_ex += 256*(n - i);
      }
      k_left -= mag;
      prev = i;
      if (k_left == 0) break;
    }
  }
  if (k > 0) {
    curr[OD_ADAPT_COUNT_Q8] = 256*sum_c;
    curr[OD_ADAPT_COUNT_EX_Q8] = sum_ex;
  }
  else {
    curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
    curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
  }
  curr[OD_ADAPT_K_Q8] = 0;
  curr[OD_ADAPT_SUM_EX_Q8] = 0;
}

/** Encodes a vector of integers assumed to come from rounding a sequence of
 * Laplace-distributed real values in decreasing order of variance.
 *
 * @param [in,out] enc range encoder
 * @param [in]     y     vector to encode
 * @param [in]     N     dimension of the vector
 * @param [in]     K     sum of the absolute value of components of y
 * @param [out]    curr  Adaptation context output, may alias means.
 * @param [in]     means Adaptation context input.
 */
void laplace_encode_vector(od_ec_enc *enc, const od_coeff *y, int n, int k,
                           int32_t *curr, const int32_t *means) {
  int i;
  int sum_ex;
  int kn;
  int exp_q8;
  int mean_k_q8;
  int mean_sum_ex_q8;
  int ran_delta;
  ran_delta = 0;
  if (k <= 1) {
    laplace_encode_vector_delta(enc, y, n, k, curr, means);
    return;
  }
  sum_ex = 0;
  kn = k;
  /* Estimates the factor relating pulses_left and positions_left to E(|x|) */
  mean_k_q8 = means[OD_ADAPT_K_Q8];
  mean_sum_ex_q8 = means[OD_ADAPT_SUM_EX_Q8];
  if (mean_k_q8 < 1 << 23) exp_q8 = 256*mean_k_q8/(1 + mean_sum_ex_q8);
  else exp_q8 = mean_k_q8/(1 + (mean_sum_ex_q8 >> 8));
  for (i = 0; i < n; i++) {
    int ex;
    int x;
    if (kn == 0) break;
    if (kn <= 1 && i != n - 1) {
      laplace_encode_vector_delta(enc, y + i, n - i, kn, curr, means);
      ran_delta = 1;
      break;
    }
    x = abs(y[i]);
    /* Expected value of x (round-to-nearest) is
       expQ8*pulses_left/positions_left */
    ex = (2*exp_q8*kn + (n - i))/(2*(n - i));
    if (ex > kn*256) ex = kn*256;
    sum_ex += (2*256*kn + (n - i))/(2*(n - i));
    /* No need to encode the magnitude for the last bin. */
    if (i != n - 1) laplace_encode(enc, x, ex, kn);
    if (x != 0) od_ec_enc_bits(enc, y[i] < 0, 1);
    kn -= x;
  }
  /* Adapting the estimates for expQ8 */
  if (!ran_delta) {
    curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
    curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
  }
  curr[OD_ADAPT_K_Q8] = k - kn;
  curr[OD_ADAPT_SUM_EX_Q8] = sum_ex;
}
