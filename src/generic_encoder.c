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
#include "generic_code.h"
#include "entenc.h"
#include "entdec.h"
#include "logging.h"
#include "odintrin.h"

/** Encodes a value from 0 to N-1 (with N up to 16) based on a cdf and adapts
 * the cdf accordingly.
 *
 * @param [in,out] enc   range encoder
 * @param [in]     val   variable being encoded
 * @param [in]     cdf   CDF of the variable (Q15)
 * @param [in]     n     number of values possible
 * @param [in]     increment adaptation speed (Q15)
 */
void od_encode_cdf_adapt(od_ec_enc *ec, int val, uint16_t *cdf, int n,
 int increment) {
  int i;
  od_ec_encode_cdf_unscaled(ec, val, cdf, n);
  if (cdf[n-1] + increment > 32767) {
    for (i = 0; i < n; i++) {
      /* Second term ensures that the pdf is non-null */
      cdf[i] = (cdf[i] >> 1) + i + 1;
    }
  }
  for (i = val; i < n; i++) cdf[i] += increment;
}

/** Encodes a random variable using a "generic" model, assuming that the
 * distribution is one-sided (zero and up), has a single mode, and decays
 * exponentially past the model.
 *
 * @param [in,out] enc   range encoder
 * @param [in,out] model generic probability model
 * @param [in]     x     variable being encoded
 * @param [in]     max   largest value possible
 * @param [in,out] ExQ16 expectation of x (adapted)
 * @param [in]     integration integration period of ExQ16 (leaky average over
 * 1<<integration samples)
 */
void generic_encode(od_ec_enc *enc, generic_encoder *model, int x, int max,
 int *ex_q16, int integration) {
  int lg_q1;
  int shift;
  int id;
  uint16_t *cdf;
  int xs;
  int ms;
  if (max == 0) return;
  lg_q1 = log_ex(*ex_q16);
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG,
   "%d %d", *ex_q16, lg_q1));
  /* If expectation is too large, shift x to ensure that
     all we have past xs=15 is the exponentially decaying tail
     of the distribution */
  shift = OD_MAXI(0, (lg_q1 - 5) >> 1);
  /* Choose the cdf to use: we have two per "octave" of ExQ16 */
  id = OD_MINI(GENERIC_TABLES - 1, lg_q1);
  cdf = model->cdf[id];
  xs = (x + (1 << shift >> 1)) >> shift;
  ms = (max + (1 << shift >> 1)) >> shift;
  OD_ASSERT(max == -1 || xs <= ms);
  if (max == -1) od_ec_encode_cdf_unscaled(enc, OD_MINI(15, xs), cdf, 16);
  else {
    od_ec_encode_cdf_unscaled(enc, OD_MINI(15, xs), cdf, OD_MINI(ms + 1, 16));
  }
  if (xs >= 15) {
    int e;
    unsigned decay;
    /* Estimate decay based on the assumption that the distribution is close
       to Laplacian for large values. We should probably have an adaptive
       estimate instead. Note: The 2* is a kludge that's not fully understood
       yet. */
    e = ((2**ex_q16 >> 8) + (1 << shift >> 1)) >> shift;
    decay = OD_MAXI(2, OD_MINI(254, 256*e/(e + 256)));
    /* Encode the tail of the distribution assuming exponential decay. */
    laplace_encode_special(enc, xs - 15, decay, (max == -1) ? -1 : ms - 15);
  }
  if (shift != 0) {
    int special;
    /* Because of the rounding, there's only half the number of possibilities
       for xs=0. */
    special = xs == 0;
    if (shift - special > 0) {
      od_ec_enc_bits(enc, x - (xs << shift) + (!special << (shift - 1)),
       shift - special);
    }
  }
  generic_model_update(model, ex_q16, x, xs, id, integration);
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG,
   "enc: %d %d %d %d %d %x", *ex_q16, x, shift, id, xs, enc->rng));
}

/** Estimates the cost of encoding a value with generic_encode().
 *
 * @param [in,out] model generic probability model
 * @param [in]     x     variable being encoded
 * @param [in]     max   largest value possible
 * @param [in,out] ExQ16 expectation of x (adapted)
 * @return number of bits (approximation)
 */
double generic_encode_cost(generic_encoder *model, int x, int max,
 int *ex_q16) {
  int lg_q1;
  int shift;
  int id;
  uint16_t *cdf;
  int xs;
  int ms;
  int extra;
  if (max == 0) return 0;
  lg_q1 = log_ex(*ex_q16);
  /* If expectation is too large, shift x to ensure that
       all we have past xs=15 is the exponentially decaying tail
       of the distribution */
  shift = OD_MAXI(0, (lg_q1 - 5) >> 1);
  /* Choose the cdf to use: we have two per "octave" of ExQ16 */
  id = OD_MINI(GENERIC_TABLES - 1, lg_q1);
  cdf = model->cdf[id];
  xs = (x + (1 << shift >> 1)) >> shift;
  ms = (max + (1 << shift >> 1)) >> shift;
  OD_ASSERT(max == -1 || xs <= ms);
  extra = 0;
  if (shift) extra = shift - (xs == 0);
  xs = OD_MINI(15, xs);
  /* Shortcut: assume it's going to cost 2 bits for the Laplace coder. */
  if (xs == 15) extra += 2;
  if (max == -1) {
    return extra - OD_LOG2((double)(cdf[xs] - (xs == 0 ? 0 : cdf[xs - 1]))/
     cdf[15]);
  }
  else {
    return extra - OD_LOG2((double)(cdf[xs] - (xs == 0 ? 0 : cdf[xs - 1]))/
     cdf[OD_MINI(ms, 15)]);
  }
}
