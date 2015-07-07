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
#include "entdec.h"
#include "logging.h"
#include "odintrin.h"

/** Decodes a value from 0 to N-1 (with N up to 16) based on a cdf and adapts
 * the cdf accordingly.
 *
 * @param [in,out] enc   range encoder
 * @param [in]     cdf   CDF of the variable (Q15)
 * @param [in]     n     number of values possible
 * @param [in]     increment adaptation speed (Q15)
 *
 * @retval decoded variable
 */
int od_decode_cdf_adapt(od_ec_dec *ec, uint16_t *cdf, int n,
 int increment) {
  int i;
  int val;
  val = od_ec_decode_cdf_unscaled(ec, cdf, n);
  if (cdf[n-1] + increment > 32767) {
    for (i = 0; i < n; i++) {
      /* Second term ensures that the pdf is non-null */
      cdf[i] = (cdf[i] >> 1) + i + 1;
    }
  }
  for (i = val; i < n; i++) cdf[i] += increment;
  return val;
}

/** Encodes a random variable using a "generic" model, assuming that the
 * distribution is one-sided (zero and up), has a single mode, and decays
 * exponentially past the model.
 *
 * @param [in,out] dec   range decoder
 * @param [in,out] model generic probability model
 * @param [in]     x     variable being encoded
 * @param [in,out] ExQ16 expectation of x (adapted)
 * @param [in]     integration integration period of ExQ16 (leaky average over
 * 1<<integration samples)
 *
 * @retval decoded variable x
 */
int generic_decode(od_ec_dec *dec, generic_encoder *model, int max,
 int *ex_q16, int integration) {
  int lg_q1;
  int shift;
  int id;
  uint16_t *cdf;
  int xs;
  int lsb;
  int x;
  int ms;
  lsb = 0;
  if (max == 0) return 0;
  lg_q1 = log_ex(*ex_q16);
  /* If expectation is too large, shift x to ensure that
     all we have past xs=15 is the exponentially decaying tail
     of the distribution. */
  shift = OD_MAXI(0, (lg_q1 - 5) >> 1);
  /* Choose the cdf to use: we have two per "octave" of ExQ16. */
  id = OD_MINI(GENERIC_TABLES - 1, lg_q1);
  cdf = model->cdf[id];
  ms = (max + (1 << shift >> 1)) >> shift;
  if (max == -1) xs = od_ec_decode_cdf_unscaled(dec, cdf, 16);
  else xs = od_ec_decode_cdf_unscaled(dec, cdf, OD_MINI(ms + 1, 16));
  if (xs == 15) {
    int e;
    unsigned decay;
    /* Estimate decay based on the assumption that the distribution is close
       to Laplacian for large values. We should probably have an adaptive
       estimate instead. Note: The 2* is a kludge that's not fully understood
       yet. */
    e = ((2**ex_q16 >> 8) + (1 << shift >> 1)) >> shift;
    decay = OD_MAXI(2, OD_MINI(254, 256*e/(e + 256)));
    xs += laplace_decode_special(dec, decay, (max == -1) ? -1 : ms - 15);
  }
  if (shift != 0) {
    int special;
    /* Because of the rounding, there's only half the number of possibilities
       for xs=0 */
    special = xs == 0;
    if (shift - special > 0) lsb = od_ec_dec_bits(dec, shift - special);
    lsb -= !special << (shift - 1);
  }
  x = (xs << shift) + lsb;
  generic_model_update(model, ex_q16, x, xs, id, integration);
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG,
   "dec: %d %d %d %d %d %x", *ex_q16, x, shift, id, xs, dec->rng));
  return x;
}
