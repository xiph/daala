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

#include "generic_code.h"

void od_cdf_init(uint16_t *cdf, int ncdfs, int nsyms, int val, int first) {
  int i;
  int j;
  for (i = 0; i < ncdfs; i++) {
    cdf[i*nsyms] = first;
    for (j = 0; j < nsyms; j++) {
      cdf[i*nsyms + j] = val*j + first;
    }
  }
}

/** Initializes the cdfs and freq counts for a model.
 *
 * @param [out] model model being initialized
 */
void generic_model_init(generic_encoder *model) {
  int i;
  int j;
  model->increment = 64;
  for (i = 0; i < GENERIC_TABLES; i++) {
    for (j = 0; j < 16; j++) {
      /* Do flat initialization equivalent to a single symbol in each bin. */
      model->cdf[i][j] = (j + 1) * model->increment;
    }
  }
}

/** Takes the base-2 log of E(x) in Q1.
 *
 * @param [in] ExQ16 expectation of x in Q16
 *
 * @retval 2*log2(ExQ16/2^16)
 */
int log_ex(int ex_q16) {
  int lg;
  int lg_q1;
  int odd;
  lg = od_ilog(ex_q16);
  if (lg < 15) {
    odd = ex_q16*ex_q16 > 2 << 2*lg;
  }
  else {
    int tmp;
    tmp = ex_q16 >> (lg - 8);
    odd = tmp*tmp > (1 << 15);
  }
  lg_q1 = OD_MAXI(0, 2*lg - 33 + odd);
  return lg_q1;
}

/** Updates the probability model based on the encoded/decoded value
 *
 * @param [in,out] model generic prob model
 * @param [in,out] ExQ16 expectation of x
 * @param [in]     x     variable encoded/decoded (used for ExQ16)
 * @param [in]     xs    variable x after shift (used for the model)
 * @param [in]     id    id of the icdf to adapt
 * @param [in]     integration integration period of ExQ16 (leaky average over
 * 1<<integration samples)
 */
void generic_model_update(generic_encoder *model, int *ex_q16, int x, int xs,
 int id, int integration) {
  int i;
  int xenc;
  uint16_t *cdf;
  cdf = model->cdf[id];
  /* Renormalize if we cannot add increment */
  if (cdf[15] + model->increment > 32767) {
    for (i = 0; i < 16; i++) {
      /* Second term ensures that the pdf is non-null */
      cdf[i] = (cdf[i] >> 1) + i + 1;
    }
  }
  /* Update freq count */
  xenc = OD_MINI(15, xs);
  /* This can be easily vectorized */
  for (i = xenc; i < 16; i++) cdf[i] += model->increment;
  /* We could have saturated ExQ16 directly, but this is safe and simpler */
  x = OD_MINI(x, 32767);
  OD_IIR_DIADIC(*ex_q16, x << 16, integration);
}
