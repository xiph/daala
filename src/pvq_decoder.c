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

#include <stdlib.h>
#include <stdio.h>
#include "internal.h"
#include "decint.h"
#include "entdec.h"
#include "entcode.h"
#include "logging.h"
#include "laplace_code.h"
#include "pvq_code.h"

/** Decodes a single vector of integers (eg, a partition within a
 *  coefficient block) encoded using PVQ
 *
 * @param [in,out] ec      range encoder
 * @param [in]     q       scale/quantizer
 * @param [in]     n       number of coefficients in partition
 * @param [in,out] model   entropy decoder state
 * @param [in,out] adapt   adaptation context
 * @param [in,out] exg     ExQ16 expectation of decoded gain value
 * @param [in,out] ext     ExQ16 expectation of decoded theta value
 * @param [in]     ref     'reference' (prediction) vector
 * @param [out]    out     decoded partition
 * @param [in]     noref   boolean indicating absence of reference
 */
static void pvq_decode_partition(od_ec_dec *ec,
                                 int q,
                                 int n,
                                 generic_encoder *model,
                                 int *adapt,
                                 int *exg,
                                 int *ext,
                                 od_coeff *ref,
                                 od_coeff *out,
                                 int noref) {
  int adapt_curr[OD_NSB_ADAPT_CTXS] = {0};
  int speed;
  int k;
  double qcg;
  int max_theta;
  int itheta;
  double theta;
  double gr;
  double gain_offset;
  od_coeff y[1024];
  double r[1024];
  int qg;

  speed = 5;
  theta = 0;
  gr = 0;
  gain_offset = 0;

  /* read quantized gain */
  qg = generic_decode(ec, model, exg, 2);

  if(!noref){
    /* we have a reference; compute its gain */
    double cgr;
    int icgr;
    int i;
    cgr = pvq_compute_gain(ref, n, q, &gr);
    icgr = floor(.5+cgr);
    /* quantized gain is interleave encoded when there's a reference;
       deinterleave it now */
    qg = neg_deinterleave(qg, icgr);
    gain_offset = cgr-icgr;
    qcg = qg + gain_offset;
    /* read and decode first-stage PVQ error theta */
    max_theta = pvq_compute_max_theta(qcg);
    itheta = generic_decode(ec, model, ext, 2);
    theta = pvq_compute_theta(itheta, max_theta);
    for (i = 0; i < n; i++) r[i] = ref[i];
  }
  else{
    qcg=qg;
  }

  k = pvq_compute_k(qcg, theta, noref, n);
  /* when noref==0, y is actually size n-1 */
  laplace_decode_vector(ec, y, n-(!noref), k, adapt_curr, adapt);
  pvq_synthesis(out, y, r, n, gr, noref, qg, gain_offset, theta, q);

  if (adapt_curr[OD_ADAPT_K_Q8] > 0) {
    adapt[OD_ADAPT_K_Q8]
      += 256*adapt_curr[OD_ADAPT_K_Q8]-adapt[OD_ADAPT_K_Q8]>>speed;
    adapt[OD_ADAPT_SUM_EX_Q8]
      += adapt_curr[OD_ADAPT_SUM_EX_Q8]-adapt[OD_ADAPT_SUM_EX_Q8]>>speed;
  }
  if (adapt_curr[OD_ADAPT_COUNT_Q8] > 0) {
    adapt[OD_ADAPT_COUNT_Q8]
      += adapt_curr[OD_ADAPT_COUNT_Q8]-adapt[OD_ADAPT_COUNT_Q8]>>speed;
    adapt[OD_ADAPT_COUNT_EX_Q8]
      += adapt_curr[OD_ADAPT_COUNT_EX_Q8]-adapt[OD_ADAPT_COUNT_EX_Q8]>>speed;
  }
}

/** Decodes a coefficient block (except for DC) encoded using PVQ
 *
 * @param [in,out] dec     daala decoder context
 * @param [in]     ref     'reference' (prediction) vector
 * @param [out]    out     decoded partition
 * @param [in]     q       quantizer
 * @param [in]     n       number of coefficients on one side of block
 */
void pvq_decode(daala_dec_ctx *dec,
                od_coeff *ref,
                od_coeff *out,
                int q,
                int n){

  int noref[PVQ_MAX_PARTITIONS];
  int *adapt;
  int *exg;
  int *ext;
  int predflags8;
  int predflags16;
  generic_encoder *model;
  adapt = dec->state.pvq_adapt;
  exg = dec->state.pvq_exg;
  ext = dec->state.pvq_ext;
  model = &dec->state.pvq_gain_model;

  if (n == 4) {
    noref[0] = !od_ec_decode_bool_q15(&dec->ec, PRED4_PROB);
    pvq_decode_partition(&dec->ec, q, 15, model, adapt, exg, ext, ref+1,
                   out+1, noref[0]);
  }
  else {
    predflags8 = od_ec_decode_cdf_q15(&dec->ec, pred8_cdf, 16);
    noref[0] = !(predflags8>>3);
    noref[1] = !((predflags8>>2) & 0x1);
    noref[2] = !((predflags8>>1) & 0x1);
    noref[3] = !(predflags8 & 0x1);

    if(n >= 16) {
      predflags16 = od_ec_decode_cdf_q15(&dec->ec, pred16_cdf[predflags8], 8);
      noref[4] = !((predflags16>>2) & 0x1);
      noref[5] = !((predflags16>>1) & 0x1);
      noref[6] = !(predflags16 & 0x1);
    }

    pvq_decode_partition(&dec->ec, q, 15, model, adapt, exg, ext, ref+1,
                   out+1, noref[0]);
    pvq_decode_partition(&dec->ec, q, 8, model, adapt, exg+1, ext+1, ref+16,
                   out+16, noref[1]);
    pvq_decode_partition(&dec->ec, q, 8, model, adapt, exg+2, ext+2, ref+24,
                   out+24, noref[2]);
    pvq_decode_partition(&dec->ec, q, 32, model, adapt, exg+3, ext+3, ref+32,
                   out+32, noref[3]);

    if(n >= 16) {
      pvq_decode_partition(&dec->ec, q, 32, model, adapt, exg+4, ext+4, ref+64,
                     out+64, noref[4]);
      pvq_decode_partition(&dec->ec, q, 32, model, adapt, exg+5, ext+5, ref+96,
                     out+96, noref[5]);
      pvq_decode_partition(&dec->ec, q, 128, model, adapt, exg+6, ext+6, ref+128,
                     out+128, noref[6]);
    }
  }
}
