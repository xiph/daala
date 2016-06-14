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

#if !defined(_pvq_encoder_H)
# define _pvq_encoder_H (1)
# include "internal.h"
# include "filter.h"
# include "pvq.h"
# include "entenc.h"
# include "encint.h"

void od_encode_band_pvq_splits(od_ec_enc *ec, od_pvq_codeword_ctx *adapt,
 const int *y, int n, int k, int level);

void laplace_encode_special(od_ec_enc *enc, int x, unsigned decay, int max);
void laplace_encode(od_ec_enc *enc, int x, int ex_q8, int k);
void laplace_encode_vector(od_ec_enc *enc, const od_coeff *y, int n, int k,
                                  int32_t *curr, const int32_t *means);

#if OD_SIGNAL_Q_SCALING
void od_encode_quantizer_scaling(daala_enc_ctx *enc, int q_scaling, int bx,
 int by, int skip);
#endif

int od_pvq_encode(daala_enc_ctx *enc, od_coeff *ref, const od_coeff *in,
 od_coeff *out, int q0, int pli, int bs, const double *beta, int robust,
 int is_keyframe, int q_scaling, int bx, int by, const int16_t *qm,
 const int16_t *qm_inv, int speed);

#endif
