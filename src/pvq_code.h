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

#if !defined(_pvq_code_H)
# define _pvq_code_H

# include "entenc.h"
# include "entdec.h"
# include "adapt.h"
# include "filter.h"

extern const ogg_uint16_t EXP_CDF_TABLE[][16];
extern const ogg_uint16_t LAPLACE_OFFSET[];

void laplace_encode_special(od_ec_enc *enc, int pos, unsigned decay, int max);
int laplace_decode_special(od_ec_dec *dec, unsigned decay, int max);

void laplace_encode(od_ec_enc *enc, int x, int ex_q8, int k);
int laplace_decode(od_ec_dec *dec, int ex_q8, int k);

void pvq_encoder(od_ec_enc *enc, const od_coeff *y, int n, int k,
 ogg_int32_t *curr, const ogg_int32_t *means);
void pvq_decoder(od_ec_dec *dec, od_coeff *y, int n, int k,
 ogg_int32_t *curr, const ogg_int32_t *means);

void pvq_encode_delta(od_ec_enc *enc, const od_coeff *y, int n, int k,
 ogg_int32_t *curr, const ogg_int32_t *means);
void pvq_decode_delta(od_ec_dec *dec, od_coeff *y, int n, int k,
 ogg_int32_t *curr, const ogg_int32_t *means);

#endif
