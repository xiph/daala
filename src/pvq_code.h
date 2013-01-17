/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifndef PVQ_CODE_H_
#define PVQ_CODE_H_

#include "entenc.h"
#include "entdec.h"

extern const ogg_uint16_t cdf_table[][16];
extern const unsigned char decayE[];

typedef struct od_pvq_adapt_ctx od_pvq_adapt_ctx;

# define OD_DISABLE_PVQ_CODE1 (1)

/*TODO: These should be adapted based on target quality/bitrate.*/
# define OD_K_INIT_Q7 (2031)

# define OD_COUNT_INIT_Q7 (104)
# define OD_COUNT_EX_INIT_Q7 (128)

# if !defined(OD_DISABLE_PVQ_CODE1)
#  define OD_SUM_EX_INIT_Q7 (216)
#  define OD_POS_INIT_Q3 (13)

#  define OD_K_ADAPT_SPEED (2)
#  define OD_SUM_EX_ADAPT_SPEED (2)
#  define OD_POS_ADAPT_SPEED (1)
# else
#  define OD_SUM_EX_INIT_Q7 (216)

#  define OD_K_ADAPT_SPEED (2)
#  define OD_SUM_EX_ADAPT_SPEED (2)
# endif

# define OD_DELTA_ADAPT_SPEED (1)

/*The scaling for the row contexts is different.*/
# define OD_K_ROW_INIT_Q8 (2*OD_K_INIT_Q7-(OD_K_INIT_Q7>>OD_K_ADAPT_SPEED))
# define OD_SUM_EX_ROW_INIT_Q8 \
 (2*OD_SUM_EX_INIT_Q7-(OD_SUM_EX_INIT_Q7>>OD_SUM_EX_ADAPT_SPEED))
# define OD_COUNT_ROW_INIT_Q8 \
 (2*OD_COUNT_INIT_Q7-(OD_COUNT_INIT_Q7>>OD_DELTA_ADAPT_SPEED))
# define OD_COUNT_EX_ROW_INIT_Q8 \
 (2*OD_COUNT_EX_INIT_Q7-(OD_COUNT_EX_INIT_Q7>>OD_DELTA_ADAPT_SPEED))
# define OD_POS_ROW_INIT_Q4 \
 (2*OD_POS_INIT_Q3-(OD_POS_INIT_Q3>>OD_POS_ADAPT_SPEED+1))

struct od_pvq_adapt_ctx{
  /*Mean value of K.*/
  int mean_k_q8;
  /*Mean value of the sum of (pulses left)/(positions left) for all
     positions.*/
  int mean_sum_ex_q8;
  /*For the delta version.*/
  int mean_count_q8;
  int mean_count_ex_q8;
# if !defined(OD_DISABLE_PVQ_CODE1)
  /*Mean pulse position for the case of a single pulse (K==1).*/
  int mean_pos_q4;
# endif
  /*Actual value of K for this context.*/
  int k;
  /*Actual sum of (pulses left)/(positions left) for all positions for this
     context.*/
  int sum_ex_q8;
  int count_q8;
  int count_ex_q8;
# if !defined(OD_DISABLE_PVQ_CODE1)
  /*Actual pulse position for the case of a single pulse for this context.*/
  int pos;
# endif
};

void laplace_encode_special(od_ec_enc *enc,int pos,unsigned decay,int max);
int laplace_decode_special(od_ec_dec *dec,unsigned decay,int max);

void pvq_encoder(od_ec_enc *enc,const int *y,int N,int K,
 od_pvq_adapt_ctx *_adapt);
void pvq_decoder(od_ec_dec *dec,int *y,int N,int K,od_pvq_adapt_ctx *_adapt);

void pvq_encode_delta(od_ec_enc *enc, const int *y,int N,int K,
 od_pvq_adapt_ctx *_adapt);
void pvq_decode_delta(od_ec_dec *dec, int *y,int N,int K,
 od_pvq_adapt_ctx *_adapt);

#endif
