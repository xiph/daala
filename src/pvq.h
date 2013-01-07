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

#if !defined(_pvq_H)
# define _pvq_H (1)
# include "internal.h"

typedef struct od_pvq_adapt_ctx od_pvq_adapt_ctx;

# define OD_DISABLE_PVQ_CODE1 (1)


/*TODO: These should be adapted based on target quality/bitrate.*/
# define OD_K_INIT_Q7 (2031)
# if !defined(OD_DISABLE_PVQ_CODE1)
#  define OD_SUM_EX_INIT_Q7 (216)
#  define OD_POS_INIT_Q3 (13)

#  define OD_K_ADAPT_SPEED (5)
#  define OD_SUM_EX_ADAPT_SPEED (5)
#  define OD_POS_ADAPT_SPEED (1)
# else
#  define OD_SUM_EX_INIT_Q7 (194)

#  define OD_K_ADAPT_SPEED (4)
#  define OD_SUM_EX_ADAPT_SPEED (4)
# endif

/*The scaling for the row contexts is different.*/
# define OD_K_ROW_INIT_Q8 (2*OD_K_INIT_Q7-(OD_K_INIT_Q7>>OD_K_ADAPT_SPEED))
# define OD_SUM_EX_ROW_INIT_Q8 \
 (2*OD_SUM_EX_INIT_Q7-(OD_SUM_EX_INIT_Q7>>OD_SUM_EX_ADAPT_SPEED))
# define OD_POS_ROW_INIT_Q4 \
 (2*OD_POS_INIT_Q3-(OD_POS_INIT_Q3>>OD_POS_ADAPT_SPEED+1))

struct od_pvq_adapt_ctx{
  /*Mean value of K.*/
  int mean_k_q8;
  /*Mean value of the sum of (pulses left)/(positions left) for all
     positions.*/
  int mean_sum_ex_q8;
# if !defined(OD_DISABLE_PVQ_CODE1)
  /*Mean pulse position for the case of a single pulse (K==1).*/
  int mean_pos_q4;
# endif
  /*Actual value of K for this context.*/
  int k;
  /*Actual sum of (pulses left)/(positions left) for all positions for this
     context.*/
  int sum_ex_q8;
# if !defined(OD_DISABLE_PVQ_CODE1)
  /*Actual pulse position for the case of a single pulse for this context.*/
  int pos;
# endif
};

int quant_pvq_theta(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q);

int quant_pvq(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q,int *qg);

int quant_scalar(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q);

int quant_pvq_noref(ogg_int32_t *_x,float gr,
    ogg_int16_t *_scale,int *y,int N,int Q);

#endif
