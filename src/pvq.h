#if !defined(_pvq_H)
# define _pvq_H (1)
# include "internal.h"

int quant_pvq(ogg_int16_t *_x,const ogg_int16_t *_r,const int *_q,
    ogg_int16_t *_scale,int *out,int N,int K,int paper_scaling);

#endif
