#if !defined(_pvq_H)
# define _pvq_H (1)
# include "internal.h"

int quant_pvq(ogg_int16_t *_x,const ogg_int16_t *_r,const int *_q,
     int *out,int stride,int N,int K);

#endif
