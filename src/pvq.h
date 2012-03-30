#if !defined(_pvq_H)
# define _pvq_H (1)
# include "internal.h"

int quant_pvq_theta(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q);

int quant_pvq(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q);

int quant_scalar(ogg_int32_t *_x,const ogg_int32_t *_r,
    ogg_int16_t *_scale,int *y,int N,int Q);

#endif
