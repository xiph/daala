#if !defined(_filter_H)
# define _filter_H (1)
# include "internal.h"

/*These are the pre/post filtering functions used by Daala.
  The idea is to pre/post filter in the spatial domain (the time domain in
   signal processing terms) to improve the energy compaction as well as reduce
   or eliminate blocking artifacts.*/

void od_pre_filter_4(ogg_int16_t _y[4],const ogg_int16_t _x[4]);
void od_post_filter_4(ogg_int16_t _x[4],const ogg_int16_t _y[4]);
void od_pre_filter_8(ogg_int16_t _y[8],const ogg_int16_t _x[8]);
void od_post_filter_8(ogg_int16_t _x[8],const ogg_int16_t _y[8]);
void od_pre_filter_16(ogg_int16_t _y[16],const ogg_int16_t _x[16]);
void od_post_filter_16(ogg_int16_t _x[16],const ogg_int16_t _y[16]);

#endif
