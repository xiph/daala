#if !defined(_x86_x86int_H)
# define _x86_x86int_H (1)
# include "../state.h"

void od_state_opt_vtbl_init_x86(od_state *_state);

void od_mc_predict1imv8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],const int _m[4],int _r,int _log_xblk_sz,
 int _log_yblk_sz);
void od_mc_predict1fmv8_sse2(unsigned char *_dst,const unsigned char *_src,
 int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz);
void od_mc_blend_full8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz);
void od_mc_blend_full_split8_sse2(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz);

#endif
