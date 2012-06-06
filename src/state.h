/*
    Daala video codec
    Copyright (C) 2006-2010 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/


#if !defined(_state_H)
# define _state_H (1)

typedef struct od_state_opt_vtbl od_state_opt_vtbl;
typedef struct od_state          od_state;

# include "internal.h"
# include "mc.h"

/*The golden reference frame.*/
#define OD_FRAME_GOLD (0)
/*The previous reference frame.*/
#define OD_FRAME_PREV (1)
/*The next reference frame.*/
#define OD_FRAME_NEXT (2)
/*The current frame.*/
#define OD_FRAME_SELF (3)



extern const int        OD_VERT_D[];
/*The vector offsets in the X direction for each motion comepnsation block
   vertex from the upper-left.*/
#define OD_VERT_DX (OD_VERT_D+1)
/*The vector offsets in the Y direction for each motion compensation block
   vertex from the upper-left.*/
#define OD_VERT_DY (OD_VERT_D+0)
extern const int *const OD_VERT_SETUP_DX[4][4];
extern const int *const OD_VERT_SETUP_DY[4][4];



/*This should be a power of 2, and at least 8.*/
#define OD_UMV_PADDING (16)



/*The shared (encoder and decoder) functions that have accelerated variants.*/
struct od_state_opt_vtbl{
  void (*mc_predict1imv8)(unsigned char *_dst,int _dystride,
   const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
   const ogg_int32_t _mvy[4],const int _m[4],int _r,int _log_xblk_sz,
   int _log_yblk_sz);
  void (*mc_predict1fmv8)(unsigned char *_dst,const unsigned char *_src,
   int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
   int _log_xblk_sz,int _log_yblk_sz);
  void (*mc_blend_full8)(unsigned char *_dst,int _dystride,
   const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz);
  void (*mc_blend_full_split8)(unsigned char *_dst,int _dystride,
   const unsigned char *_src[4],int _c,int _s,
   int _log_xblk_sz,int _log_yblk_sz);
};



struct od_state{
  daala_info          info;
  od_state_opt_vtbl   opt_vtbl;
  ogg_uint32_t        cpu_flags;
  od_img              input;
  int                 ref_imgi[4];
  od_img              ref_imgs[4];
  od_img              rec_img;
  unsigned char      *ref_line_buf[8];
  unsigned char      *ref_img_data;
  ogg_int64_t         cur_time;
  od_mv_grid_pt     **mv_grid;
  int                 nhmbs;
  int                 nvmbs;
#if defined(OD_DUMP_IMAGES)
  od_img              vis_img;
#if defined(OD_ANIMATE)
  int                 ani_iter;
#endif
#endif
};


int  od_state_init(od_state *_state,const daala_info *_info);
void od_state_clear(od_state *_state);

void od_state_pred_block_from_setup(od_state *_state,unsigned char *_buf,
 int _ystride,int _ref,int _pli,int _vx,int _vy,int _c,int _s,int _log_mvb_sz);
void od_state_pred_block(od_state *_state,unsigned char *_buf,int _ystride,
 int _ref,int _pli,int _vx,int _vy,int _log_mvb_sz);
void od_state_mc_predict(od_state *_state,int _ref);
void od_state_upsample8(od_state *_state,od_img *_dst,const od_img *_src);
int od_state_dump_yuv(od_state *_state,od_img *_img,const char *_suf);
#if defined(OD_DUMP_IMAGES)
int od_state_dump_img(od_state *_state,od_img *_img,const char *_suf);
void od_img_draw_point(od_img *_img,int _x,int _y,
 const unsigned char _ycbcr[3]);
void od_img_draw_line(od_img *_img,int _x0,int _y0,int _x1,int _y1,
 const unsigned char _ycbcr[3]);
void od_state_draw_mv_grid(od_state *_state);
void od_state_draw_mvs(od_state *_state);
void od_state_fill_vis(od_state *_state);
#endif

/*Shared accelerated functions.*/

/*Default pure-C implementations.*/
void od_mc_predict1imv8_c(unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],const int _m[4],int _r,int _log_xblk_sz,
 int _log_yblk_sz);
void od_mc_predict1fmv8_c(unsigned char *_dst,const unsigned char *_src,
 int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz);
void od_mc_blend_full8_c(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz);
void od_mc_blend_full_split8_c(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz);

void od_state_opt_vtbl_init_c(od_state *_state);

#endif
