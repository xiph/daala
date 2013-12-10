/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

#if !defined(_state_H)
# define _state_H (1)

typedef struct od_state_opt_vtbl od_state_opt_vtbl;
typedef struct od_state          od_state;

# include "internal.h"
# include "mc.h"
# include "pvq_code.h"
# include "adapt.h"
# include "generic_code.h"

/*The golden reference frame.*/
# define OD_FRAME_GOLD (0)
/*The previous reference frame.*/
# define OD_FRAME_PREV (1)
/*The next reference frame.*/
# define OD_FRAME_NEXT (2)
/*The current frame.*/
# define OD_FRAME_SELF (3)

/*The reconstructed I/O frame.*/
# define OD_FRAME_REC   (0)
/*The input I/O frame.*/
# define OD_FRAME_INPUT (1)

/*Constants for the packet state machine common between encoder and decoder.*/

/*Next packet to emit/read: Codec info header.*/
# define OD_PACKET_INFO_HDR    (-3)
/*Next packet to emit/read: Comment header.*/
# define OD_PACKET_COMMENT_HDR (-2)
/*Next packet to emit/read: Codec setup header.*/
# define OD_PACKET_SETUP_HDR   (-1)
/*Next more packets to emit/read.*/
# define OD_PACKET_DONE        (INT_MAX)

extern const int        OD_VERT_D[];
/*The vector offsets in the X direction for each motion comepnsation block
   vertex from the upper-left.*/
# define OD_VERT_DX (OD_VERT_D+1)
/*The vector offsets in the Y direction for each motion compensation block
   vertex from the upper-left.*/
# define OD_VERT_DY (OD_VERT_D+0)
extern const int *const OD_VERT_SETUP_DX[4][4];
extern const int *const OD_VERT_SETUP_DY[4][4];



/*This should be a power of 2, and at least 8.*/
# define OD_UMV_PADDING (32)

# define OD_SUPERBLOCK_SIZE (32)

/*These are buckets for performance and bitrate measurements.*/
# define OD_METRIC_TOTAL 0
# define OD_METRIC_MV 1
# define OD_METRIC_BLOCK_SWITCHING 2
# define OD_METRIC_INTRA 3
# define OD_METRIC_PVQ 4
# define OD_METRIC_DC 5

# define OD_METRIC_COUNT (OD_METRIC_DC + 1)

/*The shared (encoder and decoder) functions that have accelerated variants.*/
struct od_state_opt_vtbl{
  void (*mc_predict1imv8)(unsigned char *_dst, int _dystride,
   const unsigned char *_src, int _systride, const ogg_int32_t _mvx[4],
   const ogg_int32_t _mvy[4], const int _m[4], int _r, int _log_xblk_sz,
   int _log_yblk_sz);
  void (*mc_predict1fmv8)(unsigned char *_dst, const unsigned char *_src,
   int _systride, ogg_int32_t _mvx, ogg_int32_t _mvy,
   int _log_xblk_sz, int _log_yblk_sz);
  void (*mc_blend_full8)(unsigned char *_dst, int _dystride,
   const unsigned char *_src[4], int _log_xblk_sz, int _log_yblk_sz);
  void (*mc_blend_full_split8)(unsigned char *_dst, int _dystride,
   const unsigned char *_src[4], int _c, int _s,
   int _log_xblk_sz, int _log_yblk_sz);
};



struct od_state{
  daala_info          info;
  od_state_opt_vtbl   opt_vtbl;
  ogg_uint32_t        cpu_flags;
  ogg_int32_t         frame_width;
  ogg_int32_t         frame_height;
  od_img              input;
  /** Buffer for the 4 ref images. */
  int                 ref_imgi[4];
  /** Pointers to the ref images so one can move them around without coping
      them. */
  od_img              ref_imgs[4];
  /** Pointer to input and output image. */
  od_img              io_imgs[2];
  unsigned char *ref_line_buf[8];
  unsigned char *ref_img_data;
  /** Increments by 1 for each frame. */
  ogg_int64_t         cur_time;
  od_mv_grid_pt **mv_grid;
  od_adapt_ctx        adapt_sb[OD_NPLANES_MAX];
  int                 pvq_adapt[OD_NSB_ADAPT_CTXS];
  generic_encoder     pvq_gain_model;
  int                 pvq_ext;
  int                 pvq_exg;
  /** number of horizontal macro blocks. */
  int                 nhmbs;
  /** number of vertical macro blocks. */
  int                 nvmbs;
  int                 nhsb;
  int                 nvsb;
  /** Each 8x8 block of pixels in the image (+ one superblock of
      padding on each side) has a corresponding byte in this array, and
      every 32x32 superblock is represented by 16 (4 by 4) entries
      ((4 * 8) * (4 * 8) == 32 * 32) that encode the block size decisions
      for the superblock. The entry format is:
      - 0 means the 8x8 block has been split into 4x4 blocks
      - 1 means the 8x8 block is an 8x8 block
      - 2 means the 8x8 block is part of a 16x16 block
      - 3 means the 8x8 block is part of a 32x32 block.
      The padding is filled as though it consisted of 32x32 blocks.

      E.g., `state->bsize[j * state->bstride + i]` accesses the i'th 8x8
      block in the j'th row of 8x8 blocks.

      The `bstride` member has the distance between vertically adjacent
      entries (horizontally adjacent entries are adjacent in memory). */
  unsigned char *bsize;
  int                 bstride;
  int                 mv_res;
# if defined(OD_DUMP_IMAGES)
  od_img              vis_img;
#  if defined(OD_ANIMATE)
  int                 ani_iter;
#  endif
# endif
# if defined(OD_METRICS)
  ogg_int64_t bit_metrics[OD_METRIC_COUNT];
# endif
};


int od_state_init(od_state *_state, const daala_info *_info);
void od_state_clear(od_state *_state);

void od_state_pred_block_from_setup(od_state *_state, unsigned char *_buf,
 int _ystride, int _ref, int _pli, int _vx, int _vy, int _c, int _s,
 int _log_mvb_sz);
void od_state_pred_block(od_state *_state, unsigned char *_buf, int _ystride,
 int _ref, int _pli, int _vx, int _vy, int _log_mvb_sz);
void od_state_mc_predict(od_state *_state, int _ref);
void od_state_init_border_as_32x32(od_state *_state);
void od_state_upsample8(od_state *_state, od_img *_dst, const od_img *_src);
int od_state_dump_yuv(od_state *_state, od_img *_img, const char *_suf);
# if defined(OD_DUMP_IMAGES)
int od_state_dump_img(od_state *_state, od_img *_img, const char *_suf);
void od_img_draw_point(od_img *_img, int _x, int _y,
 const unsigned char _ycbcr[3]);
void od_img_draw_line(od_img *_img, int _x0, int _y0, int _x1, int _y1,
 const unsigned char _ycbcr[3]);
void od_state_draw_mv_grid(od_state *_state);
void od_state_draw_mvs(od_state *_state);
void od_state_fill_vis(od_state *_state);
# endif

/*Shared accelerated functions.*/

/*Default pure-C implementations.*/
void od_mc_predict1imv8_c(unsigned char *_dst, int _dystride,
 const unsigned char *_src, int _systride, const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4], const int _m[4], int _r, int _log_xblk_sz,
 int _log_yblk_sz);
void od_mc_predict1fmv8_c(unsigned char *_dst, const unsigned char *_src,
 int _systride, ogg_int32_t _mvx, ogg_int32_t _mvy,
 int _log_xblk_sz, int _log_yblk_sz);
void od_mc_blend_full8_c(unsigned char *_dst, int _dystride,
 const unsigned char *_src[4], int _log_xblk_sz, int _log_yblk_sz);
void od_mc_blend_full_split8_c(unsigned char *_dst, int _dystride,
 const unsigned char *_src[4], int _c, int _s, int _log_xblk_sz,
 int _log_yblk_sz);

void od_state_opt_vtbl_init_c(od_state *_state);

#endif
