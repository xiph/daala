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
typedef struct od_yuv_dumpfile   od_yuv_dumpfile;
typedef struct od_adapt_ctx      od_adapt_ctx;

# include <stdio.h>
# include "internal.h"
# include "dct.h"
# include "mc.h"
# include "pvq.h"
# include "adapt.h"
# include "generic_code.h"
#include "intra.h"

/* Normalized RDO lambda used for the block size decision
   (different from quantization lambda). Tuned on subset1 and ntt-short
   although it's hard to tune since FASTSSIM and PSNRHVS go in opposite
   directions (did some visual inspection). */
#define OD_BS_RDO_LAMBDA (.08)

extern const od_coeff OD_DC_RES[3];

extern const od_coeff OD_DC_QM[2][OD_NBSIZES - 1][2];

extern const int OD_HAAR_QM[2][OD_LOG_BSIZE_MAX];

/*Adaptation speed of scalar Laplace encoding.*/
# define OD_SCALAR_ADAPT_SPEED (4)

#define OD_MAX_CODED_REFS (2)

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
/*The vector offsets in the X direction for each motion compensation block
   vertex from the upper-left.*/
# define OD_VERT_DX (OD_VERT_D+1)
/*The vector offsets in the Y direction for each motion compensation block
   vertex from the upper-left.*/
# define OD_VERT_DY (OD_VERT_D+0)
extern const int *const OD_VERT_SETUP_DX[4][4];
extern const int *const OD_VERT_SETUP_DY[4][4];

/*This should be a power of 2, and at least 8.*/
# define OD_UMV_PADDING (32)

/*The shared (encoder and decoder) functions that have accelerated variants.*/
struct od_state_opt_vtbl{
  void (*mc_predict1fmv8)(unsigned char *_dst, const unsigned char *_src,
   int _systride, int32_t _mvx, int32_t _mvy,
   int _log_xblk_sz, int _log_yblk_sz);
  void (*mc_blend_full8)(unsigned char *_dst, int _dystride,
   const unsigned char *_src[4], int _log_xblk_sz, int _log_yblk_sz);
  void (*mc_blend_full_split8)(unsigned char *_dst, int _dystride,
   const unsigned char *_src[4], int _c, int _s,
   int _log_xblk_sz, int _log_yblk_sz);
  void (*restore_fpu)(void);
  od_dct_func_2d fdct_2d[OD_NBSIZES + 1];
  od_dct_func_2d idct_2d[OD_NBSIZES + 1];
};

# if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
struct od_yuv_dumpfile{
  char tag[16];
  FILE *fd;
};
# endif

struct od_adapt_ctx {
  /* Support for PVQ encode/decode */
  int                 pvq_adapt[2*OD_NBSIZES*OD_NSB_ADAPT_CTXS];
  int                 pvq_k1_increment;
  /* CDFs are size 16 despite the fact that we're using less than that. */
  uint16_t        pvq_k1_cdf[4][16];
  generic_encoder     pvq_param_model[3];
  int                 pvq_ext[OD_NBSIZES*PVQ_MAX_PARTITIONS];
  int                 pvq_exg[OD_NPLANES_MAX][OD_NBSIZES][PVQ_MAX_PARTITIONS];
  int                 pvq_gaintheta_increment;
  uint16_t        pvq_gaintheta_cdf[2*OD_NBSIZES*PVQ_MAX_PARTITIONS][16];
  int                 pvq_skip_dir_increment;
  uint16_t        pvq_skip_dir_cdf[2*(OD_NBSIZES-1)][7];
  /* Motion vectors */
  generic_encoder     mv_model;
  uint16_t        mv_ref_cdf[5][16];
  int                 mv_ex[OD_MC_NLEVELS];
  int                 mv_ey[OD_MC_NLEVELS];
  uint16_t        mv_small_cdf[5][16];
  int                 mv_small_increment;
  uint16_t        split_flag_cdf[OD_MC_LEVEL_MAX][9][2];
  int                 split_flag_increment;

  generic_encoder model_dc[OD_NPLANES_MAX];
  generic_encoder model_g[OD_NPLANES_MAX];

  int ex_sb_dc[OD_NPLANES_MAX];
  int ex_dc[OD_NPLANES_MAX][OD_NBSIZES][3];
  int ex_g[OD_NPLANES_MAX][OD_NBSIZES];

  /* Joint skip flag for DC and AC */
  uint16_t skip_cdf[OD_NBSIZES*2][5];
  int skip_increment;
  uint16_t haar_coeff_cdf[15*3*(OD_NBSIZES + 1)][16];
  int haar_coeff_increment;
  uint16_t haar_split_cdf[15*2*5][16];
  int haar_split_increment;
  uint16_t haar_bits_cdf[3][16];
  int haar_bits_increment;
  uint16_t clpf_cdf[4][2];
  int clpf_increment;
};

struct od_state{
  od_adapt_ctx        adapt;
  daala_info          info;
  OD_ALIGN16(unsigned char mc_buf[5][OD_MVBSIZE_MAX*OD_MVBSIZE_MAX]);
  od_state_opt_vtbl   opt_vtbl;
  uint32_t        cpu_flags;
  int32_t         frame_width;
  int32_t         frame_height;
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
  int64_t         cur_time;
  od_mv_grid_pt **mv_grid;

  /** Number of horizontal motion-vector blocks. */
  int                 nhmvbs;
  /** Number of vertical motion-vector blocks. */
  int                 nvmvbs;
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
  od_coeff           *(sb_dc_mem[OD_NPLANES_MAX]);
  int                 mv_res; /* 0: 1/8, 1:1/4, 2: 1/2 pel */
# if defined(OD_DUMP_IMAGES) || defined(OD_DUMP_RECONS)
  int                 dump_tags;
  od_yuv_dumpfile    *dump_files;
# endif
# if defined(OD_DUMP_IMAGES)
  od_img              vis_img;
  od_img              tmp_vis_img;
#  if defined(OD_ANIMATE)
  int                 ani_iter;
#  endif
# endif
  od_coeff *ctmp[OD_NPLANES_MAX];
  od_coeff *dtmp[OD_NPLANES_MAX];
  od_coeff *mctmp[OD_NPLANES_MAX];
  od_coeff *mdtmp[OD_NPLANES_MAX];
  od_coeff *ltmp[OD_NPLANES_MAX];
  od_coeff *lbuf[OD_NPLANES_MAX];
  unsigned char pvq_qm_q4[OD_NPLANES_MAX][OD_QM_SIZE];
  /*These flags provide context for the CLP filter.*/
  unsigned char *clpf_flags;
  unsigned char *sb_skip_flags;
};

int od_state_init(od_state *_state, const daala_info *_info);
void od_state_clear(od_state *_state);

void od_img_copy(od_img* dest, od_img* src);
void od_adapt_ctx_reset(od_adapt_ctx *state, int is_keyframe);
void od_state_set_mv_res(od_state *state, int mv_res);
void od_state_pred_block_from_setup(od_state *_state, unsigned char *_buf,
 int _ystride, int _pli, int _vx, int _vy, int _c, int _s,
 int _log_mvb_sz);
void od_state_pred_block(od_state *_state, unsigned char *_buf, int _ystride,
 int _pli, int _vx, int _vy, int _log_mvb_sz);
void od_state_mc_predict(od_state *_state);
void od_state_init_border(od_state *_state);
void od_state_upsample8(od_state *_state, od_img *_dst, const od_img *_src);
int od_state_dump_yuv(od_state *_state, od_img *_img, const char *_tag);
void od_img_edge_ext(od_img* src);
# if defined(OD_DUMP_IMAGES)
int od_state_dump_img(od_state *_state, od_img *_img, const char *_tag);
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
void od_mc_predict1fmv8_c(unsigned char *_dst, const unsigned char *_src,
 int _systride, int32_t _mvx, int32_t _mvy,
 int _log_xblk_sz, int _log_yblk_sz);
void od_mc_blend_full8_c(unsigned char *_dst, int _dystride,
 const unsigned char *_src[4], int _log_xblk_sz, int _log_yblk_sz);
void od_mc_blend_full_split8_c(unsigned char *_dst, int _dystride,
 const unsigned char *_src[4], int _c, int _s, int _log_xblk_sz,
 int _log_yblk_sz);
void od_restore_fpu(od_state *state);

void od_state_opt_vtbl_init_c(od_state *_state);

#endif
