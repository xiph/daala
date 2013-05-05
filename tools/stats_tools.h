/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#if !defined(_tools_stats_tools_H)
# define _tools_stats_tools_H (1)

#include "intra_fit_tools.h"
#include "image.h"
#include "../src/intra.h"
#include "od_covmat.h"

#define INPUT_SCALE (16)

typedef struct mode_data mode_data;

struct mode_data{
  int       n;
  double    satd_avg[B_SZ*B_SZ];
  double    mean;
  double    var;
  od_covmat ref;
  od_covmat pred;
};

void mode_data_init(mode_data *_md);
void mode_data_clear(mode_data *_md);
void mode_data_add_input(mode_data *_md,const unsigned char *_data,int _stride);
void mode_data_add_block(mode_data *_md,const od_coeff *_block,int _stride,
 int _ref);
void mode_data_correct(mode_data *_md);
void mode_data_print(mode_data *_md,const char *_label,double *_scale);
void mode_data_combine(mode_data *_a,const mode_data *_b);
void mode_data_params(mode_data *_this,double _b[B_SZ*B_SZ],double *_scale);

typedef struct intra_stats intra_stats;

struct intra_stats{
  mode_data fr;
  mode_data md[OD_INTRA_NMODES];
};

void intra_stats_init(intra_stats *_this);
void intra_stats_clear(intra_stats *_this);
void intra_stats_reset(intra_stats *_this);
void intra_stats_update(intra_stats *_this,const unsigned char *_data,
 int _stride,int _mode,const od_coeff *_ref,int _ref_stride,
 const od_coeff *_pred,int _pred_stride);
void intra_stats_correct(intra_stats *_this);
void intra_stats_print(intra_stats *_this,const char *_label,double *_scale);
void intra_stats_combine(intra_stats *_this,const intra_stats *_that);

/* Is there a good way to have a static initializer in C? */
extern double VP8_SCALE[B_SZ];
extern double OD_SCALE[B_SZ];

void vp8_scale_init(double _vp8_scale[B_SZ]);
void od_scale_init(double _od_scale[B_SZ]);

int vp8_select_mode(const unsigned char *_data,int _stride,double *_weight);
int od_select_mode_satd(const od_coeff *_block,int _stride,double *_weight); 
int od_select_mode_bits(const od_coeff *_block,int _stride,double *_weight,
 double _b[OD_INTRA_NMODES][B_SZ*B_SZ]);

extern od_rgba16_pixel COLORS[OD_INTRA_NMODES];

void image_draw_block(od_rgba16_image *_image,int _x,int _y,
 const unsigned char *_block,int _stride);
int image_write_png(od_rgba16_image *_image,const char *_name);

typedef struct image_files image_files;

struct image_files{
  od_rgba16_image raw;
  od_rgba16_image map;
  od_rgba16_image pred;
  od_rgba16_image res;
};

void image_files_init(image_files *_this,int _nxblocks,int _nyblocks);
void image_files_clear(image_files *_this);
void image_files_write(image_files *_this,const char *_name,const char *_suf);

typedef struct image_data image_data;

struct image_data{
  const char      *name;
  int              nxblocks;
  int              nyblocks;
  unsigned char   *mode;
  double          *weight;
  od_coeff        *pre;
  int              pre_stride;
  od_coeff        *fdct;
  int              fdct_stride;
  od_coeff        *pred;
  int              pred_stride;
  od_coeff        *idct;
  int              idct_stride;
  od_coeff        *post;
  int              post_stride;
};

void image_data_init(image_data *_this,const char *_name,int _nxblocks,
 int _nyblocks);
void image_data_clear(image_data *_this);
void image_data_pre_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj);
void image_data_fdct_block(image_data *_this,int _bi,int _bj);
void image_data_print_block(image_data *_this,int _bi,int _bj,FILE *_fp);
void image_data_pred_block(image_data *_this,int _bi,int _bj);
void image_data_stats_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,intra_stats *_stats);
void image_data_idct_block(image_data *_this,int _bi,int _bj);
void image_data_post_block(image_data *_this,int _bi,int _bj);
void image_data_files_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,image_files *_files);

int image_data_save_map(image_data *_this);
int image_data_load_map(image_data *_this);

int ne_apply_to_blocks(void *_ctx,int _ctx_sz,int _plmask,int _padding,
 plane_start_func _start,int _nfuncs,const block_func *_funcs,
 plane_finish_func _finish,int _argc,const char *_argv[]);

#endif
