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

#if !defined(_tools_image_tools_H)
# define _tools_image_tools_H (0)

#include <stdio.h>
#include "od_defs.h"
#include "image.h"
#include "stats_tools.h"

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
  int              b_sz_log;
  unsigned char   *mask;
  unsigned char   *mode;
  double          *weight;
  od_coeff        *pre;
  int              pre_stride;
  od_coeff        *fdct;
  int              fdct_stride;
#if TF_BLOCKS
  od_coeff        *tf;
#endif
  double          *pred;
  int              pred_stride;
  od_coeff        *idct;
  int              idct_stride;
  od_coeff        *post;
  int              post_stride;
};

void image_data_init(image_data *_this,const char *_name,int _b_sz_log,
 int _nxblocks,int _nyblocks);
void image_data_clear(image_data *_this);
void image_data_mask(image_data *_this,const unsigned char *_data,int _stride);
void image_data_pre_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj);
void image_data_fdct_block(image_data *_this,int _bi,int _bj);
#if TF_BLOCKS
void image_data_tf_block(image_data *_this,int _bi,int _bj);
#endif
void image_data_print_block(image_data *_this,int _bi,int _bj,FILE *_fp);
void image_data_load_block(image_data *_this,int _bi,int _bj,
 od_coeff _coeffs[5*B_SZ*B_SZ]);
void image_data_pred_block(image_data *_this,int _bi,int _bj);
void image_data_stats_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,intra_stats *_stats);
void image_data_idct_block(image_data *_this,int _bi,int _bj);
void image_data_post_block(image_data *_this,int _bi,int _bj);
void image_data_files_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,image_files *_files);

int image_data_save_map(image_data *_this);
int image_data_load_map(image_data *_this);

#endif
