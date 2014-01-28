/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

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

#if !defined(_tools_init_intra_tools_H)
# define _tools_init_intra_tools_H (1)
#include "vidinput.h"
#include "../src/odintrin.h"

#define B_SZ_LOG (2)
#define B_SZ     (1<<B_SZ_LOG)

#define B_SZ_LOG_MAX (OD_LOG_BSIZE0+OD_NBSIZES-1)
#define B_SZ_MAX (OD_BSIZE_MAX)

#define OD_INTRA_DC (0)
#define OD_INTRA_TM (1)
#define OD_INTRA_HU (2)
#define OD_INTRA_HE (3)
#define OD_INTRA_HD (4)
#define OD_INTRA_RD (5)
#define OD_INTRA_VR (6)
#define OD_INTRA_VE (7)
#define OD_INTRA_VL (8)
#define OD_INTRA_LD (9)
#define OD_INTRA_NMODES (10)



typedef int (*plane_start_func)(void *_ctx,const char *_name,
 const video_input_info *_info,int _pli,int _nxblocks,int _nyblocks);
typedef void (*block_func)(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj);
typedef int (*plane_finish_func)(void *_ctx);



/*Computes the starting offset and number of blocks which can be intra
   predicted with full context (i.e., all of their neighbors) available.*/
void get_intra_dims(const video_input_info *_info,int _pli,int _padding,
 int *_x0,int *_y0,int *_nxblocks,int *_nyblocks);

char *get_map_filename(const char *_name,int _pli,int _nxblocks,int _nyblocks);
char *get_weights_filename(const char *_name,
 int _pli,int _nxblocks,int _nyblocks);

int apply_to_blocks(void *_ctx,plane_start_func _start,block_func _block,
 plane_finish_func _finish,int _argc,const char **_argv);
int apply_to_blocks2(void *_ctx,int _padding,plane_start_func _start,
 const block_func *_blocks,int _nfuncs,plane_finish_func _finish,int _plmask,
 int _argc,const char **_argv);

void vp8_intra_predict(unsigned char *_dst,int _dst_stride,
 const unsigned char *_src,int _src_stride,int _mode);

#endif
