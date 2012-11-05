/*Daala video codec
Copyright (c) 2002-2007 Daala project contributors.  All rights reserved.

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

#if !defined(_vidinput_H)
# define _vidinput_H (1)
# if !defined(_LARGEFILE_SOURCE)
#  define _LARGEFILE_SOURCE
# endif
# if !defined(_LARGEFILE64_SOURCE)
#  define _LARGEFILE64_SOURCE
# endif
# if !defined(_FILE_OFFSET_BITS)
#  define _FILE_OFFSET_BITS 64
# endif
# include <stdio.h>
# include "theora/theoradec.h"

# if defined(__cplusplus)
extern "C" {
# endif



typedef struct video_input      video_input;
typedef struct video_input_vtbl video_input_vtbl;

typedef void (*video_input_get_info_func)(void *_ctx,th_info *_ti);
typedef int (*video_input_fetch_frame_func)(void *_ctx,FILE *_fin,
 th_ycbcr_buffer _ycbcr,char _tag[5]);
typedef void (*video_input_close_func)(void *_ctx);



struct video_input_vtbl{
  video_input_get_info_func     get_info;
  video_input_fetch_frame_func  fetch_frame;
  video_input_close_func        close;
};

struct video_input{
  const video_input_vtbl *vtbl;
  void                   *ctx;
  FILE                   *fin;
};



int video_input_open(video_input *_vid,FILE *_fin);
void video_input_close(video_input *_vid);

void video_input_get_info(video_input *_vid,th_info *_ti);
int video_input_fetch_frame(video_input *_vid,
 th_ycbcr_buffer _ycbcr,char _tag[5]);

# if defined(__cplusplus)
}
# endif

#endif
