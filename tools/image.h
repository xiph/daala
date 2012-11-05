/*Daala video codec
Copyright (c) 2004-2012 Daala project contributors.  All rights reserved.

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

#if !defined(_image_H)
# define _image_H (1)
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

typedef unsigned short od_rgba16_pixel[4];

typedef struct od_rgba16_image od_rgba16_image;



struct od_rgba16_image{
  od_rgba16_pixel *data;
  int              stride;
  int              width;
  int              height;
};


int od_rgba16_image_read_png(od_rgba16_image *_this,FILE *_fin);
void od_rgba16_image_init(od_rgba16_image *_this,int _width,int _height);

void od_rgba16_image_draw_point(od_rgba16_image *_this,int _x,int _y,
 const od_rgba16_pixel _color);
void od_rgba16_image_draw_line(od_rgba16_image *_this,
 int _x0,int _y0,int _x1,int _y1,const od_rgba16_pixel _color);

int od_rgba16_image_write_png(const od_rgba16_image *_this,FILE *_fout);

void od_rgba16_image_clear(od_rgba16_image *_this);

#endif
