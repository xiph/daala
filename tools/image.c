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

#include "image.h"
#include <stdlib.h>
#include <errno.h>
#include <png.h>
#include <zlib.h>
#include <ogg/os_types.h>
#include <limits.h>


static void **od_malloc_2d(size_t _height,size_t _width,size_t _sz){
  size_t  rowsz;
  size_t  colsz;
  size_t  datsz;
  char   *ret;
  colsz=_height*sizeof(void *);
  rowsz=_sz*_width;
  datsz=rowsz*_height;
  /*Alloc array and row pointers.*/
  ret=(char *)_ogg_malloc(datsz+colsz);
  /*Initialize the array.*/
  if(ret!=NULL){
    size_t   i;
    void   **p;
    char    *datptr;
    p=(void **)ret;
    i=_height;
    for(datptr=ret+colsz;i-->0;p++,datptr+=rowsz)*p=(void *)datptr;
  }
  return (void **)ret;
}

static void od_free_2d(void *_ptr){
  _ogg_free(_ptr);
}


static void png_read(png_structp _png,png_bytep _data,png_size_t _sz){
  size_t ret;
  ret=fread(_data,_sz,1,(FILE *)png_get_io_ptr(_png));
  if(ret!=1)png_error(_png,"Read Error");
}

int od_rgba16_image_read_png(od_rgba16_image *_this,FILE *_fin){
  unsigned char    header[8];
  png_structp      png;
  png_infop        info;
  png_infop        end;
  od_rgba16_pixel *data;
  png_bytep       *rows;
  png_color_16p    bkgd;
  png_uint_32      width;
  png_uint_32      height;
  int              bit_depth;
  int              color_type;
  int              interlace_type;
  int              compression_type;
  int              filter_method;
  png_uint_32      i;
  png_uint_32      j;
  if(fread(header,8,1,_fin)<1){
    fprintf(stderr,"Error reading from file.\n");
    return -EINVAL;
  }
  if(png_sig_cmp(header,0,8)){
    fprintf(stderr,"Error: Not a PNG.\n");
    return -EINVAL;
  }
  png=png_create_read_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  if(png==NULL){
    fprintf(stderr,"Error: %s\n",strerror(ENOMEM));
    return -ENOMEM;
  }
  info=png_create_info_struct(png);
  if(info==NULL){
    fprintf(stderr,"Error: %s\n",strerror(ENOMEM));
    png_destroy_read_struct(&png,NULL,NULL);
    return -ENOMEM;
  }
  end=png_create_info_struct(png);
  if(end==NULL){
    fprintf(stderr,"Error: %s\n",strerror(ENOMEM));
    png_destroy_read_struct(&png,&info,NULL);
    return -EINVAL;
  }
  rows=NULL;
  data=NULL;
  if(setjmp(png_jmpbuf(png))){
    png_free(png,rows);
    free(data);
    png_destroy_read_struct(&png,&info,&end);
    return -EINVAL;
  }
  png_set_read_fn(png,_fin,png_read);
  png_set_sig_bytes(png,8);
  png_read_info(png,info);
  png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,
   &interlace_type,&compression_type,&filter_method);
  if(width<=0||height<=0||width>INT_MAX||height>INT_MAX||
    width*(png_size_t)height*8>INT_MAX||
    width*(png_size_t)height*8/height/8!=width){
    png_destroy_read_struct(&png,&info,&end);
    return -EINVAL;
  }
  png_set_expand(png);
  if(bit_depth<8)png_set_packing(png);
  if(!(color_type&PNG_COLOR_MASK_COLOR))png_set_gray_to_rgb(png);
  if(png_get_bKGD(png,info,&bkgd)){
    png_set_background(png,bkgd,PNG_BACKGROUND_GAMMA_FILE,1,1.0);
  }
  /*Add an alpha channel if the image doesn't contain one.*/
  png_set_add_alpha(png,0xFFFF,PNG_FILLER_AFTER);
  data=(od_rgba16_pixel *)_ogg_malloc(height*width*sizeof(*data));
  rows=(png_bytep *)png_malloc(png,height*sizeof(*rows));
  for(j=0;j<height;j++)rows[j]=(png_bytep)(data+j*width);
  png_read_image(png,rows);
  png_read_end(png,end);
  if(bit_depth<16){
    /*If the image wasn't 16-bit, expand it so it is.
      We do this by duplicating the high byte, so that, e.g., 0 maps to 0 and
       255 maps to 65535.
      This is also nicely endian-independent.*/
    for(j=0;j<height;j++){
      for(i=4*width;i-->0;)rows[j][2*i]=rows[j][2*i+1]=rows[j][i];
    }
  }
  else{
    /*Otherwise, convert from big-endian to host-endian byte order.
      Instead of bothering to try to detect endianess, we just always convert.
      Since almost all hosts will be little-endian in practice, there's little
       to be gained, and this is guaranteed to work portably.*/
    for(j=0;j<height;j++){
      for(i=width;i-->0;){
        unsigned short r;
        unsigned short g;
        unsigned short b;
        unsigned short a;
        r=(unsigned short)(rows[j][8*i+0]<<8|rows[j][8*i+1]);
        g=(unsigned short)(rows[j][8*i+2]<<8|rows[j][8*i+3]);
        b=(unsigned short)(rows[j][8*i+4]<<8|rows[j][8*i+5]);
        a=(unsigned short)(rows[j][8*i+6]<<8|rows[j][8*i+7]);
        *(data+j*width+i)[0]=r;
        *(data+j*width+i)[1]=g;
        *(data+j*width+i)[2]=b;
        *(data+j*width+i)[3]=a;
      }
    }
  }
  png_free(png,rows);
  png_destroy_read_struct(&png,&info,&end);
  _this->data=data;
  _this->stride=(int)width;
  _this->width=(int)width;
  _this->height=(int)height;
  return 0;
}

void od_rgba16_image_init(od_rgba16_image *_this,int _width,int _height){
  od_rgba16_pixel *data;
  int              i;
  int              j;
  data=(od_rgba16_pixel *)_ogg_malloc(_height*_width*sizeof(*data));
  for(j=0;j<_height;j++){
    for(i=0;i<_width;i++){
      (*(data+j*_width+i))[0]=0;
      (*(data+j*_width+i))[1]=0;
      (*(data+j*_width+i))[2]=0;
      (*(data+j*_width+i))[3]=(unsigned short)0xFFFFU;
    }
  }
  _this->data=data;
  _this->stride=(int)_width;
  _this->width=(int)_width;
  _this->height=(int)_height;
}

static void png_write(png_structp _png,png_bytep _data,png_size_t _sz){
  size_t ret;
  ret=fwrite(_data,_sz,1,(FILE *)png_get_io_ptr(_png));
  if(ret!=1)png_error(_png,"Write Error");
}

static void png_flush(png_structp _png){
  fflush((FILE *)png_get_io_ptr(_png));
}

int od_rgba16_image_write_png(const od_rgba16_image *_this,FILE *_fout){
  png_structp      png;
  png_infop        info;
  png_bytep       *image;
  od_rgba16_pixel *data;
  int              width;
  int              height;
  int              stride;
  int              i;
  int              j;
  width=_this->width;
  height=_this->height;
  image=(png_bytep *)od_malloc_2d(height,8*width,sizeof(**image));
  if(image==NULL)return -EFAULT;
  png=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  if(png==NULL){
    od_free_2d(image);
    return -EFAULT;
  }
  info=png_create_info_struct(png);
  if(info==NULL){
    png_destroy_write_struct(&png,NULL);
    od_free_2d(image);
    return -EFAULT;
  }
  if(setjmp(png_jmpbuf(png))){
    png_destroy_write_struct(&png,&info);
    od_free_2d(image);
    return -EFAULT;
  }
  data=_this->data;
  stride=_this->stride;
  for(j=0;j<height;j++){
    for(i=0;i<width;i++){
      unsigned short r;
      unsigned short g;
      unsigned short b;
      unsigned short a;
      r=(*(data+j*stride+i))[0];
      g=(*(data+j*stride+i))[1];
      b=(*(data+j*stride+i))[2];
      a=(*(data+j*stride+i))[3];
      image[j][i*8+0]=(png_byte)(r>>8);
      image[j][i*8+1]=(png_byte)(r&0xFF);
      image[j][i*8+2]=(png_byte)(g>>8);
      image[j][i*8+3]=(png_byte)(g&0xFF);
      image[j][i*8+4]=(png_byte)(b>>8);
      image[j][i*8+5]=(png_byte)(b&0xFF);
      image[j][i*8+6]=(png_byte)(a>>8);
      image[j][i*8+7]=(png_byte)(a&0xFF);
    }
  }
  png_set_write_fn(png,_fout,png_write,png_flush);
  png_set_compression_level(png,Z_BEST_COMPRESSION);
  png_set_IHDR(png,info,width,height,16,PNG_COLOR_TYPE_RGBA,
   PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  png_set_rows(png,info,image);
  png_write_png(png,info,PNG_TRANSFORM_IDENTITY,NULL);
  png_write_end(png,info);
  png_destroy_write_struct(&png,&info);
  od_free_2d(image);
  return 0;
}

void od_rgba16_image_draw_point(od_rgba16_image *_this,int _x,int _y,
 const od_rgba16_pixel _color){
  od_rgba16_pixel *data;
  int              stride;
  int              width;
  int              height;
  unsigned short   r;
  unsigned short   g;
  unsigned short   b;
  unsigned short   a;
  width=_this->width;
  height=_this->height;
  if(_x<0||_y<0||_x>=width||_y>=height)return;
  data=_this->data;
  stride=_this->stride;
  r=_color[0];
  g=_color[1];
  b=_color[2];
  a=_color[3];
  (*(data+stride*_y+_x))[0]=r;
  (*(data+stride*_y+_x))[1]=g;
  (*(data+stride*_y+_x))[2]=b;
  (*(data+stride*_y+_x))[3]=a;
}

void od_rgba16_image_draw_line(od_rgba16_image *_this,
 int _x0,int _y0,int _x1,int _y1,const od_rgba16_pixel _color){
  int x0[2];
  int x1[2];
  int dx[2];
  int step[2];
  int steep;
  int err;
  int derr;
  steep=abs(_y1-_y0)>abs(_x1-_x0);
  x0[0]=_x0;
  x0[1]=_y0;
  x1[0]=_x1;
  x1[1]=_y1;
  dx[0]=abs(x1[0]-x0[0]);
  dx[1]=abs(x1[1]-x0[1]);
  err=0;
  derr=dx[1-steep];
  step[0]=((x0[0]<x1[0])<<1)-1;
  step[1]=((x0[1]<x1[1])<<1)-1;
  od_rgba16_image_draw_point(_this,x0[0],x0[1],_color);
  while(x0[steep]!=x1[steep]){
    x0[steep]+=step[steep];
    err+=derr;
    if(err<<1>dx[steep]){
      x0[1-steep]+=step[1-steep];
      err-=dx[steep];
    }
    od_rgba16_image_draw_point(_this,x0[0],x0[1],_color);
  }
}

void od_rgba16_image_clear(od_rgba16_image *_this){
  _ogg_free(_this->data);
}
