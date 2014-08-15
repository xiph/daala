/*Daala video codec
Copyright (c) 2002-2012 Daala project contributors.  All rights reserved.

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if !defined(_LARGEFILE_SOURCE)
# define _LARGEFILE_SOURCE
#endif
#if !defined(_LARGEFILE64_SOURCE)
# define _LARGEFILE64_SOURCE
#endif
#if !defined(_FILE_OFFSET_BITS)
# define _FILE_OFFSET_BITS 64
#endif
#if !defined(_BSD_SOURCE)
# define _BSD_SOURCE
#endif
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#if !defined(_WIN32)
# include <getopt.h>
#else
# include <fcntl.h>
# include <io.h>
# include "getopt.h"
#endif
#include <png.h>
#include <zlib.h>
#include <ogg/os_types.h>
#include "vidinput.h"

#define OD_MINI(_a,_b)      ((_a)<(_b)?(_a):(_b))
#define OD_MAXI(_a,_b)      ((_a)>(_b)?(_a):(_b))
#define OD_CLAMPI(_a,_b,_c) (OD_MAXI(_a,OD_MINI(_b,_c)))
#define OD_SIGNMASK(_a)     (-((_a)<0))
#define OD_FLIPSIGNI(_a,_b) ((_a)+OD_SIGNMASK(_b)^OD_SIGNMASK(_b))
#define OD_DIV_ROUND(_x,_y) (((_x)+OD_FLIPSIGNI((_y)>>1,_x))/(_y))
#define OD_CLAMP255(_x)     ((unsigned char)((((_x)<0)-1)&((_x)|-((_x)>255))))

enum{
  PIXEL_FMT_420,
  PIXEL_FMT_RESERVED,
  PIXEL_FMT_422,
  PIXEL_FMT_444
};

static const char *output_filename=NULL;

const char *OPTSTRING="ho:";
struct option OPTIONS[]={
  {"help",no_argument,NULL,'h'},
  {"output",required_argument,NULL,'o'},
  {NULL,0,NULL,0}
};

void **od_malloc_2d(size_t _height,size_t _width,size_t _sz){
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

void od_free_2d(void *_ptr){
  _ogg_free(_ptr);
}



static void usage(const char *_argv0){
  fprintf(stderr,
   "Usage: %s [options] <input>\n\n"
   "The <input> argument uses C printf format to represent a list of files,\n"
   "  i.e. file-%%05d.png to look for files file00001.png to file99999.png.\n\n"
   "Options: \n\n"
   "  -h --help                       Display this help and exit.\n"
   "  -o --output <filename.png>      Output file name (required).\n"
   "                                  This uses a C printf format string to\n"
   "                                  represent a list of files, i.e.\n"
   "                                  file%%05d to write files file00001.png\n"
   "                                  to file99999.png.\n",
   /*TODO: Currently only rec709 is supported.
     Also: RGB primaries and gamma (currently ignored).
   "  --rec601                        Use ITU-R BT.601 matrix.\n"
   "  --rec709                        Use ITU-R BT.709 matrix.\n"*/
   _argv0);
}


static void ycbcr_to_rgb(png_bytep *_image,const video_input_info *_info,
 video_input_ycbcr _ycbcr){
  unsigned char *y_row;
  unsigned char *cb_row;
  unsigned char *cr_row;
  unsigned char *y;
  unsigned char *cb;
  unsigned char *cr;
  int            y_stride;
  int            cb_stride;
  int            cr_stride;
  int            width;
  int            height;
  int            hshift;
  int            vshift;
  int            pic_x;
  int            pic_y;
  int            i;
  int            j;
  width=_info->pic_w;
  height=_info->pic_h;
  hshift=!(_info->pixel_fmt&1);
  vshift=!(_info->pixel_fmt&2);
  y_stride=_ycbcr[0].stride;
  cb_stride=_ycbcr[1].stride;
  cr_stride=_ycbcr[2].stride;
  pic_x=_info->pic_x;
  pic_y=_info->pic_y;
  y_row=_ycbcr[0].data+pic_y*y_stride;
  cb_row=_ycbcr[1].data+(pic_y>>vshift)*cb_stride;
  cr_row=_ycbcr[2].data+(pic_y>>vshift)*cr_stride;
  /*Chroma up-sampling is just done with a box filter.
    This is very likely what will actually be used in practice on a real
     display, and also removes one more layer to search in for the source of
     artifacts.
    As an added bonus, it's dead simple.*/
  for(j=0;j<height;j++){
    int dc;
    y=y_row+_info->pic_x;
    cb=cb_row+(pic_x>>hshift);
    cr=cr_row+(pic_x>>hshift);
    for(i=0;i<6*width;){
      int64_t  yval;
      int64_t  cbval;
      int64_t  crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      yval=*y-16;
      cbval=*cb-128;
      crval=*cr-128;
      /*This is intentionally slow and very accurate.*/
      rval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval+4490222169144LL*crval,9745792000LL),65535);
      gval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval-534117096223LL*cbval-1334761232047LL*crval,
       9745792000LL),65535);
      bval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval+5290866304968LL*cbval,9745792000LL),65535);
      _image[j][i++]=(unsigned char)(rval>>8);
      _image[j][i++]=(unsigned char)(rval&0xFF);
      _image[j][i++]=(unsigned char)(gval>>8);
      _image[j][i++]=(unsigned char)(gval&0xFF);
      _image[j][i++]=(unsigned char)(bval>>8);
      _image[j][i++]=(unsigned char)(bval&0xFF);
      dc=y-y_row&1|1-hshift;
      y++;
      cb+=dc;
      cr+=dc;
    }
    dc=-(pic_y+j&1|1-vshift);
    y_row+=y_stride;
    cb_row+=dc&cb_stride;
    cr_row+=dc&cr_stride;
  }
}

static void png_write(png_structp _png,png_bytep _data,png_size_t _sz){
  size_t ret;
  ret=fwrite(_data,_sz,1,(FILE *)png_get_io_ptr(_png));
  if(ret!=1)png_error(_png,"Write Error");
}

static void png_flush(png_structp _png){
  fflush((FILE *)png_get_io_ptr(_png));
}

static int write_png(FILE *_fout,const video_input_info *_info,video_input_ycbcr _ycbcr){
  /*Dump a PNG of the reconstructed image.*/
  png_structp    png;
  png_infop      info;
  png_bytep     *image;
  int            width;
  int            height;
  width=_info->pic_w;
  height=_info->pic_h;
  image=(png_bytep *)od_malloc_2d(height,6*width,sizeof(**image));
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
  ycbcr_to_rgb(image,_info,_ycbcr);
  png_set_write_fn(png,_fout,png_write,png_flush);
  png_set_compression_level(png,Z_BEST_COMPRESSION);
  png_set_IHDR(png,info,width,height,16,PNG_COLOR_TYPE_RGB,
   PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  /*Hard-coded to Rec 709 parameters.*/
  /*Setting gamma causes display in Firefox to be very wrong.
    Perhaps it's not doing what I think it should be doing.
    Setting cHRM doesn't seem to be harmful.
  png_set_gAMA(png,info,2.5);*/
  png_set_cHRM_fixed(png,info,31271,32902,
   64000,33000,30000,60000,15000,6000);
  /*switch(_info->colorspace){
    case TH_CS_ITU_REC_470M:{
      png_set_gAMA(png,info,2.2);
      png_set_cHRM_fixed(png,info,31006,31616,
       67000,32000,21000,71000,14000,8000);
    }break;
    case TH_CS_ITU_REC_470BG:{
      png_set_gAMA(png,info,2.67);
      png_set_cHRM_fixed(png,info,31271,32902,
       64000,33000,29000,60000,15000,6000);
    }break;
    default:break;
  }*/
  /*Dodgy hack for non-square pixels.*/
  png_set_pHYs(png,info,_info->par_n,_info->par_d,0);
  png_set_rows(png,info,image);
  png_write_png(png,info,PNG_TRANSFORM_IDENTITY,NULL);
  png_write_end(png,info);
  png_destroy_write_struct(&png,&info);
  od_free_2d(image);
  return 0;
}

static FILE *open_png_file(const char *_name,int _frameno){
  FILE *fout;
  char  output_filename[8192];
  snprintf(output_filename,8191,_name,_frameno);
  fout=fopen(output_filename,"wb");
  if(fout==NULL){
    fprintf(stderr,"Error opening output file \"%s\": %s\n",
     output_filename,strerror(errno));
  }
  fprintf(stderr,"%s\n",output_filename);
  return fout;
}

int main(int _argc,char **_argv){
  video_input vid;
  video_input_info info;
  video_input_ycbcr ycbcr;
  FILE *fout;
  FILE *fin;
  const char *input_filename;
  int i;
#ifdef _WIN32
  /*We need to set stdin/stdout to binary mode.
    Damn Windows.*/
  /*Beware the evil ifdef.
    We avoid these where we can, but this one we cannot.
    Don't add any more, you'll probably go to hell if you do.*/
  _setmode(_fileno(stdin),_O_BINARY);
  _setmode(_fileno(stdout),_O_BINARY);
#endif
  for(;;){
    int long_option_index;
    int c;
    c=getopt_long(_argc,_argv,OPTSTRING,OPTIONS,&long_option_index);
    if(c==EOF)break;
    switch(c){
      case 'h':{
        usage(_argv[0]);
        return EXIT_SUCCESS;
      }
      case 'o':{
        output_filename=optarg;
        break;
      }
      default:{
        usage(_argv[0]);
        return EXIT_FAILURE;
      }
    }
  }
  if(_argc<3){
    usage(_argv[0]);
    return EXIT_FAILURE;
  }
  input_filename=_argv[optind];
  if(input_filename==NULL){
    fprintf(stderr,"No input file specified. Run with -h for help.\n");
    return EXIT_FAILURE;
  }
  fin=strcmp(input_filename,"-")==0?stdin:fopen(input_filename,"rb");
  if(fin==NULL){
    fprintf(stderr,"Could not open input file \"%s\": %s\n",
     input_filename,strerror(errno));
    return EXIT_FAILURE;
  }
  if(video_input_open(&vid,fin)<0)return EXIT_FAILURE;
  if(output_filename==NULL){
    fprintf(stderr,"No output file specified. Run with -h for help.\n");
    return EXIT_FAILURE;
  }
  video_input_get_info(&vid,&info);
  for(i=0;video_input_fetch_frame(&vid,ycbcr,NULL)>0;i++){
    fout=strcmp(output_filename,"-")==0?
     stdout:open_png_file(output_filename,i+1);
    if(write_png(fout,&info,ycbcr)<0){
      fclose(fout);
      break;
    }
    if(fout==stdout)break;
    fclose(fout);
  }
  video_input_close(&vid);
  return EXIT_SUCCESS;
}
