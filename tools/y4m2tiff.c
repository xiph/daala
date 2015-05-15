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
# define _DEFAULT_SOURCE
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
#include <tiffio.h>
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
   "  i.e. file-%%05d.tiff to look for files file00001.tiff to file99999.tiff.\n\n"
   "Options: \n\n"
   "  -h --help                       Display this help and exit.\n"
   "  -o --output <filename.tiff>      Output file name (required).\n"
   "                                  This uses a C printf format string to\n"
   "                                  represent a list of files, i.e.\n"
   "                                  file%%05d to write files file00001.tiff\n"
   "                                  to file99999.tiff.\n",
   /*TODO: Currently only rec709 is supported.
     Also: RGB primaries and gamma (currently ignored).
   "  --rec601                        Use ITU-R BT.601 matrix.\n"
   "  --rec709                        Use ITU-R BT.709 matrix.\n"*/
   _argv0);
}

static void ycbcr_to_rgb(uint8_t* *_image,const video_input_info *_info,
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
  int            xstride;
  width=_info->pic_w;
  height=_info->pic_h;
  xstride=(_info->depth>8)?2:1;
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
      int      extrabits;
      extrabits=_info->depth-8;
      if(_info->depth<=8){
        yval=*y-16;
        cbval=*cb-128;
        crval=*cr-128;
      }
      else {
        yval = *y + (*(y+1)<<8);
        cbval = *cb + (*(cb+1)<<8);
        crval = *cr + (*(cr+1)<<8);
        yval=(yval-(16<<extrabits));
        cbval=(cbval-(128<<extrabits));
        crval=(crval-(128<<extrabits));
      }
      /*This is intentionally slow and very accurate.*/
      rval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval+4490222169144LL*crval,9745792000LL<<extrabits),
       65535);
      gval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval-534117096223LL*cbval-1334761232047LL*crval,
       9745792000LL<<extrabits),65535);
      bval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval+5290866304968LL*cbval,9745792000LL<<extrabits),
       65535);
      _image[j][i++]=(unsigned char)(rval&0xFF);
      _image[j][i++]=(unsigned char)(rval>>8);
      _image[j][i++]=(unsigned char)(gval&0xFF);
      _image[j][i++]=(unsigned char)(gval>>8);
      _image[j][i++]=(unsigned char)(bval&0xFF);
      _image[j][i++]=(unsigned char)(bval>>8);
      dc=y-y_row&(xstride)|1-hshift;
      y+=xstride;
      cb+=dc;
      cr+=dc;
    }
    dc=-(pic_y+j&1|1-vshift);
    y_row+=y_stride;
    cb_row+=dc&cb_stride;
    cr_row+=dc&cr_stride;
  }
}

static int write_tiff(TIFF *_fout,const video_input_info *_info,video_input_ycbcr _ycbcr){
  uint8_t      **image;
  int            width;
  int            height;
  int row;
  width=_info->pic_w;
  height=_info->pic_h;
  image=(uint8_t **)od_malloc_2d(height,6*width,sizeof(**image));
  ycbcr_to_rgb(image,_info,_ycbcr);
  for (row = 0; row < height; row++) {
    TIFFWriteScanline(_fout, image[row], row, 0);
  }
  od_free_2d(image);
  return 0;
}

int main(int _argc,char **_argv){
  video_input vid;
  video_input_info info;
  video_input_ycbcr ycbcr;
  TIFF *fout;
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
    fout=TIFFOpen(output_filename, "w");
    TIFFSetField(fout, TIFFTAG_IMAGELENGTH, info.pic_h);
    TIFFSetField(fout, TIFFTAG_IMAGEWIDTH, info.pic_w);
    TIFFSetField(fout, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(fout, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField(fout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(fout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(fout, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    write_tiff(fout,&info,ycbcr);
    TIFFClose(fout);
  }
  video_input_close(&vid);
  return EXIT_SUCCESS;
}
