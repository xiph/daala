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

/*Adapted from png2theora, part of the OggTheora software codec source code.
  The Theora source code is copyright (C) 2002-2009
  by the Xiph.Org Foundation and contributors http://www.xiph.org/
  Based on code from Vegard Nossum.*/

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>
#include <limits.h>
#if !defined(_WIN32)
# include <getopt.h>
# include <unistd.h>
#else
# include <fcntl.h>
# include <io.h>
# include "getopt.h"
#endif
#include <libgen.h>
#include <png.h>
#include "kiss99.h"

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

static const char *CHROMA_TAGS[4]={" C420jpeg",""," C422jpeg"," C444"};

static const char *output_filename=NULL;
static int fps_numerator=24;
static int fps_denominator=1;
static int aspect_numerator=0;
static int aspect_denominator=0;
static int pixel_format=PIXEL_FMT_420;

static char *input_filter;

const char *OPTSTRING="ho:s:S:f:F:";
struct option OPTIONS[]={
  {"help",no_argument,NULL,'h'},
  {"output",required_argument,NULL,'o'},
  {"chroma-444",no_argument,NULL,'\5'},
  {"chroma-422",no_argument,NULL,'\6'},
  {"aspect-numerator",required_argument,NULL,'s'},
  {"aspect-denominator",required_argument,NULL,'S'},
  {"framerate-numerator",required_argument,NULL,'f'},
  {"framerate-denominator",required_argument,NULL,'F'},
  {NULL,0,NULL,0}
};

static void usage(const char *_argv0){
  fprintf(stderr,
   "Usage: %s [options] <input>\n\n"
   "The <input> argument uses a C printf format string to represent a list of\n"
   "files, i.e., file%%05d.png to look for files file00000.png to\n"
   "file99999.png.\n\n"
   "Options: \n\n"
   "  -h --help                       Display this help and exit.\n"
   "  -o --output <filename.y4m>      Output file name (required).\n"
   "  --chroma-444                    Use 4:4:4 chroma subsampling.\n"
   "  --chroma-422                    Use 4:2:2 chroma subsampling.\n"
   "                                  (4:2:0 is default).\n\n"
   /*TODO: Currently only rec709 is supported.
     Also: RGB primaries and gamma (currently ignored).
   "  --rec601                        Use ITU-R BT.601 matrix.\n"
   "  --rec709                        Use ITU-R BT.709 matrix.\n"*/
   "  -s --aspect-numerator <n>       Aspect ratio numerator, default is 0.\n"
   "  -S --aspect-denominator <n>     Aspect ratio denominator, default is 0.\n"
   "  -f --framerate-numerator <n>    Framerate numerator (default 24).\n"
   "                                  Determines the frame rate in frames per\n"
   "                                  second when divided by the framerate\n"
   "                                  denominator.\n"
   "  -F --framerate-denominator <n>  Frame rate denominator (default 1).\n",
   _argv0);
}

#ifdef _WIN32

int alphasort(const void *_a, const void *_b){
  return strcoll((*(const struct dirent **)_a)->d_name,
   (*(const struct dirent **)_b)->d_name);
}

typedef int (*compar_func)(const void *,const void *);

int scandir(const char *_dirp,struct dirent ***_namelist,
 int (*_filter)(const struct dirent *),
 int (*_compar)(const struct dirent *,const struct dirent *)){
  DIR            *d;
  struct dirent  *entry;
  size_t          entry_sz;
  struct dirent **names;
  int             nnames;
  int             cnames;
  d=opendir(_dirp);
  if(d==NULL)return -1;
  names=NULL;
  nnames=cnames=0;
  while((entry=readdir(d))!=NULL){
    if(_filter==NULL||(*_filter)(entry)){
      if(nnames>=cnames){
        cnames=cnames<<1|1;
        names=(struct dirent **)realloc(names,cnames*sizeof(*names));
      }
      entry_sz=sizeof(*entry)-sizeof(entry->d_name)+strlen(entry->d_name)+1;
      names[nnames]=(struct dirent *)malloc(entry_sz);
      memcpy(names[nnames],entry,entry_sz);
      nnames++;
    }
  }
  if(closedir(d)<0)return -1;
  if(nnames==0)return -1;
  if(_compar!=NULL)qsort(names,nnames,sizeof(*names),(compar_func)_compar);
  *_namelist=(struct dirent **)realloc(names,nnames*sizeof(*names));
  return nnames;
}
#endif

typedef struct img_plane img_plane;

struct img_plane{
  int            width;
  int            height;
  int            stride;
  unsigned char *data;
};

/*Generate a triangular deviate with zero mean and range [-_range,_range].*/
static int32_t triangle_rand(kiss99_ctx *_kiss,int32_t _range){
  uint32_t m;
  uint32_t r1;
  uint32_t r2;
  _range++;
  m=0xFFFFFFFFU/_range*_range;
  do r1=kiss99_rand(_kiss);
  while(r1>=m);
  do r2=kiss99_rand(_kiss);
  while(r2>=m);
  return (int32_t)(r1%_range)-(int32_t)(r2%_range);
}

/*WARNING: The constants in the following code are hard-coded for the BT.709
   matrices.*/

/*Adds triangle-shaped dither to an RGB pixel, but aligned with the Y'CbCr axes
   and scaled with the quantizer in that space.
  The dither has to be added in Y'CbCr space because that's where we're
   quantizing, but projecting it back allows us to measure error in the RGB
   space, which is especially useful for handling things like subsampling.*/
static void get_dithered_pixel(int32_t *_r,int32_t *_g,int32_t *_b,
 const png_byte *_rgb,int64_t _yd,int64_t _cbd,int64_t _crd){
  int32_t r;
  int32_t g;
  int32_t b;
  r=_rgb[0]*256+_rgb[1];
  g=_rgb[2]*256+_rgb[3];
  b=_rgb[4]*256+_rgb[5];
  r+=(int32_t)OD_DIV_ROUND(2*_yd+3*_crd,8176000);
  g+=(int32_t)OD_DIV_ROUND(2384*_yd+361*_cbd+1063*_crd,9745792000LL);
  b+=(int32_t)OD_DIV_ROUND(2*_yd+3*_cbd,8176000);
  *_r=r;
  *_g=g;
  *_b=b;
}

static int calc_y(int32_t _r,int32_t _g,int32_t _b,int _cb,int _cr){
  int64_t chroma_r;
  int64_t chroma_g;
  int64_t chroma_b;
  int64_t r_res;
  int64_t g_res;
  int64_t b_res;
  int64_t yn;
  int64_t err0;
  int64_t err1;
  int32_t r;
  int32_t g;
  int32_t b;
  int     y0;
  _cb-=128;
  _cr-=128;
  chroma_r=4490222169144LL*_cr;
  chroma_g=-534117096223LL*_cb-1334761232047LL*_cr;
  chroma_b=5290866304968LL*_cb;
  r_res=_r*9745792000LL-chroma_r+4096>>13;
  g_res=_g*9745792000LL-chroma_g+4096>>13;
  b_res=_b*9745792000LL-chroma_b+4096>>13;
  /*Take the floor here instead of rounding; we'll consider both possible
     values.*/
  yn=1063*r_res+3576*g_res+361*b_res;
  y0=(int)((yn-(1780026171874LL&OD_SIGNMASK(yn)))/1780026171875LL);
  /*Clamp before adding the offset.
    We clamp to 238 instead of 239 to ensure we can always add one and stay in
     range.*/
  y0=OD_CLAMPI(-16,y0,238);
  /*Check the reconstruction error with y0 after rounding and clamping.*/
  r=OD_CLAMPI(0,
   (int32_t)OD_DIV_ROUND(2916394880000LL*y0+chroma_r,9745792000LL),65535);
  g=OD_CLAMPI(0,
   (int32_t)OD_DIV_ROUND(2916394880000LL*y0+chroma_g,9745792000LL),65535);
  b=OD_CLAMPI(0,
   (int32_t)OD_DIV_ROUND(2916394880000LL*y0+chroma_b,9745792000LL),65535);
  err0=(_r-r)*(_r-r)+(_g-g)*(_g-g)+(_b-b)*(_b-b);
  /*Check the reconstruction error with y0+1 after rounding and clamping.*/
  r=OD_CLAMPI(0,
   (int32_t)OD_DIV_ROUND(2916394880000LL*(y0+1)+chroma_r,9745792000LL),65535);
  g=OD_CLAMPI(0,
   (int32_t)OD_DIV_ROUND(2916394880000LL*(y0+1)+chroma_g,9745792000LL),65535);
  b=OD_CLAMPI(0,
   (int32_t)OD_DIV_ROUND(2916394880000LL*(y0+1)+chroma_b,9745792000LL),65535);
  err1=(_r-r)*(_r-r)+(_g-g)*(_g-g)+(_b-b)*(_b-b);
  if(err1<err0)y0++;
  /*In the unlikely event there's a tie, round to even.*/
  else if(err1==err0)y0+=y0&1;
  return y0+16;
}



static void rgb_to_ycbcr(img_plane _ycbcr[3],png_bytep *_png){
  kiss99_ctx     kiss;
  unsigned char *ydata;
  unsigned char *cbdata;
  unsigned char *crdata;
  int            ystride;
  int            cbstride;
  int            crstride;
  int            hstep;
  int            vstep;
  int            w;
  int            h;
  int            i;
  int            j;
  w=_ycbcr[0].width;
  h=_ycbcr[0].height;
  ystride=_ycbcr[0].stride;
  ydata=_ycbcr[0].data;
  cbstride=_ycbcr[1].stride;
  cbdata=_ycbcr[1].data;
  crstride=_ycbcr[2].stride;
  crdata=_ycbcr[2].data;
  hstep=pixel_format&1;
  vstep=pixel_format&2;
  kiss99_srand(&kiss,NULL,0);
  for(j=0;j<h;j+=2){
    for(i=0;i<w;i+=2){
      int32_t yd[4];
      int32_t cbd[4];
      int32_t crd[4];
      int32_t r0;
      int32_t g0;
      int32_t b0;
      int32_t r1;
      int32_t g1;
      int32_t b1;
      int32_t r2;
      int32_t g2;
      int32_t b2;
      int32_t r3;
      int32_t g3;
      int32_t b3;
      int64_t rsum;
      int64_t gsum;
      int64_t bsum;
      int     k;
      int     cb;
      int     cr;
      /*This often generates more dither values than we use, but keeps them in
         sync for the luma plane across the different pixel formats.*/
      for(k=0;k<4;k++){
        /*The size of the dither here is chosen to be the largest divisor of
           all the corresponding coefficients in the transform that still fits
           in 31 bits.*/
        yd[k]=triangle_rand(&kiss,1223320000);
        cbd[k]=triangle_rand(&kiss,1479548743);
        crd[k]=triangle_rand(&kiss,1255654969);
      }
      get_dithered_pixel(&r0,&g0,&b0,_png[j]+6*i,yd[0],cbd[0],crd[0]);
      if(i+1<w){
        get_dithered_pixel(&r1,&g1,&b1,_png[j]+6*(i+1),
         yd[1],cbd[hstep],crd[hstep]);
      }
      else{
        r1=r0;
        g1=g0;
        b1=b0;
      }
      if(j+1<h){
        get_dithered_pixel(&r2,&g2,&b2,_png[j+1]+6*i,
         yd[2],cbd[vstep],crd[vstep]);
        if(i+1<w){
          get_dithered_pixel(&r3,&g3,&b3,_png[j+1]+6*(i+1),
           yd[3],cbd[vstep+hstep],crd[vstep+hstep]);
        }
        else{
          r3=r2;
          g3=g2;
          b3=b2;
        }
      }
      else{
        r2=r0;
        g2=g0;
        b2=b0;
        r3=r1;
        g3=g1;
        b3=b1;
      }
      if(pixel_format==PIXEL_FMT_420){
        rsum=r0+r1+r2+r3;
        gsum=g0+g1+g2+g3;
        bsum=b0+b1+b2+b3;
        cb=OD_CLAMP255(
         OD_DIV_ROUND(-29764*rsum-100128*gsum+129892*bsum,304016865)+128);
        cr=OD_CLAMP255(
         OD_DIV_ROUND(110236*rsum-100128*gsum-10108*bsum,258011295)+128);
        cbdata[(j>>1)*cbstride+(i>>1)]=(unsigned char)cb;
        crdata[(j>>1)*crstride+(i>>1)]=(unsigned char)cr;
        ydata[j*ystride+i]=calc_y(r0,g0,b0,cb,cr);
        if(i+1<w)ydata[j*ystride+i+1]=calc_y(r1,g1,b1,cb,cr);
        if(j+1<h){
          ydata[(j+1)*ystride+i]=calc_y(r2,g2,b2,cb,cr);
          if(i+1<w)ydata[(j+1)*ystride+i+1]=calc_y(r3,g3,b3,cb,cr);
        }
      }
      else if(pixel_format==PIXEL_FMT_422){
        rsum=r0+r1;
        gsum=g0+g1;
        bsum=b0+b1;
        cb=OD_CLAMP255(
         OD_DIV_ROUND(-59528*rsum-200256*gsum+259784*bsum,304016865)+128);
        cr=OD_CLAMP255(
         OD_DIV_ROUND(220472*rsum-200256*gsum-20216*bsum,258011295)+128);
        cbdata[j*cbstride+(i>>1)]=(unsigned char)cb;
        crdata[j*crstride+(i>>1)]=(unsigned char)cr;
        ydata[j*ystride+i]=calc_y(r0,g0,b0,cb,cr);
        if(i+1<w)ydata[j*ystride+i+1]=calc_y(r1,g1,b1,cb,cr);
        if(j+1<h){
          rsum=r2+r3;
          gsum=g2+g3;
          bsum=b2+b3;
          cb=OD_CLAMP255(
           OD_DIV_ROUND(-59528*rsum-200256*gsum+259784*bsum,304016865)+128);
          cr=OD_CLAMP255(
           OD_DIV_ROUND(220472*rsum-200256*gsum-20216*bsum,258011295)+128);
          cbdata[(j+1)*cbstride+(i>>1)]=(unsigned char)cb;
          crdata[(j+1)*crstride+(i>>1)]=(unsigned char)cr;
          ydata[(j+1)*ystride+i]=calc_y(r2,g2,b2,cb,cr);
          if(i+1<w)ydata[(j+1)*ystride+i+1]=calc_y(r3,g3,b3,cb,cr);
        }
      }
      else{
        rsum=r0;
        gsum=g0;
        bsum=b0;
        cb=OD_CLAMP255(
         OD_DIV_ROUND(-119056*rsum-400512*gsum+519568*bsum,304016865)+128);
        cr=OD_CLAMP255(
         OD_DIV_ROUND(440944*rsum-400512*gsum-40432*bsum,258011295)+128);
        cbdata[j*cbstride+i]=(unsigned char)cb;
        crdata[j*crstride+i]=(unsigned char)cr;
        ydata[j*ystride+i]=calc_y(r0,g0,b0,cb,cr);
        if(i+1<w){
          rsum=r1;
          gsum=g1;
          bsum=b1;
          cb=OD_CLAMP255(
           OD_DIV_ROUND(-119056*rsum-400512*gsum+519568*bsum,304016865)+128);
          cr=OD_CLAMP255(
           OD_DIV_ROUND(440944*rsum-400512*gsum-40432*bsum,258011295)+128);
          cbdata[j*cbstride+i+1]=(unsigned char)cb;
          crdata[j*crstride+i+1]=(unsigned char)cr;
          ydata[j*ystride+i+1]=calc_y(r1,g1,b1,cb,cr);
        }
        if(j+1<h){
          rsum=r2;
          gsum=g2;
          bsum=b2;
          cb=OD_CLAMP255(
           OD_DIV_ROUND(-119056*rsum-400512*gsum+519568*bsum,304016865)+128);
          cr=OD_CLAMP255(
           OD_DIV_ROUND(440944*rsum-400512*gsum-40432*bsum,258011295)+128);
          cbdata[(j+1)*cbstride+i]=(unsigned char)cb;
          crdata[(j+1)*crstride+i]=(unsigned char)cr;
          ydata[(j+1)*ystride+i]=calc_y(r2,g2,b2,cb,cr);
          if(i+1<w){
            rsum=r3;
            gsum=g3;
            bsum=b3;
            cb=OD_CLAMP255(
             OD_DIV_ROUND(-119056*rsum-400512*gsum+519568*bsum,304016865)+128);
            cr=OD_CLAMP255(
             OD_DIV_ROUND(440944*rsum-400512*gsum-40432*bsum,258011295)+128);
            cbdata[(j+1)*cbstride+i+1]=(unsigned char)cb;
            crdata[(j+1)*crstride+i+1]=(unsigned char)cr;
            ydata[(j+1)*ystride+i+1]=calc_y(r3,g3,b3,cb,cr);
          }
        }
      }
    }
  }
}

static void png_read(png_structp _png,png_bytep _data,png_size_t _sz){
  size_t ret;
  ret=fread(_data,_sz,1,(FILE *)png_get_io_ptr(_png));
  if(ret!=1)png_error(_png,"Read Error");
}

static int read_png(img_plane _ycbcr[3],FILE *_fin){
  unsigned char  header[8];
  png_structp    png;
  png_infop      info;
  png_infop      end;
  png_bytep      img;
  png_bytep     *rows;
  png_color_16p  bkgd;
  png_uint_32    width;
  png_uint_32    height;
  int            bit_depth;
  int            color_type;
  int            interlace_type;
  int            compression_type;
  int            filter_method;
  int            hshift;
  int            vshift;
  png_uint_32    i;
  png_uint_32    j;
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
  img=NULL;
  if(setjmp(png_jmpbuf(png))){
    png_free(png,rows);
    free(img);
    png_destroy_read_struct(&png,&info,&end);
    return -EINVAL;
  }
  png_set_read_fn(png,_fin,png_read);
  png_set_sig_bytes(png,8);
  png_read_info(png,info);
  png_get_IHDR(png,info,&width,&height,&bit_depth,&color_type,
   &interlace_type,&compression_type,&filter_method);
  if(width>INT_MAX||height>INT_MAX||width*(png_size_t)height>INT_MAX){
    png_destroy_read_struct(&png,&info,&end);
    return -EINVAL;
  }
  png_set_expand(png);
  if(bit_depth<8)png_set_packing(png);
  if(!(color_type&PNG_COLOR_MASK_COLOR))png_set_gray_to_rgb(png);
  if(png_get_bKGD(png,info,&bkgd)){
    png_set_background(png,bkgd,PNG_BACKGROUND_GAMMA_FILE,1,1.0);
  }
  /*Note that color_type 2 and 3 can also have alpha, despite not setting the
     PNG_COLOR_MASK_ALPHA bit.
    We always strip it to prevent libpng from overrunning our buffer.*/
  png_set_strip_alpha(png);
  img=(png_bytep)malloc(height*width*6*sizeof(*img));
  rows=(png_bytep *)png_malloc(png,height*sizeof(*rows));
  for(j=0;j<height;j++)rows[j]=img+j*width*6;
  png_read_image(png,rows);
  png_read_end(png,end);
  /*If the image wasn't 16-bit, expand it so it is.
    We do this by duplicating the high byte, so that, e.g., 0 maps to 0 and
     255 maps to 65535.
    This is also nicely endian-independent.*/
  if(bit_depth<16){
    for(j=0;j<height;j++){
      for(i=3*width;i-->0;)rows[j][2*i]=rows[j][2*i+1]=rows[j][i];
    }
  }
  hshift=!(pixel_format&1);
  vshift=!(pixel_format&2);
  if(_ycbcr[0].data==NULL){
    _ycbcr[0].width=(int)width;
    _ycbcr[0].height=(int)height;
    _ycbcr[0].stride=(int)width;
    _ycbcr[1].width=(int)(width+hshift>>hshift);
    _ycbcr[1].height=(int)(height+vshift>>vshift);
    _ycbcr[1].stride=_ycbcr[1].width;
    _ycbcr[2].width=_ycbcr[1].width;
    _ycbcr[2].stride=_ycbcr[1].stride;
    _ycbcr[2].height=_ycbcr[1].height;
    _ycbcr[0].data=(unsigned char *)malloc(_ycbcr[0].stride*_ycbcr[0].height);
    _ycbcr[1].data=(unsigned char *)malloc(_ycbcr[1].stride*_ycbcr[1].height);
    _ycbcr[2].data=(unsigned char *)malloc(_ycbcr[2].stride*_ycbcr[2].height);
  }
  else if(_ycbcr[0].width!=(int)width||_ycbcr[0].height!=(int)height){
    fprintf(stderr,"Input size %ix%i does not match %ix%i\n",
     (int)width,(int)height,_ycbcr[0].width,_ycbcr[0].height);
    return -EINVAL;
  }
  rgb_to_ycbcr(_ycbcr,rows);
  png_free(png,rows);
  free(img);
  png_destroy_read_struct(&png,&info,&end);
  return 0;
}

static int include_files(const struct dirent *de){
  char name[8192];
  int  number;
  number=-1;
  sscanf(de->d_name,input_filter,&number);
  snprintf(name,8191,input_filter,number);
  return !strcmp(name,de->d_name);
}

static FILE *open_png_file(const char *_input_directory,const char *_name){
  FILE *fin;
  char  input_filename[8192];
  snprintf(input_filename,8191,"%s/%s",_input_directory,_name);
  fin=fopen(input_filename,"rb");
  if(fin==NULL){
    fprintf(stderr,"Error opening \"%s\": %s\n",input_filename,strerror(errno));
  }
  return fin;
}

int main(int _argc,char **_argv){
  img_plane       ycbcr[3];
  FILE           *fout;
  FILE           *fin;
  struct dirent **png_files;
  int             npng_files;
  char           *input_mask;
  char           *input_directory;
  char           *scratch;
  int             i;
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
      case 's':aspect_numerator=atol(optarg);break;
      case 'S':aspect_denominator=atol(optarg);break;
      case 'f':fps_numerator=atol(optarg);break;
      case 'F':fps_denominator=atol(optarg);break;
      case '\5':pixel_format=PIXEL_FMT_444;break;
      case '\6':pixel_format=PIXEL_FMT_422;break;
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
  input_mask=_argv[optind];
  if(input_mask==NULL){
    fprintf(stderr,"No input files specified. Run with -h for help.\n");
    return EXIT_FAILURE;
  }
  ycbcr[0].data=NULL;
  if(strcmp(input_mask,"-")!=0){
    /*dirname() and basename() must operate on scratch strings.*/
    scratch=strdup(input_mask);
    input_directory=strdup(dirname(scratch));
    free(scratch);
    scratch=strdup(input_mask);
    input_filter=strdup(basename(scratch));
    free(scratch);
    npng_files=scandir(input_directory,&png_files,include_files,alphasort);
    if(npng_files<=0){
      fprintf(stderr,"No input files found. Run with -h for help.\n");
      return EXIT_FAILURE;
    }
    fin=open_png_file(input_directory,png_files[0]->d_name);
    if(fin==NULL)return EXIT_FAILURE;
  }
  else{
    /*Don't set npng_files to 1 to avoid trying to free a non-existant dirent
       later.*/
    npng_files=0;
    png_files=NULL;
    input_directory=NULL;
    fin=stdin;
  }
  if(output_filename==NULL){
    fprintf(stderr,"No output file specified. Run with -h for help.\n");
    return EXIT_FAILURE;
  }
  fout=strcmp(output_filename,"-")==0?stdout:fopen(output_filename,"wb");
  if(fout==NULL){
    fprintf(stderr,"Error opening output file \"%s\".\n",output_filename);
    return 1;
  }
  if(read_png(ycbcr,fin)<0)return EXIT_FAILURE;
  if(fin!=stdin)fclose(fin);
  fprintf(stderr,"%d frames, %dx%d\n",
   OD_MAXI(npng_files,1),ycbcr[0].width,ycbcr[0].height);
  /*Write the Y4M header.*/
  fprintf(fout,"YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   ycbcr[0].width,ycbcr[0].height,fps_numerator,fps_denominator,
   aspect_numerator,aspect_denominator,CHROMA_TAGS[pixel_format]);
  i=0;
  do{
    int pli;
    int j;
    if(i>0){
      fin=open_png_file(input_directory,png_files[i]->d_name);
      if(fin==NULL)return EXIT_FAILURE;
      fprintf(stderr,"%s\n",png_files[i]->d_name);
      if(read_png(ycbcr,fin)<0)return EXIT_FAILURE;
      if(fin!=stdin)fclose(fin);
    }
    else if(npng_files>0)fprintf(stderr,"%s\n",png_files[0]->d_name);
    fprintf(fout,"FRAME\n");
    for(pli=0;pli<3;pli++){
      for(j=0;j<ycbcr[pli].height;j++){
        if(fwrite(ycbcr[pli].data+j*ycbcr[pli].stride,
         ycbcr[pli].width,1,fout)<1){
          fprintf(stderr,"Error writing to \"%s\".\n",output_filename);
          return EXIT_FAILURE;
        }
      }
    }
  }
  while(++i<npng_files);
  free(ycbcr[2].data);
  free(ycbcr[1].data);
  free(ycbcr[0].data);
  if(fout!=stdout)fclose(fout);
  else fflush(fout);
  while(npng_files-->0)free(png_files[npng_files]);
  free(png_files);
  free(input_filter);
  free(input_directory);
  return EXIT_SUCCESS;
}
