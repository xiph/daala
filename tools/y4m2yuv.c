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
#include "vidinput.h"

/*Helper macro to validate fwrite/fread operations*/
#define CHECK_IO(func,ptr,size,nmemb,stream)      \
  if(func(ptr,size,nmemb,stream)<(size_t)nmemb){  \
    fprintf(stderr,#func"() error\n");            \
    return 1;                                     \
  }

static const char *output_filename=NULL;

const char *OPTSTRING="ho:";
struct option OPTIONS[]={
  {"help",no_argument,NULL,'h'},
  {"output",required_argument,NULL,'o'},
  {NULL,0,NULL,0}
};

static void usage(const char *_argv0){
  fprintf(stderr,
   "Usage: %s [options] -o <filename.yuv> <input.y4m>\n\n"
   "Options: \n\n"
   "  -h --help                       Display this help and exit.\n"
   "  -o --output <filename.yuv>      Output file name (required).\n",
   _argv0);
}

int main(int _argc,char **_argv){
  video_input vid;
  video_input_info info;
  video_input_ycbcr ycbcr;
  FILE *fout;
  FILE *fin;
  const char *input_filename;
  int i;
  int xdec;
  int ydec;
  int pli;
  int xstride;
  int y;
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
  fout=strcmp(output_filename,"-")==0?
    stdout:fopen(output_filename,"w");
  xstride = info.depth > 8 ? 2 : 1;
  for(i=0;video_input_fetch_frame(&vid,ycbcr,NULL)>0;i++){
    for(pli=0;pli<3;pli++) {
      xdec=pli&&!(info.pixel_fmt&1);
      ydec=pli&&!(info.pixel_fmt&2);
      for (y = info.pic_y >> ydec;
       y < (info.pic_y + info.pic_h + ydec) >> ydec; y++) {
        CHECK_IO(fwrite,ycbcr[pli].data + y*ycbcr[pli].stride +
         (info.pic_x >> xdec)*xstride, xstride,
         (info.pic_w + xdec) >> xdec, fout);
      }
    }
  }
  if(fout!=stdout){
    fclose(fout);
  }
  video_input_close(&vid);
  return EXIT_SUCCESS;
}
