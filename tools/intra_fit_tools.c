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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "intra_fit_tools.h"

/*Computes the starting offset and number of blocks which can be intra
   predicted with full context (i.e., all of their neighbors) available.*/
void get_intra_dims(const video_input_info *_info,int _pli,int _padding,
 int *_x0,int *_y0,int *_nxblocks,int *_nyblocks){
  int xshift;
  int yshift;
  xshift=_pli!=0&&!(_info->pixel_fmt&1);
  yshift=_pli!=0&&!(_info->pixel_fmt&2);
  /*An offset of 1 would be fine to provide enough context for VP8-style intra
     prediction, but for frequency-domain prediction, we'll want a full block,
     plus overlap.*/
  *_x0=(_info->pic_x>>xshift)+(_padding>>1);
  *_y0=(_info->pic_y>>yshift)+(_padding>>1);
  /*We take an extra block off the end to give enough context for above-right
     intra prediction.*/
  *_nxblocks=_info->pic_w-(_padding<<xshift)>>B_SZ_LOG+xshift;
  *_nyblocks=_info->pic_h-(_padding<<yshift)>>B_SZ_LOG+yshift;
}

char *get_map_filename(const char *_name,int _pli,int _nxblocks,int _nyblocks){
  char  *ret;
  char   fname[8192];
  size_t fname_len;
  sprintf(fname,"%s-pli%i-%i-%ix%i.intra.map",
   _name,_pli,B_SZ,_nxblocks,_nyblocks);
  fname_len=strlen(fname);
  ret=(char *)malloc(fname_len+1);
  memcpy(ret,fname,fname_len+1);
  return ret;
}

char *get_weights_filename(const char *_name,
 int _pli,int _nxblocks,int _nyblocks){
  char  *ret;
  char   fname[8192];
  size_t fname_len;
  sprintf(fname,"%s-pli%i-%i-%ix%i.intra.weights",
   _name,_pli,B_SZ,_nxblocks,_nyblocks);
  fname_len=strlen(fname);
  ret=(char *)malloc(fname_len+1);
  memcpy(ret,fname,fname_len+1);
  return ret;
}

int apply_to_blocks(void *_ctx,plane_start_func _start,block_func _block,
 plane_finish_func _finish,int _argc,const char **_argv){
  block_func blocks[1];
  blocks[0]=_block;
  return apply_to_blocks2(_ctx,4*B_SZ,_start,blocks,1,_finish,0x7,_argc,_argv);
}

int apply_to_blocks2(void *_ctx,int _padding,plane_start_func _start,
 const block_func *_blocks,int _nfuncs,plane_finish_func _finish,int _plmask,
 int _argc,const char **_argv){
  int ai;
  for(ai=1;ai<_argc;ai++){
    video_input vid;
    video_input_info info;
    video_input_ycbcr ycbcr;
    FILE *fin;
    int pli;
    fin=fopen(_argv[ai],"rb");
    if(fin==NULL){
      fprintf(stderr,"Could not open '%s' for reading.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    if(video_input_open(&vid,fin)<0){
      fprintf(stderr,"Error reading video info from '%s'.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    video_input_get_info(&vid,&info);
    if(video_input_fetch_frame(&vid,ycbcr,NULL)<0){
      fprintf(stderr,"Error reading first frame from '%s'.\n",_argv[ai]);
      return EXIT_FAILURE;
    }
    for(pli=0;pli<3;pli++){
      if(_plmask&1<<pli){
        int x0;
        int y0;
        int nxblocks;
        int nyblocks;
        int ret;
        get_intra_dims(&info,pli,_padding,&x0,&y0,&nxblocks,&nyblocks);
        if(_start!=NULL){
          ret=(*_start)(_ctx,_argv[ai],&info,pli,nxblocks,nyblocks);
          if(ret)return ret;
        }
        if(_blocks!=NULL){
          int f;
          for(f=0;f<_nfuncs;f++){
            if(_blocks[f]!=NULL){
              const unsigned char *data;
              int                  stride;
              int                  bj;
              int                  bi;
              data=ycbcr[pli].data;
              stride=ycbcr[pli].stride;
                for(bj=0;bj<nyblocks;bj++){
                  int y;
                  y=y0+B_SZ*bj;
                  for(bi=0;bi<nxblocks;bi++){
                    int x;
                    x=x0+B_SZ*bi;
                    (*_blocks[f])(_ctx,&data[stride*y+x],stride,bi,bj);
                  }
                }
              }
           }
        }
        if(_finish!=NULL){
          ret=(*_finish)(_ctx);
          if(ret)return ret;
        }
      }
    }
    video_input_close(&vid);
  }
  return EXIT_SUCCESS;
}

void vp8_intra_predict(unsigned char *_dst,int _dst_stride,
 const unsigned char *_src,int _src_stride,int _mode){
  const unsigned char *above;
  unsigned char       *left;
  unsigned char        p[4*B_SZ];
  int                  x;
  int                  y;
  above=_src-_src_stride;
  p[2*B_SZ-1]=*(above-1);
  left=p+2*B_SZ;
  for(y=0;y<B_SZ;y++)left[y]=*(_src+_src_stride*y-1);
  switch(_mode){
    case OD_INTRA_DC:{
      int dc;
      dc=0;
      for(x=0;x<B_SZ;x++)dc+=above[x];
      for(y=0;y<B_SZ;y++)dc+=left[y];
      dc=dc+B_SZ>>B_SZ_LOG+1;
      for(y=0;y<B_SZ;y++){
        for(x=0;x<B_SZ;x++){
          *(_dst+y*_dst_stride+x)=(unsigned char)dc;
        }
      }
    }break;
    case OD_INTRA_TM:{
      for(y=0;y<B_SZ;y++){
        for(x=0;x<B_SZ;x++){
          *(_dst+y*_dst_stride+x)=OD_CLAMP255(above[x]+left[y]-*(above-1));
        }
      }
    }break;
    case OD_INTRA_HU:{
      for(y=B_SZ;y<B_SZ+(B_SZ>>1)+1;y++)left[y]=left[B_SZ-1];
      for(y=0;y<B_SZ;y++){
        for(x=0;x<(B_SZ>>1);x++){
          *(_dst+y*_dst_stride+2*x)=
           (unsigned char)(left[y+x]+left[y+x+1]+1>>1);
          *(_dst+y*_dst_stride+2*x+1)=
           (unsigned char)(left[y+x]+2*left[y+x+1]+left[x+y+2]+2>>2);
        }
      }
    }break;
    case OD_INTRA_HE:{
      unsigned char q[B_SZ];
      p[2*B_SZ-1]=*(above-1);
      left[B_SZ]=left[B_SZ-1];
      for(y=0;y<B_SZ;y++){
        q[y]=(unsigned char)(p[2*B_SZ+y-1]+2*p[2*B_SZ+y]+p[2*B_SZ+y+1]+2>>2);
      }
      for(y=0;y<B_SZ;y++)for(x=0;x<B_SZ;x++)*(_dst+y*_dst_stride+x)=q[y];
    }break;
    case OD_INTRA_HD:{
      for(x=0;x<B_SZ;x++)p[2*B_SZ-x-1]=*(above+x-1);
      for(y=0;y<B_SZ;y++){
        for(x=0;x<(B_SZ>>1);x++){
          if(y<x){
            *(_dst+y*_dst_stride+2*x)=(unsigned char)(
             p[2*B_SZ+y-2*x+1]+2*p[2*B_SZ+y-2*x]+p[2*B_SZ+y-2*x-1]+2>>2);
            *(_dst+y*_dst_stride+2*x+1)=(unsigned char)(
             p[2*B_SZ+y-2*x]+2*p[2*B_SZ+y-2*x-1]+p[2*B_SZ+y-2*x-2]+2>>2);
          }
          else{
            *(_dst+y*_dst_stride+2*x)=
             (unsigned char)(p[2*B_SZ+y-x-1]+p[2*B_SZ+y-x]+1>>1);
            *(_dst+y*_dst_stride+2*x+1)=(unsigned char)(
             p[2*B_SZ+y-x]+2*p[2*B_SZ+y-x-1]+p[2*B_SZ+y-x-2]+2>>2);
          }
        }
      }
    }break;
    case OD_INTRA_RD:{
      for(x=0;x<=B_SZ;x++)p[2*B_SZ-x-1]=*(above+x-1);
      for(y=0;y<B_SZ;y++){
        for(x=0;x<B_SZ;x++){
          *(_dst+y*_dst_stride+x)=(unsigned char)(
           p[2*B_SZ+y-x]+2*p[2*B_SZ+y-x-1]+p[2*B_SZ+y-x-2]+2>>2);
        }
      }
    }break;
    case OD_INTRA_VR:{
      for(x=0;x<=B_SZ;x++)p[2*B_SZ-x-1]=*(above+x-1);
      for(y=0;y<(B_SZ>>1);y++){
        for(x=0;x<B_SZ;x++){
          if(x<y){
            *(_dst+2*y*_dst_stride+x)=(unsigned char)(
             p[2*B_SZ+2*y-x-1]+2*p[2*B_SZ+2*y-x-2]+p[2*B_SZ+2*y-x-3]+2>>2);
            *(_dst+(2*y+1)*_dst_stride+x)=(unsigned char)(
             p[2*B_SZ+2*y-x]+2*p[2*B_SZ+2*y-x-1]+p[2*B_SZ+2*y-x-2]+2>>2);
          }
          else{
            *(_dst+2*y*_dst_stride+x)=
             (unsigned char)(p[2*B_SZ+y-x-1]+p[2*B_SZ+y-x-2]+1>>1);
            *(_dst+(2*y+1)*_dst_stride+x)=(unsigned char)(
             p[2*B_SZ+y-x]+2*p[2*B_SZ+y-x-1]+p[2*B_SZ+y-x-2]+2>>2);
          }
        }
      }
    }break;
    case OD_INTRA_VE:{
      unsigned char q[B_SZ];
      for(x=0;x<B_SZ;x++){
        q[x]=(unsigned char)(*(above+x-1)+2*above[x]+above[x+1]+2>>2);
      }
      for(y=0;y<B_SZ;y++)for(x=0;x<B_SZ;x++)*(_dst+y*_dst_stride+x)=q[x];
    }break;
    case OD_INTRA_VL:{
      for(y=0;y<(B_SZ>>1);y++){
        /*This ignores the pattern mis-match for the last few values of the
           VP8 predictor.*/
        for(x=0;x<B_SZ;x++){
          *(_dst+2*y*_dst_stride+x)=
           (unsigned char)(above[y+x]+above[y+x+1]+1>>1);
          *(_dst+(2*y+1)*_dst_stride+x)=
           (unsigned char)(above[y+x]+2*above[y+x+1]+above[y+x+2]+2>>2);
        }
      }
    }break;
    case OD_INTRA_LD:{
      for(y=0;y<B_SZ;y++){
        for(x=0;x<B_SZ;x++){
          *(_dst+y*_dst_stride+x)=(unsigned char)(above[y+x]+2*above[y+x+1]
           +above[OD_MINI(y+x+2,2*B_SZ-1)]+2>>2);
        }
      }
    }break;
  }
}
