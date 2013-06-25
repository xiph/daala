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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include "od_defs.h"
#include "od_filter.h"
#include "od_intra.h"
#include "image_tools.h"
#include "../src/dct.h"
#include <stdlib.h>

od_rgba16_pixel COLORS[OD_INTRA_NMODES];

void image_draw_block(od_rgba16_image *_image,int _x,int _y,
 const unsigned char *_block,int _stride){
  od_rgba16_pixel color;
  int             i;
  int             j;
  color[3]=(unsigned short)0xFFFFU;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      color[0]=color[1]=color[2]=_block[_stride*i+j]*0x101;
      od_rgba16_image_draw_point(_image,_x+j,_y+i,color);
    }
  }
}

int image_write_png(od_rgba16_image *_image,const char *_name){
  char  fout_name[8192];
  FILE *fout;
  sprintf(fout_name,"%s.png",_name);
  fout=fopen(fout_name,"wb");
  if(fout==NULL){
    fprintf(stderr,"Could not open '%s' for reading.\n",fout_name);
    return EXIT_FAILURE;
  }
  od_rgba16_image_write_png(_image,fout);
  fclose(fout);
  return EXIT_SUCCESS;
}

void image_files_init(image_files *_this,int _nxblocks,int _nyblocks){
  od_rgba16_image_init(&_this->raw,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->map,_nxblocks,_nyblocks);
  od_rgba16_image_init(&_this->pred,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->res,B_SZ*_nxblocks,B_SZ*_nyblocks);
}

void image_files_clear(image_files *_this){
  od_rgba16_image_clear(&_this->raw);
  od_rgba16_image_clear(&_this->map);
  od_rgba16_image_clear(&_this->pred);
  od_rgba16_image_clear(&_this->res);
}

void image_files_write(image_files *_this,const char *_name,const char *_suf){
  char name[8192];
  sprintf(name,"%s-raw%s",_name,_suf==NULL?"":_suf);
  image_write_png(&_this->raw,name);
  sprintf(name,"%s-map%s",_name,_suf==NULL?"":_suf);
  image_write_png(&_this->map,name);
  sprintf(name,"%s-pred%s",_name,_suf==NULL?"":_suf);
  image_write_png(&_this->pred,name);
  sprintf(name,"%s-res%s",_name,_suf==NULL?"":_suf);
  image_write_png(&_this->res,name);
}

static void od_pre_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int by;
  int bx;
  int j;
  int i;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
      for(i=0;i<B_SZ;i++){
        od_coeff col[B_SZ];
#if APPLY_FILTER
        for(j=0;j<B_SZ;j++){
          col[j]=_in[_in_stride*(y+j)+x+i];
        }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*NE_PRE_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(col,col);
#else
# error "Need a prefilter implementation for this block size."
#endif
        for(j=0;j<B_SZ;j++){
          _out[_out_stride*(y+j)+x+i]=col[j];
        }
#else
        for(j=0;j<B_SZ;j++){
          _out[_out_stride*(y+j)+x+i]=_in[_in_stride*(y+j)+x+i];
        }
#endif
      }
#if APPLY_FILTER
      for(j=0;j<B_SZ;j++){
        od_coeff *row;
        row=&_out[_out_stride*(y+j)+x];
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*NE_PRE_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(row,row);
#else
# error "Need a prefilter implementation for this block size."
#endif
      }
#endif
    }
  }
}

static void od_post_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  int j;
  int i;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
      for(j=0;j<B_SZ;j++){
        od_coeff *row;
        for(i=0;i<B_SZ;i++){
          _out[_out_stride*(y+j)+x+i]=_in[_in_stride*(y+j)+x+i];
        }
#if APPLY_FILTER
        row=&_out[_out_stride*(y+j)+x];
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*NE_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(row,row);
#else
# error "Need a postfilter implementation for this block size."
#endif
#endif
      }
#if APPLY_FILTER
      for(i=0;i<B_SZ;i++){
        od_coeff col[B_SZ];
        for(j=0;j<B_SZ;j++){
          col[j]=_out[_out_stride*(y+j)+x+i];
        }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*NE_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(col,col);
#else
# error "Need a postfilter implementation for this block size."
#endif
        for(j=0;j<B_SZ;j++){
          _out[_out_stride*(j+y)+x+i]=col[j];
        }
      }
#endif
    }
  }
}

static void od_fdct_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
      (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(&_out[_out_stride*y+x],_out_stride,
       &_in[_in_stride*y+x],_in_stride);
#else
# error "Need an fDCT implementation for this block size."
#endif
    }
  }
}

static void od_idct_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
      (*OD_IDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(&_out[_out_stride*y+x],_out_stride,
       &_in[_in_stride*y+x],_in_stride);
#else
# error "Need an iDCT implementation for this block size."
#endif
    }
  }
}

void image_data_init(image_data *_this,const char *_name,int _nxblocks,
 int _nyblocks){
  int w;
  int h;
  _this->name=_name;
  _this->nxblocks=_nxblocks;
  _this->nyblocks=_nyblocks;
  _this->mode=(unsigned char *)malloc(sizeof(*_this->mode)*_nxblocks*_nyblocks);
  _this->weight=(double *)malloc(sizeof(*_this->weight)*_nxblocks*_nyblocks);
  w=B_SZ*(_nxblocks+3);
  h=B_SZ*(_nyblocks+3);
  _this->pre=(od_coeff *)malloc(sizeof(*_this->pre)*w*h);
  _this->pre_stride=w;
  w=B_SZ*(_nxblocks+2);
  h=B_SZ*(_nyblocks+2);
  _this->fdct=(od_coeff *)malloc(sizeof(*_this->fdct)*w*h);
  _this->fdct_stride=w;
  w=B_SZ*(_nxblocks+0);
  h=B_SZ*(_nyblocks+0);
  _this->pred=(double *)malloc(sizeof(*_this->pred)*w*h);
  _this->pred_stride=w;
  w=B_SZ*(_nxblocks+2);
  h=B_SZ*(_nyblocks+2);
  _this->idct=(od_coeff *)malloc(sizeof(*_this->idct)*w*h);
  _this->idct_stride=w;
  w=B_SZ*(_nxblocks+1);
  h=B_SZ*(_nyblocks+1);
  _this->post=(od_coeff *)malloc(sizeof(*_this->post)*w*h);
  _this->post_stride=w;
}

void image_data_clear(image_data *_this){
  free(_this->mode);
  free(_this->weight);
  free(_this->pre);
  free(_this->fdct);
  free(_this->pred);
  free(_this->idct);
  free(_this->post);
}

void image_data_pre_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  int x0;
  int y0;
  int bx;
  int by;
  int x;
  int y;
  int bi;
  int bj;
  int i;
  int j;
  od_coeff buf[B_SZ*B_SZ];
  x0=-(B_SZ>>1);
  y0=-(B_SZ>>1);
  bx=by=1;
  if(_bi==0){
    x0-=B_SZ;
    bx++;
  }
  if(_bj==0){
    y0-=B_SZ;
    by++;
  }
  if(_bi==_this->nxblocks-1){
    bx+=2;
  }
  if(_bj==_this->nyblocks-1){
    by+=2;
  }
  x=x0+_bi*B_SZ+(3*B_SZ>>1);
  y=y0+_bj*B_SZ+(3*B_SZ>>1);
  for(bj=0;bj<by;bj++){
    for(bi=0;bi<bx;bi++){
      for(j=0;j<B_SZ;j++){
        for(i=0;i<B_SZ;i++){
          buf[B_SZ*j+i]=
           (_data[_stride*(y0+B_SZ*bj+j)+x0+B_SZ*bi+i]-128)*INPUT_SCALE;
        }
      }
      od_pre_blocks(&_this->pre[_this->pre_stride*(y+B_SZ*bj)+x+B_SZ*bi],
       _this->pre_stride,buf,B_SZ,1,1);
    }
  }
}

void image_data_fdct_block(image_data *_this,int _bi,int _bj){
  int x0;
  int y0;
  int bx;
  int by;
  int x;
  int y;
  x0=_bi*B_SZ+(3*B_SZ>>1);
  y0=_bj*B_SZ+(3*B_SZ>>1);
  bx=by=1;
  if(_bi==0){
    x0-=B_SZ;
    bx++;
  }
  if(_bj==0){
    y0-=B_SZ;
    by++;
  }
  if(_bi==_this->nxblocks-1){
    bx++;
  }
  if(_bj==_this->nyblocks-1){
    by++;
  }
  x=x0-(B_SZ>>1);
  y=y0-(B_SZ>>1);
  od_fdct_blocks(&_this->fdct[_this->fdct_stride*y+x],_this->fdct_stride,
   &_this->pre[_this->pre_stride*y0+x0],_this->pre_stride,bx,by);
}

void image_data_print_block(image_data *_this,int _bi,int _bj,FILE *_fp){
  int by;
  int bx;
  int j;
  int i;
  fprintf(_fp,"%i",_this->mode[_this->nxblocks*_bj+_bi]);
  for(by=0;by<=1;by++){
    for(bx=0;bx<=2-by;bx++){
      od_coeff *block;
      block=&_this->fdct[_this->fdct_stride*B_SZ*(_bj+by)+B_SZ*(_bi+bx)];
      for(j=0;j<B_SZ;j++){
        for(i=0;i<B_SZ;i++){
          fprintf(_fp," %i",block[_this->fdct_stride*j+i]);
        }
      }
    }
  }
  fprintf(_fp,"\n");
  fflush(_fp);
}

void image_data_load_block(image_data *_this,int _bi,int _bj,
 od_coeff _coeffs[5*B_SZ*B_SZ]){
  od_coeff *fdct;
  int       by;
  int       bx;
  int       y;
  int       x;
  fdct=&_this->fdct[_this->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];
  for(by=0;by<=1;by++){
    for(bx=0;bx<=2-by;bx++){
      for(y=0;y<B_SZ;y++){
        for(x=0;x<B_SZ;x++){
          (*_coeffs)=fdct[_this->fdct_stride*(B_SZ*(by-1)+y)+B_SZ*(bx-1)+x];
          _coeffs++;
        }
      }
    }
  }
}

void image_data_pred_block(image_data *_this,int _bi,int _bj){
  double   *pred;
  int       mode;
  od_coeff  coeffs[5*B_SZ*B_SZ];
  pred=&_this->pred[_this->pred_stride*B_SZ*_bj+B_SZ*_bi];
  mode=_this->mode[_this->nxblocks*_bj+_bi];
  image_data_load_block(_this,_bi,_bj,coeffs);
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*NE_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(pred,_this->pred_stride,coeffs,
   mode);
#else
# error "Need a predictor implementation for this block size."
#endif
}

void image_data_stats_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,intra_stats *_stats){
  int       mode;
  od_coeff *ref;
  double   *pred;
  int       j;
  int       i;
  double    buf[B_SZ*B_SZ];
  mode=_this->mode[_this->nxblocks*_bj+_bi];
  ref=&_this->fdct[_this->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];
  pred=&_this->pred[_this->pred_stride*B_SZ*_bj+B_SZ*_bi];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=ref[_this->fdct_stride*j+i]-pred[_this->pred_stride*j+i];
    }
  }
  intra_stats_update(_stats,_data,_stride,mode,ref,_this->fdct_stride,buf,B_SZ);
}

void image_data_idct_block(image_data *_this,int _bi,int _bj){
  int      x0;
  int      y0;
  int      x;
  int      y;
  int      bx;
  int      by;
  int      j;
  int      i;
  double  *p;
  od_coeff buf[B_SZ*B_SZ];
  x0=B_SZ*_bi;
  y0=B_SZ*_bj;
  x=x0+B_SZ;
  y=y0+B_SZ;
  bx=by=1;
  if(_bi==0){
    x-=B_SZ;
    bx++;
  }
  if(_bj==0){
    y-=B_SZ;
    by++;
  }
  if(_bi==_this->nxblocks-1){
    bx++;
  }
  if(_bj==_this->nyblocks-1){
    by++;
  }
  /* TODO remove redundant computations here */
  if(bx!=1||by!=1){
    od_idct_blocks(&_this->idct[_this->idct_stride*y+x],_this->idct_stride,
     &_this->fdct[_this->fdct_stride*y+x],_this->fdct_stride,bx,by);
  }
  p=&_this->pred[_this->pred_stride*y0+x0];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[j*B_SZ+i]=(od_coeff)floor(p[_this->pred_stride*j+i]+0.5);
    }
  }
  x=x0+B_SZ;
  y=y0+B_SZ;
  od_idct_blocks(&_this->idct[_this->idct_stride*y+x],_this->idct_stride,
   buf,B_SZ,1,1);
}

void image_data_post_block(image_data *_this,int _bi,int _bj){
  int x0;
  int y0;
  int x;
  int y;
  int bx;
  int by;
  x=B_SZ*_bi;
  y=B_SZ*_bj;
  x0=x+(B_SZ>>1);
  y0=y+(B_SZ>>1);
  bx=by=1;
  if(_bi==_this->nxblocks-1){
    bx++;
  }
  if(_bj==_this->nyblocks-1){
    by++;
  }
  od_post_blocks(&_this->post[_this->post_stride*y+x],_this->post_stride,
   &_this->idct[_this->idct_stride*y0+x0],_this->idct_stride,bx,by);
}

void image_data_files_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,image_files *_files){
  int            mode;
  od_coeff      *p;
  int            j;
  int            i;
  od_coeff       v;
  unsigned char  buf[B_SZ*B_SZ];

  mode=_this->mode[_bj*_this->nxblocks+_bi];

  od_rgba16_image_draw_point(&_files->map,_bi,_bj,COLORS[mode]);

  p=&_this->pre[_this->pre_stride*(B_SZ*_bj+(3*B_SZ>>1))+B_SZ*_bi+(3*B_SZ>>1)];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      v=(p[_this->pre_stride*j+i]+INPUT_SCALE*128+INPUT_SCALE/2)/INPUT_SCALE;
      buf[B_SZ*j+i]=OD_CLAMPI(0,v,255);
    }
  }
  image_draw_block(&_files->raw,B_SZ*_bi,B_SZ*_bj,buf,B_SZ);

  p=&_this->post[_this->post_stride*(B_SZ*_bj+(B_SZ>>1))+B_SZ*_bi+(B_SZ>>1)];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      v=(p[_this->post_stride*j+i]+INPUT_SCALE*128+INPUT_SCALE/2)/INPUT_SCALE;
      buf[B_SZ*j+i]=OD_CLAMPI(0,v,255);
    }
  }
  image_draw_block(&_files->pred,B_SZ*_bi,B_SZ*_bj,buf,B_SZ);

  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=OD_CLAMPI(0,_data[_stride*j+i]-buf[B_SZ*j+i]+128,255);
    }
  }
  image_draw_block(&_files->res,B_SZ*_bi,B_SZ*_bj,buf,B_SZ);
}

int image_data_save_map(image_data *_this){
  char  name[8192];
  char *pos;
  int   eos;
  FILE *fout;
  strcpy(name,_this->name);
  pos=strrchr(name,'.');
  if(!pos){
    eos=strlen(name);
  }
  else{
    eos=pos-name;
  }
  sprintf(&name[eos],".map");
  fout=fopen(name,"wb");
  if(fout==NULL){
    fprintf(stderr,"Error opening output file '%s'.\n",name);
    return EXIT_FAILURE;
  }
  if(fwrite(_this->mode,
   _this->nxblocks*(size_t)_this->nyblocks*sizeof(*_this->mode),1,fout)<1){
    fprintf(stderr,"Error writing to output file '%s'.\n",name);
    return EXIT_FAILURE;
  }
  if(fwrite(_this->weight,
   _this->nxblocks*(size_t)_this->nyblocks*sizeof(*_this->weight),1,fout)<1){
    fprintf(stderr,"Error writing to output file '%s'.\n",name);
    return EXIT_FAILURE;
  }
  fclose(fout);
  return EXIT_SUCCESS;
}

int image_data_load_map(image_data *_this){
  char  name[8192];
  char *pos;
  int   eos;
  FILE *fin;
  strcpy(name,_this->name);
  pos=strrchr(name,'.');
  if(!pos){
    eos=strlen(name);
  }
  else{
    eos=pos-name;
  }
  sprintf(&name[eos],".map");
  fin=fopen(name,"rb");
  if(fin==NULL){
    fprintf(stderr,"Error opening input file '%s'.\n",name);
    return EXIT_FAILURE;
  }
  if(fread(_this->mode,
   _this->nxblocks*(size_t)_this->nyblocks*sizeof(*_this->mode),1,fin)<1) {
    fprintf(stderr,"Error reading from input file '%s'.\n",name);
    return EXIT_FAILURE;
  }
  if(fread(_this->weight,
   _this->nxblocks*(size_t)_this->nyblocks*sizeof(*_this->weight),1,fin)<1) {
    fprintf(stderr,"Error reading from input file '%s'.\n",name);
    return EXIT_FAILURE;
  }
  fclose(fin);
  return EXIT_SUCCESS;
}
