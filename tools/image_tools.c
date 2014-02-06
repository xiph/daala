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
#include "../src/block_size_enc.h"
#include "../src/dct.h"
#include "../src/tf.h"
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
 int _in_stride,int _bx,int _by,int _b_sz_log){
  int b_sz;
  int by;
  int bx;
  int j;
  int i;
  b_sz=1<<_b_sz_log;
  for(by=0;by<_by;by++){
    int y;
    y=by<<_b_sz_log;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=bx<<_b_sz_log;
      for(i=0;i<b_sz;i++){
        od_coeff col[B_SZ_MAX];
#if APPLY_FILTER
        for(j=0;j<b_sz;j++){
          col[j]=_in[_in_stride*(y+j)+x+i];
        }
        OD_ASSERT(_b_sz_log>=OD_LOG_BSIZE0&&_b_sz_log<=B_SZ_LOG_MAX);
        (*NE_PRE_FILTER[_b_sz_log-OD_LOG_BSIZE0])(col,col);
        for(j=0;j<b_sz;j++){
          _out[_out_stride*(y+j)+x+i]=col[j];
        }
#else
        for(j=0;j<b_sz;j++){
          _out[_out_stride*(y+j)+x+i]=_in[_in_stride*(y+j)+x+i];
        }
#endif
      }
#if APPLY_FILTER
      for(j=0;j<b_sz;j++){
        od_coeff *row;
        row=&_out[_out_stride*(y+j)+x];
        OD_ASSERT(_b_sz_log>=OD_LOG_BSIZE0&&_b_sz_log<=B_SZ_LOG_MAX);
        (*NE_PRE_FILTER[_b_sz_log-OD_LOG_BSIZE0])(row,row);
      }
#endif
    }
  }
}

static void od_post_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by,int _b_sz_log){
  int b_sz;
  int bx;
  int by;
  int j;
  int i;
  b_sz=1<<_b_sz_log;
  for(by=0;by<_by;by++){
    int y;
    y=by<<_b_sz_log;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=bx<<_b_sz_log;
      for(j=0;j<b_sz;j++){
        od_coeff *row;
        for(i=0;i<b_sz;i++){
          _out[_out_stride*(y+j)+x+i]=_in[_in_stride*(y+j)+x+i];
        }
#if APPLY_FILTER
        row=&_out[_out_stride*(y+j)+x];
        OD_ASSERT(_b_sz_log>=OD_LOG_BSIZE0&&_b_sz_log<=B_SZ_LOG_MAX);
        (*NE_POST_FILTER[_b_sz_log-OD_LOG_BSIZE0])(row,row);
#endif
      }
#if APPLY_FILTER
      for(i=0;i<b_sz;i++){
        od_coeff col[B_SZ_MAX];
        for(j=0;j<b_sz;j++){
          col[j]=_out[_out_stride*(y+j)+x+i];
        }
        OD_ASSERT(_b_sz_log>=OD_LOG_BSIZE0&&_b_sz_log<=B_SZ_LOG_MAX);
        (*NE_POST_FILTER[_b_sz_log-OD_LOG_BSIZE0])(col,col);
        for(j=0;j<b_sz;j++){
          _out[_out_stride*(j+y)+x+i]=col[j];
        }
      }
#endif
    }
  }
}

static void od_fdct_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by,int _b_sz_log){
  int bx;
  int by;
  for(by=0;by<_by;by++){
    int y;
    y=by<<_b_sz_log;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=bx<<_b_sz_log;
      OD_ASSERT(_b_sz_log>=OD_LOG_BSIZE0&&_b_sz_log<=B_SZ_LOG_MAX);
      (*OD_FDCT_2D[_b_sz_log-OD_LOG_BSIZE0])(&_out[_out_stride*y+x],_out_stride,
       &_in[_in_stride*y+x],_in_stride);
    }
  }
}

#if TF_BLOCKS
static void od_tf_blocks_down(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by,int _b_sz_log){
  int b_sz;
  int by;
  int bx;
  b_sz=1<<_b_sz_log;
  for(by=0;by<_by;by++){
    int y;
    y=by<<_b_sz_log;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=bx<<_b_sz_log;
      if(_b_sz_log==2){
        int j;
        int i;
        for(j=0;j<b_sz;j++){
          for(i=0;i<b_sz;i++){
            _out[_out_stride*(y+j)+x+i]=_in[_in_stride*(y+j)+x+i];
          }
        }
      }
      else {
        od_convert_block_down(&_out[_out_stride*y+x],_out_stride,
         &_in[_in_stride*y+x],_in_stride,_b_sz_log-2,0,0);
      }
    }
  }
}
#endif

static void od_idct_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by,int _b_sz_log){
  int bx;
  int by;
  for(by=0;by<_by;by++){
    int y;
    y=by<<_b_sz_log;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=bx<<_b_sz_log;
      OD_ASSERT(_b_sz_log>=OD_LOG_BSIZE0&&_b_sz_log<=B_SZ_LOG_MAX);
      (*OD_IDCT_2D[_b_sz_log-OD_LOG_BSIZE0])(&_out[_out_stride*y+x],_out_stride,
       &_in[_in_stride*y+x],_in_stride);
    }
  }
}

void image_data_init(image_data *_this,const char *_name,int _b_sz_log,
 int _nxblocks,int _nyblocks){
  int w;
  int h;
  _this->name=_name;
  _this->b_sz_log=_b_sz_log;
  _this->nxblocks=_nxblocks;
  _this->nyblocks=_nyblocks;
  _this->mask=(unsigned char *)malloc(sizeof(*_this->mask)*_nxblocks*_nyblocks);
  _this->mode=(unsigned char *)malloc(sizeof(*_this->mode)*_nxblocks*_nyblocks);
  _this->weight=(double *)malloc(sizeof(*_this->weight)*_nxblocks*_nyblocks);
  w=(_nxblocks+3)<<_b_sz_log;
  h=(_nyblocks+3)<<_b_sz_log;
  _this->pre=(od_coeff *)malloc(sizeof(*_this->pre)*w*h);
  _this->pre_stride=w;
  w=(_nxblocks+2)<<_b_sz_log;
  h=(_nyblocks+2)<<_b_sz_log;
  _this->fdct=(od_coeff *)malloc(sizeof(*_this->fdct)*w*h);
  _this->fdct_stride=w;
#if TF_BLOCKS
  _this->tf=(od_coeff *)malloc(sizeof(*_this->fdct)*w*h);
#endif
  w=(_nxblocks+0)<<_b_sz_log;
  h=(_nyblocks+0)<<_b_sz_log;
  _this->pred=(double *)malloc(sizeof(*_this->pred)*w*h);
  _this->pred_stride=w;
  w=(_nxblocks+2)<<_b_sz_log;
  h=(_nyblocks+2)<<_b_sz_log;
  _this->idct=(od_coeff *)malloc(sizeof(*_this->idct)*w*h);
  _this->idct_stride=w;
  w=(_nxblocks+1)<<_b_sz_log;
  h=(_nyblocks+1)<<_b_sz_log;
  _this->post=(od_coeff *)malloc(sizeof(*_this->post)*w*h);
  _this->post_stride=w;
}

void image_data_clear(image_data *_this){
  free(_this->mask);
  free(_this->mode);
  free(_this->weight);
  free(_this->pre);
  free(_this->fdct);
#if TF_BLOCKS
  free(_this->tf);
#endif
  free(_this->pred);
  free(_this->idct);
  free(_this->post);
}

void image_data_mask(image_data *_this,const unsigned char *_data,int _stride){
  int j;
  int i;
  memset(_this->mask,0,_this->nyblocks*_this->nxblocks);
  /* process_block_size32 needs 32x32 image data with 6 pixel padding */
  OD_ASSERT(((_stride-_this->nxblocks*B_SZ)/2)>=6);
  for(j=0;j<_this->nyblocks*B_SZ/32;j++){
    for(i=0;i<_this->nxblocks*B_SZ/32;i++){
      const unsigned char *b;
      BlockSizeComp        bs;
      int                  dec[4][4];
      int                  k;
      int                  l;
      b=&_data[_stride*32*j+32*i];
      process_block_size32(&bs,b,_stride,b,_stride,dec);
      for(l=0;l<32/B_SZ;l++){
        for(k=0;k<32/B_SZ;k++){
          /*printf("i=%i j=%i k=%i l=%i\n",i,j,k,l);
          printf("bx=%i by=%i\n",i*32/B_SZ+k*B_SZ,j*32/B_SZ+l*B_SZ);
          fflush(stdout);*/
#if B_SZ==4
          _this->mask[_this->nxblocks*(j*32/B_SZ+l)+i*32/B_SZ+k]=
           (dec[l>>1][k>>1]==0);
#elif B_SZ==8
          _this->mask[_this->nxblocks*(j*32/B_SZ+l)+i*32/B_SZ+k]=
           (dec[l][1]==1);
#elif B_SZ==16
          _this->mask[_this->nxblocks*(j*32/B_SZ+l)+i*32/B_SZ+k]=
           (dec[l<<1][k<<1]==2||dec[l<<1][k<<1]==3);
#else
# error "Invalid block size."
#endif
        }
      }
    }
  }
}

void image_data_pre_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  int b_sz;
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
  od_coeff buf[B_SZ_MAX*B_SZ_MAX];
  b_sz=1<<_this->b_sz_log;
  x0=-(b_sz>>1);
  y0=-(b_sz>>1);
  bx=by=1;
  if(_bi==0){
    x0-=b_sz;
    bx++;
  }
  if(_bj==0){
    y0-=b_sz;
    by++;
  }
  if(_bi==_this->nxblocks-1){
    bx+=2;
  }
  if(_bj==_this->nyblocks-1){
    by+=2;
  }
  x=x0+_bi*b_sz+(3*b_sz>>1);
  y=y0+_bj*b_sz+(3*b_sz>>1);
  for(bj=0;bj<by;bj++){
    for(bi=0;bi<bx;bi++){
      for(j=0;j<b_sz;j++){
        for(i=0;i<b_sz;i++){
          buf[b_sz*j+i]=
           (_data[_stride*(y0+b_sz*bj+j)+x0+b_sz*bi+i]-128)*INPUT_SCALE;
        }
      }
      od_pre_blocks(&_this->pre[_this->pre_stride*(y+b_sz*bj)+x+b_sz*bi],
       _this->pre_stride,buf,b_sz,1,1,_this->b_sz_log);
    }
  }
}

void image_data_fdct_block(image_data *_this,int _bi,int _bj){
  int b_sz;
  int x0;
  int y0;
  int bx;
  int by;
  int x;
  int y;
  b_sz=1<<_this->b_sz_log;
  x0=_bi*b_sz+(3*b_sz>>1);
  y0=_bj*b_sz+(3*b_sz>>1);
  bx=by=1;
  if(_bi==0){
    x0-=b_sz;
    bx++;
  }
  if(_bj==0){
    y0-=b_sz;
    by++;
  }
  if(_bi==_this->nxblocks-1){
    bx++;
  }
  if(_bj==_this->nyblocks-1){
    by++;
  }
  x=x0-(b_sz>>1);
  y=y0-(b_sz>>1);
  od_fdct_blocks(&_this->fdct[_this->fdct_stride*y+x],_this->fdct_stride,
   &_this->pre[_this->pre_stride*y0+x0],_this->pre_stride,bx,by,
    _this->b_sz_log);
}

#if TF_BLOCKS
void image_data_tf_block(image_data *_this,int _bi,int _bj){
  int b_sz;
  int x;
  int y;
  int bx;
  int by;
  b_sz=1<<_this->b_sz_log;
  x=_bi*b_sz+b_sz;
  y=_bj*b_sz+b_sz;
  bx=by=1;
  if(_bi==0){
    x-=b_sz;
    bx++;
  }
  if(_bj==0){
    y-=b_sz;
    by++;
  }
  if(_bi==_this->nxblocks-1){
    bx++;
  }
  if(_bj==_this->nyblocks-1){
    by++;
  }
  od_tf_blocks_down(&_this->tf[_this->fdct_stride*y+x],_this->fdct_stride,
   &_this->fdct[_this->fdct_stride*y+x],_this->fdct_stride,bx,by,
   _this->b_sz_log);
}
#endif

void image_data_print_block(image_data *_this,int _bi,int _bj,FILE *_fp){
  int b_sz;
  int by;
  int bx;
  od_coeff *block;
  int j;
  int i;
  b_sz=1<<_this->b_sz_log;
#if MASK_BLOCKS
  if(!_this->mask[_this->nxblocks*_bj+_bi]){
    return;
  }
#endif
  fprintf(_fp,"%i",_this->mode[_this->nxblocks*_bj+_bi]);
  for(by=0;by<=1;by++){
    for(bx=0;bx<=(1-by)<<1;bx++){
#if TF_BLOCKS
      block=&_this->tf[_this->fdct_stride*b_sz*(_bj+by)+b_sz*(_bi+bx)];
#else
      block=&_this->fdct[_this->fdct_stride*b_sz*(_bj+by)+b_sz*(_bi+bx)];
#endif
      for(j=0;j<b_sz;j++){
        for(i=0;i<b_sz;i++){
          fprintf(_fp," %i",block[_this->fdct_stride*j+i]);
        }
      }
    }
  }
  block=&_this->fdct[_this->fdct_stride*b_sz*(_bj+1)+b_sz*(_bi+1)];
  for(j=0;j<b_sz;j++){
    for(i=0;i<b_sz;i++){
      fprintf(_fp," %i",block[_this->fdct_stride*j+i]);
    }
  }
  fprintf(_fp,"\n");
  fflush(_fp);
}

void image_data_load_block(image_data *_this,int _bi,int _bj,
 od_coeff *_coeffs){
  int       b_sz;
  od_coeff *block;
  int       by;
  int       bx;
  int       y;
  int       x;
  b_sz=1<<_this->b_sz_log;
#if TF_BLOCKS
  block=&_this->tf[_this->fdct_stride*b_sz*_bj+b_sz*_bi];
#else
  block=&_this->fdct[_this->fdct_stride*b_sz*_bj+b_sz*_bi];
#endif
  for(by=0;by<=1;by++){
    for(bx=0;bx<=(1-by)<<1;bx++){
      for(y=0;y<b_sz;y++){
        for(x=0;x<b_sz;x++){
          (*_coeffs)=block[_this->fdct_stride*(b_sz*by+y)+b_sz*bx+x];
          _coeffs++;
        }
      }
    }
  }
  block=&_this->fdct[_this->fdct_stride*b_sz*(_bj+1)+b_sz*(_bi+1)];
  for(y=0;y<b_sz;y++){
    for(x=0;x<b_sz;x++){
      (*_coeffs)=block[_this->fdct_stride*y+x];
      _coeffs++;
    }
  }
}

void image_data_pred_block(image_data *_this,int _bi,int _bj){
  int       b_sz;
  double   *pred;
  int       mode;
  od_coeff  coeffs[5*B_SZ_MAX*B_SZ_MAX];
  b_sz=1<<_this->b_sz_log;
  pred=&_this->pred[_this->pred_stride*b_sz*_bj+b_sz*_bi];
#if MASK_BLOCKS
  if(!_this->mask[_this->nxblocks*_bj+_bi]){
    od_coeff *fdct;
    int       j;
    int       i;
    fdct=&_this->fdct[_this->fdct_stride*b_sz*(_bj+1)+b_sz*(_bi+1)];
    for(j=0;j<b_sz;j++){
      for(i=0;i<b_sz;i++){
        pred[_this->pred_stride*j+i]=fdct[_this->fdct_stride*j+i];
      }
    }
    return;
  }
#endif
  mode=_this->mode[_this->nxblocks*_bj+_bi];
  image_data_load_block(_this,_bi,_bj,coeffs);
  OD_ASSERT(_this->b_sz_log>=OD_LOG_BSIZE0&&_this->b_sz_log<=B_SZ_LOG_MAX);
  (*NE_INTRA_MULT[_this->b_sz_log-OD_LOG_BSIZE0])(pred,_this->pred_stride,
   coeffs,mode);
}

void image_data_stats_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,intra_stats *_stats){
  int       b_sz;
  int       mode;
  od_coeff *ref;
  double   *pred;
  int       j;
  int       i;
  double    buf[B_SZ_MAX*B_SZ_MAX];
  b_sz=1<<_this->b_sz_log;
#if MASK_BLOCKS
  if(!_this->mask[_this->nxblocks*_bj+_bi]){
    return;
  }
#endif
  mode=_this->mode[_this->nxblocks*_bj+_bi];
  ref=&_this->fdct[_this->fdct_stride*b_sz*(_bj+1)+b_sz*(_bi+1)];
  pred=&_this->pred[_this->pred_stride*b_sz*_bj+b_sz*_bi];
  for(j=0;j<b_sz;j++){
    for(i=0;i<b_sz;i++){
      buf[b_sz*j+i]=ref[_this->fdct_stride*j+i]-pred[_this->pred_stride*j+i];
    }
  }
  intra_stats_update(_stats,_data,_stride,mode,ref,_this->fdct_stride,buf,b_sz);
}

void image_data_idct_block(image_data *_this,int _bi,int _bj){
  int      b_sz;
  int      x0;
  int      y0;
  int      x;
  int      y;
  int      bx;
  int      by;
  int      j;
  int      i;
  double  *p;
  od_coeff buf[B_SZ_MAX*B_SZ_MAX];
  b_sz=1<<_this->b_sz_log;
  x0=b_sz*_bi;
  y0=b_sz*_bj;
  x=x0+b_sz;
  y=y0+b_sz;
  bx=by=1;
  if(_bi==0){
    x-=b_sz;
    bx++;
  }
  if(_bj==0){
    y-=b_sz;
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
     &_this->fdct[_this->fdct_stride*y+x],_this->fdct_stride,bx,by,
     _this->b_sz_log);
  }
  p=&_this->pred[_this->pred_stride*y0+x0];
  for(j=0;j<b_sz;j++){
    for(i=0;i<b_sz;i++){
      buf[j*b_sz+i]=(od_coeff)floor(p[_this->pred_stride*j+i]+0.5);
    }
  }
  x=x0+b_sz;
  y=y0+b_sz;
  od_idct_blocks(&_this->idct[_this->idct_stride*y+x],_this->idct_stride,
   buf,b_sz,1,1,_this->b_sz_log);
}

void image_data_post_block(image_data *_this,int _bi,int _bj){
  int b_sz;
  int x0;
  int y0;
  int x;
  int y;
  int bx;
  int by;
  b_sz=1<<_this->b_sz_log;
  x=b_sz*_bi;
  y=b_sz*_bj;
  x0=x+(b_sz>>1);
  y0=y+(b_sz>>1);
  bx=by=1;
  if(_bi==_this->nxblocks-1){
    bx++;
  }
  if(_bj==_this->nyblocks-1){
    by++;
  }
  od_post_blocks(&_this->post[_this->post_stride*y+x],_this->post_stride,
   &_this->idct[_this->idct_stride*y0+x0],_this->idct_stride,bx,by,
   _this->b_sz_log);
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
