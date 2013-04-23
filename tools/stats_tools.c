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

#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "stats_tools.h"
#include "od_filter.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define PRINT_SCALE (0)

void mode_data_init(mode_data *_md){
  int i;
  _md->n=0;
  _md->mean=0;
  _md->var=0;
  for(i=0;i<B_SZ*B_SZ;i++){
    _md->satd_avg[i]=0;
  }
  od_covmat_init(&_md->ref,B_SZ*B_SZ);
  od_covmat_init(&_md->pred,B_SZ*B_SZ);
}

void mode_data_clear(mode_data *_md){
  od_covmat_clear(&_md->ref);
  od_covmat_clear(&_md->pred);
}

void mode_data_reset(mode_data *_md){
  int i;
  _md->n=0;
  _md->mean=0;
  _md->var=0;
  for(i=0;i<B_SZ*B_SZ;i++){
    _md->satd_avg[i]=0;
  }
  od_covmat_reset(&_md->ref);
  od_covmat_reset(&_md->pred);
}

/* update the input mean and variance */
void mode_data_add_input(mode_data *_md,const unsigned char *_data,int _stride){
  int n;
  int i;
  int j;
  n=_md->n*B_SZ*B_SZ;
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      double delta;
      double s;
      n++;
      s=1.0/n;
      delta=_data[_stride*j+i]*INPUT_SCALE-_md->mean;
      _md->mean+=delta*s;
      _md->var+=delta*delta*(n-1)*s;
    }
  }
  _md->n++;
}

void mode_data_add_block(mode_data *_md,const od_coeff *_block,int _stride,
 int _ref){
  int    j;
  int    i;
  double buf[B_SZ*B_SZ];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=_block[_stride*j+i];
    }
  }
  if(_ref){
    od_covmat_add(&_md->ref,buf,1);
  }
  else{
    od_covmat_add(&_md->pred,buf,1);
  }
}

void mode_data_combine(mode_data *_a,const mode_data *_b){
  double s;
  double delta;
  int    i;
  if(_b->n==0){
    return;
  }
  s=((double)_b->n)/(_a->n+_b->n);
  delta=_b->mean-_a->mean;
  _a->mean+=delta*s;
  for(i=0;i<B_SZ*B_SZ;i++){
    _a->satd_avg[i]+=(_b->satd_avg[i]-_a->satd_avg[i])*s;
  }
  s*=_a->n;
  _a->var+=_b->var+delta*s;
  od_covmat_combine(&_a->ref,&_b->ref);
  od_covmat_combine(&_a->pred,&_b->pred);
  _a->n+=_b->n;
}

void mode_data_correct(mode_data *_md){
  _md->var/=_md->n*B_SZ*B_SZ;
  od_covmat_correct(&_md->ref);
  od_covmat_correct(&_md->pred);
}

void mode_data_print(mode_data *_md,const char *_label,double *_scale){
  double     cg_ref;
  double     cg_pred;
  double     pg;
  int        v;
  int        u;
  double     satd_avg;
  double     bits_avg;
  cg_ref=10*log10(_md->var);
  cg_pred=10*log10(_md->var);
  pg=0;
  satd_avg=0;
  bits_avg=0;
  for(v=0;v<B_SZ;v++){
    for(u=0;u<B_SZ;u++){
      int    i;
      int    ii;
      double b;
      i=B_SZ*v+u;
      ii=B_SZ*B_SZ*i+i;
      cg_ref-=10*log10(_md->ref.cov[ii]*_scale[v]*_scale[u])/(B_SZ*B_SZ);
      cg_pred-=10*log10(_md->pred.cov[ii]*_scale[v]*_scale[u])/(B_SZ*B_SZ);
      /* compute prediction gain 10*log10(prod_j(ref_cov[j]/pred_cov[j])) */
      pg+=10*log10(_md->ref.cov[ii]/_md->pred.cov[ii])/(B_SZ*B_SZ);
      satd_avg+=sqrt(_scale[v]*_scale[u])*_md->satd_avg[i];
      b=sqrt(_scale[v]*_scale[u]*_md->pred.cov[ii]/2);
      bits_avg+=1+OD_LOG2(b)+M_LOG2E/b*_md->satd_avg[i];
    }
  }
  printf("%s Blocks %5i SATD %G Bits %G Mean %G Var %G CgRef %G CgPred %G Pg %G\n",
   _label,_md->n,satd_avg,bits_avg,_md->mean,_md->var,cg_ref,cg_pred,pg);
}

void mode_data_params(mode_data *_this,double _b[B_SZ*B_SZ],double *_scale){
  int v;
  int u;
  int i;
  int ii;
  for(v=0;v<B_SZ;v++){
    for(u=0;u<B_SZ;u++){
      i=(v*B_SZ+u);
      ii=B_SZ*B_SZ*i+i;
      _b[i]=sqrt(_scale[v]*_scale[u]*_this->pred.cov[ii]/2);
    }
  }
}

void intra_stats_init(intra_stats *_this){
  int mode;
  mode_data_init(&_this->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_init(&_this->md[mode]);
  }
}

void intra_stats_clear(intra_stats *_this){
  int i;
  mode_data_clear(&_this->fr);
  for(i=0;i<OD_INTRA_NMODES;i++){
    mode_data_clear(&_this->md[i]);
  }
}

void intra_stats_reset(intra_stats *_this){
  int i;
  mode_data_reset(&_this->fr);
  for(i=0;i<OD_INTRA_NMODES;i++){
    mode_data_reset(&_this->md[i]);
  }
}

void intra_stats_update(intra_stats *_this,const unsigned char *_data,
 int _stride,int _mode,const od_coeff *_ref,int _ref_stride,
 const od_coeff *_pred,int _pred_stride){
  mode_data *fr;
  mode_data *md;
  int        j;
  int        i;
  double     buf[B_SZ*B_SZ];

  fr=&_this->fr;
  md=&_this->md[_mode];

  /* update the input mean and variance */
  mode_data_add_input(fr,_data,_stride);
  mode_data_add_input(md,_data,_stride);

  /* update the reference mean and covariance */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=_ref[_ref_stride*j+i];
    }
  }
  od_covmat_add(&fr->ref,buf,1);
  od_covmat_add(&md->ref,buf,1);

  /* update the prediction mean and covariance */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=_pred[_pred_stride*j+i];
    }
  }
  od_covmat_add(&fr->pred,buf,1);
  od_covmat_add(&md->pred,buf,1);

  /* update the average satd */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      double satd;
      satd=abs(buf[B_SZ*j+i]);
      fr->satd_avg[B_SZ*j+i]+=(satd-fr->satd_avg[B_SZ*j+i])/fr->n;
      md->satd_avg[B_SZ*j+i]+=(satd-md->satd_avg[B_SZ*j+i])/md->n;
    }
  }
}

void intra_stats_correct(intra_stats *_this){
  int mode;
  mode_data_correct(&_this->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_correct(&_this->md[mode]);
  }
}

void intra_stats_print(intra_stats *_this,const char *_label,
 double *_scale){
  int mode;
  printf("%s\n",_label);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    char label[16];
    sprintf(label,"Mode %i",mode);
    mode_data_print(&_this->md[mode],label,_scale);
  }
  mode_data_print(&_this->fr,"Pooled",_scale);
}

void intra_stats_combine(intra_stats *_this,const intra_stats *_that){
  int mode;
  mode_data_combine(&_this->fr,&_that->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_combine(&_this->md[mode],&_that->md[mode]);
  }
}

/* compute the scale factors for the DCT and TDLT transforms */
double VP8_SCALE[B_SZ];
double OD_SCALE[B_SZ];

#define SCALE_BITS (14)

void vp8_scale_init(double _vp8_scale[B_SZ]){
  int j;
  int i;
  od_coeff buf[B_SZ];
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      buf[j]=i!=j?0:(1<<SCALE_BITS);
    }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_IDCT_1D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,1,buf);
#else
# error "Need an iDCT implementation for this block size."
#endif
    _vp8_scale[i]=0;
    for(j=0;j<B_SZ;j++){
      double c=((double)buf[j])/(1<<SCALE_BITS);
      _vp8_scale[i]+=c*c;
    }
#if PRINT_SCALE
    printf("%s%- 24.18G",i==0?"":" ",_vp8_scale[i]);
#endif
  }
#if PRINT_SCALE
  printf("\n");
#endif
}

#define APPLY_PREFILTER (1)
#define APPLY_POSTFILTER (1)

void od_scale_init(double _od_scale[B_SZ]){
  int i;
  int j;
  od_coeff buf[2*B_SZ];
  for(i=0;i<B_SZ;i++){
    for(j=0;j<2*B_SZ;j++){
      buf[j]=(B_SZ>>1)+i!=j?0:(1<<SCALE_BITS);
    }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_IDCT_1D[B_SZ_LOG-OD_LOG_BSIZE0])(&buf[B_SZ>>1],1,&buf[B_SZ>>1]);
#else
# error "Need an iDCT implementation for this block size."
#endif
#if APPLY_POSTFILTER
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*NE_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(buf,buf);
    (*NE_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(&buf[B_SZ],&buf[B_SZ]);
#else
# error "Need a postfilter implementation for this block size."
#endif
#endif
    _od_scale[i]=0;
    for(j=0;j<2*B_SZ;j++){
      double c=((double)buf[j])/(1<<SCALE_BITS);
      _od_scale[i]+=c*c;
    }
#if PRINT_SCALE
    printf("%s%- 24.18G",i==0?"":" ",_od_scale[i]);
#endif
  }
#if PRINT_SCALE
  printf("\n");
#endif
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
#if APPLY_PREFILTER
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
#if APPLY_PREFILTER
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
#if APPLY_POSTFILTER
        row=&_out[_out_stride*(y+j)+x];
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*NE_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(row,row);
#else
# error "Need a postfilter implementation for this block size."
#endif
#endif
      }
#if APPLY_POSTFILTER
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

#define SCALE_SATD (1)

/* find the best vp8 mode */
int vp8_select_mode(const unsigned char *_data,int _stride,double *_weight){
  double best_satd;
  double next_best_satd;
  int    mode;
  int    best_mode;
  int    next_best_mode;
  best_mode=0;
  best_satd=UINT_MAX;
  next_best_mode=best_mode;
  next_best_satd=best_satd;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    unsigned char block[B_SZ*B_SZ];
    od_coeff      buf[B_SZ*B_SZ];
    int           j;
    int           i;
    double        satd;

    memset(block,0,B_SZ*B_SZ);
    vp8_intra_predict(block,B_SZ,_data,_stride,mode);
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        buf[B_SZ*j+i]=block[B_SZ*j+i]-_data[_stride*j+i];
      }
    }

#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,B_SZ,buf,B_SZ);
#else
# error "Need an fDCT implementation for this block size."
#endif

    satd=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
#if SCALE_SATD
        satd+=sqrt(VP8_SCALE[j]*VP8_SCALE[i])*abs(buf[B_SZ*j+i]);
#else
        satd+=abs(buf[B_SZ*j+i]);
#endif
      }
    }
    if(satd<best_satd){
      next_best_satd=best_satd;
      next_best_mode=best_mode;
      best_satd=satd;
      best_mode=mode;
    }
    else{
      if(satd<next_best_satd){
        next_best_satd=satd;
        next_best_mode=mode;
      }
    }
  }

  if(_weight!=NULL){
    *_weight=best_mode!=0?next_best_satd-best_satd:1;
  }
  return best_mode;
}

int od_select_mode_bits(const od_coeff *_block,int _stride,double *_weight,
 double _b[OD_INTRA_NMODES][B_SZ*B_SZ]){
  int    best_mode;
  double best_bits;
  int    next_best_mode;
  double next_best_bits;
  int    mode;
  best_mode=0;
  best_bits=UINT_MAX;
  next_best_mode=best_mode;
  next_best_bits=best_bits;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double p[B_SZ*B_SZ];
    double bits;
    int    j;
    int    i;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
#if 0
    (*OD_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,_block,_stride,mode);
#else
    (*NE_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,_block,_stride,mode);
#endif
#else
# error "Need a predictor implementation for this block size."
#endif
    bits=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        bits+=1+OD_LOG2(_b[mode][j*B_SZ+i])+M_LOG2E/_b[mode][j*B_SZ+i]*
         abs(_block[_stride*j+i]-(od_coeff)floor(p[B_SZ*j+i]+0.5));
      }
    }
    if(bits<best_bits){
      next_best_mode=best_mode;
      next_best_bits=best_bits;
      best_mode=mode;
      best_bits=bits;
    }
    else{
      if(bits<next_best_bits){
        next_best_mode=mode;
        next_best_bits=bits;
      }
    }
  }
  if(_weight!=NULL){
    *_weight=best_mode!=0?next_best_bits-best_bits:1;
  }
  return best_mode;
}

int od_select_mode_satd(const od_coeff *_block,int _stride,double *_weight){
  int    best_mode;
  double best_satd;
  int    next_best_mode;
  double next_best_satd;
  int    mode;
  best_mode=0;
  best_satd=UINT_MAX;
  next_best_mode=best_mode;
  next_best_satd=best_satd;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double p[B_SZ*B_SZ];
    double satd;
    int    j;
    int    i;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
#if 0
    (*OD_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,_block,_stride,mode);
#else
    (*NE_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,_block,_stride,mode);
#endif
#else
# error "Need a predictor implementation for this block size."
#endif
    satd=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
#if SCALE_SATD
        satd+=sqrt(OD_SCALE[j]*OD_SCALE[i])*
         abs(_block[_stride*j+i]-(od_coeff)floor(p[B_SZ*j+i]+0.5));
#else
        satd+=abs(_block[_stride*j+i]-(od_coeff)floor(p[B_SZ*j+i]+0.5));
#endif
      }
    }
    if(satd<best_satd){
      next_best_mode=best_mode;
      next_best_satd=best_satd;
      best_mode=mode;
      best_satd=satd;
    }
    else{
      if(satd<next_best_satd){
        next_best_mode=mode;
        next_best_satd=satd;
      }
    }
  }
  if(_weight!=NULL){
    *_weight=best_mode!=0?next_best_satd-best_satd:1;
  }
  return best_mode;
}

od_rgba16_pixel COLORS[OD_INTRA_NMODES];

void image_draw_block(od_rgba16_image *_image,int _x,int _y,
 const unsigned char *_block,int _stride){
  od_rgba16_pixel color;
  int             i;
  int             j;
  color[3]=(unsigned short)0xFFFFU;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      color[0]=color[1]=color[2]=_block[_stride*i+j]<<8;
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
  _this->pred=(od_coeff *)malloc(sizeof(*_this->pred)*w*h);
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

void image_data_pred_block(image_data *_this,int _bi,int _bj){
  int       mode;
  od_coeff *block;
  double    p[B_SZ*B_SZ];
  od_coeff *pred;
  int       j;
  int       i;
  mode=_this->mode[_this->nxblocks*_bj+_bi];
  block=&_this->fdct[_this->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*NE_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,block,_this->fdct_stride,mode);
#else
# error "Need a predictor implementation for this block size."
#endif
  pred=&_this->pred[_this->pred_stride*B_SZ*_bj+B_SZ*_bi];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      pred[_this->pred_stride*j+i]=(od_coeff)floor(p[B_SZ*j+i]+0.5);
    }
  }
}

void image_data_stats_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,intra_stats *_stats){
  int       mode;
  od_coeff *ref;
  od_coeff *pred;
  int       j;
  int       i;
  od_coeff  buf[B_SZ*B_SZ];
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
  int x0;
  int y0;
  int x;
  int y;
  int bx;
  int by;
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
  x=x0+B_SZ;
  y=y0+B_SZ;
  od_idct_blocks(&_this->idct[_this->idct_stride*y+x],_this->idct_stride,
   &_this->pred[_this->pred_stride*y0+x0],_this->pred_stride,1,1);
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
  char  name[8196];
  int   eos;
  FILE *fout;
  strcpy(name,_this->name);
  eos=strlen(name)-4;
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
  char  name[8196];
  int   eos;
  FILE *fin;
  strcpy(name,_this->name);
  eos=strlen(name)-4;
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

int ne_apply_to_blocks(void *_ctx,int _ctx_sz,int _plmask,int _padding,
 plane_start_func _start,int _nfuncs,const block_func *_funcs,
 plane_finish_func _finish,int _argc,const char *_argv[]){
  int ai;
#pragma omp parallel for schedule(dynamic)
  for(ai=1;ai<_argc;ai++){
    FILE            *fin;
    video_input      vid;
    th_info          ti;
    th_ycbcr_buffer  ycbcr;
    int              pli;
    int              tid;
    unsigned char   *ctx;
    fin=fopen(_argv[ai],"rb");
    if(fin==NULL){
      fprintf(stderr,"Could not open '%s' for reading.\n",_argv[ai]);
      continue;
    }
    if(video_input_open(&vid,fin)<0){
      fprintf(stderr,"Error reading video info from '%s'.\n",_argv[ai]);
      continue;
    }
    video_input_get_info(&vid,&ti);
    if(video_input_fetch_frame(&vid,ycbcr,NULL)<0){
      fprintf(stderr,"Error reading first frame from '%s'.\n",_argv[ai]);
      continue;
    }
    tid=omp_get_thread_num();
    ctx=((unsigned char *)_ctx)+tid*_ctx_sz;
    for(pli=0;pli<3;pli++){
      if(_plmask&1<<pli){
        int x0;
        int y0;
        int nxblocks;
        int nyblocks;
        get_intra_dims(&ti,pli,_padding,&x0,&y0,&nxblocks,&nyblocks);
        if(_start!=NULL){
          (*_start)(ctx,_argv[ai],&ti,pli,nxblocks,nyblocks);
        }
        if(_funcs!=NULL){
          int f;
          for(f=0;f<_nfuncs;f++){
            if(_funcs[f]!=NULL){
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
                  (*_funcs[f])(ctx,&data[stride*y+x],stride,bi,bj);
                }
              }
            }
          }
        }
        if(_finish!=NULL){
          (*_finish)(ctx);
        }
      }
    }
    video_input_close(&vid);
  }
  return EXIT_SUCCESS;
}

static void ne_intra_pred4x4_mult(double *_p,const od_coeff *_c,int _stride,
 int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
      _p[4*j+i]=NE_PRED_OFFSETS_4x4[_mode][j][i];
      for(k=0;k<NE_PRED_MULTS_4x4[_mode][j][i];k++){
        x=NE_PRED_PARAMX_4x4[_mode][j][i][k];
        y=NE_PRED_PARAMY_4x4[_mode][j][i][k];
        _p[4*j+i]+=_c[_stride*(y-4)+(x-4)]*NE_PRED_WEIGHTS_4x4[_mode][j][i][k];
      }
    }
  }
}

static void ne_intra_pred8x8_mult(double *_p,const od_coeff *_c,int _stride,
 int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<8;j++){
    for(i=0;i<8;i++){
      _p[8*j+i]=NE_PRED_OFFSETS_8x8[_mode][j][i];
      for(k=0;k<NE_PRED_MULTS_8x8[_mode][j][i];k++){
        x=NE_PRED_PARAMX_8x8[_mode][j][i][k];
        y=NE_PRED_PARAMY_8x8[_mode][j][i][k];
        _p[8*j+i]+=_c[_stride*(y-8)+(x-8)]*NE_PRED_WEIGHTS_8x8[_mode][j][i][k];
      }
    }
  }
}

static void ne_intra_pred16x16_mult(double *_p,const od_coeff *_c,int _stride,
 int _mode){
  int j;
  int i;
  int k;
  int x;
  int y;
  for(j=0;j<16;j++){
    for(i=0;i<16;i++){
      _p[16*j+i]=NE_PRED_OFFSETS_16x16[_mode][j][i];
      for(k=0;k<NE_PRED_MULTS_16x16[_mode][j][i];k++){
        x=NE_PRED_PARAMX_16x16[_mode][j][i][k];
        y=NE_PRED_PARAMY_16x16[_mode][j][i][k];
        _p[16*j+i]+=
         _c[_stride*(y-16)+(x-16)]*NE_PRED_WEIGHTS_16x16[_mode][j][i][k];
      }
    }
  }
}

const od_intra_mult_func NE_INTRA_MULT[OD_NBSIZES]={
  ne_intra_pred4x4_mult,
  ne_intra_pred8x8_mult,
  ne_intra_pred16x16_mult
};

void print_betas(FILE *_fp){
  int m;
  int j;
  int i;
  int k;
#if B_SZ==4
  fprintf(_fp,"double NE_PRED_OFFSETS_4x4[OD_INTRA_NMODES][4][4]={\n");
#elif B_SZ==8
  fprintf(_fp,"double NE_PRED_OFFSETS_8x8[OD_INTRA_NMODES][8][8]={\n");
#elif B_SZ==16
  fprintf(_fp,"double NE_PRED_OFFSETS_16x16[OD_INTRA_NMODES][16][16]={\n");
#else
# error "Unsupported block size."
#endif
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {");
      for(i=0;i<B_SZ;i++){
#if B_SZ==4
        fprintf(_fp,"%s%- 24.18G",i>0?",":" ",NE_PRED_OFFSETS_4x4[m][j][i]);
#elif B_SZ==8
        fprintf(_fp,"%s%- 24.18G",i>0?",":"",NE_PRED_OFFSETS_8x8[m][j][i]);
#elif B_SZ==16
        fprintf(_fp,"%s%- 24.18G",i>0?",":"",NE_PRED_OFFSETS_16x16[m][j][i]);
#else
# error "Unsupported block size."
#endif
      }
      fprintf(_fp," }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
#if B_SZ==4
  fprintf(_fp,"int NE_PRED_MULTS_4x4[OD_INTRA_NMODES][4][4]={\n");
#elif B_SZ==8
  fprintf(_fp,"int NE_PRED_MULTS_8x8[OD_INTRA_NMODES][8][8]={\n");
#elif B_SZ==16
  fprintf(_fp,"int NE_PRED_MULTS_16x16[OD_INTRA_NMODES][16][16]={\n");
#else
# error "Unsupported block size."
#endif
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {");
      for(i=0;i<B_SZ;i++){
#if B_SZ==4
        fprintf(_fp,"%s%4i",i>0?",":"",NE_PRED_MULTS_4x4[m][j][i]);
#elif B_SZ==8
        fprintf(_fp,"%s%4i",i>0?",":"",NE_PRED_MULTS_8x8[m][j][i]);
#elif B_SZ==16
        fprintf(_fp,"%s%4i",i>0?",":"",NE_PRED_MULTS_16x16[m][j][i]);
#else
# error "Unsupported block size."
#endif
      }
      fprintf(_fp,"   }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
#if B_SZ==4
  fprintf(_fp,"double NE_PRED_WEIGHTS_4x4[OD_INTRA_NMODES][4][4][4*4*4]={\n");
#elif B_SZ==8
  fprintf(_fp,"double NE_PRED_WEIGHTS_8x8[OD_INTRA_NMODES][8][8][4*8*8]={\n");
#elif B_SZ==16
  fprintf(_fp,
   "double NE_PRED_WEIGHTS_16x16[OD_INTRA_NMODES][16][16][4*16*16]={\n");
#else
# error "Unsupported block size."
#endif
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"      {");
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          fprintf(_fp,"%s%s%- 24.18G",k>0?",":" ",k>0&&k%8==0?"\n        ":"",
           NE_PRED_WEIGHTS_4x4[m][j][i][k]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          fprintf(_fp,"%s%s%- 24.18G",k>0?",":" ",k>0&&k%8==0?"\n        ":"",
           NE_PRED_WEIGHTS_8x8[m][j][i][k]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          fprintf(_fp,"%s%s%- 24.18G",k>0?",":" ",k>0&&k%8==0?"\n        ":"",
           NE_PRED_WEIGHTS_16x16[m][j][i][k]);
        }
#else
# error "Unsupported block size."
#endif
        fprintf(_fp,"%s}%s\n",k==0?"  0  ":"",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
#if B_SZ==4
  fprintf(_fp,"int NE_PRED_PARAMX_4x4[OD_INTRA_NMODES][4][4][4*4*4]={\n");
#elif B_SZ==8
  fprintf(_fp,"int NE_PRED_PARAMX_8x8[OD_INTRA_NMODES][8][8][4*8*8]={\n");
#elif B_SZ==16
  fprintf(_fp,"int NE_PRED_PARAMX_16x16[OD_INTRA_NMODES][16][16][4*16*16]={\n");
#else
# error "Unsupported block size."
#endif
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"      {");
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMX_4x4[m][j][i][k]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMX_8x8[m][j][i][k]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMX_16x16[m][j][i][k]);
        }
#else
# error "Unsupported block size."
#endif
        fprintf(_fp,"%s   }%s\n",k==0?"  0":"",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
#if B_SZ==4
  fprintf(_fp,"int NE_PRED_PARAMY_4x4[OD_INTRA_NMODES][4][4][4*4*4]={\n");
#elif B_SZ==8
  fprintf(_fp,"int NE_PRED_PARAMY_8x8[OD_INTRA_NMODES][8][8][4*8*8]={\n");
#elif B_SZ==16
  fprintf(_fp,"int NE_PRED_PARAMY_16x16[OD_INTRA_NMODES][16][16][4*16*16]={\n");
#else
# error "Unsupported block size."
#endif
  for(m=0;m<OD_INTRA_NMODES;m++){
    fprintf(_fp,"/* Mode %i */\n  {\n",m);
    for(j=0;j<B_SZ;j++){
      fprintf(_fp,"    {\n");
      for(i=0;i<B_SZ;i++){
        fprintf(_fp,"      {");
#if B_SZ==4
        for(k=0;k<NE_PRED_MULTS_4x4[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMY_4x4[m][j][i][k]);
        }
#elif B_SZ==8
        for(k=0;k<NE_PRED_MULTS_8x8[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMY_8x8[m][j][i][k]);
        }
#elif B_SZ==16
        for(k=0;k<NE_PRED_MULTS_16x16[m][j][i];k++){
          fprintf(_fp,"%s%3i",k>0?",":"",NE_PRED_PARAMY_16x16[m][j][i][k]);
        }
#else
# error "Unsupported block size."
#endif
        fprintf(_fp,"%s   }%s\n",k==0?"  0":"",i<B_SZ-1?",":"");
      }
      fprintf(_fp,"    }%s\n",j<B_SZ-1?",":"");
    }
    fprintf(_fp,"  }%s\n",m<OD_INTRA_NMODES-1?",":"");
  }
  fprintf(_fp,"};\n");
  fflush(_fp);
}
