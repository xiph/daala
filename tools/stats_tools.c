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

#include <stdlib.h>
#include <string.h>
#include "stats_tools.h"
#include "od_defs.h"
#include "od_filter.h"
#include "od_intra.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define PRINT_SCALE (0)

void mode_data_init(mode_data *_md,int _b_sz){
  int i;
  _md->n=0;
  _md->mean=0;
  _md->var=0;
  for(i=0;i<B_SZ_MAX*B_SZ_MAX;i++){
    _md->satd_avg[i]=0;
  }
  od_covmat_init(&_md->ref,_b_sz*_b_sz);
  od_covmat_init(&_md->res,_b_sz*_b_sz);
}

void mode_data_clear(mode_data *_md){
  od_covmat_clear(&_md->ref);
  od_covmat_clear(&_md->res);
}

void mode_data_reset(mode_data *_md){
  int i;
  _md->n=0;
  _md->mean=0;
  _md->var=0;
  for(i=0;i<B_SZ_MAX*B_SZ_MAX;i++){
    _md->satd_avg[i]=0;
  }
  od_covmat_reset(&_md->ref);
  od_covmat_reset(&_md->res);
}

/* update the input mean and variance */
void mode_data_add_input(mode_data *_md,const unsigned char *_data,int _stride,
 int _b_sz){
  int n;
  int i;
  int j;
  n=_md->n*_b_sz*_b_sz;
  for(j=0;j<_b_sz;j++){
    for(i=0;i<_b_sz;i++){
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
    od_covmat_add(&_md->res,buf,1);
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
  for(i=0;i<B_SZ_MAX*B_SZ_MAX;i++){
    _a->satd_avg[i]+=(_b->satd_avg[i]-_a->satd_avg[i])*s;
  }
  s*=_a->n;
  _a->var+=_b->var+delta*s;
  od_covmat_combine(&_a->ref,&_b->ref);
  od_covmat_combine(&_a->res,&_b->res);
  _a->n+=_b->n;
}

void mode_data_correct(mode_data *_md,int _b_sz){
  _md->var/=_md->n*_b_sz*_b_sz;
  od_covmat_correct(&_md->ref);
  od_covmat_correct(&_md->res);
}

void mode_data_print(mode_data *_md,const char *_label,double *_scale,
 int _b_sz){
  double     cg_ref;
  double     cg_res;
  int        v;
  int        u;
  double     satd_avg;
  double     bits_avg;
  cg_ref=10*log10(_md->var);
  cg_res=10*log10(_md->var);
  satd_avg=0;
  bits_avg=0;
  for(v=0;v<_b_sz;v++){
    for(u=0;u<_b_sz;u++){
      int    i;
      int    ii;
      double b;
      i=_b_sz*v+u;
      ii=_b_sz*_b_sz*i+i;
      cg_ref-=10*log10(_md->ref.cov[ii]*_scale[v]*_scale[u])/(_b_sz*_b_sz);
      cg_res-=10*log10(_md->res.cov[ii]*_scale[v]*_scale[u])/(_b_sz*_b_sz);
      satd_avg+=sqrt(_scale[v]*_scale[u])*_md->satd_avg[i];
      b=sqrt(_scale[v]*_scale[u]*_md->res.cov[ii]/2);
      bits_avg+=1+OD_LOG2(b)+M_LOG2E/b*_md->satd_avg[i];
    }
  }
  printf("%s Blocks %5i SATD %G Bits %G Mean %G Var %G CgRef %G CgRes %G Pg %G\n",
   _label,_md->n,satd_avg,bits_avg,_md->mean,_md->var,cg_ref,cg_res,cg_res-cg_ref);
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
      _b[i]=sqrt(_scale[v]*_scale[u]*_this->res.cov[ii]/2);
    }
  }
}

void intra_stats_init(intra_stats *_this,int _b_sz_log){
  int mode;
  _this->b_sz_log=_b_sz_log;
  mode_data_init(&_this->fr,1<<_b_sz_log);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_init(&_this->md[mode],1<<_b_sz_log);
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
 const double *_res,int _res_stride){
  int        b_sz;
  mode_data *fr;
  mode_data *md;
  int        j;
  int        i;
  double     buf[B_SZ_MAX*B_SZ_MAX];

  b_sz=1<<_this->b_sz_log;

  fr=&_this->fr;
  md=&_this->md[_mode];

  /* Update the input mean and variance. */
  mode_data_add_input(fr,_data,_stride,b_sz);
  mode_data_add_input(md,_data,_stride,b_sz);

  /* Update the reference mean and covariance. */
  for(j=0;j<b_sz;j++){
    for(i=0;i<b_sz;i++){
      buf[b_sz*j+i]=_ref[_ref_stride*j+i];
    }
  }
  od_covmat_add(&fr->ref,buf,1);
  od_covmat_add(&md->ref,buf,1);

  /* Update the residual mean and covariance. */
  for(j=0;j<b_sz;j++){
    for(i=0;i<b_sz;i++){
      buf[b_sz*j+i]=_res[_res_stride*j+i];
    }
  }
  od_covmat_add(&fr->res,buf,1);
  od_covmat_add(&md->res,buf,1);

  /* Update the average SATD. */
  for(j=0;j<b_sz;j++){
    for(i=0;i<b_sz;i++){
      double satd;
      satd=abs(buf[b_sz*j+i]);
      fr->satd_avg[b_sz*j+i]+=(satd-fr->satd_avg[b_sz*j+i])/fr->n;
      md->satd_avg[b_sz*j+i]+=(satd-md->satd_avg[b_sz*j+i])/md->n;
    }
  }
}

void intra_stats_correct(intra_stats *_this){
  int mode;
  mode_data_correct(&_this->fr,1<<_this->b_sz_log);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_correct(&_this->md[mode],1<<_this->b_sz_log);
  }
}

void intra_stats_print(intra_stats *_this,const char *_label,
 double *_scale){
  int mode;
  printf("%s\n",_label);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    char label[16];
    sprintf(label,"Mode %i",mode);
    mode_data_print(&_this->md[mode],label,_scale,1<<_this->b_sz_log);
  }
  mode_data_print(&_this->fr,"Pooled",_scale,1<<_this->b_sz_log);
}

void intra_stats_combine(intra_stats *_this,const intra_stats *_that){
  int mode;
  mode_data_combine(&_this->fr,&_that->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_combine(&_this->md[mode],&_that->md[mode]);
  }
}

/* compute the scale factors for the DCT and TDLT transforms */
double VP8_SCALE[OD_NBSIZES][B_SZ_MAX];
double OD_SCALE[OD_NBSIZES][B_SZ_MAX];

#define SCALE_BITS (14)

void vp8_scale_init(double *_vp8_scale,int _b_sz_log){
  int b_sz;
  int j;
  int i;
  od_coeff buf[B_SZ_MAX];
  b_sz=1<<_b_sz_log;
  for(i=0;i<b_sz;i++){
    for(j=0;j<b_sz;j++){
      buf[j]=i!=j?0:(1<<SCALE_BITS);
    }
    (*OD_IDCT_1D[_b_sz_log-OD_LOG_BSIZE0])(buf,1,buf);
    _vp8_scale[i]=0;
    for(j=0;j<b_sz;j++){
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

void od_scale_init(double *_od_scale,int _b_sz_log){
  int b_sz;
  int i;
  int j;
  od_coeff buf[2*B_SZ_MAX];
  b_sz=1<<_b_sz_log;
  for(i=0;i<b_sz;i++){
    for(j=0;j<2*b_sz;j++){
      buf[j]=(b_sz>>1)+i!=j?0:(1<<SCALE_BITS);
    }
    (*OD_IDCT_1D[_b_sz_log-OD_LOG_BSIZE0])(&buf[b_sz>>1],1,&buf[b_sz>>1]);
#if APPLY_POSTFILTER
    (*NE_POST_FILTER[_b_sz_log-OD_LOG_BSIZE0])(buf,buf);
    (*NE_POST_FILTER[_b_sz_log-OD_LOG_BSIZE0])(&buf[b_sz],&buf[b_sz]);
#endif
    _od_scale[i]=0;
    for(j=0;j<2*b_sz;j++){
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

#define SCALE_SATD (1)

/* find the best vp8 mode */
int vp8_select_mode(const unsigned char *_data,int _stride,double *_weight){
  int     best_mode;
  double  best_satd;
  double  next_best_satd;
  double *vp8_scale;
  int     mode;
  best_mode=0;
  best_satd=UINT_MAX;
  next_best_satd=best_satd;
  vp8_scale=VP8_SCALE[B_SZ_LOG-OD_LOG_BSIZE0];
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
        satd+=sqrt(vp8_scale[j]*vp8_scale[i])*abs(buf[B_SZ*j+i]);
#else
        satd+=abs(buf[B_SZ*j+i]);
#endif
      }
    }
    if(satd<best_satd){
      next_best_satd=best_satd;
      best_satd=satd;
      best_mode=mode;
    }
    else{
      if(satd<next_best_satd){
        next_best_satd=satd;
      }
    }
  }

  if(_weight!=NULL){
    *_weight=best_mode!=0?next_best_satd-best_satd:1;
  }
  return best_mode;
}

int od_select_mode_bits(const od_coeff *_block,double *_weight,
 double _b[OD_INTRA_NMODES][B_SZ*B_SZ]){
  const od_coeff *c;
  int             best_mode;
  double          best_bits;
  double          next_best_bits;
  double         *od_scale;
  int             mode;
  c=_block+4*B_SZ*B_SZ;
  best_mode=0;
  best_bits=UINT_MAX;
  next_best_bits=best_bits;
  od_scale=OD_SCALE[B_SZ_LOG-OD_LOG_BSIZE0];
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double p[B_SZ*B_SZ];
    double bits;
    int    j;
    int    i;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
#if 0
    (*OD_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,_block,_stride,mode);
#else
    (*NE_INTRA_MULT[B_SZ_LOG-OD_LOG_BSIZE0])(p,B_SZ,_block,mode);
#endif
#else
# error "Need a predictor implementation for this block size."
#endif
    bits=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        double res;
        res=sqrt(od_scale[j]*od_scale[i])*
         abs(c[B_SZ*j+i]-(od_coeff)floor(p[B_SZ*j+i]+0.5));
        bits+=1+OD_LOG2(_b[mode][j*B_SZ+i])+M_LOG2E/_b[mode][j*B_SZ+i]*res;
      }
    }
    if(bits<best_bits){
      next_best_bits=best_bits;
      best_bits=bits;
      best_mode=mode;
    }
    else{
      if(bits<next_best_bits){
        next_best_bits=bits;
      }
    }
  }
  if(_weight!=NULL){
    *_weight=best_mode!=0?next_best_bits-best_bits:1;
  }
  return best_mode;
}

int od_select_mode_satd(const od_coeff *_block,double *_weight,int _b_sz_log){
  int             b_sz;
  const od_coeff *c;
  int             best_mode;
  double          best_satd;
  double          next_best_satd;
  double         *od_scale;
  int             mode;
  b_sz=1<<_b_sz_log;
  c=_block+4*b_sz*b_sz;
  best_mode=0;
  best_satd=UINT_MAX;
  next_best_satd=best_satd;
  od_scale=OD_SCALE[_b_sz_log-OD_LOG_BSIZE0];
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double p[B_SZ_MAX*B_SZ_MAX];
    double satd;
    int    j;
    int    i;
    (*NE_INTRA_MULT[_b_sz_log-OD_LOG_BSIZE0])(p,b_sz,_block,mode);
    satd=0;
    for(j=0;j<b_sz;j++){
      for(i=0;i<b_sz;i++){
#if SCALE_SATD
        satd+=sqrt(od_scale[j]*od_scale[i])*
         abs(c[b_sz*j+i]-(od_coeff)floor(p[b_sz*j+i]+0.5));
#else
        satd+=abs(c[b_sz*j+i]-(od_coeff)floor(p[b_sz*j+i]+0.5));
#endif
      }
    }
    if(satd<best_satd){
      next_best_satd=best_satd;
      best_mode=mode;
      best_satd=satd;
    }
    else{
      if(satd<next_best_satd){
        next_best_satd=satd;
      }
    }
  }
  if(_weight!=NULL){
    *_weight=best_mode!=0?next_best_satd-best_satd:1;
  }
  return best_mode;
}

int ne_apply_to_blocks(void *_ctx,int _ctx_sz,int _plmask,int _padding,
 plane_start_func _start,int _nfuncs,const block_func *_funcs,
 plane_finish_func _finish,int _argc,const char *_argv[]){
  int ai;
#pragma omp parallel for schedule(dynamic)
  for(ai=1;ai<_argc;ai++){
    FILE *fin;
    video_input vid;
    video_input_info info;
    video_input_ycbcr ycbcr;
    int pli;
    int tid;
    unsigned char *ctx;
    fin=fopen(_argv[ai],"rb");
    if(fin==NULL){
      fprintf(stderr,"Could not open '%s' for reading.\n",_argv[ai]);
      continue;
    }
    if(video_input_open(&vid,fin)<0){
      fprintf(stderr,"Error reading video info from '%s'.\n",_argv[ai]);
      continue;
    }
    video_input_get_info(&vid,&info);
    if(video_input_fetch_frame(&vid,ycbcr,NULL)<0){
      fprintf(stderr,"Error reading first frame from '%s'.\n",_argv[ai]);
      continue;
    }
    tid=OD_OMP_GET_THREAD;
    ctx=((unsigned char *)_ctx)+tid*_ctx_sz;
    for(pli=0;pli<3;pli++){
      if(_plmask&1<<pli){
        int x0;
        int y0;
        int nxblocks;
        int nyblocks;
        get_intra_dims(&info,pli,_padding,&x0,&y0,&nxblocks,&nyblocks);
        if(_start!=NULL){
          (*_start)(ctx,_argv[ai],&info,pli,nxblocks,nyblocks);
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
