/* Daala video codec
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

/* 1D coding gain (dB)

   AR p=.95           4x4       8x8     16x16
   ------------------------------------------
   KLT             7.5825    8.8462    9.4781

   DCT             7.5701    8.8259    9.4555

   CDF(9/7)        8.4687    9.4592    9.7866

   LappedKLT       8.5633    9.4908    9.8951

   LappedDCT       8.5523    9.4871    9.8929



   Subset 1           4x4       8x8     16x16
   ------------------------------------------
   KLT unlord      8.7563
       collapsed   8.7563
       monty       8.7654   10.2628   11.0292

   DCT             8.7469
                   8.7470
                   8.7561   10.2467   11.0115

   CDF(9/7)        9.3429
                   9.3480
                   9.4155   10.6576   11.1965

   LappedKLT       9.6182
                   9.6183
                   9.6295   10.8056   11.3722

   LappedDCT       9.6119
                   9.6120
                   9.6232   10.8028   11.3698
*/


/* 2D coding gain figures */
/* AR p=.95        4x4       8x8     16x16
   KLT         15.1649   17.6924   18.9562
   DCT         15.1403   17.6518   18.9109
   CDF(9/7)    16.9374   18.9184   19.5731
   LappedKLT   17.1266   18.9816   19.7902
   LappedDCT   17.1047   18.9741   19.7858
*/
/* subset 1        4x4       8x8     16x16
   KLT          
   DCT          
   CDF(9/7)     
   LappedKLT    
   LappedDCT    
*/

#include <stdlib.h>
#include "od_defs.h"
#include "od_filter.h"
#include "stats_tools.h"
#include "trans_tools.h"
#include "int_search.h"
#include "kiss99.h"

#define B_SZ_LOG       (2)
#define B_KLT          (0)
#define B_DCT          (1)
#define B_LT           (1)
#define B_WAVELET      (0)

#define USE_2D         (0)
#define USE_FILES      (1)
#define USE_AR95       (0)
#define COMPUTE_NATHAN (0)
#define PRINT_COV      (0)

#define B_SZ           (1<<B_SZ_LOG)

#if B_WAVELET
#if B_SZ_LOG==1
#  define B_SUPPORT (20)
#else
#  if B_SZ_LOG==2
#    define B_SUPPORT (40)
#  else
#    if B_SZ_LOG==3
#      define B_SUPPORT (80)
#    else
#      if B_SZ_LOG==4
#        define B_SUPPORT (160)
#      else
#        error "no support configuration for transform size"
#      endif
#    endif
#  endif
#endif
#else
#if B_LT
#define B_SUPPORT (B_SZ*2)
#else
#define B_SUPPORT (B_SZ)
#endif
#endif

const int    *f;

typedef void (*ne_fdct_func_1d)(double *_out,const double *_in,int _in_stride);
typedef void (*ne_idct_func_1d)(double *_out,int _out_stride,const double *_in);
extern const ne_idct_func_1d OD_IDCT_1D_DOUBLE[OD_NBSIZES];
extern const ne_fdct_func_1d OD_FDCT_1D_DOUBLE[OD_NBSIZES];

#if USE_FILES
static void print_matrix(FILE *_fp,const char *_label,const double *_mat,
 int _stride,int _n,int _m){
  int i;
  int j;
  fprintf(_fp,"%s=\n",_label);
  for(j=0;j<_n;j++){
    for(i=0;i<_m;i++){
      fprintf(_fp,"%s  %- 12.6G",i>0?",":"",_mat[j*_m+i]);
    }
    fprintf(_fp,"\n");
  }
}

typedef struct {
  int sz;
  int *n;
  double *mean_i;
  double *mean_j;
  double *cov;
} cov_state;

static void cov_init(cov_state *_this, int _sz){
  _this->sz     = _sz;
  _this->n      = calloc(_sz,sizeof(*_this->n));
  _this->mean_i = calloc(_sz,sizeof(*_this->mean_i));
  _this->mean_j = calloc(_sz,sizeof(*_this->mean_j));
  _this->cov    = calloc(_sz,sizeof(*_this->cov));
}

static void cov_clear(cov_state *_this){
  if(_this){
    if(_this->n) free(_this->n);
    if(_this->mean_i) free(_this->mean_i);
    if(_this->mean_j) free(_this->mean_j);
    if(_this->cov) free(_this->cov);
  }
}

static void cov_accumulate_1d(cov_state *_this, const unsigned char *_data, int _stride, int _n){
  int i,j;
  for(i=0;i<_this->sz;i++){
    const unsigned char *di=_data;
    const unsigned char *dj=_data+i*_stride;
    for(j=0;j<_n-i;j++){
      double s = 1.0/++_this->n[i];
      double vi = di[j*_stride]-_this->mean_i[i];
      double vj = dj[j*_stride]-_this->mean_j[i];

      _this->mean_i[i] += vi*s;
      _this->mean_j[i] += vj*s;

      _this->cov[i] += vi*vj*(1.0-s);
    }
  }
}

static void cov_accumulate_2d(cov_state *_this, const unsigned char *_data, int _stride, int _w, int _h){
  int x,y,i,j;
  int sz = sqrt(_this->sz);
  for(i=0;i<sz;i++){
    for(j=0;j<sz;j++){
      int ij = i*sz+j;
      for(y=0;y<_h-i;y++){
        const unsigned char *di=_data+y*_stride;
        const unsigned char *dj=_data+(y+i)*_stride+j;
        for(x=0;x<_w-j;x++){
          double s = 1.0/(++_this->n[ij]);
          double vi = di[x]-_this->mean_i[ij];
          double vj = dj[x]-_this->mean_j[ij];

          _this->mean_i[ij] += vi*s;
          _this->mean_j[ij] += vj*s;

          _this->cov[ij] += vi*vj*(1.0-s);
        }
      }
    }
  }
}

static void cov_combine(cov_state *_a,const cov_state *_b){
  int i;
  for(i=0;i<_a->sz;i++){
    double s  = _b->n[i]/(double)(_a->n[i]+_b->n[i]);
    double mi = _b->mean_i[i]-_a->mean_i[i];
    double mj = _b->mean_j[i]-_a->mean_j[i];

    _a->mean_i[i] += mi*s;
    _a->mean_j[i] += mj*s;

    _a->cov[i] += _b->cov[i]+mi*mj*s*_a->n[i];
    _a->n[i]   += _b->n[i];
  }
}

static void cov_normalize(cov_state *_this){
  int i;
  _this->cov[0] /= _this->n[0];
  for(i=1;i<_this->sz;i++)
    _this->cov[i] /= _this->cov[0]*_this->n[i];
  _this->cov[0]=1.;
}

int apply(trans_ctx *_ctx,cov_state *_cov,int _argc,const char *_argv[]){

  int ai;
#pragma omp parallel for schedule(dynamic)
  for(ai=1;ai<_argc;ai++){
    FILE            *fin;
    video_input      vid;
    th_info          ti;
    th_ycbcr_buffer  ycbcr;
    int              tid;
    trans_ctx       *ctx;
    cov_state       *cov;
    int              x0,y0,x1,y1;

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
    tid=OD_OMP_GET_THREAD;
    ctx=_ctx+tid;
    cov=_cov+tid;
    x0 = ti.pic_x;
    y0 = ti.pic_y;
    x1 = x0 + ti.pic_width;
    y1 = y0 + ti.pic_height;

    /* start */
    fprintf(stderr,"%s\n",_argv[ai]);

    /* map */
    {
      int                  stride=ycbcr[0].stride;
      const unsigned char *data=ycbcr[0].data;

#if USE_2D
      unsigned char        buf[B_SUPPORT][B_SUPPORT];
      int                  x,y,i,j;
      for(y=y0;y<y1-B_SUPPORT_MAX+1;y++){
        for(x=x0;x<x1-B_SUPPORT_MAX+1;x++){
          for(j=0;j<B_SUPPORT;j++){
            for(i=0;i<B_SUPPORT;i++){
              buf[j][i]=data[(y+j+B_PAD)*stride+(x+i+B_PAD)];
            }
          }
          trans_data_add(&ctx->td,(unsigned char *)buf);
        }
      }
#else
      unsigned char        buf[B_SUPPORT];
      int                  x,y,z;
      int nxblocks=ti.pic_width>>B_SZ_LOG;
      int nyblocks=ti.pic_height>>B_SZ_LOG;

      /* Direct computation of collapsed covariance matrix (monty style) */
      for(y=y0;y<y1;y++)
        cov_accumulate_1d(cov,data+y*stride+x0,1,x1-x0);
      for(x=x0;x<x1;x++)
        cov_accumulate_1d(cov,data+y0*stride+x,stride,y1-y0);

      /* block-based full covariance computation (unlord style) */
      image_ctx_init(&ctx->img,_argv[ai],nxblocks,nyblocks);
      /* add the rows */
      for(y=0;y<nyblocks*B_SZ;y++){
        for(x=0;x<nxblocks*B_SZ-B_SUPPORT+1;x++){
          for(z=0;z<B_SUPPORT;z++){
            buf[z]=data[(y+y0)*stride+x+x0+z];
          }
          trans_data_add(&ctx->td,buf);
        }
      }
      /* add the columns */
      for(y=0;y<nyblocks*B_SZ-B_SUPPORT+1;y++){
        for(x=0;x<nxblocks*B_SZ;x++){
          for(z=0;z<B_SUPPORT;z++){
            buf[z]=data[(y0+y+z)*stride+x+x0];
          }
          trans_data_add(&ctx->td,buf);
        }
      }
#endif
    }

    video_input_close(&vid);
  }
  return EXIT_SUCCESS;
}
#endif

#if B_WAVELET

/* some lifting CDF (9/7) wavelet code from Google Code's axonlib */
/* http://code.google.com/p/axonlib/source/browse/trunk/extern/dwt97.c?spec=svn19&r=19 */
/* single stage of decomposition */
static void fwt97_i(double* x,int n){
  double temp[B_SUPPORT];
  double a;
  int i;
  /* Predict 1 */
  a=-1.586134342;
  for (i=1;i<n-2;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[n-1]+=2*a*x[n-2];

  /* Update 1 */
  a=-0.05298011854;
  for (i=2;i<n;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[0]+=2*a*x[1];

  /* Predict 2 */
  a=0.8829110762;
  for (i=1;i<n-2;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[n-1]+=2*a*x[n-2];

  /* Update 2 */
  a=0.4435068522;
  for (i=2;i<n;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[0]+=2*a*x[1];

  /* Scale */
  a=1/1.149604398;
  for (i=0;i<n;i++)
  {
    if (i%2) x[i]*=a;
    else x[i]/=a;
  }
  /* Pack */
  for (i=0;i<n;i++){
    if (i%2==0)
      temp[i/2]=x[i];
    else
      temp[n/2+i/2]=x[i];
  }
  for (i=0;i<n;i++) x[i]=temp[i];
}

/* single stage of reconstruction */
void iwt97_i(double* x,int n){
  double temp[B_SUPPORT];
  double a;
  int i;
  /* Unpack */
  for (i=0;i<n/2;i++){
    temp[i*2]=x[i];
    temp[i*2+1]=x[i+n/2];
  }
  for (i=0;i<n;i++) x[i]=temp[i];

  /* Undo scale */
  a=1.149604398;
  for (i=0;i<n;i++)
  {
    if (i%2) x[i]*=a;
    else x[i]/=a;
  }

  /* Undo update 2 */
  a=-0.4435068522;
  for (i=2;i<n;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[0]+=2*a*x[1];

  /* Undo predict 2 */
  a=-0.8829110762;
  for (i=1;i<n-2;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[n-1]+=2*a*x[n-2];

  /* Undo update 1 */
  a=0.05298011854;
  for (i=2;i<n;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[0]+=2*a*x[1];

  /* Undo predict 1 */
  a=1.586134342;
  for (i=1;i<n-2;i+=2)
    x[i]+=a*(x[i-1]+x[i+1]);
  x[n-1]+=2*a*x[n-2];
}

/* multistage decomposition */
void fwt97(double *out, int n, double *in, int support){
  int i=n,j=support,k;
  while((i&1)==0){
    fwt97_i(in,j);
    i>>=1;
    for(k=0;k<i;k++)
      out[i+k] = in[((((j*3)>>1)-i)>>1) + k];
    j>>=1;
  }
  for(k=0;k<i;k++)
    out[k] = in[((j-i)>>1) + k];
}

/* multistage reconstruction */
void iwt97(double *out, int support, double *in, int n){
  int i=n,j=support,k;
  for(k=0;k<support;k++)
    out[k]=0;
  while((i&1)==0){
    i>>=1;
    for(k=0;k<i;k++)
      out[((((j*3)>>1)-i)>>1) + k]=in[i+k];
    j>>=1;
  }
  for(k=0;k<i;k++)
    out[((j-i)>>1) + k]=in[k];

  i<<=1;
  j<<=1;
  while(j<=support){
    iwt97_i(out,j);
    i<<=1;
    j<<=1;
  }
}
#endif

#if B_KLT
void symeigen(double *out,
              double *cov,
              int support){
  int    i;
  int    j;
  int    k;

  for(i=0;i<support;i++)
    for(j=0;j<support;j++)
      out[i*support+j]=i==j;

  for(;;){
    double mod=0.;
    for(i=0,j=0,k=0;k<support;k++){
      int m;
      for(m=k+1;m<support;m++){
        double q;
        q=fabs(cov[k*support+m]);
        if(q>mod){
          mod=q;
          i=k;
          j=m;
        }
      }
    }
    if(mod<1E-11)break;

    {
      double th=0.5*atan2(2*cov[i*support+j],cov[i*support+i]-cov[j*support+j]);
      double c=cos(th);
      double s=sin(th);

      for(k=0;k<support;k++){
        double t;
        t=c*cov[k*support+i]+s*cov[k*support+j];
        cov[k*support+j]=-s*cov[k*support+i]+c*cov[k*support+j];
        cov[k*support+i]=t;
      }

      for(k=0;k<support;k++){
        double t;
        t=c*cov[i*support+k]+s*cov[j*support+k];
        cov[j*support+k]=-s*cov[i*support+k]+c*cov[j*support+k];
        cov[i*support+k]=t;
      }

      for(k=0;k<support;k++){
        double t;
        t=c*out[i*support+k]+s*out[j*support+k];
        out[j*support+k]=-s*out[i*support+k]+c*out[j*support+k];
        out[i*support+k]=t;
      }
    }
  }
  /* for(j=0;j<B_SZ;j++)eigenvalue[j]=cov[j][j]; don't need eigenvalues */

#if PRINT_CG_MATH
  print_matrix(stdout,"klt",out,support,support,support);
#endif
}

void flap_2d(double out[B_SZ][B_SZ],
            double in[B_SUPPORT][B_SUPPORT],
            const int _f[]){
  int i,j;

#if (B_LT)
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES

  /* columns */
  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    double work[B_SZ*2];
    for(j=0;j<B_SZ*2;j++)
      work[j]=in[j+B_SUPPORT/2-B_SZ][i];
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&work[0],&work[0],_f);
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&work[B_SZ],&work[B_SZ],_f);
    for(j=0;j<B_SZ*2;j++)
      in[j+B_SUPPORT/2-B_SZ][i]=work[j];
  }

  /* rows */
  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&in[i][B_SUPPORT/2-B_SZ],&in[i][B_SUPPORT/2-B_SZ],_f);
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&in[i][B_SUPPORT/2],&in[i][B_SUPPORT/2],_f);
  }

#else
#  error "Need a prefilter implementation for this block size."
#endif
#endif

  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      out[i][j]=in[i+B_SUPPORT/2-B_SZ/2][j+B_SUPPORT/2-B_SZ/2];
}

void ilap_2d(double out[B_SUPPORT][B_SUPPORT],
              double in[B_SZ][B_SZ],
              const int _f[]){
  int i,j;

  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      out[i][j]=0;

  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      out[i+B_SUPPORT/2-B_SZ/2][j+B_SUPPORT/2-B_SZ/2]=in[i][j];

#if (B_LT)
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES

  /* columns */
  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    double work[B_SZ*2];
    for(j=0;j<B_SZ*2;j++)
      work[j]=out[j+B_SUPPORT/2-B_SZ][i];
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&work[0],&work[0],_f);
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&work[B_SZ],&work[B_SZ],_f);
    for(j=0;j<B_SZ*2;j++)
      out[j+B_SUPPORT/2-B_SZ][i]=work[j];
  }

  /* rows */
  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&out[i][B_SUPPORT/2-B_SZ],&out[i][B_SUPPORT/2-B_SZ],_f);
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
      (&out[i][B_SUPPORT/2],&out[i][B_SUPPORT/2],_f);
  }

#else
#  error "Need a prefilter implementation for this block size."
#endif
#endif

}

void flap_4d(double out[B_SZ*B_SZ][B_SZ*B_SZ],
            double in[B_SUPPORT][B_SUPPORT][B_SUPPORT][B_SUPPORT],
            const int _f[]){
  int i,j,k,l;

#if (B_LT)
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES

  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    for(j=B_SUPPORT/2-B_SZ;j<B_SUPPORT/2+B_SZ;j++){
      for(k=B_SUPPORT/2-B_SZ;k<B_SUPPORT/2+B_SZ;k++){
        double work[B_SZ*2];

        /* [ ][i][j][k] */
        for(l=0;l<B_SZ*2;l++)
          work[l]=in[l+B_SUPPORT/2-B_SZ][i][j][k];
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&work[0],&work[0],_f);
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&work[B_SZ],&work[B_SZ],_f);
        for(l=0;l<B_SZ*2;l++)
          in[l+B_SUPPORT/2-B_SZ][i][j][k]=work[l];
      }
    }
  }

  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    for(j=B_SUPPORT/2-B_SZ;j<B_SUPPORT/2+B_SZ;j++){
      for(k=B_SUPPORT/2-B_SZ;k<B_SUPPORT/2+B_SZ;k++){
        double work[B_SZ*2];
        /* [i][ ][j][k] */
        for(l=0;l<B_SZ*2;l++)
          work[l]=in[i][l+B_SUPPORT/2-B_SZ][j][k];
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&work[0],&work[0],_f);
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&work[B_SZ],&work[B_SZ],_f);
        for(l=0;l<B_SZ*2;l++)
          in[i][l+B_SUPPORT/2-B_SZ][j][k]=work[l];
      }
    }
  }


  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    for(j=B_SUPPORT/2-B_SZ;j<B_SUPPORT/2+B_SZ;j++){
      for(k=B_SUPPORT/2-B_SZ;k<B_SUPPORT/2+B_SZ;k++){
        double work[B_SZ*2];
        /* [i][j][ ][k] */
        for(l=0;l<B_SZ*2;l++)
          work[l]=in[i][j][l+B_SUPPORT/2-B_SZ][k];
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&work[0],&work[0],_f);
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&work[B_SZ],&work[B_SZ],_f);
        for(l=0;l<B_SZ*2;l++)
          in[i][j][l+B_SUPPORT/2-B_SZ][k]=work[l];
      }
    }
  }

  for(i=B_SUPPORT/2-B_SZ;i<B_SUPPORT/2+B_SZ;i++){
    for(j=B_SUPPORT/2-B_SZ;j<B_SUPPORT/2+B_SZ;j++){
      for(k=B_SUPPORT/2-B_SZ;k<B_SUPPORT/2+B_SZ;k++){
        /* [i][j][k][ ] */
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&in[i][j][k][B_SUPPORT/2-B_SZ],&in[i][j][k][B_SUPPORT/2-B_SZ],_f);
        (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
          (&in[i][j][k][B_SUPPORT/2],&in[i][j][k][B_SUPPORT/2],_f);

      }
    }
  }

#else
#  error "Need a prefilter implementation for this block size."
#endif
#endif

  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      for(k=0;k<B_SZ;k++)
        for(l=0;l<B_SZ;l++)
          out[i*B_SZ+j][k*B_SZ+l]=in
            [i+B_SUPPORT/2-B_SZ/2]
            [j+B_SUPPORT/2-B_SZ/2]
            [k+B_SUPPORT/2-B_SZ/2]
            [l+B_SUPPORT/2-B_SZ/2];
}

void gklt_1d(double klt[B_SZ][B_SZ],
             double cov[B_SUPPORT][B_SUPPORT],
             const int *_f){
  static double workA[B_SUPPORT][B_SUPPORT];
  static double workB[B_SZ][B_SZ];
  int i,j;
  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      workA[i][j]=cov[i][j];
  flap_2d(workB,workA,f);
  symeigen(&klt[0][0],&workB[0][0],B_SZ);
}

void gklt_2d(double klt[B_SZ*B_SZ][B_SZ*B_SZ],
             double cov[B_SUPPORT][B_SUPPORT][B_SUPPORT][B_SUPPORT],
             const int *_f){
  static double workA[B_SUPPORT][B_SUPPORT][B_SUPPORT][B_SUPPORT];
  static double workB[B_SZ*B_SZ][B_SZ*B_SZ];
  int i,j,k,l;
  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      for(k=0;k<B_SUPPORT;k++)
        for(l=0;l<B_SUPPORT;l++)
          workA[i][j][k][l]=cov[i][j][k][l];
  flap_4d(workB,workA,f);
  symeigen(&klt[0][0],&workB[0][0],B_SZ*B_SZ);
}

void gklt_1d_collapsed(double klt[B_SZ][B_SZ],
                       double cov[B_SUPPORT],
                       const int *_f){
  static double workA[B_SUPPORT][B_SUPPORT];
  static double workB[B_SZ][B_SZ];
  int i,j;
  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      workA[i][j]=cov[abs(i-j)];
  flap_2d(workB,workA,f);
  symeigen(&klt[0][0],&workB[0][0],B_SZ);
}

void gklt_2d_collapsed(double klt[B_SZ*B_SZ][B_SZ*B_SZ],
                       double cov[B_SUPPORT][B_SUPPORT],
                       const int *_f){
  static double workA[B_SUPPORT][B_SUPPORT][B_SUPPORT][B_SUPPORT];
  static double workB[B_SZ*B_SZ][B_SZ*B_SZ];
  int i,j,k,l;
  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      for(k=0;k<B_SUPPORT;k++)
        for(l=0;l<B_SUPPORT;l++)
          workA[i][j][k][l]=cov[abs(i-k)][abs(j-l)];
  flap_4d(workB,workA,f);
  symeigen(&klt[0][0],&workB[0][0],B_SZ*B_SZ);
}

void fklt(double *out,
             double *in,
             double *klt,
             int support){
  int i,j;
  for(i=0;i<support;i++){
    double acc=0.;
    for(j=0;j<support;j++)
      acc += klt[i*support+j]*in[j];
    out[i]=acc;
  }
}

void iklt(double *out,
          double *in,
          double *klt,
          int support){
  int i,j;
  for(i=0;i<support;i++){
    double acc=0.;
    for(j=0;j<support;j++)
      acc+=klt[j*support+i]*in[j];
    out[i]=acc;
  }
}

#endif

void b_analysis_1d(double *_out,int _out_stride,const double *_in,int _in_stride,
                   const int *_f, double _klt[B_SZ][B_SZ]){
  int    j;
  double t[B_SUPPORT];
  double w[B_SZ];

  for(j=0;j<B_SUPPORT;j++)
    t[j]=_in[j*_in_stride];

#if (B_WAVELET)

  fwt97(w,B_SZ,t,B_SUPPORT);
  for(j=0;j<B_SZ;j++){
    _out[j*_out_stride]=w[j];
  }

#else
#  if (B_LT)
#    if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
    (&t[B_SUPPORT/2-B_SZ],&t[B_SUPPORT/2-B_SZ],_f);
  (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
    (&t[B_SUPPORT/2],&t[B_SUPPORT/2],_f);
#    else
#      error "Need a prefilter implementation for this block size."
#    endif
#  endif
#  if B_KLT
  fklt(&w[0],&t[B_SUPPORT/2-B_SZ/2],&_klt[0][0],B_SZ);
#  else
#    if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_1D_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
    (w,&t[B_SUPPORT/2-B_SZ/2],1);
#    else
#      error "Need an fDCT implementation for this block size."
#    endif
#  endif

  for(j=0;j<B_SZ;j++)
    _out[j*_out_stride]=w[j];

#endif
}

void b_analysis_2d(double *_out,int _out_stride_i,int _out_stride_j,
                   const double *_in,int _in_stride_i,int _in_stride_j,
                   const int *_f, double _klt[B_SZ*B_SZ][B_SZ*B_SZ]){
#if !(B_KLT)
  double work[B_SUPPORT][B_SZ];
  int i;

  /* DCT and DWT are separable 1D transforms */
  /* lapping performed inside b_analysis */
  for(i=0;i<B_SUPPORT;i++)
    b_analysis_1d(&work[i][0],1,_in+i*_in_stride_i,_in_stride_j,_f,NULL);
  for(i=0;i<B_SZ;i++)
    b_analysis_1d(_out+_out_stride_i*i,_out_stride_j,&work[0][i],B_SZ,_f,NULL);

#else
  /* KLT is a non-separable 2D transform */
  double lap[B_SUPPORT][B_SUPPORT];
  double work[B_SZ][B_SZ];
  double temp[B_SZ][B_SZ];
  int i,j;
  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      lap[i][j]=*(_in+i*_in_stride_i+j*_in_stride_j);

  flap_2d(work,lap,f);

  fklt(&temp[0][0],&work[0][0],&_klt[0][0],B_SZ*B_SZ);

  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      *(_out+i*_out_stride_i+j*_out_stride_j)=temp[i][j];

#endif
}


void b_synthesis_1d(double *_out,int _out_stride,const double *_in,int _in_stride,
                    const int *_f, double _klt[B_SZ][B_SZ]){
  int    j;
  double w[B_SUPPORT];
  double t[B_SUPPORT];

  for(j=0;j<B_SUPPORT;j++){
    t[j]=0;
    w[j]=0;
  }

#if (B_WAVELET)
  for(j=0;j<B_SZ;j++)
    w[j]=_in[j*_in_stride];
  iwt97(t,B_SUPPORT,w,B_SZ);
#else
  for(j=0;j<B_SZ;j++){
    w[B_SUPPORT/2-B_SZ/2+j]=_in[j*_in_stride];
  }
#  if B_KLT
  iklt(&t[B_SUPPORT/2-B_SZ/2],&w[B_SUPPORT/2-B_SZ/2],&_klt[0][0],B_SZ);
#  else
#    if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_IDCT_1D_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
    (&t[B_SUPPORT/2-B_SZ/2],1,&w[B_SUPPORT/2-B_SZ/2]);
#    else
#      error "Need an iDCT implementation for this block size."
#    endif
#  endif

#  if B_LT
#    if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
    (&t[B_SUPPORT/2-B_SZ],&t[B_SUPPORT/2-B_SZ],_f);
  (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])
    (&t[B_SUPPORT/2],&t[B_SUPPORT/2],_f);
#    else
# error "Need a postfilter implementation for this block size."
#    endif
#  endif

#endif
  for(j=0;j<B_SUPPORT;j++)
    _out[j*_out_stride]=t[j];
}

void b_synthesis_2d(double *_out,int _out_stride_i,int _out_stride_j,
                    const double *_in,int _in_stride_i,int _in_stride_j,
                    const int *_f, double _klt[B_SZ*B_SZ][B_SZ*B_SZ]){
#if !(B_KLT)
  double work[B_SUPPORT][B_SZ];
  int i;

  /* DCT and DWT are separable 1D transforms */
  /* lapping performed inside b_analysis */
  for(i=0;i<B_SZ;i++)
    b_synthesis_1d(&work[0][i],B_SZ,_in+i*_in_stride_i,_in_stride_j,_f,NULL);
  for(i=0;i<B_SUPPORT;i++)
    b_synthesis_1d(_out+_out_stride_i*i,_out_stride_j,&work[i][0],1,_f,NULL);

#else

  /* KLT is a non-separable 2D transform */
  double temp[B_SZ][B_SZ];
  double work[B_SZ][B_SZ];
  double lap[B_SUPPORT][B_SUPPORT];
  int i,j;

  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      temp[i][j]=*(_in+i*_in_stride_i+j*_in_stride_j);

  iklt(&work[0][0],&temp[0][0],&_klt[0][0],B_SZ*B_SZ);

  ilap_2d(lap,work,_f);

  for(i=0;i<B_SUPPORT;i++)
    for(j=0;j<B_SUPPORT;j++)
      *(_out+i*_out_stride_i+j*_out_stride_j)=lap[i][j];

#endif
}

#if USE_2D
static double cg_2d_i(double rggt[B_SUPPORT][B_SUPPORT][B_SZ][B_SZ], const int *_f, double _klt[B_SZ*B_SZ][B_SZ*B_SZ]){
  double r[B_SZ][B_SZ][B_SZ][B_SZ];
  double s[B_SUPPORT][B_SZ];
  double ggrggt[B_SZ][B_SZ][B_SZ][B_SZ];
  double cg;
  int i;
  int j;
  int v;
  int u;
  int k;
  int l;

#if PRINT_CG_MATH
  print_matrix(stdout,"rggt",&rggt[0][0][0][0],B_SZ*B_SZ,B_SUPPORT*B_SZ,B_SZ*B_SZ);
#endif

  /* G1*P*G2*R*(G2*P*G1)^T */
  for(v=0;v<B_SZ;v++)
    for(j=0;j<B_SZ;j++)
      b_analysis_2d(&ggrggt[v][j][0][0],1,B_SZ,&rggt[0][0][v][j],B_SZ*B_SZ*B_SUPPORT,B_SZ*B_SZ,f,_klt);

#if PRINT_CG_MATH
  print_matrix(stdout,"ggrggt",&ggrggt[0][0][0][0],B_SZ*B_SZ,B_SZ*B_SZ,B_SZ*B_SZ);
#endif

  /* H1*P*H2 */
  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      for(k=0;k<B_SZ;k++)
        for(l=0;l<B_SZ;l++)
          r[i][j][k][l] = (i*B_SZ+j==k*B_SZ+l)?1:0;

  for(i=0;i<B_SZ;i++)
    for(j=0;j<B_SZ;j++)
      b_synthesis_2d(&rggt[0][0][i][j],B_SZ*B_SZ,B_SUPPORT*B_SZ*B_SZ,
                     &r[i][j][0][0],B_SZ,1,f,_klt);

#if PRINT_CG_MATH
  print_matrix(stdout,"hhi",&rggt[0][0][0][0],B_SZ*B_SZ,B_SUPPORT*B_SUPPORT,B_SZ*B_SZ);
  fprintf(stdout,"G,H=\n");
#endif

  /* ((H1*P*H2)^T*H1*P*H2)_ii */
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      s[i][j]=0;
      for(u=0;u<B_SUPPORT;u++){
        for(v=0;v<B_SUPPORT;v++){
          s[i][j]+=rggt[u][v][i][j]*rggt[u][v][i][j];
        }
      }
    }
  }

  /* (G1*P*G2*R*(G1*P*G2)^T)_ii * ((H1*P*H2)^T*H1*P*H2)_ii */
  cg=0;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
#if PRINT_CG_MATH
      fprintf(stdout,"  %- 12.6G,  %- 12.6G\n",ggrggt[i][j][i][j],s[i][j]);
#endif
      cg-=10*log10(ggrggt[i][j][i][j]*s[i][j]);
    }
  }
  return cg/(B_SZ*B_SZ);
}

double cg_2d(double _in[B_SUPPORT][B_SUPPORT][B_SUPPORT][B_SUPPORT],
             const int *_f){
  int    v;
  int    j;
  double (*rggt)[B_SUPPORT][B_SZ][B_SZ]=malloc(B_SUPPORT*B_SUPPORT*B_SZ*B_SZ*sizeof(****rggt));
  double klt[B_SZ*B_SZ][B_SZ*B_SZ];

#if B_KLT
  gklt_2d(klt,_in,_f);
#endif

  /* R*(G2*P*G1)^T */
  for(v=0;v<B_SUPPORT;v++)
    for(j=0;j<B_SUPPORT;j++)
      b_analysis_2d(&rggt[v][j][0][0],1,B_SZ,&_in[0][0][v][j],B_SUPPORT*B_SUPPORT*B_SUPPORT,B_SUPPORT*B_SUPPORT,_f,klt);

  return cg_2d_i(rggt,f,klt);
}

double cg_2d_collapsed(double _in[B_SUPPORT][B_SUPPORT],const int *_f){
  int    v;
  int    u;
  int    j;
  int    i;
  double r[B_SUPPORT][B_SUPPORT];
  double (*rggt)[B_SUPPORT][B_SZ][B_SZ]=malloc(B_SUPPORT*B_SUPPORT*B_SZ*B_SZ*sizeof(****rggt));
  double klt[B_SZ*B_SZ][B_SZ*B_SZ];

#if B_KLT
  gklt_2d_collapsed(klt,_in,_f);
#endif

  /* R*(G2*P*G1)^T */
  for(v=0;v<B_SUPPORT;v++){
    for(j=0;j<B_SUPPORT;j++){
      for(u=0;u<B_SUPPORT;u++)
        for(i=0;i<B_SUPPORT;i++)
          r[u][i]=_in[abs(u-v)][abs(i-j)];

      b_analysis_2d(&rggt[v][j][0][0],1,B_SZ,&r[0][0],B_SUPPORT,1,_f,klt);
    }
  }
  return cg_2d_i(rggt,f,klt);
}

#else

static double cg_1d_i(double rgt[B_SUPPORT][B_SZ],const int *_f, double klt[B_SZ][B_SZ]){
  int    j;
  int    i;
  double r[B_SZ];
  double grgt[B_SZ][B_SZ];
  double cg;

#if PRINT_CG_MATH
  print_matrix(stdout,"rgt",&rgt[0][0],B_SZ,B_SUPPORT,B_SZ);
#endif

  /* G*R*G^T */
  for(i=0;i<B_SZ;i++)
    b_analysis_1d(&grgt[0][i],B_SZ,&rgt[0][i],B_SZ,_f,klt);

#if PRINT_CG_MATH
  print_matrix(stdout,"grgt",&grgt[0][0],B_SZ,B_SZ,B_SZ);
#endif

  /* H */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      r[i]=i==j?1:0;
    }
    b_synthesis_1d(&rgt[0][j],B_SZ,r,1,_f,klt);
  }

#if PRINT_CG_MATH
  print_matrix(stdout,"hi",&rgt[0][0],B_SZ,B_SUPPORT,B_SZ);
  fprintf(stdout,"G,H=\n");
#endif

  /* (G*R*G^T)_ii * (H^T*H)_ii */
  cg=0;
  for(j=0;j<B_SZ;j++){
    double h;
    h=0;
    for(i=0;i<B_SUPPORT;i++){
      h+=rgt[i][j]*rgt[i][j];
    }
#if PRINT_CG_MATH
    fprintf(stdout,"  %- 12.6G,  %- 12.6G\n",grgt[j][j],h);
#endif
    cg-=10*log10(grgt[j][j]*h);
  }
  return cg/B_SZ;
}

double cg_1d(double in[B_SUPPORT][B_SUPPORT],const int *_f){
  int    j;
  double rgt[B_SUPPORT][B_SZ];
  double klt[B_SZ][B_SZ];

#if B_KLT
  gklt_1d(klt,in,_f);
#endif

  /* R*G^T */
  for(j=0;j<B_SUPPORT;j++){
    b_analysis_1d(&rgt[j][0],1,in[j],1,_f,klt);
  }

  return cg_1d_i(rgt,f,klt);
}

double cg_1d_collapsed(double in[B_SUPPORT],const int *_f){
  int    j;
  int    i;
  double r[B_SUPPORT];
  double rgt[B_SUPPORT][B_SZ];
  double klt[B_SZ][B_SZ];

#if B_KLT
  gklt_1d_collapsed(klt,in,_f);
#endif

  /* R*G^T */
  for(j=0;j<B_SUPPORT;j++){
    for(i=0;i<B_SUPPORT;i++){
      r[i]=in[abs(i-j)];
    }
    b_analysis_1d(&rgt[j][0],1,r,1,_f,klt);
  }

  return cg_1d_i(rgt,f,klt);
}
#endif


int main(int _argc,const char *_argv[]){

#if USE_FILES
  trans_ctx     ctx[NUM_PROCS];
  cov_state     cvs[NUM_PROCS];
  int           i;
#endif

#if USE_2D
  double        cov[B_SUPPORT][B_SUPPORT];
  double        *r=&cov[0][0];
#else
  double        cov[B_SUPPORT];
  double        *r=&cov[0];
#endif

#if B_SZ==4
  f=OD_FILTER_PARAMS4;
#elif B_SZ==8
  f=OD_FILTER_PARAMS8;
#elif B_SZ==16
  f=OD_FILTER_PARAMS16;
#else
# error "Need filter params for this block size."
#endif

#if USE_FILES
  for(i=0;i<NUM_PROCS;i++){
#if USE_2D
    trans_data_init(&ctx[i].td,B_SUPPORT*B_SUPPORT);
    cov_init(&cvs[i],B_SUPPORT*B_SUPPORT);
#else
    trans_data_init(&ctx[i].td,B_SUPPORT);
    cov_init(&cvs[i],B_SUPPORT);
#endif
  }

  OD_OMP_SET_THREADS(NUM_PROCS);
  apply(ctx,cvs,_argc,_argv);
  print_matrix(stdout,"cvs",cvs[0].cov,1,B_SUPPORT,1);
  for(i=1;i<NUM_PROCS;i++){
    trans_data_combine(&ctx[0].td,&ctx[i].td);
    print_matrix(stdout,"cvs",cvs[i].cov,1,B_SUPPORT,1);
    cov_combine(&cvs[0],&cvs[i]);
  }
  print_matrix(stdout,"cvs_comb",cvs[0].cov,1,B_SUPPORT,1);
  trans_data_normalize(&ctx[0].td);
  cov_normalize(&cvs[0]);
# if PRINT_COV
#if USE_2D
  print_matrix(stdout,"cvs",cvs[0].cov,B_SUPPORT,B_SUPPORT*B_SUPPORT,B_SUPPORT);
#else
  fprintf(stdout,"mean_i={ ");
  for(i=0;i<B_SUPPORT;i++){
    fprintf(stdout,"%s  %- 24.18G",i>0?",":"",cvs[0].mean_i[i]);
  }
  fprintf(stdout,"};\nmean_j={ ");
  for(i=0;i<B_SUPPORT;i++){
    fprintf(stdout,"%s  %- 24.18G",i>0?",":"",cvs[0].mean_j[i]);
  }
  fprintf(stdout,"}\n");
  print_matrix(stdout,"cvs",cvs[0].cov,1,B_SUPPORT,1);
#endif
  trans_data_print(&ctx[0].td,stderr);
# endif

# if USE_2D
  fprintf(stdout,"original cg=%- 24.16G\n",
          cg_2d((double(*)[B_SUPPORT][B_SUPPORT][B_SUPPORT])ctx[0].td.cov,f));
  trans_data_collapse(&ctx[0].td,B_SUPPORT,r);
  fprintf(stdout,"collapse cg=%- 24.16G\n",
          cg_2d_collapsed((double(*)[B_SUPPORT])r,f));
  trans_data_expand(&ctx[0].td,B_SUPPORT,r);
  fprintf(stdout,"expanded cg=%- 24.16G\n",
          cg_2d((double(*)[B_SUPPORT][B_SUPPORT][B_SUPPORT])ctx[0].td.cov,f));
# else
  fprintf(stdout,"original cg=%- 24.16G\n",
          cg_1d((double (*)[B_SUPPORT])ctx[0].td.cov,f));
  trans_data_collapse(&ctx[0].td,1,r);
  fprintf(stdout,"collapse cg=%- 24.16G\n",
          cg_1d_collapsed(r,f));
  trans_data_expand(&ctx[0].td,1,r);
  fprintf(stdout,"expanded cg=%- 24.16G\n",
          cg_1d((double (*)[B_SUPPORT])ctx[0].td.cov,f));
  fprintf(stdout,"monty cg=%- 24.16G\n",
          cg_1d_collapsed(cvs[0].cov,f));
# endif

#elif USE_AR95
# if USE_2D
  auto_regressive_collapsed(r,B_SUPPORT*B_SUPPORT,B_SUPPORT,0.95);
# else
  auto_regressive_collapsed(r,B_SUPPORT,1,0.95);
# endif
#endif

#if USE_2D
  fprintf(stdout,"cg=%-24.18G\n",cg_2d_collapsed(cov,f));
#else
  fprintf(stdout,"cg=%-24.18G\n",cg_1d_collapsed(cov,f));
#endif

#if USE_FILES
  for(i=0;i<NUM_PROCS;i++){
    trans_data_clear(&ctx[i].td);
    cov_clear(&cvs[i]);
  }
#endif
  return EXIT_SUCCESS;
}




