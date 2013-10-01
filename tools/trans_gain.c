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

/* 1D coding gain (dB) **********************

   AR p=.95           4x4       8x8     16x16
   ------------------------------------------
   KLT             7.5825    8.8462    9.4781

   DCT             7.5701    8.8259    9.4555

   CDF(9/7)        8.4687    9.4592    9.7866

   LappedKLT       8.5633    9.4908    9.8951

   LappedDCT       8.5523    9.4871    9.8929



   Subset 1           4x4       8x8     16x16
   ------------------------------------------
   KLT original    8.7714   10.2588   11.0039
       collapsed   8.7714   10.2588   11.0039
       monty       8.7654   10.2628   11.0292

   DCT             8.7620   10.2427   10.9861
                   8.7620   10.2427   10.9861
                   8.7561   10.2467   11.0115

   CDF(9/7)        9.3794   10.5932   11.0685
                   9.3845   10.5957   11.0825
                   9.4155   10.6576   11.1965

   LappedKLT       9.6276   10.7860   11.3254
                   9.6277   10.7867   11.3296
                   9.6295   10.8056   11.3722

   LappedDCT       9.6213   10.7832   11.3230
                   9.6214   10.7839   11.3272
                   9.6232   10.8028   11.3698


   Subset 3           4x4       8x8     16x16
   ------------------------------------------
   KLT original   10.5669   12.3711   13.2694
       collapsed  10.5669   12.3711   13.2694
       monty      10.5495   12.3573   13.2729

   DCT            10.5546   12.3532   13.2535
                  10.5547   12.3532   13.2535
                  10.5373   12.3395   13.2572

   CDF(9/7)       11.3102   12.6838   13.1845
                  11.3106   12.6871   13.2009
                  11.3389   12.7764   13.4084

   LappedKLT      11.6048   13.0138   13.6488
                  11.6046   13.0136   13.6491
                  11.5922   13.0126   13.6790

   LappedDCT      11.5970   13.0111   13.6464
                  11.5968   13.0110   13.6467
                  11.5844   13.0099   13.6766

*/

/* 2D coding gain (dB) **********************

   AR p=.95           4x4       8x8     16x16
   ------------------------------------------
   KLT            15.1649   17.6924   18.9562

   DCT            15.1403   17.6518   18.9109

   CDF(9/7)       16.9374   18.9183   19.5731

   LappedKLT      17.1265   18.9816   19.7902

   LappedDCT      17.1047   18.9741   19.7858



   Subset 1           4x4       8x8     16x16
   ------------------------------------------
   KLT original   12.4432   -------   -------
       collapsed  12.4428   -------   -------
       monty      12.4732   13.6167   14.1170

   DCT            12.3695   -------   -------
                  12.3698   -------   -------
                  12.4182   13.5473   14.0536

   CDF(9/7)       -------   -------   -------
                  -------   -------   -------
                  13.1425   13.8184   14.0110

   LappedKLT      13.2807   -------   -------
                  13.2808   -------   -------
                  13.3452   14.1273   14.4041

   LappedDCT      13.2682   -------   -------
                  13.2685   -------   -------
                  13.3330   14.1215   14.3981



   Subset 3           4x4       8x8     16x16
   ------------------------------------------
   KLT monty      14.9078   16.2416   16.7839

   DCT            14.8313   16.1578   16.7221

   CDF(9/7)       15.7553   16.4760   16.6656

   LappedKLT      15.9763   16.8549   17.1181

   LappedDCT      15.9627   16.8507   17.1152


*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include "od_defs.h"
#include "od_filter.h"
#include "stats_tools.h"
#include "trans_tools.h"

#define BLOCKSIZE_LOG  (4)

#define USE_LAPPING    (1)
#define USE_KLT        (1)
#define USE_DCT        (0)
#define USE_WAVELET    (0)

#define USE_2D         (1)
#define USE_FILES      (1)
#define USE_AR95       (0)
#define COMPUTE_NATHAN (1)
#define PRINT_COV      (0)


#define BLOCKSIZE      (1<<BLOCKSIZE_LOG)
#if USE_WAVELET
#if BLOCKSIZE_LOG==1
#  define SUPPORT (20)
#else
#  if BLOCKSIZE_LOG==2
#    define SUPPORT (40)
#  else
#    if BLOCKSIZE_LOG==3
#      define SUPPORT (80)
#    else
#      if BLOCKSIZE_LOG==4
#        define SUPPORT (160)
#      else
#        error "no support configuration for transform size"
#      endif
#    endif
#  endif
#endif
#else
#if USE_LAPPING||COMPUTE_NATHAN
/* larger than needed for 'new' covariance code, but it won't alter
   the answer, just produce a larger than needed covariance matrix.
   It is needed to make the boundary conditions of the 'old'
   covariance code match the trans and trans2d utils */
#define SUPPORT (BLOCKSIZE*2)
#else
#define SUPPORT (BLOCKSIZE)
#endif
#endif

const int    *f;

typedef void (*ne_fdct_func_1d)(double *_out,const double *_in,int _in_stride);
typedef void (*ne_idct_func_1d)(double *_out,int _out_stride,const double *_in);
extern const ne_idct_func_1d OD_IDCT_1D_DOUBLE[OD_NBSIZES];
extern const ne_fdct_func_1d OD_FDCT_1D_DOUBLE[OD_NBSIZES];

#if USE_FILES

typedef struct {
  int sz;
  u_int64_t *n;
  u_int64_t *acc_i;
  u_int64_t *acc_j;
  u_int64_t *acc_ij;
  double *cov;
} cov_state;

static void cov_init(cov_state *_this, int _sz){
  _this->sz    = _sz;
  _this->n     = calloc(_sz,sizeof(*_this->n));
  _this->acc_i = calloc(_sz,sizeof(*_this->acc_i));
  _this->acc_j = calloc(_sz,sizeof(*_this->acc_j));
  _this->acc_ij= calloc(_sz,sizeof(*_this->acc_ij));
  _this->cov   = calloc(_sz,sizeof(*_this->cov));
}

static void cov_clear(cov_state *_this){
  if(_this){
    if(_this->n) free(_this->n);
    if(_this->acc_i) free(_this->acc_i);
    if(_this->acc_j) free(_this->acc_j);
    if(_this->acc_ij) free(_this->acc_ij);
    if(_this->cov) free(_this->cov);
  }
}

#if USE_2D
/* 1D and 2D could both use the same generalized code, but it would be
   harder to read */
static void cov_accumulate_2d(cov_state *_this,
                              const unsigned char *_data,
                              int _stride, int _w, int _h){
  int x,y,i,j;
  int sz = sqrt(_this->sz);
  for(i=0;i<sz;i++){
    for(j=0;j<sz;j++){
      int ij = i*sz+j;
      for(y=0;y<_h-i;y++){
        const unsigned char *di=_data+y*_stride;
        const unsigned char *dj=_data+(y+i)*_stride+j;
        for(x=0;x<_w-j;x++){
          ++_this->n[ij];
          _this->acc_i[ij]  += di[x];
          _this->acc_j[ij]  += dj[x];
          _this->acc_ij[ij] += di[x]*dj[x];
        }
      }
    }
  }
}
#else
static void cov_accumulate_1d(cov_state *_this,
                              const unsigned char *_data,
                              int _stride, int _n){
  int i,j;
  for(i=0;i<_this->sz;i++){
    const unsigned char *di=_data;
    const unsigned char *dj=_data+i*_stride;
    for(j=0;j<_n-i;j++){
      ++_this->n[i];
      _this->acc_i[i] += di[j*_stride];
      _this->acc_j[i] += dj[j*_stride];
      _this->acc_ij[i] += di[j*_stride]*dj[j*_stride];
    }
  }
}
#endif

static void cov_combine(cov_state *_a,const cov_state *_b){
  int i;
  for(i=0;i<_a->sz;i++){
    _a->acc_i[i]  += _b->acc_i[i];
    _a->acc_j[i]  += _b->acc_j[i];
    _a->acc_ij[i] += _b->acc_ij[i];
    _a->n[i]      += _b->n[i];
  }
}

static void cov_compute(cov_state *_this){
  int i;

  for(i=0;i<_this->sz;i++)
    _this->cov[i] =
      ((double)_this->acc_ij[i] -
       (double)_this->acc_i[i]*
       _this->acc_j[i]/_this->n[i])/_this->n[i];
  for(i=1;i<_this->sz;i++)
    _this->cov[i] /= _this->cov[0];

  _this->cov[0]=1.;
}

static void process_files(trans_ctx *_ctx,
                          cov_state *_cov,
                          int _argc,
                          const char *_argv[]){

  int ai;
#pragma omp parallel for schedule(dynamic)
  for(ai=1;ai<_argc;ai++){
    FILE *fin;
    video_input vid;
    video_input_info info;
    video_input_ycbcr ycbcr;
    int tid;
    cov_state *cov;
    int x0,y0,x1,y1;

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
    cov=_cov+tid;
    x0 = info.pic_x;
    y0 = info.pic_y;
    x1 = x0 + info.pic_w;
    y1 = y0 + info.pic_h;

    fprintf(stderr,"%s\n",_argv[ai]);

    /* map */
    {
      int                  stride=ycbcr[0].stride;
      const unsigned char *data=ycbcr[0].data;

#if COMPUTE_NATHAN
      /* block-based full covariance computation (unlord style) */

      int nxblocks=info.pic_w>>BLOCKSIZE_LOG;
      int nyblocks=info.pic_h>>BLOCKSIZE_LOG;
      trans_ctx *ctx=_ctx+tid;

# if USE_2D
      unsigned char        buf[SUPPORT][SUPPORT];
      int                  x,y,i,j;
      image_ctx_init(&ctx->img,_argv[ai],nxblocks,nyblocks);
      for(y=0;y<nyblocks*BLOCKSIZE-SUPPORT+1;y++){
        for(x=0;x<nxblocks*BLOCKSIZE-SUPPORT+1;x++){
          for(j=0;j<SUPPORT;j++){
            for(i=0;i<SUPPORT;i++){
              buf[j][i]=data[(y0+y+j)*stride+(x0+x+i)];
            }
          }
          trans_data_add(&ctx->td,(unsigned char *)buf);
        }
      }
# else
      unsigned char        buf[SUPPORT];
      int                  x,y,z;
      image_ctx_init(&ctx->img,_argv[ai],nxblocks,nyblocks);
      /* add the rows */
      for(y=0;y<nyblocks*BLOCKSIZE;y++){
        for(x=0;x<nxblocks*BLOCKSIZE-SUPPORT+1;x++){
          for(z=0;z<SUPPORT;z++){
            buf[z]=data[(y+y0)*stride+x+x0+z];
          }
          trans_data_add(&ctx->td,buf);
        }
      }
      /* add the columns */
      for(y=0;y<nyblocks*BLOCKSIZE-SUPPORT+1;y++){
        for(x=0;x<nxblocks*BLOCKSIZE;x++){
          for(z=0;z<SUPPORT;z++){
            buf[z]=data[(y0+y+z)*stride+x+x0];
          }
          trans_data_add(&ctx->td,buf);
        }
      }

# endif
#endif

      /* Direct computation of collapsed covariance matrix (monty style) */
#if USE_2D
      cov_accumulate_2d(cov,data+y0*stride+x0,stride,x1-x0,y1-y0);
#else
      {
        int x,y;

        for(y=y0;y<y1;y++)
          cov_accumulate_1d(cov,data+y*stride+x0,1,x1-x0);
        for(x=x0;x<x1;x++)
          cov_accumulate_1d(cov,data+y0*stride+x,stride,y1-y0);
      }
#endif

    }

    video_input_close(&vid);
  }
}
#endif

#if USE_WAVELET

/* some lifting CDF (9/7) wavelet code from Google Code's axonlib */
/* http://code.google.com/p/axonlib/source/browse/trunk/extern/dwt97.c?spec=svn19&r=19 */
/* single stage of decomposition */

static void fwt97_i(double* x,int n){
  double temp[SUPPORT];
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
  double temp[SUPPORT];
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

#if USE_KLT
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
  /* for(j=0;j<BLOCKSIZE;j++)eigenvalue[j]=cov[j][j]; don't need eigenvalues */
}

void flap_2d(double out[BLOCKSIZE][BLOCKSIZE],
            double in[SUPPORT][SUPPORT],
            const int _f[]){
  int i,j;

#if USE_LAPPING
#if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES

  /* columns */
  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    double work[BLOCKSIZE*2];
    for(j=0;j<BLOCKSIZE*2;j++)
      work[j]=in[j+SUPPORT/2-BLOCKSIZE][i];
    (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&work[0],&work[0],_f);
    (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&work[BLOCKSIZE],&work[BLOCKSIZE],_f);
    for(j=0;j<BLOCKSIZE*2;j++)
      in[j+SUPPORT/2-BLOCKSIZE][i]=work[j];
  }

  /* rows */
  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&in[i][SUPPORT/2-BLOCKSIZE],&in[i][SUPPORT/2-BLOCKSIZE],_f);
    (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&in[i][SUPPORT/2],&in[i][SUPPORT/2],_f);
  }

#else
#  error "Need a prefilter implementation for this block size."
#endif
#endif

  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      out[i][j]=in[i+SUPPORT/2-BLOCKSIZE/2][j+SUPPORT/2-BLOCKSIZE/2];
}

void ilap_2d(double out[SUPPORT][SUPPORT],
              double in[BLOCKSIZE][BLOCKSIZE],
              const int _f[]){
  int i,j;

  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      out[i][j]=0;

  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      out[i+SUPPORT/2-BLOCKSIZE/2][j+SUPPORT/2-BLOCKSIZE/2]=in[i][j];

#if USE_LAPPING
#if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES

  /* columns */
  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    double work[BLOCKSIZE*2];
    for(j=0;j<BLOCKSIZE*2;j++)
      work[j]=out[j+SUPPORT/2-BLOCKSIZE][i];
    (*NE_POST_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&work[0],&work[0],_f);
    (*NE_POST_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&work[BLOCKSIZE],&work[BLOCKSIZE],_f);
    for(j=0;j<BLOCKSIZE*2;j++)
      out[j+SUPPORT/2-BLOCKSIZE][i]=work[j];
  }

  /* rows */
  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    (*NE_POST_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&out[i][SUPPORT/2-BLOCKSIZE],&out[i][SUPPORT/2-BLOCKSIZE],_f);
    (*NE_POST_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
      (&out[i][SUPPORT/2],&out[i][SUPPORT/2],_f);
  }

#else
#  error "Need a prefilter implementation for this block size."
#endif
#endif

}

void flap_4d(double out[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE],
            double in[SUPPORT][SUPPORT][SUPPORT][SUPPORT],
            const int _f[]){
  int i,j,k,l;

#if USE_LAPPING
#if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES

  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    for(j=SUPPORT/2-BLOCKSIZE;j<SUPPORT/2+BLOCKSIZE;j++){
      for(k=SUPPORT/2-BLOCKSIZE;k<SUPPORT/2+BLOCKSIZE;k++){
        double work[BLOCKSIZE*2];

        /* [ ][i][j][k] */
        for(l=0;l<BLOCKSIZE*2;l++)
          work[l]=in[l+SUPPORT/2-BLOCKSIZE][i][j][k];
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&work[0],&work[0],_f);
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&work[BLOCKSIZE],&work[BLOCKSIZE],_f);
        for(l=0;l<BLOCKSIZE*2;l++)
          in[l+SUPPORT/2-BLOCKSIZE][i][j][k]=work[l];
      }
    }
  }

  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    for(j=SUPPORT/2-BLOCKSIZE;j<SUPPORT/2+BLOCKSIZE;j++){
      for(k=SUPPORT/2-BLOCKSIZE;k<SUPPORT/2+BLOCKSIZE;k++){
        double work[BLOCKSIZE*2];
        /* [i][ ][j][k] */
        for(l=0;l<BLOCKSIZE*2;l++)
          work[l]=in[i][l+SUPPORT/2-BLOCKSIZE][j][k];
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&work[0],&work[0],_f);
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&work[BLOCKSIZE],&work[BLOCKSIZE],_f);
        for(l=0;l<BLOCKSIZE*2;l++)
          in[i][l+SUPPORT/2-BLOCKSIZE][j][k]=work[l];
      }
    }
  }


  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    for(j=SUPPORT/2-BLOCKSIZE;j<SUPPORT/2+BLOCKSIZE;j++){
      for(k=SUPPORT/2-BLOCKSIZE;k<SUPPORT/2+BLOCKSIZE;k++){
        double work[BLOCKSIZE*2];
        /* [i][j][ ][k] */
        for(l=0;l<BLOCKSIZE*2;l++)
          work[l]=in[i][j][l+SUPPORT/2-BLOCKSIZE][k];
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&work[0],&work[0],_f);
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&work[BLOCKSIZE],&work[BLOCKSIZE],_f);
        for(l=0;l<BLOCKSIZE*2;l++)
          in[i][j][l+SUPPORT/2-BLOCKSIZE][k]=work[l];
      }
    }
  }

  for(i=SUPPORT/2-BLOCKSIZE;i<SUPPORT/2+BLOCKSIZE;i++){
    for(j=SUPPORT/2-BLOCKSIZE;j<SUPPORT/2+BLOCKSIZE;j++){
      for(k=SUPPORT/2-BLOCKSIZE;k<SUPPORT/2+BLOCKSIZE;k++){
        /* [i][j][k][ ] */
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&in[i][j][k][SUPPORT/2-BLOCKSIZE],&in[i][j][k][SUPPORT/2-BLOCKSIZE],_f);
        (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
          (&in[i][j][k][SUPPORT/2],&in[i][j][k][SUPPORT/2],_f);

      }
    }
  }

#else
#  error "Need a prefilter implementation for this block size."
#endif
#endif

  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      for(k=0;k<BLOCKSIZE;k++)
        for(l=0;l<BLOCKSIZE;l++)
          out[i*BLOCKSIZE+j][k*BLOCKSIZE+l]=in
            [i+SUPPORT/2-BLOCKSIZE/2]
            [j+SUPPORT/2-BLOCKSIZE/2]
            [k+SUPPORT/2-BLOCKSIZE/2]
            [l+SUPPORT/2-BLOCKSIZE/2];
}

void gklt_1d(double klt[BLOCKSIZE][BLOCKSIZE],
             double cov[SUPPORT][SUPPORT],
             const int *_f){
  static double workA[SUPPORT][SUPPORT];
  static double workB[BLOCKSIZE][BLOCKSIZE];
  int i,j;
  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      workA[i][j]=cov[i][j];
  flap_2d(workB,workA,_f);
  symeigen(&klt[0][0],&workB[0][0],BLOCKSIZE);
}

void gklt_2d(double klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE],
             double cov[SUPPORT][SUPPORT][SUPPORT][SUPPORT],
             const int *_f){
  static double workA[SUPPORT][SUPPORT][SUPPORT][SUPPORT];
  static double workB[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE];
  int i,j,k,l;
  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      for(k=0;k<SUPPORT;k++)
        for(l=0;l<SUPPORT;l++)
          workA[i][j][k][l]=cov[i][j][k][l];
  flap_4d(workB,workA,_f);
  symeigen(&klt[0][0],&workB[0][0],BLOCKSIZE*BLOCKSIZE);
}

void gklt_1d_collapsed(double klt[BLOCKSIZE][BLOCKSIZE],
                       double cov[SUPPORT],
                       const int *_f){
  static double workA[SUPPORT][SUPPORT];
  static double workB[BLOCKSIZE][BLOCKSIZE];
  int i,j;
  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      workA[i][j]=cov[abs(i-j)];
  flap_2d(workB,workA,_f);
  symeigen(&klt[0][0],&workB[0][0],BLOCKSIZE);
}

void gklt_2d_collapsed(double klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE],
                       double cov[SUPPORT][SUPPORT],
                       const int *_f){
  static double workA[SUPPORT][SUPPORT][SUPPORT][SUPPORT];
  static double workB[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE];
  int i,j,k,l;
  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      for(k=0;k<SUPPORT;k++)
        for(l=0;l<SUPPORT;l++)
          workA[i][j][k][l]=cov[abs(i-k)][abs(j-l)];
  flap_4d(workB,workA,_f);
  symeigen(&klt[0][0],&workB[0][0],BLOCKSIZE*BLOCKSIZE);
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
                   const int *_f, double _klt[BLOCKSIZE][BLOCKSIZE]){
  int    j;
  double t[SUPPORT];
  double w[BLOCKSIZE];

  for(j=0;j<SUPPORT;j++)
    t[j]=_in[j*_in_stride];

#if USE_WAVELET

  fwt97(w,BLOCKSIZE,t,SUPPORT);
  for(j=0;j<BLOCKSIZE;j++){
    _out[j*_out_stride]=w[j];
  }

#else
#  if USE_LAPPING
#    if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
    (&t[SUPPORT/2-BLOCKSIZE],&t[SUPPORT/2-BLOCKSIZE],_f);
  (*NE_PRE_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
    (&t[SUPPORT/2],&t[SUPPORT/2],_f);
#    else
#      error "Need a prefilter implementation for this block size."
#    endif
#  endif
#  if USE_KLT
  fklt(&w[0],&t[SUPPORT/2-BLOCKSIZE/2],&_klt[0][0],BLOCKSIZE);
#  elif USE_DCT
#    if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_1D_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
    (w,&t[SUPPORT/2-BLOCKSIZE/2],1);
#    else
#      error "Need an fDCT implementation for this block size."
#    endif
#  else
    for(j=0;j<BLOCKSIZE;j++)
      w[j]=t[j+SUPPORT/2-BLOCKSIZE/2];
#  endif

  for(j=0;j<BLOCKSIZE;j++)
    _out[j*_out_stride]=w[j];

#endif
}

void b_analysis_2d(double *_out,int _out_stride_i,int _out_stride_j,
                   const double *_in,int _in_stride_i,int _in_stride_j,
                   const int *_f, double _klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE]){
#if USE_KLT

  /* KLT is a non-separable 2D transform */
  double lap[SUPPORT][SUPPORT];
  double work[BLOCKSIZE][BLOCKSIZE];
  double temp[BLOCKSIZE][BLOCKSIZE];
  int i,j;
  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      lap[i][j]=*(_in+i*_in_stride_i+j*_in_stride_j);

  flap_2d(work,lap,_f);

  fklt(&temp[0][0],&work[0][0],&_klt[0][0],BLOCKSIZE*BLOCKSIZE);

  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      *(_out+i*_out_stride_i+j*_out_stride_j)=temp[i][j];

#else

  double work[SUPPORT][BLOCKSIZE];
  int i;

  /* DCT and DWT are separable 1D transforms */
  /* lapping performed inside b_analysis */
  for(i=0;i<SUPPORT;i++)
    b_analysis_1d(&work[i][0],1,_in+i*_in_stride_i,_in_stride_j,_f,NULL);
  for(i=0;i<BLOCKSIZE;i++)
    b_analysis_1d(_out+_out_stride_i*i,_out_stride_j,&work[0][i],BLOCKSIZE,_f,NULL);

#endif

}


void b_synthesis_1d(double *_out,int _out_stride,const double *_in,int _in_stride,
                    const int *_f, double _klt[BLOCKSIZE][BLOCKSIZE]){
  int    j;
  double w[SUPPORT];
  double t[SUPPORT];

  for(j=0;j<SUPPORT;j++){
    t[j]=0;
    w[j]=0;
  }

#if USE_WAVELET
  for(j=0;j<BLOCKSIZE;j++)
    w[j]=_in[j*_in_stride];
  iwt97(t,SUPPORT,w,BLOCKSIZE);
#else
  for(j=0;j<BLOCKSIZE;j++){
    w[SUPPORT/2-BLOCKSIZE/2+j]=_in[j*_in_stride];
  }
#  if USE_KLT
  iklt(&t[SUPPORT/2-BLOCKSIZE/2],&w[SUPPORT/2-BLOCKSIZE/2],&_klt[0][0],BLOCKSIZE);
#  elif USE_DCT
#    if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_IDCT_1D_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
    (&t[SUPPORT/2-BLOCKSIZE/2],1,&w[SUPPORT/2-BLOCKSIZE/2]);
#    else
#      error "Need an iDCT implementation for this block size."
#    endif
#  else
    for(j=0;j<SUPPORT;j++)
      t[j]=w[j];
#  endif

#  if USE_LAPPING
#    if BLOCKSIZE_LOG>=OD_LOG_BSIZE0&&BLOCKSIZE_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*NE_POST_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
    (&t[SUPPORT/2-BLOCKSIZE],&t[SUPPORT/2-BLOCKSIZE],_f);
  (*NE_POST_FILTER_DOUBLE[BLOCKSIZE_LOG-OD_LOG_BSIZE0])
    (&t[SUPPORT/2],&t[SUPPORT/2],_f);
#    else
# error "Need a postfilter implementation for this block size."
#    endif
#  endif

#endif
  for(j=0;j<SUPPORT;j++)
    _out[j*_out_stride]=t[j];
}

void b_synthesis_2d(double *_out,int _out_stride_i,int _out_stride_j,
                    const double *_in,int _in_stride_i,int _in_stride_j,
                    const int *_f, double _klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE]){
#if USE_KLT

  /* KLT is a non-separable 2D transform */
  double temp[BLOCKSIZE][BLOCKSIZE];
  double work[BLOCKSIZE][BLOCKSIZE];
  double lap[SUPPORT][SUPPORT];
  int i,j;

  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      temp[i][j]=*(_in+i*_in_stride_i+j*_in_stride_j);

  iklt(&work[0][0],&temp[0][0],&_klt[0][0],BLOCKSIZE*BLOCKSIZE);

  ilap_2d(lap,work,_f);

  for(i=0;i<SUPPORT;i++)
    for(j=0;j<SUPPORT;j++)
      *(_out+i*_out_stride_i+j*_out_stride_j)=lap[i][j];

#else

  double work[SUPPORT][BLOCKSIZE];
  int i;

  /* DCT and DWT are separable 1D transforms */
  /* lapping performed inside b_analysis */
  for(i=0;i<BLOCKSIZE;i++)
    b_synthesis_1d(&work[0][i],BLOCKSIZE,_in+i*_in_stride_i,_in_stride_j,_f,NULL);
  for(i=0;i<SUPPORT;i++)
    b_synthesis_1d(_out+_out_stride_i*i,_out_stride_j,&work[i][0],1,_f,NULL);

#endif

}

#if USE_2D
static double cg_2d_i(double rggt[SUPPORT][SUPPORT][BLOCKSIZE][BLOCKSIZE],
                      const int *_f,
                      double _klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE]){
  double r[BLOCKSIZE][BLOCKSIZE][BLOCKSIZE][BLOCKSIZE];
  double s[SUPPORT][BLOCKSIZE];
  double ggrggt[BLOCKSIZE][BLOCKSIZE][BLOCKSIZE][BLOCKSIZE];
  double cg=0;
  int i;
  int j;
  int v;
  int u;
  int k;
  int l;

  /* G1*P*G2*R*(G2*P*G1)^T */
  for(v=0;v<BLOCKSIZE;v++)
    for(j=0;j<BLOCKSIZE;j++)
      b_analysis_2d(&ggrggt[v][j][0][0],
                    1,BLOCKSIZE,
                    &rggt[0][0][v][j],
                    BLOCKSIZE*BLOCKSIZE*SUPPORT,
                    BLOCKSIZE*BLOCKSIZE,
                    f,_klt);

  /* H1*P*H2 */
  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      for(k=0;k<BLOCKSIZE;k++)
        for(l=0;l<BLOCKSIZE;l++)
          r[i][j][k][l] = (i*BLOCKSIZE+j==k*BLOCKSIZE+l)?1:0;

  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      b_synthesis_2d(&rggt[0][0][i][j],
                     BLOCKSIZE*BLOCKSIZE,
                     SUPPORT*BLOCKSIZE*BLOCKSIZE,
                     &r[i][j][0][0],
                     BLOCKSIZE,1,
                     _f,_klt);

  /* ((H1*P*H2)^T*H1*P*H2)_ii */
  for(i=0;i<BLOCKSIZE;i++){
    for(j=0;j<BLOCKSIZE;j++){
      s[i][j]=0;
      for(u=0;u<SUPPORT;u++){
        for(v=0;v<SUPPORT;v++){
          s[i][j]+=rggt[u][v][i][j]*rggt[u][v][i][j];
        }
      }
    }
  }

  /* (G1*P*G2*R*(G1*P*G2)^T)_ii * ((H1*P*H2)^T*H1*P*H2)_ii */
  for(i=0;i<BLOCKSIZE;i++)
    for(j=0;j<BLOCKSIZE;j++)
      cg-=10*log10(ggrggt[i][j][i][j]*s[i][j]);

  return cg/(BLOCKSIZE*BLOCKSIZE);
}

double cg_2d(double _in[SUPPORT][SUPPORT][SUPPORT][SUPPORT],
             const int *_f){
  int    v;
  int    j;
  double ret;
  double (*rggt)[SUPPORT][BLOCKSIZE][BLOCKSIZE] =
    malloc(SUPPORT*SUPPORT*BLOCKSIZE*BLOCKSIZE*sizeof(****rggt));
  double klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE];

#if USE_KLT
  gklt_2d(klt,_in,_f);
#endif

  /* R*(G2*P*G1)^T */
  for(v=0;v<SUPPORT;v++)
    for(j=0;j<SUPPORT;j++)
      b_analysis_2d(&rggt[v][j][0][0],
                    1,BLOCKSIZE,
                    &_in[0][0][v][j],
                    SUPPORT*SUPPORT*SUPPORT,
                    SUPPORT*SUPPORT,
                    _f,klt);

  ret = cg_2d_i(rggt,f,klt);
  free(rggt);
  return ret;
}

double cg_2d_collapsed(double _in[SUPPORT][SUPPORT],const int *_f){
  int    v;
  int    u;
  int    j;
  int    i;
  double ret;
  double r[SUPPORT][SUPPORT];
  double (*rggt)[SUPPORT][BLOCKSIZE][BLOCKSIZE] =
    malloc(SUPPORT*SUPPORT*BLOCKSIZE*BLOCKSIZE*sizeof(****rggt));
  double klt[BLOCKSIZE*BLOCKSIZE][BLOCKSIZE*BLOCKSIZE];

#if USE_KLT
  gklt_2d_collapsed(klt,_in,_f);
#endif

  /* R*(G2*P*G1)^T */
  for(v=0;v<SUPPORT;v++){
    for(j=0;j<SUPPORT;j++){
      for(u=0;u<SUPPORT;u++)
        for(i=0;i<SUPPORT;i++)
          r[u][i]=_in[abs(u-v)][abs(i-j)];

      b_analysis_2d(&rggt[v][j][0][0],
                    1,BLOCKSIZE,
                    &r[0][0],SUPPORT,1,_f,klt);
    }
  }
  ret = cg_2d_i(rggt,f,klt);
  free(rggt);
  return ret;
}

#else

static double cg_1d_i(double rgt[SUPPORT][BLOCKSIZE],
                      const int *_f,
                      double klt[BLOCKSIZE][BLOCKSIZE]){
  int    j;
  int    i;
  double r[BLOCKSIZE];
  double grgt[BLOCKSIZE][BLOCKSIZE];
  double cg=0;

  /* G*R*G^T */
  for(i=0;i<BLOCKSIZE;i++)
    b_analysis_1d(&grgt[0][i],BLOCKSIZE,&rgt[0][i],BLOCKSIZE,_f,klt);

  /* H */
  for(j=0;j<BLOCKSIZE;j++){
    for(i=0;i<BLOCKSIZE;i++){
      r[i]=i==j?1:0;
    }
    b_synthesis_1d(&rgt[0][j],BLOCKSIZE,r,1,_f,klt);
  }

  /* (G*R*G^T)_ii * (H^T*H)_ii */
  for(j=0;j<BLOCKSIZE;j++){
    double h=0;
    for(i=0;i<SUPPORT;i++){
      h+=rgt[i][j]*rgt[i][j];
    }
    cg-=10*log10(grgt[j][j]*h);
  }
  return cg/BLOCKSIZE;
}

double cg_1d(double in[SUPPORT][SUPPORT],const int *_f){
  int    j;
  double rgt[SUPPORT][BLOCKSIZE];
  double klt[BLOCKSIZE][BLOCKSIZE];

#if USE_KLT
  gklt_1d(klt,in,_f);
#endif

  /* R*G^T */
  for(j=0;j<SUPPORT;j++){
    b_analysis_1d(&rgt[j][0],1,in[j],1,_f,klt);
  }

  return cg_1d_i(rgt,f,klt);
}

double cg_1d_collapsed(double in[SUPPORT],const int *_f){
  int    j;
  int    i;
  double r[SUPPORT];
  double rgt[SUPPORT][BLOCKSIZE];
  double klt[BLOCKSIZE][BLOCKSIZE];

#if USE_KLT
  gklt_1d_collapsed(klt,in,_f);
#endif

  /* R*G^T */
  for(j=0;j<SUPPORT;j++){
    for(i=0;i<SUPPORT;i++){
      r[i]=in[abs(i-j)];
    }
    b_analysis_1d(&rgt[j][0],1,r,1,_f,klt);
  }

  return cg_1d_i(rgt,f,klt);
}
#endif

#if USE_FILES
int main(int _argc,const char *_argv[]){

  cov_state     cvs[NUM_PROCS];
#if COMPUTE_NATHAN
  trans_ctx     ctx[NUM_PROCS];
  double        r[SUPPORT*SUPPORT]; /* maximum for 2d */
#else
  trans_ctx     *ctx=NULL;
#endif
  int           i;

#if BLOCKSIZE==4
  f=OD_FILTER_PARAMS4;
#elif BLOCKSIZE==8
  f=OD_FILTER_PARAMS8;
#elif BLOCKSIZE==16
  f=OD_FILTER_PARAMS16;
#else
# error "Need filter params for this block size."
#endif

  for(i=0;i<NUM_PROCS;i++){
#if USE_2D
    cov_init(&cvs[i],SUPPORT*SUPPORT);
#else
    cov_init(&cvs[i],SUPPORT);
#endif
  }

#if COMPUTE_NATHAN
  for(i=0;i<NUM_PROCS;i++){
#if USE_2D
    trans_data_init(&ctx[i].td,SUPPORT*SUPPORT);
#else
    trans_data_init(&ctx[i].td,SUPPORT);
#endif
  }
#endif

  OD_OMP_SET_THREADS(NUM_PROCS);
  process_files(ctx,cvs,_argc,_argv);

  for(i=1;i<NUM_PROCS;i++)
      cov_combine(&cvs[0],&cvs[i]);
  cov_compute(&cvs[0]);

#if COMPUTE_NATHAN
  for(i=1;i<NUM_PROCS;i++)
    trans_data_combine(&ctx[0].td,&ctx[i].td);
  trans_data_normalize(&ctx[0].td);
#endif

#if PRINT_COV
  {
    int i,j;
    fprintf(stdout,"collapsed_cov=\n");
    for(j=0;j<cvs[0].sz/SUPPORT;j++){
      for(i=0;i<SUPPORT;i++){
        fprintf(stdout,"%s  %- 12.6G",i>0?",":"",cvs[0].cov[j*SUPPORT+i]);
      }
      fprintf(stdout,"\n");
    }
  }
#endif

#if USE_2D
#if COMPUTE_NATHAN
  fprintf(stdout,"original cg=%-24.16G\n",
          cg_2d((double(*)[SUPPORT][SUPPORT][SUPPORT])ctx[0].td.cov,f));
  trans_data_collapse(&ctx[0].td,SUPPORT,r);
  fprintf(stdout,"collapse cg=%-24.16G\n",
          cg_2d_collapsed((double(*)[SUPPORT])r,f));
#endif
  fprintf(stdout,"monty cg=%-24.16G\n",
          cg_2d_collapsed((double(*)[SUPPORT])cvs[0].cov,f));
#else
#if COMPUTE_NATHAN
  fprintf(stdout,"original cg=%-24.16G\n",
          cg_1d((double (*)[SUPPORT])ctx[0].td.cov,f));
  trans_data_collapse(&ctx[0].td,1,r);
  fprintf(stdout,"collapse cg=%-24.16G\n",
          cg_1d_collapsed(r,f));
#endif
  fprintf(stdout,"monty cg=%-24.16G\n",
          cg_1d_collapsed(cvs[0].cov,f));
#endif

  for(i=0;i<NUM_PROCS;i++)
    cov_clear(&cvs[i]);

#if COMPUTE_NATHAN
  for(i=0;i<NUM_PROCS;i++)
    trans_data_clear(&ctx[i].td);
#endif

  return EXIT_SUCCESS;
}

#else

int main(int _argc,const char *_argv[]){

#if USE_2D
  double        cov[SUPPORT][SUPPORT];
  double        *r=&cov[0][0];
#else
  double        cov[SUPPORT];
  double        *r=&cov[0];
#endif

#if BLOCKSIZE==4
  f=OD_FILTER_PARAMS4;
#elif BLOCKSIZE==8
  f=OD_FILTER_PARAMS8;
#elif BLOCKSIZE==16
  f=OD_FILTER_PARAMS16;
#else
# error "Need filter params for this block size."
#endif

# if USE_2D
  auto_regressive_collapsed(r,SUPPORT*SUPPORT,SUPPORT,0.95);
  fprintf(stdout,"AR p=.95 cg=%-24.18G\n",cg_2d_collapsed(cov,f));
# else
  auto_regressive_collapsed(r,SUPPORT,1,0.95);
  fprintf(stdout,"AR p=.95 cg=%-24.18G\n",cg_1d_collapsed(cov,f));
# endif

  return EXIT_SUCCESS;
}

#endif
