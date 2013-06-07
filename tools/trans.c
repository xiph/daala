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

#include <stdlib.h>
#include "od_defs.h"
#include "od_filter.h"
#include "stats_tools.h"
#include "trans_tools.h"
#include "int_search.h"
#include "kiss99.h"

#define USE_FILES   (0)
#define USE_AR95    (1)
#define USE_SUBSET1 (0)
#define USE_SUBSET3 (0)

#define PRINT_COV   (0)

#define CG_SEARCH   (0)
#define USE_SIMPLEX (1)
#define RAMP_DYADIC (0)

#if CG_SEARCH
# if USE_TYPE3 && RAMP_DYADIC
#  error "Dyadic ramp constraint not supported for Type-III transform."
# endif
# if USE_SIMPLEX && RAMP_DYADIC
#  error "Dyadic ramp constraint not supported with simplex search."
# endif
static void coding_gain_search(const double _r[2*B_SZ]){
# if !USE_SIMPLEX
#  if B_SZ==4
  {
    int    f[4];
    int    p0;
    int    q0;
    int    s0;
    int    s1;
    double cg;
    double best_cg;
    best_cg=0;
#   if RAMP_DYADIC
    for(q0=(1<<FILTER_BITS);q0>=-(1<<FILTER_BITS);q0--){
      int t0;
      f[3]=q0;
      /* S0 = 4/1*(1-q0/64)
       * S0 >= 1 -> 64-q0 >= 16
       */
      t0=(1<<FILTER_BITS)-q0;
      s0=1*t0-0;
      if(s0>=(1<<FILTER_BITS-2)){
        s0*=4;
        f[0]=s0;
        for(p0=-(1<<FILTER_BITS);p0<=(1<<FILTER_BITS);p0++){
          f[2]=p0;
          /* S1 = 4/3*(1-(1-q0/64)*p0/64)
           * S1 >= 1   -> 64^2-(64-q0)*p0 >= 64*48
           * S1 = x/64 -> 64^2-(64-q0)*p0 = 0 MOD 48
           */
          s1=(1<<2*FILTER_BITS)-t0*p0;
          if(s1>=(1<<FILTER_BITS)*(3<<FILTER_BITS-2)&&s1%(3<<FILTER_BITS-2)==0){
            s1/=(3<<FILTER_BITS-2);
            f[1]=s1;
            cg=coding_gain_1d_collapsed(_r,f);
            if(cg>best_cg){
              best_cg=cg;
              printf("%i %i %i %i %G\n",p0,q0,s0,s1,cg);
            }
          }
        }
      }
    }
#   else
    for(p0=-(1<<FILTER_BITS);p0<=(1<<FILTER_BITS);p0++){
      f[2]=p0;
      for(q0=(1<<FILTER_BITS);q0>=-(1<<FILTER_BITS);q0--){
        f[3]=q0;
        for(s0=(1<<FILTER_BITS);s0<=2*(1<<FILTER_BITS);s0++){
          f[0]=s0;
          for(s1=(1<<FILTER_BITS);s1<=2*(1<<FILTER_BITS);s1++){
            f[1]=s1;
            cg=coding_gain_1d_collapsed(_r,f);
            if(cg>best_cg){
              best_cg=cg;
              printf("%i %i %i %i %G\n",p0,q0,s0,s1,cg);
            }
          }
        }
      }
    }
#   endif
  }
#  elif B_SZ==8
  {
    int    f[10];
    int    p0;
    int    p1;
    int    p2;
    int    q0;
    int    q1;
    int    q2;
    int    s0;
    int    s1;
    int    s2;
    int    s3;
    double cg;
    double best_cg;
    best_cg=0;
#   if RAMP_DYADIC
    for(q0=(1<<FILTER_BITS);q0>=-(1<<FILTER_BITS);q0--){
      int t0;
      f[7]=q0;
      /* S0 = 8/1*(1-q0/64)
       * S0 >= 1 -> 64-q0 >= 8
       */
      t0=(1<<FILTER_BITS)-q0;
      s0=1*t0-0;
      if(s0>=(1<<FILTER_BITS-3)){
        s0*=8;
        f[0]=s0;
        for(p0=-(1<<FILTER_BITS);p0<=(1<<FILTER_BITS);p0++){
          f[4]=p0;
          for(q1=(1<<FILTER_BITS);q1>=-(1<<FILTER_BITS);q1--){
            int t1;
            f[8]=q1;
            /* S1 = 8/3*((1-q1/64)-(1-q0/64)*p0/64)
             * S1 >= 1   -> 64*t1-t0*p0 >= 64*24
             * S1 = x/64 -> 64*t1-t0*p0 = 0 MOD 24
             */
            t1=(1<<FILTER_BITS)-q1;
            s1=(1<<FILTER_BITS)*t1-t0*p0;
            if(s1>=(1<<FILTER_BITS)*(3<<FILTER_BITS-3)&&
             s1%(3<<FILTER_BITS-3)==0){
              s1/=(3<<FILTER_BITS-3);
              f[1]=s1;
              for(p1=-(1<<FILTER_BITS);p1<=(1<<FILTER_BITS);p1++){
                f[5]=p1;
                for(q2=(1<<FILTER_BITS);q2>=-(1<<FILTER_BITS);q2--){
                  int t2;
                  f[9]=q2;
                  /* S2 = 8/5*((1-q2/64)-(1-q1/64)*p1/64)
                   * S2 >= 1   -> 64*t2-t1*p1) >= 64*40
                   * S2 = x/64 -> 64*t2-t1*p1 = 0 MOD 40
                   */
                  t2=(1<<FILTER_BITS)-q2;
                  s2=(1<<FILTER_BITS)*t2-t1*p1;
                  if(s2>=(1<<FILTER_BITS)*(5<<FILTER_BITS-3)&&
                   s2%(5<<FILTER_BITS-3)==0){
                    s2/=(5<<FILTER_BITS-3);
                    f[2]=s2;
                    for(p2=-(1<<FILTER_BITS);p2<=(1<<FILTER_BITS);p2++){
                      f[6]=p2;
                      /* S3 = 8/7*(1-(1-q2/64)*p2/64)
                       * S3 >= 1   -> 64^2-t2*p2 >= 64*56
                       * S3 = x/64 -> 64^2-t2*p2 = 0 MOD 56
                       */
                      s3=(1<<2*FILTER_BITS)-t2*p2;
                      if(s3>=(1<<FILTER_BITS)*(7<<FILTER_BITS-3)&&
                       s3%(7<<FILTER_BITS-3)==0){
                        s3/=(7<<FILTER_BITS-3);
                        f[3]=s3;
                        cg=coding_gain_1d_collapsed(_r,f);
                        if(cg>best_cg){
                          best_cg=cg;
                          printf("%i %i %i %i %i %i %i %i %i %i %-24.18G\n",
                           p0,p1,p2,q0,q1,q2,s0,s1,s2,s3,cg);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
#   else
#    error "Exhaustive search for B_SZ==8 only supported using RAMP_DYADIC (1)."
#   endif
  }
#  else
#   error "Exhaustive search not supported for this block size."
#  endif
# else
  {
    int         dims;
    int         i;
    kiss99_ctx  ks[NUM_PROCS];
    int         lb[22];
    int         ub[22];
#  if B_SZ==4
    dims=4;
#  elif B_SZ==8
    dims=10;
#  elif B_SZ==16
    dims=22;
#  else
#   error "Unsupported block size."
#  endif
    for(i=0;i<dims;i++){
      lb[i]=i<(B_SZ>>1)?(1<<FILTER_BITS):-(1<<FILTER_BITS);
      ub[i]=i<(B_SZ>>1)?2*(1<<FILTER_BITS):(1<<FILTER_BITS);
    }
    for(i=0;i<NUM_PROCS;i++){
      uint32_t srand;
      srand=i*16843009; /*Broadcast char to 4xchar*/
      kiss99_srand(&ks[i],(unsigned char *)&srand,sizeof(srand));
    }
    #pragma omp parallel for schedule(dynamic)
    for(i=0;i<128;i++){
      int    tid;
      int    j;
#  if B_SZ==4
      int    f[4];
#  elif B_SZ==8
      int    f[10];
#  elif B_SZ==16
      int    f[22];
#  else
#   error "Unsupported block size."
#  endif
      double cg;
      tid=OD_OMP_GET_THREAD;
      for(j=0;j<dims;j++){
        int range;
        int mask;
        int rng;
        range=ub[j]-lb[j];
        mask=(1<<OD_ILOG_NZ(range))-1;
        do {
          rng=((int)kiss99_rand(&ks[tid]))&mask;
        }
        while(rng>range);
        f[j]=lb[j]+rng;
      }
      j=int_simplex_max(&cg,dims,coding_gain_1d_collapsed,_r,lb,ub,f);
      fprintf(stdout,"obj=%-24.18G steps=%4d params={",cg,j);
      for(j=0;j<dims;j++){
        fprintf(stdout,"%3d%c",f[j],j==dims-1?'}':',');
      }
      fprintf(stdout,"\n");
    }
  }
# endif
}
#endif

#if USE_FILES
static int t_start(void *_ctx,const char *_name,const th_info *_ti,int _pli,
 int _nxblocks,int _nyblocks){
  trans_ctx *ctx;
  fprintf(stdout,"%s %i %i\n",_name,_nxblocks,_nyblocks);
  fflush(stdout);
  ctx=(trans_ctx *)_ctx;
  image_ctx_init(&ctx->img,_name,_nxblocks,_nyblocks);
  return EXIT_SUCCESS;
}

static void t_load_data(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  trans_ctx *ctx;
  ctx=(trans_ctx *)_ctx;
  if(_bi==0&&_bj==0){
    int           x;
    int           y;
    int           z;
    unsigned char buf[2*B_SZ];
    /* add the rows */
    for(y=0;y<ctx->img.nyblocks*B_SZ;y++){
      for(x=0;x<ctx->img.nxblocks*B_SZ-(2*B_SZ-1);x++){
        for(z=0;z<2*B_SZ;z++){
          buf[z]=_data[y*_stride+(x+z)];
        }
        trans_data_add(&ctx->td,buf);
      }
    }
    /* add the columns */
    for(y=0;y<ctx->img.nyblocks*B_SZ-(2*B_SZ-1);y++){
      for(x=0;x<ctx->img.nxblocks*B_SZ;x++){
        for(z=0;z<2*B_SZ;z++){
          buf[z]=_data[(y+z)*_stride+x];
        }
        trans_data_add(&ctx->td,buf);
      }
    }
  }
}

#define PADDING (0)

const block_func BLOCKS[]={
  t_load_data
};

const int NBLOCKS=sizeof(BLOCKS)/sizeof(*BLOCKS);
#endif

int main(int _argc,const char *_argv[]){
  trans_ctx     ctx[NUM_PROCS];
  const int    *f;
  int           i;
  double        r[2*B_SZ];
  const double *cov;
  (void)_argc;
  (void)_argv;
#if B_SZ==4
  f=OD_FILTER_PARAMS4;
#elif B_SZ==8
  f=OD_FILTER_PARAMS8;
#elif B_SZ==16
  f=OD_FILTER_PARAMS16;
#else
# error "Need filter params for this block size."
#endif
  for(i=0;i<NUM_PROCS;i++){
    trans_data_init(&ctx[i].td,2*B_SZ);
  }
  cov=r;
#if USE_FILES
  OD_OMP_SET_THREADS(NUM_PROCS);
  ne_apply_to_blocks(ctx,sizeof(*ctx),0x1,PADDING,t_start,NBLOCKS,BLOCKS,NULL,
   _argc,_argv);
  for(i=1;i<NUM_PROCS;i++){
    trans_data_combine(&ctx[0].td,&ctx[i].td);
  }
  trans_data_normalize(&ctx[0].td);
# if PRINT_COV
  trans_data_print(&ctx[0].td,stderr);
# endif
  fprintf(stdout,"original cg=%- 24.16G\n",coding_gain_1d(ctx[0].td.cov,f));
  trans_data_collapse(&ctx[0].td,1,r);
  fprintf(stdout,"collapse cg=%- 24.16G\n",coding_gain_1d_collapsed(r,f));
  trans_data_expand(&ctx[0].td,1,r);
  fprintf(stdout,"expanded cg=%- 24.16G\n",coding_gain_1d(ctx[0].td.cov,f));
#elif USE_AR95
  auto_regressive_collapsed(r,2*B_SZ,1,0.95);
#elif USE_SUBSET1
# if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  cov=SUBSET1_1D[B_SZ_LOG-OD_LOG_BSIZE0];
# else
#  error "Need auto-correlation matrix for subset1 for this block size."
# endif
#elif USE_SUBSET3
# if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  cov=SUBSET3_1D[B_SZ_LOG-OD_LOG_BSIZE0];
# else
#  error "Need auto-correlation matrix for subset3 for this block size."
# endif
#endif
#if CG_SEARCH
  coding_gain_search(cov);
#else
  fprintf(stdout,"cg=%-24.18G\n",coding_gain_1d_collapsed(cov,f));
#endif
  for(i=0;i<NUM_PROCS;i++){
    trans_data_clear(&ctx[i].td);
  }
  return EXIT_SUCCESS;
}
