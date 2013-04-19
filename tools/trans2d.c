#include <omp.h>
#include <stdlib.h>
#include "intra_fit_tools.h"
#include "stats_tools.h"
#include "trans_tools.h"
#include "int_search.h"
#include "kiss99.h"

#define NUM_PROCS (4)
#define CG_SEARCH (0)
#define PRINT_COV (0)

#define RAMP_DYADIC (0)

#define USE_FILES   (1)
#define USE_AR95    (0)
#define USE_SUBSET1 (0)
#define USE_SUBSET3 (0)

static void coding_gain_search(const double _r[2*B_SZ*2*B_SZ],const int *_f){
#if CG_SEARCH
# if B_SZ==4
  {
    int    f[4];
    int    p0;
    int    q0;
    int    s0;
    int    s1;
    double cg;
    double best_cg;
    best_cg=0;
#if RAMP_DYADIC
    for(q0=(1<<NE_BITS);q0>=-(1<<NE_BITS);q0--){
      f[3]=q0;
      /* S0 = 4/1*(1-q0/64)
       * S0 >= 1 -> 64-q0 >= 16
       */
      s0=(1<<NE_BITS)-q0;
      if(s0>=(1<<NE_BITS-2)){
        f[0]=s0*4;
        for(p0=-(1<<NE_BITS);p0<=(1<<NE_BITS);p0++){
          f[2]=p0;
          /* S1 = 4/3*(1-(1-q0/64)*p0/64)
           * S1 >= 1   -> 64^2-(64-q0)*p0 >= 64*48
           * S1 = x/64 -> 64^2-(64-q0)*p0 = 0 MOD 48
           */
          s1=(1<<2*NE_BITS)-s0*p0;
          if(s1>=(1<<NE_BITS)*(3<<NE_BITS-2)&&s1%(3<<NE_BITS-2)==0){
            f[1]=s1/(3<<NE_BITS-2);
            cg=coding_gain_2d_collapsed(_r,f);
            if(cg>best_cg){
              best_cg=cg;
              printf("%i %i %i %i %G\n",f[2],f[3],f[0],f[1],cg);
            }
          }
        }
      }
    }
#else
    for(p0=-(1<<NE_BITS);p0<=(1<<NE_BITS);p0++){
      f[2]=p0;
      for(q0=(1<<NE_BITS);q0>=-(1<<NE_BITS);q0--){
        f[3]=q0;
        for(s0=(1<<NE_BITS);s0<=2*(1<<NE_BITS);s0++){
          f[0]=s0;
          for(s1=(1<<NE_BITS);s1<=2*(1<<NE_BITS);s1++){
            f[1]=s1;
            cg=coding_gain_2d_collapsed(_r,f);
            if(cg>best_cg){
              best_cg=cg;
              printf("%i %i %i %i %G\n",f[2],f[3],f[0],f[1],cg);
            }
          }
        }
      }
    }
#endif
  }
# elif B_SZ==8 && RAMP_DYADIC
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
    int    t1;
    int    t2;
    double cg;
    double best_cg;
    for(q0=(1<<NE_BITS);q0>=-(1<<NE_BITS);q0--){
      f[7]=q0;
      /* S0 = 8/1*(1-q0/64)
       * S0 >= 1 -> 64-q0 >= 8
       */
      s0=(1<<NE_BITS)-q0;
      if(s0>=(1<<NE_BITS-3)){
        f[0]=s0*8;
        for(p0=-(1<<NE_BITS);p0<=(1<<NE_BITS);p0++){
          f[4]=p0;
          for(q1=(1<<NE_BITS);q1>=-(1<<NE_BITS);q1--){
            f[8]=q1;
            /* S1 = 8/3*((1-q1/64)-(1-q0/64)*p0/64)
             * S1 >= 1   -> 64*t1-t0*p0 >= 64*24
             * S1 = x/64 -> 64*t1-t0*p0 = 0 MOD 24
             */
            t1=(1<<NE_BITS)-q1;
            s1=(1<<NE_BITS)*t1-s0*p0;
            if(s1>=(1<<NE_BITS)*(3<<NE_BITS-3)&&s1%(3<<NE_BITS-3)==0){
              f[1]=s1/(3<<NE_BITS-3);
              for(p1=-(1<<NE_BITS);p1<=(1<<NE_BITS);p1++){
                f[5]=p1;
                for(q2=(1<<NE_BITS);q2>=-(1<<NE_BITS);q2--){
                  f[9]=q2;
                  /* S2 = 8/5*((1-q2/64)-(1-q1/64)*p1/64)
                   * S2 >= 1   -> 64*t2-t1*p1) >= 64*40
                   * S2 = x/64 -> 64*t2-t1*p1 = 0 MOD 40
                   */
                  t2=(1<<NE_BITS)-q2;
                  s2=(1<<NE_BITS)*t2-t1*p1;
                  if(s2>=(1<<NE_BITS)*(5<<NE_BITS-3)&&s2%(5<<NE_BITS-3)==0){
                    f[2]=s2/(5<<NE_BITS-3);
                    for(p2=-(1<<NE_BITS);p2<=(1<<NE_BITS);p2++){
                      f[6]=p2;
                      /* S3 = 8/7*(1-(1-q2/64)*p2/64)
                       * S3 >= 1   -> 64^2-t2*p2 >= 64*56
                       * S3 = x/64 -> 64^2-t2*p2 = 0 MOD 56
                       */
                      s3=(1<<2*NE_BITS)-t2*p2;
                      if(s3>=(1<<NE_BITS)*(7<<NE_BITS-3)&&s3%(7<<NE_BITS-3)==0){
                        f[3]=s3/(7<<NE_BITS-3);
                        cg=coding_gain_2d_collapsed(_r,f);
                        if(cg>best_cg){
                          best_cg=cg;
                          printf("%i %i %i %i %i %i %i %i %i %i %-24.18G\n",
                          f[4],f[5],f[6],f[7],f[8],f[9],f[0],f[1],f[2],f[3],cg);
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
  }
# else
  {
    int         i;
    kiss99_ctx  ks[NUM_PROCS];
    int         lb[22];
    int         ub[22];
#   if B_SZ==4
    int dims=4;
#   elif B_SZ==8
    int dims=10;
#   elif B_SZ==16
    int dims=22;
#   endif
    OD_ASSERT(dims<=22);
    for(i=0;i<dims;i++){
      lb[i]=i<(B_SZ>>1)?(1<<NE_BITS):-(1<<NE_BITS);
      ub[i]=i<(B_SZ>>1)?2*(1<<NE_BITS):(1<<NE_BITS);
    }
    for(i=0;i<NUM_PROCS;i++){
      uint32_t srand;
      srand=i*16843009; /*Broadcast char to 4xchar*/
      kiss99_srand(&ks[i],(unsigned char*)&srand,sizeof(int));
    }
    #pragma omp parallel for schedule(dynamic)
    for(i=0;i<128;i++){
      int rsol[22];
      double cg;
      int j;
      int tid;
      tid=omp_get_thread_num();
      for(j=0;j<dims;j++){
        int range;
        int rng;
        range=ub[j]-lb[j];
        if(!range||i==0){
          rsol[j]=_f[j];
        }else{
          int mask;
          mask=(1<<OD_ILOG_NZ(range))-1;
          do{rng=((int)kiss99_rand(&ks[tid]))&mask;}while(rng>range);
          rsol[j]=lb[j]+rng;
        }
      }
      j=int_simplex_max(&cg,dims,coding_gain_2d_collapsed,_r,lb,ub,rsol);
      fprintf(stdout,"obj=%-24.18G steps=%4d params={",cg,j);
      for(j=0;j<dims;j++)fprintf(stdout,"%3d%c",rsol[j],j==dims-1?'}':',');
      fprintf(stdout,"\n");
    }
  }
# endif
#else
  fprintf(stdout,"cg=%-24.18G\n",coding_gain_2d_collapsed(_r,_f));
#endif
}

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
    int           y;
    int           x;
    int           j;
    int           i;
    unsigned char buf[2*B_SZ*2*B_SZ];
    for(y=0;y<ctx->img.nyblocks*B_SZ-(2*B_SZ-1);y++){
      for(x=0;x<ctx->img.nxblocks*B_SZ-(2*B_SZ-1);x++){
        for(j=0;j<2*B_SZ;j++){
          for(i=0;i<2*B_SZ;i++){
            buf[j*2*B_SZ+i]=_data[(y+j)*_stride+(x+i)];
          }
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

int main(int _argc,const char *_argv[]){
  trans_ctx     ctx[NUM_PROCS];
  const int    *filter;
  int           i;
  double        r[2*B_SZ*2*B_SZ];
  const double *cov;
#if B_SZ==4
  filter=OD_FILTER_PARAMS4;
#elif B_SZ==8
  filter=OD_FILTER_PARAMS8;
#elif B_SZ==16
  filter=OD_FILTER_PARAMS16;
#else
# error "Need filter params for this block size."
#endif
  for(i=0;i<NUM_PROCS;i++){
    trans_data_init(&ctx[i].td,2*B_SZ*2*B_SZ);
  }
  cov=r;
#if USE_FILES
  omp_set_num_threads(NUM_PROCS);
  ne_apply_to_blocks(&ctx,sizeof(*ctx),0x1,PADDING,t_start,NBLOCKS,BLOCKS,NULL,
   _argc,_argv);
  for(i=1;i<NUM_PROCS;i++){
    trans_data_combine(&ctx[0].td,&ctx[i].td);
  }
  trans_data_normalize(&ctx[0].td);
#if PRINT_COV
  trans_data_print(&ctx[0].td,stderr);
#endif
  fprintf(stdout,"original cg=%- 24.16G\n",coding_gain_2d(ctx[0].td.cov,filter));
  trans_data_collapse(&ctx[0].td,2*B_SZ,r);
  fprintf(stdout,"collapse cg=%- 24.16G\n",coding_gain_2d_collapsed(r,filter));
  trans_data_expand(&ctx[0].td,2*B_SZ,r);
  /*fprintf(stdout,"expanded cg=%- 24.16G\n",coding_gain_2d(ctx[1].td.cov,filter));*/
#elif USE_AR95
  auto_regressive_collapsed(r,2*B_SZ*2*B_SZ,2*B_SZ,0.95);
#elif USE_SUBSET1
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  cov=SUBSET1_2D[B_SZ_LOG-OD_LOG_BSIZE0];
#else
# error "Need auto-correlation matrix for subset1 for this block size."
#endif
#elif USE_SUBSET3
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  cov=SUBSET3_2D[B_SZ_LOG-OD_LOG_BSIZE0];
#else
# error "Need auto-correlation matrix for subset3 for this block size."
#endif
#endif
  coding_gain_search(cov,filter);
  for(i=0;i<NUM_PROCS;i++){
    trans_data_clear(&ctx[i].td);
  }
  return EXIT_SUCCESS;
}
