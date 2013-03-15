#include <omp.h>
#include <stdlib.h>
#include "intra_fit_tools.h"
#include "stats_tools.h"
#include "trans_tools.h"

#define NUM_PROCS (4)
#define CG_SEARCH (0)
#define PRINT_COV (0)

#define USE_FILES   (0)
#define USE_AR95    (1)
#define USE_SUBSET1 (0)
#define USE_SUBSET3 (0)

static void coding_gain_search(const double _r[2*B_SZ],const int *_f){
#if CG_SEARCH
#if B_SZ==4
  {
    int    filt[4];
    int    p0;
    int    q0;
    int    s0;
    int    s1;
    double cg;
    double best_cg;
    best_cg=0;
    for(p0=-(1<<NE_BITS);p0<=(1<<NE_BITS);p0++){
      filt[2]=p0;
      for(q0=(1<<NE_BITS);q0>=-(1<<NE_BITS);q0--){
        filt[3]=q0;
        for(s0=(1<<NE_BITS);s0<=2*(1<<NE_BITS);s0++){
          filt[0]=s0;
          for(s1=(1<<NE_BITS);s1<=2*(1<<NE_BITS);s1++){
            filt[1]=s1;
            cg=coding_gain_1d_collapsed(_r,filt);
            if(cg>best_cg){
              best_cg=cg;
              fprintf(stdout,"%i %i %i %i %G\n",p0,q0,s0,s1,cg);
            }
          }
        }
      }
    }
  }
#endif
#else
  fprintf(stdout,"cg=%-24.18G\n",coding_gain_1d_collapsed(_r,_f));
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

int main(int _argc,const char *_argv[]){
  trans_ctx     ctx[NUM_PROCS];
  const int    *filter;
  int           i;
  double        r[2*B_SZ];
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
    trans_data_init(&ctx[i].td,2*B_SZ);
 }
  cov=r;
#if USE_FILES
  omp_set_num_threads(NUM_PROCS);
  ne_apply_to_blocks(ctx,sizeof(*ctx),0x1,PADDING,t_start,NBLOCKS,BLOCKS,NULL,
   _argc,_argv);
  for(i=1;i<NUM_PROCS;i++){
    trans_data_combine(&ctx[0].td,&ctx[i].td);
  }
  trans_data_normalize(&ctx[0].td);
#if PRINT_COV
  trans_data_print(&ctx[0].td,stderr);
#endif
  fprintf(stdout,"original cg=%- 24.16G\n",coding_gain_1d(ctx[0].td.cov,filter));
  trans_data_collapse(&ctx[0].td,1,r);
  fprintf(stdout,"collapse cg=%- 24.16G\n",coding_gain_1d_collapsed(r,filter));
  trans_data_expand(&ctx[0].td,1,r);
  fprintf(stdout,"expanded cg=%- 24.16G\n",coding_gain_1d(ctx[0].td.cov,filter));
#elif USE_AR95
  auto_regressive_collapsed(r,2*B_SZ,1,0.95);
#elif USE_SUBSET1
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  cov=SUBSET1_1D[B_SZ_LOG-OD_LOG_BSIZE0];
#else
# error "Need auto-correlation matrix for subset1 for this block size."
#endif
#elif USE_SUBSET3
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  cov=SUBSET3_1D[B_SZ_LOG-OD_LOG_BSIZE0];
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
