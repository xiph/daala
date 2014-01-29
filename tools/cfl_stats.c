#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include "od_defs.h"
#include "od_filter.h"
#include "od_intra.h"
#include "image_tools.h"
#include "stats_tools.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define PRINT_PROGRESS (0)
#define DEBUG_CFL_PRED (0)
#define OUTPUT_PLANES  (0)

#if B_SZ_LOG!=3
# error "Only B_SZ=8 currently supported."
#endif

typedef struct cfl_stats_ctx cfl_stats_ctx;

struct cfl_stats_ctx {
  int n;
  int curr_pli;
  image_data img[3];
  intra_stats gb_fdp[2];
  intra_stats gb_cfl[2];
  intra_stats st_fdp[2];
  intra_stats st_cfl[2];
};

static void cfl_stats_ctx_init(cfl_stats_ctx *_this){
  int i;
  _this->n=0;
  for (i=0;i<2;i++) {
    intra_stats_init(&_this->gb_fdp[i],B_SZ_LOG-1);
    intra_stats_init(&_this->gb_cfl[i],B_SZ_LOG-1);
    intra_stats_init(&_this->st_fdp[i],B_SZ_LOG-1);
    intra_stats_init(&_this->st_cfl[i],B_SZ_LOG-1);
  }
}

static void cfl_stats_ctx_clear(cfl_stats_ctx *_this){
  int i;
  for (i=0;i<2;i++) {
    intra_stats_clear(&_this->gb_fdp[i]);
    intra_stats_clear(&_this->gb_cfl[i]);
    intra_stats_clear(&_this->st_fdp[i]);
    intra_stats_clear(&_this->st_cfl[i]);
  }
}

static void cfl_stats_ctx_combine(cfl_stats_ctx *_a,cfl_stats_ctx *_b) {
  int i;
  if (_b->n==0) {
    return;
  }
  for (i=0;i<2;i++) {
    intra_stats_combine(&_a->gb_fdp[i],&_b->gb_fdp[i]);
    intra_stats_combine(&_a->gb_cfl[i],&_b->gb_cfl[i]);
  }
  _a->n+=_b->n;
}

static void od_pre_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  cfl_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_pre_block\n");
  }
#endif
  ctx=(cfl_stats_ctx *)_ctx;
  image_data_pre_block(&ctx->img[ctx->curr_pli],_data,_stride,_bi,_bj);
}

static void od_fdct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_fdct_block\n");
  }
#endif
  ctx=(cfl_stats_ctx *)_ctx;
  image_data_fdct_block(&ctx->img[ctx->curr_pli],_bi,_bj);
}

#if TF_BLOCKS
static void od_tf_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_tf_block\n");
  }
#endif
  ctx=(cfl_stats_ctx *)_ctx;
  image_data_tf_block(&ctx->img[ctx->curr_pli],_bi,_bj);
}
#endif

static void od_mode_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
  ctx=(cfl_stats_ctx *)_ctx;
  if (ctx->curr_pli==0) {
    od_coeff block[5*B_SZ_MAX*B_SZ_MAX];
#if PRINT_PROGRESS
    if(_bi==0&&_bj==0){
      fprintf(stdout,"in od_mode_block\n");
    }
#endif
    image_data_load_block(&ctx->img[0],_bi,_bj,block);
    ctx->img[0].mode[ctx->img[0].nxblocks*_bj+_bi]=
     od_select_mode_satd(block,NULL,ctx->img[0].b_sz_log);
  }
}

static void od_fdp_pred_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
  ctx=(cfl_stats_ctx *)_ctx;
  if (ctx->curr_pli>0) {
#if PRINT_PROGRESS
    if(_bi==0&&_bj==0){
      fprintf(stdout,"in od_fdp_pred_block\n");
    }
#endif
    ctx->img[ctx->curr_pli].mode[ctx->img[ctx->curr_pli].nxblocks*_bj+_bi]=
     ctx->img[0].mode[ctx->img[0].nxblocks*_bj+_bi];
    image_data_pred_block(&ctx->img[ctx->curr_pli],_bi,_bj);
  }
}

static void od_fdp_stats_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
  ctx=(cfl_stats_ctx *)_ctx;
  if (ctx->curr_pli>0) {
#if PRINT_PROGRESS
    if(_bi==0&&_bj==0){
      fprintf(stdout,"in od_fdp_stats_block\n");
    }
#endif
    image_data_stats_block(&ctx->img[ctx->curr_pli],_data,_stride,_bi,_bj,
     &ctx->st_fdp[ctx->curr_pli-1]);
  }
}

static void od_cfl_pred_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
  ctx=(cfl_stats_ctx *)_ctx;
  if (ctx->curr_pli>0) {
    image_data   *img;
    image_data   *img0;
    int           b_sz;
    int           mode;
    od_coeff      p[B_SZ_MAX*B_SZ_MAX];
    od_coeff      c[2*B_SZ_MAX*2*B_SZ_MAX];
    od_coeff      l[2*B_SZ_MAX*2*B_SZ_MAX];
    unsigned char bsize[8*8];
    int           bx;
    int           by;
    int           i;
    int           j;
    double       *pred;
#if PRINT_PROGRESS
    if(_bi==0&&_bj==0){
      fprintf(stdout,"in od_cfl_pred_block\n");
    }
#endif
    img=&ctx->img[ctx->curr_pli];
    img0=&ctx->img[0];
    b_sz=1<<img->b_sz_log;
    for (i=0;i<8*8;i++) {
      bsize[i]=B_SZ_LOG-OD_LOG_BSIZE0;
    }
    mode=img->mode[img->nxblocks*_bj+_bi]=
     ctx->img[0].mode[ctx->img[0].nxblocks*_bj+_bi];
    for (by=0;by<2;by++) {
      for (bx=0;bx<2;bx++) {
        od_coeff *fdct;
        fdct=&img->fdct[img->fdct_stride*b_sz*(_bj+by)+b_sz*(_bi+bx)];
        for (j=0;j<b_sz;j++) {
          for (i=0;i<b_sz;i++) {
            c[(by*b_sz+j)*2*b_sz+bx*b_sz+i]=fdct[img->fdct_stride*j+i];
          }
        }
        fdct=&img0->fdct[img0->fdct_stride*2*b_sz*(_bj+by)+2*b_sz*(_bi+bx)];
        od_resample_luma_coeffs(&l[by*b_sz*2*b_sz+bx*b_sz],2*b_sz,fdct,img0->fdct_stride,1,1,0,1);
      }
    }
    od_chroma_pred(p,c,l,2*b_sz,1,1,
     img->b_sz_log-OD_LOG_BSIZE0,1,1,bsize,8,OD_INTRA_CHROMA_WEIGHTS_Q8[mode]);
    pred=&img->pred[img->pred_stride*b_sz*_bj+b_sz*_bi];
    for (j=0;j<b_sz;j++) {
      for (i=0;i<b_sz;i++) {
        pred[img->pred_stride*j+i]=p[b_sz*j+i];
      }
    }
#if DEBUG_CFL_PRED
    if (ctx->curr_pli>0) {
      if (_bi==0&&_bj==0) {
        od_coeff *block;
        fprintf(stderr,"%i %i\n",_bi,_bj);
        fprintf(stderr,"luma (time)\n");
        block=&img0->pre[img0->pre_stride*b_sz*(2*_bj+1)+b_sz*(2*_bi+1)];
        for (j=0;j<4*b_sz;j++) {
          if (j==2*b_sz) {
            for (i=0;i<4*b_sz;i++) {
              fprintf(stderr,"%s-----",i==2*b_sz?"+":"");
            }
            fprintf(stderr,"\n");
          }
          for (i=0;i<4*b_sz;i++) {
            fprintf(stderr,"%s%4i ",i==2*b_sz?"|":"",block[j*img0->pre_stride+i]);
          }
          fprintf(stderr,"\n");
        }
        fprintf(stderr,"luma (freq)\n");
        block=&img0->fdct[img0->fdct_stride*2*b_sz*_bj+2*b_sz*_bi];
        for (j=0;j<4*b_sz;j++) {
          if (j==2*b_sz) {
            for (i=0;i<4*b_sz;i++) {
              fprintf(stderr,"%s-----",i==2*b_sz?"+":"");
            }
            fprintf(stderr,"\n");
          }
          for (i=0;i<4*b_sz;i++) {
            fprintf(stderr,"%s%4i ",i==2*b_sz?"|":"",block[j*img0->fdct_stride+i]);
          }
          fprintf(stderr,"\n");
        }
        fprintf(stderr,"luma (resampled)\n");
        for (j=0;j<2*b_sz;j++) {
          if (j==b_sz) {
            for (i=0;i<2*b_sz;i++) {
              fprintf(stderr,"%s-----",i==b_sz?"+":"");
            }
            fprintf(stderr,"\n");
          }
          for (i=0;i<2*b_sz;i++) {
            fprintf(stderr,"%s%4i ",i==b_sz?"|":"",l[j*2*b_sz+i]);
          }
          fprintf(stderr,"\n");
        }
        fprintf(stderr,"chroma %i\n",ctx->curr_pli);
        for (j=0;j<2*b_sz;j++) {
          if (j==b_sz) {
            for (i=0;i<2*b_sz;i++) {
              fprintf(stderr,"%s-----",i==b_sz?"+":"");
            }
            fprintf(stderr,"\n");
          }
          for (i=0;i<2*b_sz;i++) {
            fprintf(stderr,"%s%4i ",i==b_sz?"|":"",c[j*2*b_sz+i]);
          }
          fprintf(stderr,"\n");
        }
        fprintf(stderr,"pred\n");
        for (j=0;j<b_sz;j++) {
          for (i=0;i<b_sz;i++) {
            fprintf(stderr,"%4i ",p[b_sz*j+i]);
          }
          fprintf(stderr,"\n");
        }
      }
    }
#endif
  }
}

static void od_cfl_stats_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  cfl_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
  ctx=(cfl_stats_ctx *)_ctx;
  if (ctx->curr_pli>0) {
#if PRINT_PROGRESS
    if(_bi==0&&_bj==0){
      fprintf(stdout,"in od_cfl_stats_block\n");
    }
#endif
    image_data_stats_block(&ctx->img[ctx->curr_pli],_data,_stride,_bi,_bj,
     &ctx->st_cfl[ctx->curr_pli-1]);
  }
}

static int stats_start(void *_ctx,const char *_name,
 const video_input_info *_info,int _pli,int _nxblocks,int _nyblocks){
  cfl_stats_ctx *ctx;
  ctx=(cfl_stats_ctx *)_ctx;
  ctx->curr_pli=_pli;
  image_data_init(&ctx->img[_pli],_name,B_SZ_LOG-(_pli!=0),_nxblocks,_nyblocks);
  if (_pli==0) {
    fprintf(stdout,"%s\n",_name);
    ctx->n++;
  }
  if (_pli>0) {
    intra_stats_reset(&ctx->st_fdp[_pli-1]);
    intra_stats_reset(&ctx->st_cfl[_pli-1]);
  }
  return EXIT_SUCCESS;
}

static int stats_finish(void *_ctx){
  int i;
  cfl_stats_ctx *ctx;
  ctx=(cfl_stats_ctx *)_ctx;
  if (ctx->curr_pli==2) {
    double *od_scale;
    od_scale=OD_SCALE[B_SZ_LOG-1-OD_LOG_BSIZE0];
    for (i=0;i<2;i++) {
      char label[128];
      sprintf(label,"Frequency-Domain Predictors (plane %i)",i+1);
      intra_stats_combine(&ctx->gb_fdp[i],&ctx->st_fdp[i]);
      intra_stats_correct(&ctx->st_fdp[i]);
      intra_stats_print(&ctx->st_fdp[i],label,od_scale);
      sprintf(label,"Chroma-from-Luma Predictors (plane %i)",i+1);
      intra_stats_combine(&ctx->gb_cfl[i],&ctx->st_cfl[i]);
      intra_stats_correct(&ctx->st_cfl[i]);
      intra_stats_print(&ctx->st_cfl[i],label,od_scale);
    }
    for (i=0;i<3;i++) {
      image_data_clear(&ctx->img[i]);
    }
  }
  return EXIT_SUCCESS;
}

const block_func BLOCKS[]={
  od_pre_block,
  od_fdct_block,
#if TF_BLOCKS
  od_tf_block,
#endif
  od_mode_block,
  od_fdp_pred_block,
  od_fdp_stats_block,
  od_cfl_pred_block,
  od_cfl_stats_block,
};

const int NBLOCKS=sizeof(BLOCKS)/sizeof(*BLOCKS);

/* There is a bug in get_intra_dims where chroma blocks do not line up under
    luma blocks. */
static void cfl_get_intra_dims(const video_input_info *_info,int _pli,
 int _padding,int *_x0,int *_y0,int *_nxblocks,int *_nyblocks){
  int xshift;
  int yshift;
  int padding;
  xshift=_pli!=0&&!(_info->pixel_fmt&1);
  yshift=_pli!=0&&!(_info->pixel_fmt&2);
  OD_ASSERT(xshift==yshift);
  padding=_padding<<B_SZ_LOG;
  /*An offset of 1 would be fine to provide enough context for VP8-style intra
     prediction, but for frequency-domain prediction, we'll want a full block,
     plus overlap.*/
  *_x0=(_info->pic_x>>xshift)+(padding>>1+xshift);
  *_y0=(_info->pic_y>>yshift)+(padding>>1+yshift);
  /*We take an extra block off the end to give enough context for above-right
     intra prediction.*/
  *_nxblocks=_info->pic_w-padding>>B_SZ_LOG;
  *_nyblocks=_info->pic_h-padding>>B_SZ_LOG;
}

static int cfl_apply_to_blocks(void *_ctx,int _ctx_sz,int _plmask,int _padding,
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
#if OUTPUT_PLANES
    {
      char *names[]={"y","cb","cr"};
      for (pli=0;pli<3;pli++) {
        int w;
        int h;
        int u;
        int v;
        int x0;
        int y0;
        const unsigned char *data;
        w=info.pic_w>>(pli!=0);
        h=info.pic_h>>(pli!=0);
        x0=info.pic_x>>(pli!=0);
        y0=info.pic_y>>(pli!=0);
        data=ycbcr[pli].data;
        fprintf(stderr,"%s=[",names[pli]);
        for (v=0;v<h;v++) {
          for (u=0;u<w;u++) {
            if (u>0) {
              fprintf(stderr,",");
            }
            fprintf(stderr,"%i",data[(v+y0)*w+u+x0]);
          }
          if (v<h) {
            fprintf(stderr,";");
          }
        }
        fprintf(stderr,"];\n");
      }
    }
#endif
    for(pli=0;pli<3;pli++){
      if(_plmask&1<<pli){
        int x0;
        int y0;
        int nxblocks;
        int nyblocks;
        cfl_get_intra_dims(&info,pli,_padding,&x0,&y0,&nxblocks,&nyblocks);
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
                y=y0+(bj<<B_SZ_LOG-(pli!=0));
                for(bi=0;bi<nxblocks;bi++){
                  int x;
                  x=x0+(bi<<B_SZ_LOG-(pli!=0));
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

#define PADDING (4)
#if PADDING<3
# error "PADDING must be at least 3 luma blocks"
#endif

int main(int _argc,const char *_argv[]) {
  cfl_stats_ctx ctx[NUM_PROCS];
  int i;
  ne_filter_params_init();
  for (i=0;i<OD_NBSIZES;i++) {
    od_scale_init(OD_SCALE[i],OD_LOG_BSIZE0+i);
  }
  for (i=0;i<NUM_PROCS;i++) {
    cfl_stats_ctx_init(&ctx[i]);
  }
  od_intra_init();
  OD_OMP_SET_THREADS(NUM_PROCS);
  cfl_apply_to_blocks(ctx,sizeof(*ctx),0x7,PADDING,stats_start,NBLOCKS,BLOCKS,
   stats_finish,_argc,_argv);
  for (i=1;i<NUM_PROCS;i++) {
    cfl_stats_ctx_combine(&ctx[0],&ctx[i]);
  }
  printf("Processed %i image(s)\n",ctx[0].n);
  if (ctx[0].n>0) {
    double *od_scale;
    od_scale=OD_SCALE[B_SZ_LOG-1-OD_LOG_BSIZE0];
    for (i=0;i<2;i++) {
      char label[128];
      sprintf(label,"Frequency-Domain Predictors (plane %i)",i+1);
      intra_stats_correct(&ctx[0].gb_fdp[i]);
      intra_stats_print(&ctx[0].gb_fdp[i],label,od_scale);
      sprintf(label,"Chroma-from-Luma Predictors (plane %i)",i+1);
      intra_stats_correct(&ctx[0].gb_cfl[i]);
      intra_stats_print(&ctx[0].gb_cfl[i],label,od_scale);
    }
  }
  for (i=0;i<NUM_PROCS;i++) {
    cfl_stats_ctx_clear(&ctx[i]);
  }
  od_intra_clear();
  return EXIT_SUCCESS;
}
