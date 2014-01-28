#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "od_defs.h"
#include "od_filter.h"
#include "od_intra.h"
#include "image_tools.h"
#include "stats_tools.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define WRITE_IMAGES   (0)
#define PRINT_PROGRESS (0)
#define PRINT_BLOCKS   (0)

typedef struct intra_stats_ctx intra_stats_ctx;

struct intra_stats_ctx{
  int         n;
  intra_stats gb_vp8;
  intra_stats gb_od;
  intra_stats st_vp8;
  intra_stats st_od;
  image_data  img;
#if WRITE_IMAGES
  image_files files_vp8;
  image_files files_od;
#endif
};

static void intra_stats_ctx_init(intra_stats_ctx *_this){
  _this->n=0;
  intra_stats_init(&_this->gb_vp8,B_SZ_LOG);
  intra_stats_init(&_this->gb_od,B_SZ_LOG);
  intra_stats_init(&_this->st_vp8,B_SZ_LOG);
  intra_stats_init(&_this->st_od,B_SZ_LOG);
}

static void intra_stats_ctx_clear(intra_stats_ctx *_this){
  intra_stats_clear(&_this->gb_vp8);
  intra_stats_clear(&_this->gb_od);
  intra_stats_clear(&_this->st_vp8);
  intra_stats_clear(&_this->st_od);
}

static void intra_stats_ctx_set_image(intra_stats_ctx *_this,const char *_name,
 int _nxblocks,int _nyblocks){
  _this->n++;
  intra_stats_reset(&_this->st_vp8);
  intra_stats_reset(&_this->st_od);
  image_data_init(&_this->img,_name,B_SZ_LOG,_nxblocks,_nyblocks);
#if WRITE_IMAGES
  image_files_init(&_this->files_vp8,_nxblocks,_nyblocks);
  image_files_init(&_this->files_od,_nxblocks,_nyblocks);
#endif
}

static void intra_stats_ctx_clear_image(intra_stats_ctx *_this){
  image_data_clear(&_this->img);
#if WRITE_IMAGES
  image_files_clear(&_this->files_vp8);
  image_files_clear(&_this->files_od);
#endif
}

static void intra_stats_ctx_combine(intra_stats_ctx *_a,intra_stats_ctx *_b){
  if(_b->n==0){
    return;
  }
  intra_stats_combine(&_a->gb_vp8,&_b->gb_vp8);
  intra_stats_combine(&_a->gb_od,&_b->gb_od);
  _a->n+=_b->n;
}

static void vp8_stats_block(intra_stats_ctx *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj,int _mode,const unsigned char *_pred){
  int      j;
  int      i;
  od_coeff ref[B_SZ*B_SZ];
  od_coeff buf[B_SZ*B_SZ];
  double   res[B_SZ*B_SZ];
  (void)_bi;
  (void)_bj;
  /* Compute reference transform coefficients. */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      ref[B_SZ*j+i]=(_data[_stride*j+i]-128)*INPUT_SCALE;
    }
  }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(ref,B_SZ,ref,B_SZ);
#else
# error "Need an fDCT implementation for this block size."
#endif
  /* Compute residual transform coefficients. */
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=(_data[_stride*j+i]-_pred[B_SZ*j+i])*INPUT_SCALE;
    }
  }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,B_SZ,buf,B_SZ);
#else
# error "Need an fDCT implementation for this block size."
#endif
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      res[B_SZ*j+i]=buf[B_SZ*j+i];
    }
  }
  intra_stats_update(&_ctx->st_vp8,_data,_stride,_mode,ref,B_SZ,res,B_SZ);
}

#if WRITE_IMAGES
static void vp8_files_block(intra_stats_ctx *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj,int _mode,const unsigned char *_pred){
  int           j;
  int           i;
  unsigned char res[B_SZ*B_SZ];
  image_draw_block(&_ctx->files_vp8.raw,B_SZ*_bi,B_SZ*_bj,_data,_stride);
  od_rgba16_image_draw_point(&_ctx->files_vp8.map,_bi,_bj,COLORS[_mode]);
  image_draw_block(&_ctx->files_vp8.pred,B_SZ*_bi,B_SZ*_bj,_pred,B_SZ);
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      res[j*B_SZ+i]=OD_CLAMP255(_data[_stride*j+i]-_pred[B_SZ*j+i]+128);
    }
  }
  image_draw_block(&_ctx->files_vp8.res,B_SZ*_bi,B_SZ*_bj,res,B_SZ);
}
#endif

#if MASK_BLOCKS
static void od_mask_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    print_progress(stdout,"od_mask_block");
  }
#endif
  if(_bi==0&&_bj==0){
    intra_stats_ctx *ctx;
    ctx=(intra_stats_ctx *)_ctx;
    image_data_mask(&ctx->img,_data,_stride);
  }
}
#endif

static void vp8_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  int              mode;
  unsigned char    pred[B_SZ*B_SZ];
  ctx=(intra_stats_ctx *)_ctx;
  mode=vp8_select_mode(_data,_stride,NULL);
  memset(pred,0,B_SZ*B_SZ);
  vp8_intra_predict(pred,B_SZ,_data,_stride,mode);
#if MASK_BLOCKS
  if(!ctx->img.mask[ctx->img.nxblocks*_bj+_bi]){
    int j;
    int i;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        pred[j*B_SZ+i]=_data[_stride*j+i];
      }
    }
  }
#endif
#if WRITE_IMAGES
  vp8_files_block(ctx,_data,_stride,_bi,_bj,mode,pred);
#endif
  ctx->img.mode[ctx->img.nxblocks*_bj+_bi]=mode;
#if MASK_BLOCKS
  if(!ctx->img.mask[ctx->img.nxblocks*_bj+_bi]){
    return;
  }
#endif
  vp8_stats_block(ctx,_data,_stride,_bi,_bj,mode,pred);
}

static void od_pre_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_pre_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_pre_block(&ctx->img,_data,_stride,_bi,_bj);
}

static void od_fdct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_fdct_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_fdct_block(&ctx->img,_bi,_bj);
}

#if TF_BLOCKS
static void od_tf_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_tf_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_tf_block(&ctx->img,_bi,_bj);
}
#endif

static void od_mode_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  od_coeff         block[5*B_SZ*B_SZ];
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_mode_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_load_block(&ctx->img,_bi,_bj,block);
  ctx->img.mode[ctx->img.nxblocks*_bj+_bi]=
   od_select_mode_satd(block,NULL,ctx->img.b_sz_log);
}

#if PRINT_BLOCKS
static void od_print_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    print_progress(stdout,"od_print_block");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_print_block(&ctx->img,_bi,_bj,stderr);
}
#endif

static void od_pred_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  (void)_data;
  (void)_stride;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_pred_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_pred_block(&ctx->img,_bi,_bj);
}

#if WRITE_IMAGES
static void od_idct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_idct_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_idct_block(&ctx->img,_bi,_bj);
}

static void od_post_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_post_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_post_block(&ctx->img,_bi,_bj);
}
#endif

static void od_stats_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_stats_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_stats_block(&ctx->img,_data,_stride,_bi,_bj,&ctx->st_od);
}

#if WRITE_IMAGES
static void od_image_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    fprintf(stdout,"in od_files_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_files_block(&ctx->img,_data,_stride,_bi,_bj,&ctx->files_od);
}
#endif

static int stats_start(void *_ctx,const char *_name,
 const video_input_info *_info,int _pli, int _nxblocks,int _nyblocks){
  intra_stats_ctx *ctx;
  (void)_info;
  (void)_pli;
  fprintf(stdout,"%s\n",_name);
  ctx=(intra_stats_ctx *)_ctx;
  intra_stats_ctx_set_image(ctx,_name,_nxblocks,_nyblocks);
  return EXIT_SUCCESS;
}

static int stats_finish(void *_ctx){
  intra_stats_ctx *ctx;
#if WRITE_IMAGES
  char             name[8192];
  int              eos;
#endif
  ctx=(intra_stats_ctx *)_ctx;
#if WRITE_IMAGES
  strcpy(name,ctx->img.name);
  eos=strlen(name)-4;
  sprintf(&name[eos],"%s","-vp8");
  image_files_write(&ctx->files_vp8,name,NULL);
  sprintf(&name[eos],"%s","-daala");
  image_files_write(&ctx->files_od,name,NULL);
#endif
  intra_stats_combine(&ctx->gb_vp8,&ctx->st_vp8);
  intra_stats_correct(&ctx->st_vp8);
  intra_stats_print(&ctx->st_vp8,"VP8 Intra Predictors",VP8_SCALE[B_SZ_LOG-OD_LOG_BSIZE0]);
  intra_stats_combine(&ctx->gb_od,&ctx->st_od);
  intra_stats_correct(&ctx->st_od);
  intra_stats_print(&ctx->st_od,"Daala Intra Predictors",OD_SCALE[B_SZ_LOG-OD_LOG_BSIZE0]);
  intra_stats_ctx_clear_image(ctx);
  return EXIT_SUCCESS;
}

const block_func BLOCKS[]={
#if MASK_BLOCKS
  od_mask_block,
#endif
  vp8_block,
  od_pre_block,
  od_fdct_block,
#if TF_BLOCKS
  od_tf_block,
#endif
#if PRINT_BLOCKS
  od_print_block,
#endif
  od_mode_block,
  od_pred_block,
#if WRITE_IMAGES
  od_idct_block,
  od_post_block,
#endif
  od_stats_block,
#if WRITE_IMAGES
  od_image_block,
#endif
};

const int NBLOCKS=sizeof(BLOCKS)/sizeof(*BLOCKS);

#define PADDING (4*B_SZ)
#if PADDING<3*B_SZ
# error "PADDING must be at least 3*B_SZ"
#endif

int main(int _argc,const char *_argv[]){
  intra_stats_ctx ctx[NUM_PROCS];
  int             i;
  ne_filter_params_init();
  vp8_scale_init(VP8_SCALE[B_SZ_LOG-OD_LOG_BSIZE0],B_SZ_LOG);
  od_scale_init(OD_SCALE[B_SZ_LOG-OD_LOG_BSIZE0],B_SZ_LOG);
#if WRITE_IMAGES
  intra_map_colors(COLORS,OD_INTRA_NMODES);
#endif
  for(i=0;i<NUM_PROCS;i++){
    intra_stats_ctx_init(&ctx[i]);
  }
  od_intra_init();
  OD_OMP_SET_THREADS(NUM_PROCS);
  ne_apply_to_blocks(ctx,sizeof(*ctx),0x1,PADDING,stats_start,NBLOCKS,BLOCKS,
   stats_finish,_argc,_argv);
  for(i=1;i<NUM_PROCS;i++){
    intra_stats_ctx_combine(&ctx[0],&ctx[i]);
  }
  printf("Processed %i image(s)\n",ctx[0].n);
  if(ctx[0].n>0){
    intra_stats_correct(&ctx[0].gb_vp8);
    intra_stats_print(&ctx[0].gb_vp8,"VP8 Intra Predictors",
     VP8_SCALE[B_SZ_LOG-OD_LOG_BSIZE0]);
    intra_stats_correct(&ctx[0].gb_od);
    intra_stats_print(&ctx[0].gb_od,"Daala Intra Predictors",
     OD_SCALE[B_SZ_LOG-OD_LOG_BSIZE0]);
  }
  for(i=0;i<NUM_PROCS;i++){
    intra_stats_ctx_clear(&ctx[i]);
  }
  od_intra_clear();
  return EXIT_SUCCESS;
}
