#include <stdlib.h>
#include <string.h>
#include "intra_fit_tools.h"
#include "stats_tools.h"
#include "image.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define WRITE_IMAGES (0)
#define PRINT_PROGRESS (0)

typedef struct intra_stats_image intra_stats_image;

struct intra_stats_image{
  intra_stats     vp8_stats;
  intra_stats     od_stats;
  image_data      img_data;
#if WRITE_IMAGES
  image_files     vp8_files;
  image_files     od_files;
#endif
};

static void intra_stats_image_init(intra_stats_image *_this,const char *_name,
 int _nxblocks,int _nyblocks){
  intra_stats_init(&_this->vp8_stats);
  intra_stats_init(&_this->od_stats);
  image_data_init(&_this->img_data,_name,_nxblocks,_nyblocks);
#if WRITE_IMAGES
  image_files_init(&_this->vp8_files,_nxblocks,_nyblocks);
  image_files_init(&_this->od_files,_nxblocks,_nyblocks);
#endif
}

static void intra_stats_image_clear(intra_stats_image *_this){
  image_data_clear(&_this->img_data);
#if WRITE_IMAGES
  image_files_clear(&_this->vp8_files);
  image_files_clear(&_this->od_files);
#endif
}

typedef struct intra_stats_global intra_stats_global;

struct intra_stats_global{
  intra_stats vp8_stats;
  intra_stats od_stats;
};

static void intra_stats_global_init(intra_stats_global *_this){
  intra_stats_init(&_this->vp8_stats);
  intra_stats_init(&_this->od_stats);
}

typedef struct intra_stats_ctx intra_stats_ctx;

struct intra_stats_ctx{
  int                n;
  intra_stats_global gb;
  intra_stats_image  img;
};

static void intra_stats_ctx_init(intra_stats_ctx *_this){
  _this->n=0;
  intra_stats_global_init(&_this->gb);
}

static int stats_start(void *_ctx,const char *_name,const th_info *_ti,int _pli,
 int _nxblocks,int _nyblocks){
  intra_stats_ctx *ctx;
  fprintf(stdout,"%s\n",_name);
  ctx=(intra_stats_ctx *)_ctx;
  ctx->n++;
  intra_stats_image_init(&ctx->img,_name,_nxblocks,_nyblocks);
  return EXIT_SUCCESS;
}

static void vp8_update_stats(intra_stats_image *_img,const unsigned char *_data,
 int _stride,int _bi,int _bj,int _mode){
  mode_data     *fr;
  mode_data     *md;
  int            i;
  int            j;
  od_coeff       buf[B_SZ*B_SZ];
  unsigned char  block[B_SZ*B_SZ];

#if WRITE_IMAGES
  od_rgba16_image_draw_point(&_img->vp8_files.map,_bi,_bj,COLORS[_mode]);
#endif

  fr=&_img->vp8_stats.fr;
  md=&_img->vp8_stats.md[_mode];

  fr->n++;
  md->n++;

  /* update the input mean and variance */
  mode_data_add_input(fr,_data,_stride);
  mode_data_add_input(md,_data,_stride);

  /* update the vp8 reference mean and covariance */
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      buf[B_SZ*i+j]=(_data[_stride*i+j]-128)*INPUT_SCALE;
    }
  }

#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,B_SZ,buf,B_SZ);
#else
# error "Need an fDCT implementation for this block size."
#endif

  mode_data_add_block(fr,buf,B_SZ,1);
  mode_data_add_block(md,buf,B_SZ,1);

  /* update the vp8 predicted mean and covariance */
  memset(block,0,B_SZ*B_SZ);
  vp8_intra_predict(block,B_SZ,_data,_stride,_mode);
#if WRITE_IMAGES
  image_draw_block(&_img->vp8_files.pred,B_SZ*_bi,B_SZ*_bj,block,B_SZ);
#endif
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      buf[B_SZ*i+j]=(_data[_stride*i+j]-block[B_SZ*i+j])*INPUT_SCALE;
#if WRITE_IMAGES
      block[B_SZ*i+j]=abs(_data[_stride*i+j]-block[B_SZ*i+j]);
#endif
    }
  }
#if WRITE_IMAGES
  image_draw_block(&_img->vp8_files.res,B_SZ*_bi,B_SZ*_bj,block,B_SZ);
#endif

#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,B_SZ,buf,B_SZ);
#else
# error "Need an fDCT implementation for this block size."
#endif

  mode_data_add_block(fr,buf,B_SZ,0);
  mode_data_add_block(md,buf,B_SZ,0);

  /* update the vp8 satd */
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      double satd;
      satd=abs(buf[B_SZ*j+i]);
      md->satd_avg[B_SZ*j+i]+=(satd-md->satd_avg[B_SZ*j+i])/md->n;
      fr->satd_avg[B_SZ*j+i]+=(satd-fr->satd_avg[B_SZ*j+i])/fr->n;
    }
  }
}

static void vp8_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  int              mode;
  ctx=(intra_stats_ctx *)_ctx;
#if WRITE_IMAGES
  image_draw_block(&ctx->img.vp8_files.raw,B_SZ*_bi,B_SZ*_bj,_data,_stride);
#endif
  mode=vp8_select_mode(_data,_stride,NULL);
  vp8_update_stats(&ctx->img,_data,_stride,_bi,_bj,mode);
}

static void od_pre_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_pre_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_pre_block(&ctx->img.img_data,_data,_stride,_bi,_bj);
}

static void od_fdct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_fdct_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_fdct_block(&ctx->img.img_data,_bi,_bj);
}

static void od_mode_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_mode_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_mode_block(&ctx->img.img_data,_bi,_bj);
}

static void od_pred_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_pred_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_pred_block(&ctx->img.img_data,_bi,_bj);
}

#if WRITE_IMAGES
static void od_idct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_idct_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_idct_block(&ctx->img.img_data,_bi,_bj);
}

static void od_post_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_post_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  image_data_post_block(&ctx->img.img_data,_bi,_bj);
}
#endif

static void od_stats_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  image_data      *img;
  intra_stats     *stats;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_stats_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  img=&ctx->img.img_data;
  stats=&ctx->img.od_stats;
  image_data_stats_block(img,_data,_stride,_bi,_bj,stats);
}

#if WRITE_IMAGES
static void od_image_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  image_data      *img;
  image_files     *files;
#if PRINT_PROGRESS 
  if(_bi==0&&_bj==0){
    printf("in od_image_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  img=&ctx->img.img_data;
  files=&ctx->img.od_files;
  image_data_files_block(img,_data,_stride,_bi,_bj,files);
}
#endif

static int stats_finish(void *_ctx){
  intra_stats_ctx *ctx;
#if WRITE_IMAGES
  char             name[8192];
  int              eos;
#endif
#if PRINT_PROGRESS
  printf("in stats_finish\n");
#endif
  ctx=(intra_stats_ctx *)_ctx;
#if WRITE_IMAGES
  strcpy(name,ctx->img.img_data.name);
  eos=strlen(name)-4;
  sprintf(&name[eos],"%s","-vp8");
  image_files_write(&ctx->img.vp8_files,name,NULL);
  sprintf(&name[eos],"%s","-daala");
  image_files_write(&ctx->img.od_files,name,NULL);
#endif
  intra_stats_combine(&ctx->gb.vp8_stats,&ctx->img.vp8_stats);
  intra_stats_correct(&ctx->img.vp8_stats);
  intra_stats_print(&ctx->img.vp8_stats,"VP8 Intra Predictors",VP8_SCALE);
  intra_stats_combine(&ctx->gb.od_stats,&ctx->img.od_stats);
  intra_stats_correct(&ctx->img.od_stats);
  intra_stats_print(&ctx->img.od_stats,"Daala Intra Predictors",OD_SCALE);
  intra_stats_image_clear(&ctx->img);
  return EXIT_SUCCESS;
}

const block_func BLOCKS[]={
  vp8_block,
  od_pre_block,
  od_fdct_block,
  od_mode_block,
  od_pred_block,
#if WRITE_IMAGES
  od_idct_block,
  od_post_block,
#endif
  od_stats_block,
#if WRITE_IMAGES
  od_image_block
#endif
};

const int NBLOCKS=sizeof(BLOCKS)/sizeof(*BLOCKS);

#define PADDING (4*B_SZ)
#if PADDING<3*B_SZ
# error "PADDING must be at least 3*B_SZ"
#endif

int main(int _argc,const char **_argv){
  intra_stats_ctx ctx;
  int ret;
#if WRITE_IMAGES
  intra_map_colors(COLORS,OD_INTRA_NMODES);
#endif
  vp8_scale_init(VP8_SCALE);
  od_scale_init(OD_SCALE);
  intra_map_colors(COLORS,OD_INTRA_NMODES);
  intra_stats_ctx_init(&ctx);
  ret=apply_to_blocks2(&ctx,PADDING,stats_start,BLOCKS,NBLOCKS,stats_finish,0x1,
   _argc,_argv);
  if(ret)return ret;
  printf("Processed %i image(s)\n",ctx.n);
  if(ctx.n>0){
    intra_stats_correct(&ctx.gb.vp8_stats);
    intra_stats_print(&ctx.gb.vp8_stats,"VP8 Intra Predictors",VP8_SCALE);
    intra_stats_correct(&ctx.gb.od_stats);
    intra_stats_print(&ctx.gb.od_stats,"Daala Intra Predictors",OD_SCALE);
  }
  return ret;
}
