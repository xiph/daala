#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "cholesky.h"
#include "intra_fit_tools.h"
#include "stats_tools.h"
#include "image.h"
#include "svd.h"

#define PRINT_PROGRESS (0)
#define PRINT_BLOCKS (0)
#define WRITE_IMAGES (0)
#define PRINT_BETAS (0)
#define WEIGHT_BLOCKS (0)
#define BITS_SELECT (0)
#define USE_SVD (0)

typedef struct pred_data pred_data;

struct pred_data{
  double w;
  double mean[5*B_SZ*B_SZ];
  double cov[5*B_SZ*B_SZ][5*B_SZ*B_SZ];
  int    mask[B_SZ*B_SZ][4*B_SZ*B_SZ];
};

static void pred_data_init(pred_data *_this,int _reset_mask){
  int i;
  int j;
  _this->w=0;
  for(i=0;i<5*B_SZ*B_SZ;i++){
    _this->mean[i]=0;
    for(j=0;j<5*B_SZ*B_SZ;j++){
      _this->cov[i][j]=0;
    }
  }
  if(_reset_mask){
    for(j=0;j<B_SZ*B_SZ;j++){
      for(i=0;i<4*B_SZ*B_SZ;i++){
        _this->mask[j][i]=1;
      }
    }
  }
}

static void pred_data_add_block(pred_data *_this,const od_coeff *_block,
 int _stride,double _w){
  double delta[5*B_SZ*B_SZ];
  int    by;
  int    bx;
  int    j;
  int    i;
  _this->w+=_w;
  for(by=0;by<=1;by++){
    for(bx=0;bx<=2-by;bx++){
      int             bo;
      const od_coeff *block;
      bo=B_SZ*B_SZ*(3*by+bx);
      block=&_block[_stride*B_SZ*(by-1)+B_SZ*(bx-1)];
      for(j=0;j<B_SZ;j++){
        for(i=0;i<B_SZ;i++){
          delta[bo+B_SZ*j+i]=block[_stride*j+i]-_this->mean[bo+B_SZ*j+i];
          _this->mean[bo+B_SZ*j+i]+=delta[bo+B_SZ*j+i]*_w/_this->w;
        }
      }
    }
  }
  for(j=0;j<5*B_SZ*B_SZ;j++){
    for(i=0;i<5*B_SZ*B_SZ;i++){
      _this->cov[i][j]+=delta[i]*delta[j]*_w*(_this->w-_w)/_this->w;
    }
  }
}

typedef struct solve_ctx solve_ctx;

struct solve_ctx{
  double scale[5*B_SZ*B_SZ];
  double xtx[4*B_SZ*B_SZ][4*B_SZ*B_SZ];
  double xty[4*B_SZ*B_SZ][B_SZ*B_SZ];
  double beta_0[B_SZ*B_SZ];
  double beta_1[B_SZ*B_SZ][4*B_SZ*B_SZ];
};

/* solve for beta_0[_y] and beta_1[_y] */
static void pred_data_solve(pred_data *_this,solve_ctx *_ctx,int _y){
  int nmi;
  int mi[4*B_SZ*B_SZ];
  int i;
  int j;

  nmi=0;
  for(i=0;i<4*B_SZ*B_SZ;i++){
    if(_this->mask[_y][i]){
      mi[nmi]=i;
      nmi++;
    }
    _ctx->beta_1[_y][i]=0;
  }
#if USE_SVD
  {
    static double  xtx[2*4*B_SZ*B_SZ][4*B_SZ*B_SZ];
    static double *xtxp[2*4*B_SZ*B_SZ];
    static double  s[4*B_SZ*B_SZ];

    for(j=0;j<nmi;j++){
      for(i=0;i<nmi;i++){
        xtx[j][i]=_ctx->xtx[mi[j]][mi[i]];
      }
    }

    for(i=0;i<2*nmi;i++){
      xtxp[i]=xtx[i];
    }
    svd_pseudoinverse(xtxp,s,nmi,nmi);

    /* compute beta_1 = (X^T*X)^-1 * X^T*Y and beta_0 = Ym - Xm * beta_1 */
    _ctx->beta_0[_y]=_this->mean[4*B_SZ*B_SZ+_y];
    for(j=0;j<nmi;j++){
      _ctx->beta_1[_y][mi[j]]=0;
      for(i=0;i<nmi;i++){
        _ctx->beta_1[_y][mi[j]]+=xtx[j][i]*_ctx->xty[mi[i]][_y];
      }
      _ctx->beta_1[_y][mi[j]]*=_ctx->scale[4*B_SZ*B_SZ+_y]/_ctx->scale[mi[j]];
      _ctx->beta_0[_y]-=_this->mean[mi[j]]*_ctx->beta_1[_y][mi[j]];
    }
  }
#else
  {
    static double r[UT_SZ(4*B_SZ*B_SZ,4*B_SZ*B_SZ)];
    static int    pivot[4*B_SZ*B_SZ];
    int           rank;
    static double tau[4*B_SZ*B_SZ];
    static double b[4*B_SZ*B_SZ];
    static double work[4*B_SZ*B_SZ];

    for(j=0;j<nmi;j++){
      for(i=j;i<nmi;i++){
        r[UT_IDX(j,i,nmi)]=_ctx->xtx[mi[j]][mi[i]];
      }
      b[j]=_ctx->xty[mi[j]][_y];
    }

    rank=cholesky(r,pivot,DBL_EPSILON,nmi);
    chdecomp(r,tau,rank,nmi);
    chsolve(r,pivot,tau,b,b,work,rank,nmi);

    /* compute beta_1 = (X^T*X)^-1 * X^T*Y and beta_0 = Ym - Xm * beta_1 */
    _ctx->beta_0[_y]=_this->mean[4*B_SZ*B_SZ+_y];
    for(j=0;j<nmi;j++){
      _ctx->beta_1[_y][mi[j]]=
       _ctx->scale[4*B_SZ*B_SZ+_y]/_ctx->scale[mi[j]]*b[j];
      _ctx->beta_0[_y]-=_this->mean[mi[j]]*_ctx->beta_1[_y][mi[j]];
    }
  }
#endif
}

/* compute the coefficient prediction gain based on the mask */
static double pred_data_pg(pred_data *_this,solve_ctx *_ctx,int _y){
  double var;
  double err;
  int    i;
  pred_data_solve(_this,_ctx,_y);
  var=_this->cov[4*B_SZ*B_SZ+_y][4*B_SZ*B_SZ+_y];
  /* compute the covariance of the prediction error
      err = (Y-X*beta_1)^T * (Y-X*beta_1) = Y^T*Y - Y^T*X * beta_1 */
  err=var;
  for(i=0;i<4*B_SZ*B_SZ;i++){
    if(_this->mask[_y][i]){
      err-=_this->cov[4*B_SZ*B_SZ+_y][i]*_ctx->beta_1[_y][i];
    }
  }
  return(10*log10(var/err));
}

static int pred_data_delta_pg(pred_data *_this,solve_ctx *_ctx,int _y,
 double *_delta_pg){
  double pg;
  int    i;
  int    idx=-1;
  double min=UINT_MAX;
  pg=pred_data_pg(_this,_ctx,_y);
  for(i=0;i<4*B_SZ*B_SZ;i++){
    if(_this->mask[_y][i]){
      _this->mask[_y][i]=0;
      _delta_pg[i]=pg-pred_data_pg(_this,_ctx,_y);
      _this->mask[_y][i]=1;
      if(_delta_pg[i]<min){
        idx=i;
        min=_delta_pg[i];
      }
    }
  }
  return(idx);
}

#define PRINT_COMP (0)
#define PRINT_DELTA_PG (0)
#define PRINT_DROPS (0)

static void pred_data_update(pred_data *_this,int _mode,int _drop){
  int              i;
  int              j;
  static solve_ctx ctx;

  /* compute the scale factors */
  for(i=0;i<5*B_SZ*B_SZ;i++){
    ctx.scale[i]=sqrt(_this->cov[i][i]);
  }

  /* normalize X^T*X and X^T*Y */
  for(j=0;j<4*B_SZ*B_SZ;j++){
    for(i=0;i<4*B_SZ*B_SZ;i++){
      ctx.xtx[j][i]=_this->cov[j][i]/(ctx.scale[j]*ctx.scale[i]);
    }
    for(i=0;i<B_SZ*B_SZ;i++){
      ctx.xty[j][i]=
       _this->cov[j][4*B_SZ*B_SZ+i]/(ctx.scale[j]*ctx.scale[4*B_SZ*B_SZ+i]);
    }
  }

#if PRINT_COMP
  fprintf(stderr,"xtx=[");
  for(j=0;j<4*B_SZ*B_SZ;j++){
    fprintf(stderr,"%s",j!=0?";":"");
    for(i=0;i<4*B_SZ*B_SZ;i++){
      fprintf(stderr,"%s%- 24.18G",i!=0?",":"",ctx.xtx[j][i]);
    }
  }
  fprintf(stderr,"];\n");

  fprintf(stderr,"xty=[");
  for(j=0;j<4*B_SZ*B_SZ;j++){
    fprintf(stderr,"%s",j!=0?";":"");
    for(i=0;i<B_SZ*B_SZ;i++){
      fprintf(stderr,"%s%- 24.18G",i!=0?",":"",ctx.xty[j][i]);
    }
  }
  fprintf(stderr,"];\n");
#endif

  /* if we are dropping coefficients */
  if(_drop>0){
    double delta_pg[B_SZ*B_SZ][4*B_SZ*B_SZ];
    int idx[B_SZ*B_SZ];
    for(j=0;j<B_SZ*B_SZ;j++){
      idx[j]=pred_data_delta_pg(_this,&ctx,j,delta_pg[j]);
    }
#if PRINT_DELTA_PG
    for(i=0;i<4*B_SZ*B_SZ;i++){
      fprintf(stdout,"%s%- 24.18G",i>0?",":"",delta_pg[j][i]);
    }
    fprintf(stdout,"\n");
#endif
    while(_drop-->0){
      j=-1;
      for(i=0;i<B_SZ*B_SZ;i++){
        if(idx[i]!=-1){
          if(j==-1||delta_pg[i][idx[i]]<delta_pg[j][idx[j]]){
            j=i;
          }
        }
      }
#if PRINT_DROPS
      fprintf(stdout,"Dropping (%i,%i) with pg=%g\n",j,idx[j],delta_pg[j][idx[j]]);
#endif
      _this->mask[j][idx[j]]=0;
      idx[j]=pred_data_delta_pg(_this,&ctx,j,delta_pg[j]);
    }
  }

  for(j=0;j<B_SZ*B_SZ;j++){
    pred_data_solve(_this,&ctx,j);
  }

#if PRINT_COMP
  fprintf(stderr,"beta_1=[");
  for(j=0;j<4*B_SZ*B_SZ;j++){
    fprintf(stderr,"%s",j!=0?";":"");
    for(i=0;i<B_SZ*B_SZ;i++){
      fprintf(stderr,"%s%- 24.18G",i!=0?",":"",ctx.beta_1[i][j]);
    }
  }
  fprintf(stderr,"];\n");

  fprintf(stderr,"beta_0=[");
  for(i=0;i<B_SZ*B_SZ;i++){
    fprintf(stderr,"%s%- 24.18G",i!=0?",":"",ctx.beta_0[i]);
  }
  fprintf(stderr,"];\n");
#endif

  /* update the predictor constants */
  for(i=0;i<B_SZ*B_SZ;i++){
    int y;
    int x;
    y=i/B_SZ;
    x=i%B_SZ;
#if B_SZ==4
    NE_PRED_OFFSETS_4x4[_mode][y][x]=ctx.beta_0[i];
#elif B_SZ==8
    NE_PRED_OFFSETS_8x8[_mode][y][x]=ctx.beta_0[i];
#elif B_SZ==16
    NE_PRED_OFFSETS_16x16[_mode][y][x]=ctx.beta_0[i];
#else
# error "Need predictors for this block size."
#endif
    for(j=0;j<B_SZ*B_SZ;j++){
      int v;
      int u;
      v=j/B_SZ;
      u=j%B_SZ;
#if B_SZ==4
      NE_PRED_WEIGHTS_4x4[_mode][y][x][0*B_SZ+v][0*B_SZ+u]=ctx.beta_1[i][0*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_4x4[_mode][y][x][0*B_SZ+v][1*B_SZ+u]=ctx.beta_1[i][1*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_4x4[_mode][y][x][0*B_SZ+v][2*B_SZ+u]=ctx.beta_1[i][2*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_4x4[_mode][y][x][1*B_SZ+v][0*B_SZ+u]=ctx.beta_1[i][3*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_4x4[_mode][y][x][1*B_SZ+v][1*B_SZ+u]=0;
#elif B_SZ==8
      NE_PRED_WEIGHTS_8x8[_mode][y][x][0*B_SZ+v][0*B_SZ+u]=ctx.beta_1[i][0*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_8x8[_mode][y][x][0*B_SZ+v][1*B_SZ+u]=ctx.beta_1[i][1*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_8x8[_mode][y][x][0*B_SZ+v][2*B_SZ+u]=ctx.beta_1[i][2*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_8x8[_mode][y][x][1*B_SZ+v][0*B_SZ+u]=ctx.beta_1[i][3*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_8x8[_mode][y][x][1*B_SZ+v][1*B_SZ+u]=0;
#elif B_SZ==16
      NE_PRED_WEIGHTS_16x16[_mode][y][x][0*B_SZ+v][0*B_SZ+u]=ctx.beta_1[i][0*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_16x16[_mode][y][x][0*B_SZ+v][1*B_SZ+u]=ctx.beta_1[i][1*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_16x16[_mode][y][x][0*B_SZ+v][2*B_SZ+u]=ctx.beta_1[i][2*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_16x16[_mode][y][x][1*B_SZ+v][0*B_SZ+u]=ctx.beta_1[i][3*B_SZ*B_SZ+j];
      NE_PRED_WEIGHTS_16x16[_mode][y][x][1*B_SZ+v][1*B_SZ+u]=0;
#else
# error "Need predictors for this block size."
#endif
    }
  }
}

typedef struct intra_pred_ctx intra_pred_ctx;

struct intra_pred_ctx{
  int         step;
  intra_stats stats;
  intra_stats gb;
  pred_data   pd[OD_INTRA_NMODES];
  double      b[OD_INTRA_NMODES][B_SZ*B_SZ];
  image_data  img;
#if WRITE_IMAGES
  image_files files;
#endif
};

static void intra_pred_ctx_init(intra_pred_ctx *_this,int _reset_mask){
  int i;
  intra_stats_init(&_this->stats);
  intra_stats_init(&_this->gb);
  for(i=0;i<OD_INTRA_NMODES;i++){
    pred_data_init(&_this->pd[i],_reset_mask);
  }
}

static void intra_pred_ctx_update(intra_pred_ctx *_this,int _drop){
  int i;
  for(i=0;i<OD_INTRA_NMODES;i++){
    mode_data_params(&_this->gb.md[i],_this->b[i],OD_SCALE);
    pred_data_update(&_this->pd[i],i,_drop);
  }
}

static int init_start(void *_ctx,const char *_name,const th_info *_ti,int _pli,
 int _nxblocks,int _nyblocks){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  printf("in init_start\n");
#endif
  fprintf(stdout,"%s\n",_name);
  fflush(stdout);
  ctx=(intra_pred_ctx *)_ctx;
  image_data_init(&ctx->img,_name,_nxblocks,_nyblocks);
  return EXIT_SUCCESS;
}

static int init_finish(void *_ctx){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  printf("in init_finish\n");
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_save_map(&ctx->img);
  image_data_clear(&ctx->img);
  return EXIT_SUCCESS;
}

static void ip_pre_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_pre_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_pre_block(&ctx->img,_data,_stride,_bi,_bj);
}

static void ip_fdct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_fdct_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_fdct_block(&ctx->img,_bi,_bj);
}

/* select the initial mode based on the VP8 mode classification */
static void ip_init_mode(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_init_mode\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  ctx->img.mode[ctx->img.nxblocks*_bj+_bi]=
   vp8_select_mode(_data,_stride,&ctx->img.weight[ctx->img.nxblocks*_bj+_bi]);
#if !WEIGHT_BLOCKS
  ctx->img.weight[ctx->img.nxblocks*_bj+_bi]=1;
#endif
}

static void ip_print_blocks(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
  image_data     *img;
  od_coeff       *block;
  int             mode;
  int             by;
  int             bx;
  int             j;
  int             i;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_print_blocks\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  img=&ctx->img;
  mode=img->mode[_bj*img->nxblocks+_bi];
  fprintf(stderr,"%i",mode);
  for(by=0;by<=1;by++){
    for(bx=0;bx<=2-by;bx++){
      block=&img->fdct[img->fdct_stride*B_SZ*(_bj+by)+B_SZ*(_bi+bx)];
      for(j=0;j<B_SZ;j++){
        for(i=0;i<B_SZ;i++){
          fprintf(stderr," %i",block[img->fdct_stride*j+i]);
        }
      }
    }
  }
  fprintf(stderr,"\n");
}

static void ip_update_pred(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
  image_data     *img;
  od_coeff       *block;
  int             mode;
  double          weight;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_update_pred\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  img=&ctx->img;
  block=&img->fdct[img->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];

  mode=img->mode[img->nxblocks*_bj+_bi];
  weight=img->weight[img->nxblocks*_bj+_bi];

  pred_data_add_block(&ctx->pd[mode],block,img->fdct_stride,weight);
}

static int pred_start(void *_ctx,const char *_name,const th_info *_ti,int _pli,
 int _nxblocks,int _nyblocks){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  printf("in pred_start\n");
#endif
  fprintf(stdout,"%s\n",_name);
  ctx=(intra_pred_ctx *)_ctx;
  intra_stats_init(&ctx->stats);
  image_data_init(&ctx->img,_name,_nxblocks,_nyblocks);
#if WRITE_IMAGES
  image_files_init(&ctx->files,_nxblocks,_nyblocks);
#endif
  image_data_load_map(&ctx->img);
  return EXIT_SUCCESS;
}

static int pred_finish(void *_ctx){
  intra_pred_ctx *ctx;
#if WRITE_IMAGES
  char step[16];
#endif
#if PRINT_PROGRESS
  printf("in pred_finish\n");
#endif
  ctx=(intra_pred_ctx *)_ctx;
  intra_stats_combine(&ctx->gb,&ctx->stats);
  intra_stats_correct(&ctx->stats);
  intra_stats_print(&ctx->stats,"Daala Intra Predictors",OD_SCALE);
  image_data_save_map(&ctx->img);
  image_data_clear(&ctx->img);
#if WRITE_IMAGES
  sprintf(step,"-step%02i",ctx->step);
  image_files_write(&ctx->files,ctx->img.name,step);
  image_files_clear(&ctx->files);
#endif
  return EXIT_SUCCESS;
}

static void ip_pred_mode(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
  image_data     *img;
  od_coeff       *block;
  double         *weight;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_pred_mode\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  img=&ctx->img;
  block=&img->fdct[img->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];
  weight=&img->weight[img->nxblocks*_bj+_bi];
#if BITS_SELECT
  /* use satd on first step since we don't have daala statistics yet*/
  if(ctx->step==1){
    img->mode[img->nxblocks*_bj+_bi]=
     od_select_mode_satd(block,img->fdct_stride,weight);
  }
  else{
    img->mode[img->nxblocks*_bj+_bi]=
     od_select_mode_bits(block,img->fdct_stride,weight,ctx->b);
  }
#else
  img->mode[img->nxblocks*_bj+_bi]=
   od_select_mode_satd(block,img->fdct_stride,weight);
#endif
#if !WEIGHT_BLOCKS
  img->weight[img->nxblocks*_bj+_bi]=1;
#endif
}

static void ip_pred_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_pred_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_pred_block(&ctx->img,_bi,_bj);
}

static void ip_idct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_idct_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_idct_block(&ctx->img,_bi,_bj);
}

static void ip_post_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_post_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_post_block(&ctx->img,_bi,_bj);
}

static void ip_stats_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_stats_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_stats_block(&ctx->img,_data,_stride,_bi,_bj,&ctx->stats);
}

#if WRITE_IMAGES
static void ip_files_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_pred_ctx *ctx;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in ip_files_block\n");
  }
#endif
  ctx=(intra_pred_ctx *)_ctx;
  image_data_files_block(&ctx->img,_data,_stride,_bi,_bj,&ctx->files);
}
#endif

const block_func INIT[]={
  ip_pre_block,
  ip_fdct_block,
  ip_init_mode,
#if PRINT_BLOCKS
  ip_print_blocks,
#endif
  ip_update_pred
};

const int NINIT=sizeof(INIT)/sizeof(*INIT);

const block_func KMEANS[]={
  ip_pre_block,
  ip_fdct_block,
  ip_pred_mode,
  ip_update_pred,
  ip_pred_block,
  ip_stats_block,
  ip_idct_block,
  ip_post_block,
#if WRITE_IMAGES
  ip_files_block
#endif
};

const int NKMEANS=sizeof(KMEANS)/sizeof(*KMEANS);

#define PADDING (4*B_SZ)
#if PADDING<3*B_SZ
# error "PADDING must be at least 3*B_SZ"
#endif

#define INIT_STEPS (10)
#define DROP_STEPS (60)
#define DROPS_PER_STEP (16)

int main(int _argc,const char *_argv[]){
  static intra_pred_ctx ctx;
  int                   ret;
  int                   s;
  /* initialize some constants */
#if B_SZ==4
  ne_filter_params4_init(OD_FILTER_PARAMS4);
#elif B_SZ==8
  ne_filter_params8_init(OD_FILTER_PARAMS8);
#elif B_SZ==16
  ne_filter_params16_init(OD_FILTER_PARAMS16);
#else
# error "Need filter params for this block size."
#endif
  vp8_scale_init(VP8_SCALE);
  od_scale_init(OD_SCALE);
#if WRITE_IMAGES
  intra_map_colors(COLORS,OD_INTRA_NMODES);
#endif
  /* first pass across images uses VP8 SATD mode selection */
  intra_pred_ctx_init(&ctx,1);
  ret=apply_to_blocks2(&ctx,PADDING,init_start,INIT,NINIT,init_finish,0x1,
   _argc,_argv);
  /* each k-means step uses Daala SATD mode selection */
  for(s=1;s<=INIT_STEPS;s++){
    printf("Starting Step %02i (%i mults / block)\n",s,64*B_SZ*B_SZ);
    ctx.step=s;
    /* update the intra predictors model */
    intra_pred_ctx_update(&ctx,0);
    intra_pred_ctx_init(&ctx,0);
    /* reclassify based on the new model */
    ret=apply_to_blocks2(&ctx,PADDING,pred_start,KMEANS,NKMEANS,pred_finish,0x1,
     _argc,_argv);
    printf("Finished Step %02i\n",s);
    intra_stats_correct(&ctx.gb);
    intra_stats_print(&ctx.gb,"Daala Intra Predictors",OD_SCALE);
  }
  for(s=INIT_STEPS+1;s<=INIT_STEPS+DROP_STEPS;s++){
    printf("Starting Step %02i (%i mults / block)\n",s,64*B_SZ*B_SZ-DROPS_PER_STEP*(s-INIT_STEPS));
    ctx.step=s;
    /* update the intra predictors model */
    intra_pred_ctx_update(&ctx,DROPS_PER_STEP);
    intra_pred_ctx_init(&ctx,0);
    /* reclassify based on the new model */
    ret=apply_to_blocks2(&ctx,PADDING,pred_start,KMEANS,NKMEANS,pred_finish,0x1,
     _argc,_argv);
    printf("Finished Step %02i\n",s);
    intra_stats_correct(&ctx.gb);
    intra_stats_print(&ctx.gb,"Daala Intra Predictors",OD_SCALE);
  }
#if PRINT_BETAS
  print_betas(stderr);
#endif
  return ret;
}
