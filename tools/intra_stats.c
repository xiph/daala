#include <stdlib.h>
#include <string.h>
#include "intra_fit_tools.h"
#include "image.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define WRITE_IMAGES (0)

typedef struct mode_data mode_data;

struct mode_data{
  int    n;
  double satd_avg;
  double mean;
  double var;
  double ref_mean[B_SZ*B_SZ];
  double ref_cov[B_SZ*B_SZ][B_SZ*B_SZ];
  double pred_mean[B_SZ*B_SZ];
  double pred_cov[B_SZ*B_SZ][B_SZ*B_SZ];
};

static void mode_data_init(mode_data *_md){
  int i;
  int j;
  _md->n=0;
  _md->satd_avg=0;
  _md->mean=0;
  _md->var=0;
  for(i=0;i<B_SZ*B_SZ;i++){
    _md->ref_mean[i]=0;
    _md->pred_mean[i]=0;
    for(j=0;j<B_SZ*B_SZ;j++){
      _md->ref_cov[i][j]=0;
      _md->pred_cov[i][j]=0;
    }
  }
}

/* update the input mean and variance */
static void mode_data_add_input(mode_data *_md,const unsigned char *_data,
 int _stride){
  int n;
  int i;
  int j;
  n=(_md->n-1)*B_SZ*B_SZ;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      double delta;
      n++;
      delta=_data[_stride*i+j]-_md->mean;
      _md->mean+=delta/n;
      _md->var+=delta*delta*(n-1)/n;
    }
  }
}

static void add_block(double _mean[B_SZ*B_SZ],double _cov[B_SZ*B_SZ][B_SZ*B_SZ],
 int _n,const od_coeff *_block,int _stride){
  double delta[B_SZ*B_SZ];
  int    i;
  int    j;
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      delta[B_SZ*j+i]=_block[_stride*j+i]-_mean[B_SZ*j+i];
      _mean[B_SZ*j+i]+=delta[B_SZ*j+i]/_n;
    }
  }
  for(i=0;i<B_SZ*B_SZ;i++){
    for(j=0;j<B_SZ*B_SZ;j++){
      _cov[i][j]+=delta[i]*delta[j]*(_n-1)/_n;
    }
  }
}

static void mode_data_combine(mode_data *_a,const mode_data *_b){
  double s;
  double delta;
  int    i;
  int    j;
  double delta_ref[B_SZ*B_SZ];
  double delta_pred[B_SZ*B_SZ];
  s=((double)_b->n)/(_a->n+_b->n);
  delta=_b->mean-_a->mean;
  _a->mean+=delta*s;
  _a->satd_avg+=(_b->satd_avg-_a->satd_avg)*s;
  for(i=0;i<B_SZ*B_SZ;i++){
    delta_ref[i]=_b->ref_mean[i]-_a->ref_mean[i];
    _a->ref_mean[i]+=delta_ref[i]*s;
    delta_pred[i]=_b->pred_mean[i]-_a->pred_mean[i];
    _a->pred_mean[i]+=delta_pred[i]*s;
  }
  s*=_a->n;
  _a->var+=_b->var+delta*s;
  for(i=0;i<B_SZ*B_SZ;i++){
    for(j=0;j<B_SZ*B_SZ;j++){
      _a->ref_cov[i][j]+=_b->ref_cov[i][j]+delta_ref[i]*delta_ref[j]*s;
      _a->pred_cov[i][j]+=_b->pred_cov[i][j]+delta_pred[i]*delta_pred[j]*s;
    }
  }
  _a->n+=_b->n;
}

static void mode_data_correct(mode_data *_md){
  int i;
  int j;
  _md->var/=_md->n*B_SZ*B_SZ;
  for(i=0;i<B_SZ*B_SZ;i++){
    for(j=0;j<B_SZ*B_SZ;j++){
      _md->ref_cov[i][j]/=_md->n;
      _md->pred_cov[i][j]/=_md->n;
    }
  }
}

static void mode_data_print(mode_data *_md,const char *_label,double *_scale){
  double     cg_ref;
  double     cg_pred;
  double     pg;
  int        j;
  int        i;
  /* compute prediction gain 10*log_{10}(prod_{j}(ref_cov[j]/pred_cov[j])) */
  cg_ref=10*log10(_md->var);
  cg_pred=10*log10(_md->var);
  pg=0;
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      int x;
      x=B_SZ*j+i;
      cg_ref-=10*log10(_md->ref_cov[x][x]*_scale[j]*_scale[i])/(B_SZ*B_SZ);
      cg_pred-=10*log10(_md->pred_cov[x][x]*_scale[j]*_scale[i])/(B_SZ*B_SZ);
      pg+=10*log10(_md->ref_cov[x][x]/_md->pred_cov[x][x])/(B_SZ*B_SZ);
    }
  }
  printf("%s Blocks %5i Avg. SATD %G Mean %G Var %G CgRef %G CgPred %G Pg %G\n",
   _label,_md->n,_md->satd_avg,_md->mean,_md->var,cg_ref,cg_pred,pg);
}

typedef struct intra_stats intra_stats;

struct intra_stats{
  mode_data fr;
  mode_data md[OD_INTRA_NMODES];
};

static void intra_stats_init(intra_stats *_this){
  int mode;
  mode_data_init(&_this->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_init(&_this->md[mode]);
  }
}

static void intra_stats_correct(intra_stats *_this){
  int mode;
  mode_data_correct(&_this->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_correct(&_this->md[mode]);
  }
}

static void intra_stats_print(intra_stats *_this,const char *_label,
 double *_scale){
  int mode;
  printf("%s\n",_label);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    char label[16];
    sprintf(label,"Mode %i",mode);
    mode_data_print(&_this->md[mode],label,_scale);
  }
  mode_data_print(&_this->fr,"Pooled",_scale);
}

static void intra_stats_combine(intra_stats *_this,const intra_stats *_that){
  int mode;
  mode_data_combine(&_this->fr,&_that->fr);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    mode_data_combine(&_this->md[mode],&_that->md[mode]);
  }
}

typedef struct intra_stats_od intra_stats_od;

struct intra_stats_od{
  intra_stats  super;
  int         *mode;
  od_coeff    *pre;
  int          pre_stride;
  od_coeff    *fdct;
  int          fdct_stride;
  od_coeff    *pred;
  int          pred_stride;
  od_coeff    *idct;
  int          idct_stride;
  od_coeff    *post;
  int          post_stride;
};

static void intra_stats_od_init(intra_stats_od *_this,int _nxblocks,
 int _nyblocks){
  int w;
  int h;
  intra_stats_init(&_this->super);
  _this->mode=(int *)malloc(sizeof(*_this->mode)*_nxblocks*_nyblocks);
  memset(_this->mode,0,sizeof(*_this->mode)*_nxblocks*_nyblocks);
  w=B_SZ*(_nxblocks+3);
  h=B_SZ*(_nyblocks+3);
  _this->pre=(od_coeff *)malloc(sizeof(*_this->pre)*w*h);
  _this->pre_stride=w;
  w=B_SZ*(_nxblocks+2);
  h=B_SZ*(_nyblocks+2);
  _this->fdct=(od_coeff *)malloc(sizeof(*_this->fdct)*w*h);
  _this->fdct_stride=w;
  w=B_SZ*_nxblocks;
  h=B_SZ*_nyblocks;
  _this->pred=(od_coeff *)malloc(sizeof(*_this->pred)*w*h);
  _this->pred_stride=w;
  w=B_SZ*(_nxblocks+2);
  h=B_SZ*(_nyblocks+2);
  _this->idct=(od_coeff *)malloc(sizeof(*_this->idct)*w*h);
  _this->idct_stride=w;
  w=B_SZ*(_nxblocks+1);
  h=B_SZ*(_nyblocks+1);
  _this->post=(od_coeff *)malloc(sizeof(*_this->post)*w*h);
  _this->post_stride=w;
}

static void intra_stats_od_clear(intra_stats_od *_this){
  free(_this->mode);
  free(_this->pre);
  free(_this->fdct);
  free(_this->pred);
  free(_this->idct);
  free(_this->post);
}

typedef struct intra_stats_image intra_stats_image;

struct intra_stats_image{
  const char      *name;
  int              nxblocks;
  int              nyblocks;
  intra_stats      vp8_stats;
  intra_stats_od   od_stats;
#if WRITE_IMAGES
  od_rgba16_image  luma;
  od_rgba16_image  vp8_map;
  od_rgba16_image  vp8_pred;
  od_rgba16_image  vp8_res;
  od_rgba16_image  pre;
  od_rgba16_image  od_map;
  od_rgba16_image  od_pred;
  od_rgba16_image  od_res;
#endif
};

static void intra_stats_image_init(intra_stats_image *_this,const char *_name,
 int _nxblocks,int _nyblocks){
  _this->name=_name;
  _this->nxblocks=_nxblocks;
  _this->nyblocks=_nyblocks;
  intra_stats_init(&_this->vp8_stats);
  intra_stats_od_init(&_this->od_stats,_nxblocks,_nyblocks);
#if WRITE_IMAGES
  od_rgba16_image_init(&_this->luma,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->vp8_map,_nxblocks,_nyblocks);
  od_rgba16_image_init(&_this->vp8_pred,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->vp8_res,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->pre,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->od_map,_nxblocks,_nyblocks);
  od_rgba16_image_init(&_this->od_pred,B_SZ*_nxblocks,B_SZ*_nyblocks);
  od_rgba16_image_init(&_this->od_res,B_SZ*_nxblocks,B_SZ*_nyblocks);
#endif
}

static void intra_stats_image_clear(intra_stats_image *_this){
  intra_stats_od_clear(&_this->od_stats);
#if WRITE_IMAGES
  od_rgba16_image_clear(&_this->luma);
  od_rgba16_image_clear(&_this->vp8_map);
  od_rgba16_image_clear(&_this->vp8_pred);
  od_rgba16_image_clear(&_this->vp8_res);
  od_rgba16_image_clear(&_this->pre);
  od_rgba16_image_clear(&_this->od_map);
  od_rgba16_image_clear(&_this->od_pred);
  od_rgba16_image_clear(&_this->od_res);
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

#if WRITE_IMAGES
od_rgba16_pixel colors[OD_INTRA_NMODES];

static void image_draw_block(od_rgba16_image *_image,int _x,int _y,
 const unsigned char *_block,int _stride){
  od_rgba16_pixel color;
  int             i;
  int             j;
  color[3]=(unsigned short)0xFFFFU;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      color[0]=color[1]=color[2]=_block[_stride*i+j]<<8;
      od_rgba16_image_draw_point(_image,_x+j,_y+i,color);
    }
  }
}

static int image_write_png(od_rgba16_image *_image,const char *_name){
  char  fout_name[8192];
  FILE *fout;
  sprintf(fout_name,"%s.png",_name);
  fout=fopen(fout_name,"wb");
  if(fout==NULL){
    fprintf(stderr,"Could not open '%s' for reading.\n",fout_name);
    return EXIT_FAILURE;
  }
  od_rgba16_image_write_png(_image,fout);
  return EXIT_SUCCESS;
}
#endif

static int stats_start(void *_ctx,const char *_name,const th_info *_ti,int _pli,
 int _nxblocks,int _nyblocks){
  intra_stats_ctx *ctx;
  fprintf(stdout,"%s\n",_name);
  ctx=(intra_stats_ctx *)_ctx;
  ctx->n++;
  intra_stats_image_init(&ctx->img,_name,_nxblocks,_nyblocks);
  return EXIT_SUCCESS;
}

static void od_fdct_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
      (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(&_out[_out_stride*y+x],_out_stride,
       &_in[_in_stride*y+x],_in_stride);
#else
# error "Need an fDCT implementation for this block size."
#endif
    }
  }
}

#if WRITE_IMAGES
static void od_idct_blocks(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
      (*OD_IDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(&_out[_out_stride*y+x],_out_stride,
       &_in[_in_stride*y+x],_in_stride);
#else
# error "Need an iDCT implementation for this block size."
#endif
    }
  }
}
#endif

/* compute the scale factors for the DCT and TDLT transforms */
double vp8_scale[B_SZ];

double od_scale[B_SZ];

#define SCALE_BITS (14)
#define PRINT_SCALE (0)

static void vp8_scale_init(double _vp8_scale[B_SZ]){
  int j;
  int i;
  od_coeff buf[B_SZ];
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      buf[j]=i!=j?0:(1<<SCALE_BITS);
    }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_IDCT_1D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,1,buf);
#else
# error "Need an iDCT implementation for this block size."
#endif
    _vp8_scale[i]=0;
    for(j=0;j<B_SZ;j++){
      double c=((double)buf[j])/(1<<SCALE_BITS);
      _vp8_scale[i]+=c*c;
    }
#if PRINT_SCALE
    printf("%s%G",i==0?"":" ",_vp8_scale[i]);
#endif
  }
#if PRINT_SCALE
  printf("\n");
#endif
}

static void od_scale_init(double _od_scale[B_SZ]){
  int i;
  int j;
  od_coeff buf[2*B_SZ];
  for(i=0;i<B_SZ;i++){
    for(j=0;j<2*B_SZ;j++){
      buf[j]=(B_SZ>>1)+i!=j?0:(1<<SCALE_BITS);
    }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_IDCT_1D[B_SZ_LOG-OD_LOG_BSIZE0])(&buf[B_SZ>>1],1,&buf[B_SZ>>1]);
#else
# error "Need an iDCT implementation for this block size."
#endif
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(buf,buf);
    (*OD_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(&buf[B_SZ],&buf[B_SZ]);
#else
# error "Need a postfilter implementation for this block size."
#endif
    _od_scale[i]=0;
    for(j=0;j<2*B_SZ;j++){
      double c=((double)buf[j])/(1<<SCALE_BITS);
      _od_scale[i]+=c*c;
    }
#if PRINT_SCALE
    printf("%s%G",i==0?"":" ",_od_scale[i]);
#endif
  }
#if PRINT_SCALE
  printf("\n");
#endif
}

#define VP8_SCALE_SATD (1)

/* find the best vp8 mode */
static int vp8_select_mode(const unsigned char *_data,int _stride){
  double best_satd;
  int    mode;
  int    best_mode;
  best_mode=0;
  best_satd=UINT_MAX;
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

    od_fdct_blocks(buf,B_SZ,buf,B_SZ,1,1);

    satd=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
#if VP8_SCALE_SATD
        satd+=sqrt(vp8_scale[j]*vp8_scale[i])*abs(buf[B_SZ*j+i]);
#else
        satd+=abs(buf[B_SZ*j+i]);
#endif
      }
    }
    if(satd<best_satd){
      best_satd=satd;
      best_mode=mode;
    }
  }

  return best_mode;
}

static void vp8_update_stats(intra_stats_image *_img,const unsigned char *_data,
 int _stride,int _bi,int _bj,int _mode){
  mode_data     *fr;
  mode_data     *md;
  int            i;
  int            j;
  od_coeff       buf[B_SZ*B_SZ];
  unsigned char  block[B_SZ*B_SZ];
  double         satd;

#if WRITE_IMAGES
  od_rgba16_image_draw_point(&_img->vp8_map,_bi,_bj,colors[_mode]);
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
      buf[B_SZ*i+j]=_data[_stride*i+j]-128;
    }
  }

  od_fdct_blocks(buf,B_SZ,buf,B_SZ,1,1);

  add_block(fr->ref_mean,fr->ref_cov,fr->n,buf,B_SZ);
  add_block(md->ref_mean,md->ref_cov,md->n,buf,B_SZ);

  /* update the vp8 predicted mean and covariance */
  memset(block,0,B_SZ*B_SZ);
  vp8_intra_predict(block,B_SZ,_data,_stride,_mode);
#if WRITE_IMAGES
  image_draw_block(&_img->vp8_pred,B_SZ*_bi,B_SZ*_bj,block,B_SZ);
#endif
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      buf[B_SZ*i+j]=_data[_stride*i+j]-block[B_SZ*i+j];
#if WRITE_IMAGES
      block[B_SZ*i+j]=abs(_data[_stride*i+j]-block[B_SZ*i+j]);
#endif
    }
  }
#if WRITE_IMAGES
  image_draw_block(&_img->vp8_res,B_SZ*_bi,B_SZ*_bj,block,B_SZ);
#endif

  od_fdct_blocks(buf,B_SZ,buf,B_SZ,1,1);

  add_block(fr->pred_mean,fr->pred_cov,fr->n,buf,B_SZ);
  add_block(md->pred_mean,md->pred_cov,md->n,buf,B_SZ);

  /* update the vp8 satd */
  satd=0;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
#if VP8_SCALE_SATD
        satd+=sqrt(vp8_scale[j]*vp8_scale[i])*abs(buf[B_SZ*j+i]);
#else
        satd+=abs(buf[B_SZ*j+i]);
#endif
    }
  }

  md->satd_avg+=(satd-md->satd_avg)/md->n;
  fr->satd_avg+=(satd-fr->satd_avg)/fr->n;
}

#define APPLY_PREFILTER (1)
#define APPLY_POSTFILTER (1)

static void od_prefilter(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  int j;
  int i;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
      for(j=0;j<B_SZ;j++){
        od_coeff col[B_SZ];
#if APPLY_PREFILTER
        for(i=0;i<B_SZ;i++){
          col[i]=_in[_in_stride*(y+i)+x+j];
        }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*OD_PRE_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(col,col);
#else
# error "Need a prefilter implementation for this block size."
#endif
        for(i=0;i<B_SZ;i++){
          _out[_out_stride*(y+i)+x+j]=col[i];
        }
#else
        for(i=0;i<B_SZ;i++){
          _out[_out_stride*(y+i)+x+j]=_in[_in_stride*(y+i)+x+j];
        }
#endif
      }
#if APPLY_PREFILTER
      for(j=0;j<B_SZ;j++){
        od_coeff *row;
        row=&_out[_out_stride*(y+j)+x];
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*OD_PRE_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(row,row);
#else
# error "Need a prefilter implementation for this block size."
#endif
      }
#endif
    }
  }
}

#if WRITE_BLOCKS
static void od_postfilter(od_coeff *_out,int _out_stride,od_coeff *_in,
 int _in_stride,int _bx,int _by){
  int bx;
  int by;
  int j;
  int i;
  for(by=0;by<_by;by++){
    int y;
    y=B_SZ*by;
    for(bx=0;bx<_bx;bx++){
      int x;
      x=B_SZ*bx;
      for(j=0;j<B_SZ;j++){
        od_coeff *row;
        for(i=0;i<B_SZ;i++){
          _out[_out_stride*(y+j)+x+i]=_in[_in_stride*(y+j)+x+i];
        }
#if APPLY_POSTFILTER
        row=&_out[_out_stride*(y+j)+x];
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*OD_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(row,row);
#else
# error "Need a postfilter implementation for this block size."
#endif
#endif
      }
#if APPLY_POSTFILTER
      for(i=0;i<B_SZ;i++){
        od_coeff col[B_SZ];
        for(j=0;j<B_SZ;j++){
          col[j]=_out[_out_stride*(y+j)+x+i];
        }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*OD_POST_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(col,col);
#else
# error "Need a postfilter implementation for this block size."
#endif
        for(j=0;j<B_SZ;j++){
          _out[_out_stride*(j+y)+x+i]=col[j];
        }
      }
    }
  }
#endif
}
#endif

static void vp8_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  int              mode;
  ctx=(intra_stats_ctx *)_ctx;
#if WRITE_IMAGES
  image_draw_block(&ctx->img.luma,B_SZ*_bi,B_SZ*_bj,_data,_stride);
#endif
  mode=vp8_select_mode(_data,_stride);
  vp8_update_stats(&ctx->img,_data,_stride,_bi,_bj,mode);
}

#define PRINT_PROGRESS (0)

static void od_pre_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              x0;
  int              y0;
  int              bx;
  int              by;
  int              x;
  int              y;
  int              bi;
  int              bj;
  int              i;
  int              j;
  od_coeff         buf[B_SZ*B_SZ];
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_pre_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  x0=-(B_SZ>>1);
  y0=-(B_SZ>>1);
  bx=by=1;
  if(_bi==0){
    x0-=B_SZ;
    bx++;
  }
  if(_bj==0){
    y0-=B_SZ;
    by++;
  }
  if(_bi==ctx->img.nxblocks-1){
    bx+=2;
  }
  if(_bj==ctx->img.nyblocks-1){
    by+=2;
  }
  x=x0+B_SZ*_bi+(3*B_SZ>>1);
  y=y0+B_SZ*_bj+(3*B_SZ>>1);
  for(bj=0;bj<by;bj++){
    for(bi=0;bi<bx;bi++){
      for(j=0;j<B_SZ;j++){
        for(i=0;i<B_SZ;i++){
          buf[B_SZ*j+i]=_data[_stride*(y0+B_SZ*bj+j)+x0+B_SZ*bi+i]-128;
        }
      }
      od_prefilter(&od->pre[od->pre_stride*(y+B_SZ*bj)+x+B_SZ*bi],
       od->pre_stride,buf,B_SZ,1,1);
    }
  }
}

static void od_fdct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              x0;
  int              y0;
  int              bx;
  int              by;
  int              x;
  int              y;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_fdct_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  x0=_bi*B_SZ+(3*B_SZ>>1);
  y0=_bj*B_SZ+(3*B_SZ>>1);
  bx=by=1;
  if(_bi==0){
    x0-=B_SZ;
    bx++;
  }
  if(_bj==0){
    y0-=B_SZ;
    by++;
  }
  if(_bi==ctx->img.nxblocks-1){
    bx++;
  }
  if(_bj==ctx->img.nyblocks-1){
    by++;
  }
  x=x0-(B_SZ>>1);
  y=y0-(B_SZ>>1);
  od_fdct_blocks(&od->fdct[y*od->fdct_stride+x],od->fdct_stride,
   &od->pre[y0*od->pre_stride+x0],od->pre_stride,bx,by);
}

static void od_mode_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  od_coeff        *block;
  int              best_mode;
  double           best_satd;
  int              mode;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_mode_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  block=&od->fdct[od->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];
  best_mode=0;
  best_satd=UINT_MAX;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double    p[B_SZ*B_SZ];
    double    satd;
    int       i;
    int       j;
    od_intra_pred4x4_mult(block,od->fdct_stride,mode,p);
    satd=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        satd+=sqrt(od_scale[j]*od_scale[i])*
	 abs(block[od->fdct_stride*j+i]-(od_coeff)floor(p[B_SZ*j+i]+0.5));
      }
    }
    if(satd<best_satd){
      best_mode=mode;
      best_satd=satd;
    }
  }
  od->mode[_bj*ctx->img.nxblocks+_bi]=best_mode;
}

static void od_pred_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              mode;
  od_coeff        *block;
  double           p[B_SZ*B_SZ];
  od_coeff        *pred;
  int              j;
  int              i;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_pred_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  mode=od->mode[ctx->img.nxblocks*_bj+_bi];
  block=&od->fdct[od->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];
  od_intra_pred4x4_mult(block,od->fdct_stride,mode,p);
  pred=&od->pred[od->pred_stride*B_SZ*_bj+B_SZ*_bi];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      pred[od->pred_stride*j+i]=(od_coeff)floor(p[B_SZ*j+i]+0.5);
    }
  }
}

#if WRITE_IMAGES
static void od_idct_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              x0;
  int              y0;
  int              x;
  int              y;
  int              bx;
  int              by;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_idct_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  x0=B_SZ*_bi;
  y0=B_SZ*_bj;
  x=x0+B_SZ;
  y=y0+B_SZ;
  bx=by=1;
  if(_bi==0){
    x-=B_SZ;
    bx++;
  }
  if(_bj==0){
    y-=B_SZ;
    by++;
  }
  if(_bi==ctx->img.nxblocks-1){
    bx++;
  }
  if(_bj==ctx->img.nyblocks-1){
    by++;
  }
  if(bx!=1||by!=1){
    od_idct_blocks(&od->idct[od->idct_stride*y+x],od->idct_stride,
     &od->fdct[od->fdct_stride*y+x],od->fdct_stride,bx,by);
  }
  x=x0+B_SZ;
  y=y0+B_SZ;
  od_idct_blocks(&od->idct[od->idct_stride*y+x],od->idct_stride,
   &od->pred[od->pred_stride*y0+x0],od->pred_stride,1,1);
}

static void od_post_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              x0;
  int              y0;
  int              x;
  int              y;
  int              bx;
  int              by;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_post_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  x=B_SZ*_bi;
  y=B_SZ*_bj;
  x0=x+(B_SZ>>1);
  y0=y+(B_SZ>>1);
  bx=by=1;
  if(_bi==ctx->img.nxblocks-1){
    bx++;
  }
  if(_bj==ctx->img.nyblocks-1){
    by++;
  }
  od_postfilter(&od->post[od->post_stride*y+x],od->post_stride,
   &od->idct[od->idct_stride*y0+x0],od->idct_stride,bx,by);
}
#endif

static void od_stats_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              mode;
  mode_data       *fr;
  mode_data       *md;
  od_coeff        *pred;
  od_coeff        *block;
  int              j;
  int              i;
  od_coeff         buf[B_SZ*B_SZ];
  double           satd;
#if PRINT_PROGRESS
  if(_bi==0&&_bj==0){
    printf("in od_stats_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  mode=od->mode[_bj*ctx->img.nxblocks+_bi];
  fr=&od->super.fr;
  md=&od->super.md[mode];

  fr->n++;
  md->n++;

  /* update the input mean and variance */
  mode_data_add_input(fr,_data,_stride);
  mode_data_add_input(md,_data,_stride);

  /* update the daala reference mean and covariance */
  block=&od->fdct[od->fdct_stride*B_SZ*(_bj+1)+B_SZ*(_bi+1)];

  add_block(fr->ref_mean,fr->ref_cov,fr->n,block,od->fdct_stride);
  add_block(md->ref_mean,md->ref_cov,md->n,block,od->fdct_stride);
  
  /* update the daala predicted mean and covariance */
  pred=&od->pred[od->pred_stride*B_SZ*_bj+B_SZ*_bi];

  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=block[od->fdct_stride*j+i]-pred[od->pred_stride*j+i];
    }
  }

  add_block(fr->pred_mean,fr->pred_cov,fr->n,buf,B_SZ);
  add_block(md->pred_mean,md->pred_cov,md->n,buf,B_SZ);

  /* update the daala satd */
  satd=0;
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      satd+=sqrt(od_scale[i]*od_scale[j])*abs(buf[B_SZ*j+i]);
    }
  }
  md->satd_avg+=(satd-md->satd_avg)/md->n;
  fr->satd_avg+=(satd-fr->satd_avg)/fr->n;
}

#if WRITE_IMAGES
static void od_image_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  intra_stats_ctx *ctx;
  intra_stats_od  *od;
  int              mode;
  od_coeff        *p;
  int              j;
  int              i;
  unsigned char    buf[B_SZ*B_SZ];
#if PRINT_PROGRESS 
  if(_bi==0&&_bj==0){
    printf("in od_image_block\n");
  }
#endif
  ctx=(intra_stats_ctx *)_ctx;
  od=&ctx->img.od_stats;
  mode=od->mode[_bj*ctx->img.nxblocks+_bi];

  od_rgba16_image_draw_point(&ctx->img.od_map,_bi,_bj,colors[mode]);

  p=&od->pre[od->pre_stride*(B_SZ*_bj+(3*B_SZ>>1))+B_SZ*_bi+(3*B_SZ>>1)];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=OD_CLAMPI(p[od->pre_stride*j+i]+128,0,255);
    }
  }
  image_draw_block(&ctx->img.pre,B_SZ*_bi,B_SZ*_bj,buf,B_SZ);

  p=&od->post[od->post_stride*(B_SZ*_bj+(B_SZ>>1))+B_SZ*_bi+(B_SZ>>1)];
  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=OD_CLAMPI(p[od->post_stride*j+i]+128,0,255);
    }
  }
  image_draw_block(&ctx->img.od_pred,B_SZ*_bi,B_SZ*_bj,buf,B_SZ);

  for(j=0;j<B_SZ;j++){
    for(i=0;i<B_SZ;i++){
      buf[B_SZ*j+i]=abs(_data[_stride*j+i]-buf[B_SZ*j+i]);
    }
  }
  image_draw_block(&ctx->img.od_res,B_SZ*_bi,B_SZ*_bj,buf,B_SZ);
}
#endif

static int stats_finish(void *_ctx){
  intra_stats_ctx *ctx;
#if WRITE_IMAGES
  char             name[8196];
  int              eos;
#endif
#if PRINT_PROGRESS
  printf("in stats_finish\n");
#endif
  ctx=(intra_stats_ctx *)_ctx;
#if WRITE_IMAGES
  strcpy(name,ctx->img.name);
  eos=strlen(ctx->img.name)-4;
  sprintf(&name[eos],".luma");
  image_write_png(&ctx->img.luma,name);
  sprintf(&name[eos],".luma_pre");
  image_write_png(&ctx->img.pre,name);
  sprintf(&name[eos],".luma_vp8_map");
  image_write_png(&ctx->img.vp8_map,name);
  sprintf(&name[eos],".luma_vp8_pred");
  image_write_png(&ctx->img.vp8_pred,name);
  sprintf(&name[eos],".luma_vp8_res");
  image_write_png(&ctx->img.vp8_res,name);
  sprintf(&name[eos],".luma_daala_map");
  image_write_png(&ctx->img.od_map,name);
  sprintf(&name[eos],".luma_daala_pred");
  image_write_png(&ctx->img.od_pred,name);
  sprintf(&name[eos],".luma_daala_res");
  image_write_png(&ctx->img.od_res,name);
#endif
  intra_stats_combine(&ctx->gb.vp8_stats,&ctx->img.vp8_stats);
  intra_stats_correct(&ctx->img.vp8_stats);
  intra_stats_print(&ctx->img.vp8_stats,"VP8 Intra Predictors",vp8_scale);
  intra_stats_combine(&ctx->gb.od_stats,&ctx->img.od_stats.super);
  intra_stats_correct(&ctx->img.od_stats.super);
  intra_stats_print(&ctx->img.od_stats.super,"Daala Intra Predictors",od_scale);
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

#define PADDING (4*B_SZ)
#if PADDING<3*B_SZ
# error "PADDING must be at least 3*B_SZ"
#endif

int main(int _argc,const char **_argv){
  intra_stats_ctx ctx;
  int ret;
#if WRITE_IMAGES
  intra_map_colors(colors,OD_INTRA_NMODES);
#endif
  vp8_scale_init(vp8_scale);
  od_scale_init(od_scale);
  intra_stats_ctx_init(&ctx);
  ret=apply_to_blocks2(&ctx,PADDING,stats_start,BLOCKS,
   sizeof(BLOCKS)/sizeof(*BLOCKS),stats_finish,0x1,_argc,_argv);
  if(ret)return ret;
  printf("Processed %i image(s)\n",ctx.n);
  intra_stats_correct(&ctx.gb.vp8_stats);
  intra_stats_print(&ctx.gb.vp8_stats,"VP8 Intra Predictors",vp8_scale);
  intra_stats_correct(&ctx.gb.od_stats);
  intra_stats_print(&ctx.gb.od_stats,"Daala Intra Predictors",od_scale);
  return ret;
}
