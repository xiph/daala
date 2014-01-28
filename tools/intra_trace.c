#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include "od_defs.h"
#include "od_filter.h"
#include "image_tools.h"
#include "stats_tools.h"
#include "trans_tools.h"
#include "od_intra.h"

typedef struct intra_trace_ctx intra_trace_ctx;

struct intra_trace_ctx{
  int        x;
  int        y;
  int        d;
  int        c;
  image_data img;
  od_coeff   min[OD_INTRA_NMODES];
  od_coeff   max[OD_INTRA_NMODES];
};

void intra_trace_ctx_reset(intra_trace_ctx *_this){
  _this->x=0;
  _this->y=0;
  _this->d=3;
  _this->c=0;
}

void intra_trace_ctx_init(intra_trace_ctx *_this){
  int i;
  intra_trace_ctx_reset(_this);
  image_data_init(&_this->img,NULL,B_SZ_LOG,1,1);
  for(i=0;i<OD_INTRA_NMODES;i++){
    _this->min[i]=INT_MAX;
    _this->max[i]=-INT_MAX;
  }
}

void intra_trace_ctx_clear(intra_trace_ctx *_this){
  image_data_clear(&_this->img);
}

int intra_trace_ctx_next(intra_trace_ctx *_this,int _loop){
  int m;
  /*printf("%i %i %i\n",_this->x,_this->y,_this->d);*/
  m=0;
  if(_this->d&2){
    if(_this->d&1){
      if(_this->y==B_SZ*3-1){
        _this->x+=_this->x&1?-1:1;
        _this->d=2;
        m=1;
      }
    }
    else{
      if(_this->y==(_this->x^1)){
        _this->x++;
        _this->d=1;
        m=1;
      }
    }
  }
  else{
    if(_this->d&1){
      if(_this->x==B_SZ*4-1){
        _this->y+=_this->y&1?-1:1;
        _this->d=0;
        m=1;
      }
    }
    else{
      if(_this->x==(_this->y&1?_this->y+1:_this->y-1)){
        _this->y++;
        _this->d=3;
        m=1;
      }
    }
  }
  if(!m){
    if(_this->d&2){
      _this->y+=_this->d&1?1:-1;
    }
    else{
      _this->x+=_this->d&1?1:-1;
    }
  }
  if(_loop&&_this->x==2*B_SZ&&_this->y==2*B_SZ){
    _this->x=2*B_SZ-1;
    _this->y=2*B_SZ-1;
    _this->d=3;
  }
  _this->c++;
  if(_loop){
    return(_this->x==0&&_this->y==0);
  }
  else{
    return(_this->x==B_SZ&&_this->y==B_SZ);
  }
}

#define SUPPORT_SZ 32
#define SPACE_SZ   16
#define BLOCK_SZ   128
#define IMAGE_X    (5*(SPACE_SZ+BLOCK_SZ)+SPACE_SZ)
#define IMAGE_Y    (3*(SPACE_SZ+BLOCK_SZ)+SPACE_SZ)

static void draw_pulse(od_rgba16_image *_image,int _x,int _y){
  int x0;
  int y0;
  od_rgba16_pixel c;
  int n;
  int j;
  int i;
  x0=(IMAGE_X-4*SUPPORT_SZ)/2;
  y0=SPACE_SZ+(BLOCK_SZ-3*SUPPORT_SZ)/2;
  n=SUPPORT_SZ/B_SZ;
  c[3]=(unsigned short)0xFFFFU;
  /* Clear region. */
  c[0]=c[1]=c[2]=0;
  for(j=0;j<3*SUPPORT_SZ;j++){
    for(i=0;i<4*SUPPORT_SZ;i++){
      od_rgba16_image_draw_point(_image,x0+i,y0+j,c);
    }
  }
  /* Draw blocks. */
  x0+=SUPPORT_SZ/2;
  y0+=SUPPORT_SZ/2;
  c[0]=c[1]=c[2]=0x7F00;
  for(j=0;j<SUPPORT_SZ;j++){
    for(i=0;i<SUPPORT_SZ;i++){
      od_rgba16_image_draw_point(_image,x0+SUPPORT_SZ*0+i,y0+SUPPORT_SZ*0+j,c);
      od_rgba16_image_draw_point(_image,x0+SUPPORT_SZ*1+i,y0+SUPPORT_SZ*0+j,c);
      od_rgba16_image_draw_point(_image,x0+SUPPORT_SZ*2+i,y0+SUPPORT_SZ*0+j,c);
      od_rgba16_image_draw_point(_image,x0+SUPPORT_SZ*0+i,y0+SUPPORT_SZ*1+j,c);
    }
  }
  /* Draw pulse. */
  x0-=SUPPORT_SZ/2;
  y0-=SUPPORT_SZ/2;
  c[0]=c[1]=c[2]=0xFF00;
  for(j=0;j<n;j++){
    for(i=0;i<n;i++){
      od_rgba16_image_draw_point(_image,x0+n*_x+i,y0+n*_y+j,c);
    }
  }
}

static void draw_point(od_rgba16_image *_image,int _m,int _x,int _y,
 unsigned char _v){
  int             x0;
  int             y0;
  int             n;
  od_rgba16_pixel c;
  int             j;
  int             i;
  x0=SPACE_SZ+(SPACE_SZ+BLOCK_SZ)*(_m%5);
  y0=(SPACE_SZ+BLOCK_SZ)*(_m/5+1);
  n=BLOCK_SZ/(2*B_SZ);
  c[3]=(unsigned short)0xFFFFU;
  c[0]=c[1]=c[2]=_v*0x101;
  for(j=0;j<n;j++){
    for(i=0;i<n;i++){
      od_rgba16_image_draw_point(_image,x0+n*_x+i,y0+n*_y+j,c);
    }
  }
}

static void intra_trace_ctx_step(intra_trace_ctx *_ctx,od_rgba16_image *_image,
 const int *_f){
  unsigned char data[4*B_SZ*4*B_SZ];
  int           m;
  (void)_f;
  memset(data,0,4*B_SZ*4*B_SZ);
  data[4*B_SZ*_ctx->y+_ctx->x]=255;
  draw_pulse(_image,_ctx->x,_ctx->y);
  image_data_pre_block(&_ctx->img,&data[4*B_SZ*(3*B_SZ>>1)+(3*B_SZ>>1)],4*B_SZ,
   0,0);
  image_data_fdct_block(&_ctx->img,0,0);
#if TF_BLOCKS
  image_data_tf_block(&_ctx->img,0,0);
#endif
  for(m=0;m<OD_INTRA_NMODES;m++){
    int j;
    int i;
    _ctx->img.mode[0]=m;
    image_data_pred_block(&_ctx->img,0,0);
    image_data_idct_block(&_ctx->img,0,0);
    /* Set non-predicted blocks to zero before post-filtering. */
    for(j=0;j<3*B_SZ;j++){
      for(i=0;i<3*B_SZ;i++){
        if(j/B_SZ!=1||i/B_SZ!=1){
          _ctx->img.idct[3*B_SZ*j+i]=-128*INPUT_SCALE;
        }
      }
    }
    image_data_post_block(&_ctx->img,0,0);
    for(j=0;j<2*B_SZ;j++){
      for(i=0;i<2*B_SZ;i++){
        od_coeff v;
        v=(_ctx->img.post[2*B_SZ*j+i]+INPUT_SCALE*128+INPUT_SCALE/2)/INPUT_SCALE;
        if(_ctx->min[m]>v){
          _ctx->min[m]=v;
        }
        if(_ctx->max[m]<v){
          _ctx->max[m]=v;
        }
        draw_point(_image,m,i,j,OD_CLAMPI(0,v+128,255));
      }
    }
  }
}

int main(int _argc,char *_argv[]){
  intra_trace_ctx  ctx;
  const int       *f;
  od_rgba16_image  image;
  (void)_argc;
  (void)_argv;
  ne_filter_params_init();
  od_scale_init(OD_SCALE[B_SZ_LOG-OD_LOG_BSIZE0],B_SZ_LOG);
#if B_SZ==4
    f=NE_FILTER_PARAMS4;
#elif B_SZ==8
    f=NE_FILTER_PARAMS8;
#elif B_SZ==16
    f=NE_FILTER_PARAMS16;
#else
# error "Unsupported block size."
#endif
  od_intra_init();
  intra_trace_ctx_init(&ctx);
  od_rgba16_image_init(&image,IMAGE_X,IMAGE_Y);
  do {
    char name[16];
    intra_trace_ctx_step(&ctx,&image,f);
    sprintf(name,"trace%04i",ctx.c);
    image_write_png(&image,name);
  }
  while(!intra_trace_ctx_next(&ctx,1));
  {
    int i;
    for(i=0;i<OD_INTRA_NMODES;i++){
      printf("%4i %4i\n",ctx.min[i],ctx.max[i]);
    }
  }
  intra_trace_ctx_clear(&ctx);
  od_rgba16_image_clear(&image);
  return EXIT_SUCCESS;
}
