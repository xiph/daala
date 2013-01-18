#include <stdlib.h>
#include "intra_fit_tools.h"
#include "../src/internal.h"

#define FORCE_DYADIC (0)
#define BITS 6
#define USE_AR95 (0)

typedef struct trans_data trans_data;

struct trans_data{
  int    n;
  double mean[2*B_SZ];
  double cov[2*B_SZ][2*B_SZ];
};

static void trans_data_init(trans_data *_this){
  int j;
  int i;
  _this->n=0;
  for(j=0;j<2*B_SZ;j++){
    _this->mean[j]=0;
    for(i=0;i<2*B_SZ;i++){
      _this->cov[j][i]=0;
    }
  }
}

static void trans_data_add(trans_data *_this,const unsigned char *_data){
  double delta[2*B_SZ];
  int    i;
  int    j;
  _this->n++;
  for(i=0;i<2*B_SZ;i++){
    delta[i]=_data[i]-_this->mean[i];
    _this->mean[i]+=delta[i]/_this->n;
  }
  for(j=0;j<2*B_SZ;j++){
    for(i=0;i<2*B_SZ;i++){
      _this->cov[j][i]+=delta[j]*delta[i]*(_this->n-1)/_this->n;
    }
  }
}

static void trans_data_correct(trans_data *_this){
  int j;
  int i;
  for(j=0;j<2*B_SZ;j++){
    for(i=0;i<2*B_SZ;i++){
      _this->cov[j][i]/=_this->n;
    }
  }
}

static void trans_data_normalize(trans_data *_this){
  int    i;
  int    j;
  double scale[2*B_SZ];
  for(i=0;i<2*B_SZ;i++){
    scale[i]=sqrt(_this->cov[i][i]);
  }
  for(j=0;j<2*B_SZ;j++){
    for(i=0;i<2*B_SZ;i++){
      _this->cov[j][i]/=scale[j]*scale[i];
    }
  }
}

static void trans_data_print(trans_data *_this,FILE *_fp){
  int i;
  int j;
  for(i=0;i<2*B_SZ;i++){
    fprintf(_fp,"%s  %- 24.18G%s",i>0?",":"",_this->mean[i],i==2*B_SZ-1?"\n":"");
  }
  for(j=0;j<2*B_SZ;j++){
    for(i=0;i<2*B_SZ;i++){
      fprintf(_fp,"%s  %- 24.18G%s",i>0?",":"",_this->cov[j][i],i==2*B_SZ-1?"\n":"");
    }
  }
}

typedef struct trans_ctx trans_ctx;

struct trans_ctx{
  const char *name;
  int         nxblocks;
  int         nyblocks;
  trans_data  gb;
};

static void trans_ctx_init(trans_ctx *_this){
  trans_data_init(&_this->gb);
}

static void trans_ctx_image(trans_ctx *_this,const char *_name,int _nxblocks,
 int _nyblocks){
  _this->name=_name;
  _this->nxblocks=_nxblocks;
  _this->nyblocks=_nyblocks;
}

int NE_FILTER_PARAMS4[4];

static void ne_pre_filter4(double _y[4],const double _x[4]){
  double t[4];
  t[3]=_x[0]-_x[3];
  t[2]=_x[1]-_x[2];
  t[1]=_x[1]-(t[2]/2);
  t[0]=_x[0]-(t[3]/2);
  t[2]=t[2]*NE_FILTER_PARAMS4[0]/(1<<BITS);
  t[3]=t[3]*NE_FILTER_PARAMS4[1]/(1<<BITS);
  t[3]+=t[2]*NE_FILTER_PARAMS4[2]/(1<<BITS);
  t[2]+=t[3]*NE_FILTER_PARAMS4[3]/(1<<BITS);
  t[0]+=t[3]/2;
  _y[0]=t[0];
  t[1]+=t[2]/2;
  _y[1]=t[1];
  _y[2]=(t[1]-t[2]);
  _y[3]=(t[0]-t[3]);
}

static void ne_post_filter4(double _x[4],const double _y[4]){
  double t[4];
  t[3]=_y[0]-_y[3];
  t[2]=_y[1]-_y[2];
  t[1]=_y[1]-(t[2]/2);
  t[0]=_y[0]-(t[3]/2);
  t[2]-=t[3]*NE_FILTER_PARAMS4[3]/(1<<BITS);
  t[3]-=t[2]*NE_FILTER_PARAMS4[2]/(1<<BITS);
  t[3]=(t[3]*(1<<BITS))/NE_FILTER_PARAMS4[1];
  t[2]=(t[2]*(1<<BITS))/NE_FILTER_PARAMS4[0];
  t[0]+=t[3]/2;
  _x[0]=t[0];
  t[1]+=t[2]/2;
  _x[1]=t[1];
  _x[2]=(t[1]-t[2]);
  _x[3]=(t[0]-t[3]);
}

typedef void (*ne_filter_func)(double _out[],const double _in[]);

const ne_filter_func NE_PRE_FILTER_DOUBLE[OD_NBSIZES]={
  ne_pre_filter4,
  NULL,
  NULL
};

const ne_filter_func NE_POST_FILTER_DOUBLE[OD_NBSIZES]={
  ne_post_filter4,
  NULL,
  NULL
};

/*The true forward 4-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT4_BASIS[4][4]={
  {
     0.5,                                 0.5,
     0.5,                                 0.5
  },
  {
     0.65328148243818826392832158671359,  0.27059805007309849219986160268319,
    -0.27059805007309849219986160268319, -0.65328148243818826392832158671359
  },
  {
     0.5,                                -0.5,
    -0.5,                                 0.5
  },
  {
     0.27059805007309849219986160268319, -0.65328148243818826392832158671359,
     0.65328148243818826392832158671359, -0.27059805007309849219986160268319
  },
};

void ne_fdct4_double(double _x[],const double _y[],int _in_stride){
  double t[4];
  int    i;
  int    j;
  for(j=0;j<4;j++){
    t[j]=0;
    for(i=0;i<4;i++)t[j]+=DCT4_BASIS[j][i]*_y[i*_in_stride];
  }
  for(j=0;j<4;j++)_x[j]=t[j];
}

void ne_idct4_double(double _x[],int _out_stride,const double _y[]){
  double t[4];
  int    i;
  int    j;
  for(j=0;j<4;j++){
    t[j]=0;
    for(i=0;i<4;i++)t[j]+=DCT4_BASIS[i][j]*_y[i];
  }
  for(j=0;j<4;j++)_x[j*_out_stride]=t[j];
}

typedef void (*ne_fdct_func_1d)(double *_out,const double *_in,int _in_stride);
typedef void (*ne_idct_func_1d)(double *_out,int _out_stride,const double *_in);

const ne_fdct_func_1d OD_FDCT_1D_DOUBLE[OD_NBSIZES]={
  ne_fdct4_double,
  NULL,
  NULL
};

const ne_idct_func_1d OD_IDCT_1D_DOUBLE[OD_NBSIZES]={
  ne_idct4_double,
  NULL,
  NULL
};

static void analysis(double *_y,const double *_x,int _stride){
  int    i;
  int    j;
  double t[2*B_SZ];
  for(i=0;i<_stride;i++){
    for(j=0;j<2*B_SZ;j++){
      t[j]=_x[j*_stride+i];
    }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[0],&t[0]);
    (*NE_PRE_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ],&t[B_SZ]);
#else
# error "Need a prefilter implementation for this block size."
#endif
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_FDCT_1D_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ/2],&t[B_SZ/2],1);
#else
# error "Need an fDCT implementation for this block size."
#endif
    for(j=0;j<B_SZ;j++){
      _y[j*_stride+i]=t[B_SZ/2+j];
    }
  }
}

static void synthesis(double *_y,const double *_x,int _stride){
  int    i;
  int    j;
  double t[2*B_SZ];
  for(i=0;i<_stride;i++){
    for(j=0;j<2*B_SZ;j++){
      t[j]=0;
    }
    for(j=0;j<B_SZ;j++){
      t[B_SZ/2+j]=_x[j*_stride+i];
    }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*OD_IDCT_1D_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ/2],1,&t[B_SZ/2]);
#else
# error "Need an iDCT implementation for this block size."
#endif
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[0],&t[0]);
    (*NE_POST_FILTER_DOUBLE[B_SZ_LOG-OD_LOG_BSIZE0])(&t[B_SZ],&t[B_SZ]);
#else
# error "Need a postfilter implementation for this block size."
#endif
    for(j=0;j<2*B_SZ;j++){
      _y[j*_stride+i]=t[j];
    }
  }
}

static double coding_gain(const double _r[2*B_SZ*2*B_SZ]){
  double gr[B_SZ*2*B_SZ];
  double rgt[2*B_SZ*B_SZ];
  double grgt[B_SZ*B_SZ];
  double i[B_SZ*B_SZ];
  double hi[2*B_SZ*B_SZ];
  int    j;
  int    k;
  double cg;
  analysis(gr,_r,2*B_SZ);
  for(j=0;j<2*B_SZ;j++){
    for(k=0;k<B_SZ;k++){
      rgt[j*B_SZ+k]=gr[k*2*B_SZ+j];
    }
  }
  analysis(grgt,rgt,B_SZ);
  for(j=0;j<B_SZ;j++){
    for(k=0;k<B_SZ;k++){
      i[j*B_SZ+k]=j!=k?0:1;
    }
  }
  synthesis(hi,i,B_SZ);
  cg=1;
  for(j=0;j<B_SZ;j++){
    double h;
    h=0;
    for(k=0;k<2*B_SZ;k++){
      h+=hi[k*B_SZ+j]*hi[k*B_SZ+j];
    }
    cg*=grgt[j*B_SZ+j]*h;
  }
  return(10*log10(1/cg)/B_SZ);
}

static void trans_ctx_search(trans_ctx *_this){
  int    j;
  int    i;
  double r[2*B_SZ*2*B_SZ];
#if USE_AR95
  for(j=0;j<2*B_SZ;j++){
    r[j*2*B_SZ+j]=1;
    for(i=j+1;i<2*B_SZ;i++){
      r[j*2*B_SZ+i]=r[i*2*B_SZ+j]=r[j*2*B_SZ+(i-1)]*0.95;
    }
  }
#else
  for(j=0;j<2*B_SZ;j++){
    for(i=0;i<2*B_SZ;i++){
      r[j*2*B_SZ+i]=_this->gb.cov[j][i];
    }
  }
#endif

#if B_SZ==4
  {
    int    p0;
    int    q0;
    int    s0;
    int    s1;
    double cg;
    double best_cg;
    best_cg=0;
    for(p0=-(1<<BITS);p0<=(1<<BITS);p0++){
      NE_FILTER_PARAMS4[2]=p0;
      for(q0=(1<<BITS);q0>=-(1<<BITS);q0--){
        NE_FILTER_PARAMS4[3]=q0;
        for(s0=(1<<BITS);s0<=2*(1<<BITS);s0++){
          NE_FILTER_PARAMS4[0]=s0;
          for(s1=(1<<BITS);s1<=2*(1<<BITS);s1++){
            NE_FILTER_PARAMS4[1]=s1;
            cg=coding_gain(r);
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

  fprintf(stdout,"cg=%G\n",coding_gain(r));
}

static int t_start(void *_ctx,const char *_name,const th_info *_ti,int _pli,
 int _nxblocks,int _nyblocks){
  trans_ctx *ctx;
  fprintf(stdout,"%s %i %i\n",_name,_nxblocks,_nyblocks);
  fflush(stdout);
  ctx=(trans_ctx *)_ctx;
  trans_ctx_image(ctx,_name,_nxblocks,_nyblocks);
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
    for(y=0;y<ctx->nyblocks*B_SZ;y++){
      for(x=0;x<ctx->nxblocks*B_SZ-(2*B_SZ-1);x++){
        for(z=0;z<2*B_SZ;z++){
          buf[z]=_data[y*_stride+(x+z)];
        }
        trans_data_add(&ctx->gb,buf);
      }
    }
    /* add the columns */
    for(y=0;y<ctx->nyblocks*B_SZ-(2*B_SZ-1);y++){
      for(x=0;x<ctx->nxblocks*B_SZ;x++){
        for(z=0;z<2*B_SZ;z++){
          buf[z]=_data[(y+z)*_stride+x];
        }
        trans_data_add(&ctx->gb,buf);
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
  trans_ctx ctx;
  int       ret;
  trans_ctx_init(&ctx);
  ret=apply_to_blocks2(&ctx,PADDING,t_start,BLOCKS,NBLOCKS,NULL,0x1,_argc,_argv);
  trans_data_normalize(&ctx.gb);
  trans_data_print(&ctx.gb,stdout);
  trans_ctx_search(&ctx);
  return EXIT_SUCCESS;
}
