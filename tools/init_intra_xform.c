#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "intra_fit_tools.h"
#include "../src/dct.h"
#include "svd.h"

typedef struct intra_xform_ctx intra_xform_ctx;

struct intra_xform_ctx{
  char          *map_filename;
  unsigned char *map;
  int            nxblocks;
  int            nyblocks;
  long long      r_n[OD_INTRA_NMODES];
  double         r_x[OD_INTRA_NMODES][2*B_SZ*2*B_SZ];
  double         r_xx[OD_INTRA_NMODES][2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];
  double         scale[OD_INTRA_NMODES][2*B_SZ*2*B_SZ];
  double         beta[OD_INTRA_NMODES][B_SZ*B_SZ][2*B_SZ*2*B_SZ];
  long long      n;
  double         satd_avg;
};

static int intra_xform_train_plane_start(void *_ctx,const char *_name,
 const th_info *_ti,int _pli,int _nxblocks,int _nyblocks){
  intra_xform_ctx *ctx;
  FILE            *map_file;
  char            *map_filename;
  ctx=(intra_xform_ctx *)_ctx;
  map_filename=get_map_filename(_name,_pli,_nxblocks,_nyblocks);
  map_file=fopen(map_filename,"rb");
  if(map_file==NULL){
    fprintf(stderr,"Error opening input file '%s'.\n",map_filename);
    return EXIT_FAILURE;
  }
  ctx->map_filename=map_filename;
  ctx->map=(unsigned char *)malloc(_nxblocks*(size_t)_nyblocks);
  if(fread(ctx->map,_nxblocks*(size_t)_nyblocks,1,map_file)<1){
    fprintf(stderr,"Error reading from input file '%s'.\n",map_filename);
    return EXIT_FAILURE;
  }
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  return EXIT_SUCCESS;
}

#define APPLY_PREFILTER (1)
#define APPLY_DCT (1)

static void xform_blocks(od_coeff _buf[3*B_SZ*3*B_SZ],
 const unsigned char *_data,int _stride){
  od_coeff            *buf2;
  od_coeff             col[B_SZ];
  od_coeff            *row;
  const unsigned char *origin;
  int                  bx;
  int                  by;
  int                  x;
  int                  j;
  int                  i;
  origin=_data-(3*B_SZ>>1)*_stride-(3*B_SZ>>1);
  for(by=0;by<3;by++){
    for(bx=0;bx<3;bx++){
      for(j=0;j<B_SZ;j++){
        x=B_SZ*bx+j;
#if APPLY_PREFILTER
        for(i=0;i<B_SZ;i++)col[i]=origin[_stride*(B_SZ*by+i)+x]-128;
        od_pre_filter4(col,col);
        for(i=0;i<B_SZ;i++)_buf[3*B_SZ*(B_SZ*by+i)+x]=col[i];
#else
        for(i=0;i<B_SZ;i++){
          _buf[3*B_SZ*(B_SZ*by+i)+x]=origin[_stride*(B_SZ*by+i)+x]-128;
        }
#endif
      }
    }
  }
#if APPLY_PREFILTER
  for(by=0;by<3;by++){
    for(bx=0;bx<3;bx++){
      for(i=0;i<B_SZ;i++){
        row=_buf+3*B_SZ*(B_SZ*by+i)+B_SZ*bx;
        od_pre_filter4(row,row);
      }
    }
  }
#endif
  buf2=_buf+3*B_SZ*(B_SZ>>1)+(B_SZ>>1);
#if APPLY_DCT
  for(by=0;by<2;by++){
    for(bx=0;bx<2;bx++){
      for(i=0;i<B_SZ;i++){
        row=buf2+3*B_SZ*(B_SZ*by+i)+B_SZ*bx;
        od_bin_fdct4(row,row);
      }
      for(j=0;j<B_SZ;j++){
        for(i=0;i<B_SZ;i++)col[i]=buf2[3*B_SZ*(B_SZ*by+i)+B_SZ*bx+j];
        od_bin_fdct4(col,col);
        for(i=0;i<B_SZ;i++)buf2[3*B_SZ*(B_SZ*by+i)+B_SZ*bx+j]=col[i];
      }
    }
  }
#endif
}

static void intra_xform_train_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  intra_xform_ctx     *ctx;
  double               delta[2*B_SZ*2*B_SZ];
  double               dw;
  od_coeff             buf[3*B_SZ*3*B_SZ];
  od_coeff            *buf2;
  int                  mode;
  long long            n;
  int                  j;
  int                  i;
  int                  k;
  int                  l;
  xform_blocks(buf,_data,_stride);
  buf2=buf+3*B_SZ*(B_SZ>>1)+(B_SZ>>1);
  ctx=(intra_xform_ctx *)_ctx;
  mode=ctx->map[_bj*ctx->nxblocks+_bi];
  n=ctx->r_n[mode]++;
  dw=1.0/(n+1);
  for(i=0;i<2*B_SZ;i++){
    for(j=0;j<2*B_SZ;j++){
      int    ci;
      ci=2*B_SZ*i+j;
      delta[ci]=(buf2[3*B_SZ*i+j]-ctx->r_x[mode][ci]);
      ctx->r_x[mode][ci]+=delta[ci]*dw;
    }
  }
  for(i=0;i<2*B_SZ;i++){
    for(j=0;j<2*B_SZ;j++){
      int    ci;
      ci=2*B_SZ*i+j;
      for(k=0;k<2*B_SZ;k++){
        for(l=0;l<2*B_SZ;l++){
          int cj;
          cj=2*B_SZ*k+l;
          ctx->r_xx[mode][ci][cj]+=n*dw*delta[ci]*delta[cj];
        }
      }
    }
  }
}

static int intra_xform_train_plane_finish(void *_ctx){
  intra_xform_ctx *ctx;
  ctx=(intra_xform_ctx *)_ctx;
  free(ctx->map_filename);
  free(ctx->map);
  return EXIT_SUCCESS;
}


static const char *MODE_NAME[OD_INTRA_NMODES]={
  "OD_INTRA_DC","OD_INTRA_TM","OD_INTRA_HU","OD_INTRA_HE","OD_INTRA_HD",
  "OD_INTRA_RD","OD_INTRA_VR","OD_INTRA_VE","OD_INTRA_VL","OD_INTRA_LD"
};

static void print_beta(int _mode,int _i,int _j,double *_beta){
  int i;
  int j;
  printf("      /*%s (%i,%i)*/\n",MODE_NAME[_mode],_i,_j);
  printf("      {\n");
  for(j=0;j<2*B_SZ;j++){
    if(j==B_SZ)printf("\n");
    printf("        {");
    /*printf("        {\n");
    printf("          ");*/
    for(i=0;i<2*B_SZ;i++){
      printf("%s%- 24.18G%s",i==B_SZ?"   ":"",_beta[2*B_SZ*j+i],i<2*B_SZ-1?",":"");
    }
    /*printf("        }%s\n",j<2*B_SZ-1?",":"");*/
    printf("}%s\n",j<2*B_SZ-1?",":"");
  }
  printf("      }%s\n",_j<B_SZ-1?",":"");
}

typedef double r_xx_row[2*B_SZ*2*B_SZ];

static void update_intra_xforms(intra_xform_ctx *_ctx){
  int mode;
  /*Update the model for each coefficient in each mode.*/
  printf("#include \"intra.h\"\n");
  printf("\n");
  printf("double OD_INTRA_PRED_WEIGHTS_%ix%i"
   "[OD_INTRA_NMODES][%i][%i][2*%i][2*%i]={\n",
   B_SZ,B_SZ,B_SZ,B_SZ,B_SZ,B_SZ);
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    int        xi[2*B_SZ*2*B_SZ];
    int        nxi;
    int        i;
    int        j;
    long long  r_n;
    double    *r_x;
    r_xx_row  *r_xx;
    double    *scale;
    r_n=_ctx->r_n[mode];
    r_x=_ctx->r_x[mode];
    r_xx=_ctx->r_xx[mode];
    scale=_ctx->scale[mode];
    printf("  {\n");
    for(i=0;i<2*B_SZ*2*B_SZ;i++){
      scale[i]=sqrt(r_xx[i][i]);
    }
    for(i=0;i<2*B_SZ*2*B_SZ;i++){
      for(j=0;j<2*B_SZ*2*B_SZ;j++){
        r_xx[i][j]/=scale[i]*scale[j];
      }
    }
    nxi=0;
    for(j=0;j<B_SZ;j++){
      for(i=0;i<B_SZ;i++){
        xi[nxi]=2*B_SZ*j+i;
        xi[nxi+B_SZ*B_SZ]=2*B_SZ*j+B_SZ+i;
        xi[nxi+2*B_SZ*B_SZ]=2*B_SZ*(B_SZ+j)+i;
        nxi++;
      }
    }
#if 0
    if(mode==0){
      for(i=0;i<2*B_SZ;i++){
        for(j=0;j<2*B_SZ;j++){
          int k;
          int l;
          for(k=0;k<2*B_SZ;k++){
            for(l=0;l<2*B_SZ;l++){
              printf("%0.18G%s",r_xx[2*B_SZ*i+j][2*B_SZ*k+l],2*B_SZ*k+l>=2*B_SZ*2*B_SZ-1?"\n":" ");
            }
          }
        }
      }
    }
#endif
    for(i=0;i<B_SZ;i++){
      printf("    {\n");
      for(j=0;j<B_SZ;j++){
        double  xtx[2*2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];
        double *xtxp[2*2*B_SZ*2*B_SZ];
        double  xty[2*B_SZ*2*B_SZ];
        double  s[2*B_SZ*2*B_SZ];
        double *beta;
        int     xii;
        int     xij;
        int     yi;
        int     k;
        int     l;
        nxi=3*B_SZ*B_SZ;
        for(k=0;k<=i;k++){
          for(l=0;l<=j;l++){
            xi[nxi++]=2*B_SZ*(B_SZ+k)+B_SZ+l;
          }
        }
        nxi--;
        yi=2*B_SZ*(B_SZ+i)+B_SZ+j;
        for(xii=0;xii<nxi;xii++){
          for(xij=0;xij<nxi;xij++){
            xtx[xii][xij]=r_xx[xi[xii]][xi[xij]];
          }
          xty[xii]=r_xx[xi[xii]][yi];
        }
        for(xii=0;xii<2*nxi;xii++)xtxp[xii]=xtx[xii];
        svd_pseudoinverse(xtxp,s,nxi,nxi);
        beta=_ctx->beta[mode][B_SZ*i+j];
        memset(beta,0,2*B_SZ*2*B_SZ*sizeof(*beta));
        /*beta[yi]=r_x[yi];*/
        for(xii=0;xii<nxi;xii++){
          double beta_i;
          beta_i=0;
          for(xij=0;xij<nxi;xij++)beta_i+=xtx[xij][xii]*xty[xij];
          beta[xi[xii]]=beta_i*scale[yi]/scale[xi[xii]];
          /*beta[yi]-=beta_i*r_x[xi[xii]];*/
        }
        print_beta(mode,i,j,beta);
      }
      printf("    }%s\n",i<B_SZ-1?",":"");
    }
    printf("  }%s\n",mode<OD_INTRA_NMODES-1?",":"");
  }
  printf("};\n");
}


static int intra_xform_update_plane_start(void *_ctx,const char *_name,
 const th_info *_ti,int _pli,int _nxblocks,int _nyblocks){
  intra_xform_ctx *ctx;
  ctx=(intra_xform_ctx *)_ctx;
  ctx->map_filename=get_map_filename(_name,_pli,_nxblocks,_nyblocks);
  ctx->map=(unsigned char *)malloc(_nxblocks*(size_t)_nyblocks);
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  return EXIT_SUCCESS;
}

static void intra_xform_update_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  intra_xform_ctx *ctx;
  od_coeff         buf[3*B_SZ*3*B_SZ];
  od_coeff         buf2[2*B_SZ*2*B_SZ];
  unsigned         satd;
  unsigned         best_satd;
  int              mode;
  int              best_mode;
  int              j;
  int              i;
  int              k;
  int              l;
  xform_blocks(buf,_data,_stride);
  ctx=(intra_xform_ctx *)_ctx;
  best_satd=UINT_MAX;
  best_mode=0;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    for(k=0;k<2*B_SZ;k++){
      for(l=0;l<2*B_SZ;l++){
        buf2[2*B_SZ*k+l]=buf[3*B_SZ*(k+(B_SZ>>1))+l+(B_SZ>>1)];
      }
    }
    satd=0;
    for(i=0;i<B_SZ;i++){
      for(j=0;j<B_SZ;j++){
        const double *beta;
        double        p;
        beta=ctx->beta[mode][i*B_SZ+j];
        p=0;
        for(k=0;k<2*B_SZ;k++){
          for(l=0;l<2*B_SZ;l++){
            p+=beta[2*B_SZ*k+l]*buf2[2*B_SZ*k+l];
          }
        }
        satd+=abs(buf[2*B_SZ*(i+B_SZ)+j+B_SZ]-(od_coeff)floor(p+0.5));
      }
    }
    if(satd<best_satd){
      best_satd=satd;
      best_mode=mode;
    }
  }
  ctx->satd_avg+=(best_satd-ctx->satd_avg)/++(ctx->n);
  ctx->map[_bj*ctx->nxblocks+_bi]=best_mode;
}

static int intra_xform_update_plane_finish(void *_ctx){
  intra_xform_ctx *ctx;
  FILE            *map_file;
  ctx=(intra_xform_ctx *)_ctx;
  map_file=fopen(ctx->map_filename,"wb");
  if(map_file==NULL){
    fprintf(stderr,"Error opening output file '%s'.\n",ctx->map_filename);
    return EXIT_FAILURE;
  }
  if(fwrite(ctx->map,ctx->nxblocks*(size_t)ctx->nyblocks,1,map_file)<1){
    fprintf(stderr,"Error writing to output file '%s'.\n",ctx->map_filename);
    return EXIT_FAILURE;
  }
  printf("Average SATD: %G\n",ctx->satd_avg);
  return intra_xform_train_plane_finish(_ctx);
}


int main(int _argc,const char **_argv){
  static intra_xform_ctx ctx;
  int                    ret;
  ret=apply_to_blocks(&ctx,intra_xform_train_plane_start,
   intra_xform_train_block,intra_xform_train_plane_finish,_argc,_argv);
  if(ret==EXIT_SUCCESS){
    update_intra_xforms(&ctx);
    ret=apply_to_blocks(&ctx,intra_xform_update_plane_start,
     intra_xform_update_block,intra_xform_update_plane_finish,_argc,_argv);
  }
  return ret;
}
