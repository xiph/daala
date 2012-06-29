#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "intra_fit_tools.h"
#include "../src/dct.h"
#include "svd.h"
#include <math.h>

/* #define INTRA_NO_RDO */

long long hist[B_SZ*B_SZ][8192];
#define NB_CONTEXTS 8
#define P0_COEF 0.002

/* These weight are estimated from the entropy of the residual (in one older run) on subset1-y4m */
const float satd_weights[16] =
{0.027602, 0.100551, 0.255518, 0.342399, 0.097165, 0.162053, 0.340471, 0.427589,
 0.225212, 0.320962, 0.523734, 0.608046, 0.294722, 0.389122, 0.593447, 0.651656, };

#define GET_CONTEXT(modes,pos,m,width) (((modes)[(pos)-1]==(m))*4 + ((modes)[(pos)-(width)-1]==(m))*2 + ((modes)[(pos)-(width)]==(m))*1)

typedef struct intra_xform_ctx intra_xform_ctx;

struct intra_xform_ctx{
  char          *map_filename;
  unsigned char *map;
  char          *weights_filename;
  unsigned      *weights;
  int            nxblocks;
  int            nyblocks;
  double         r_w[OD_INTRA_NMODES];
  double         r_x[OD_INTRA_NMODES][2*B_SZ*2*B_SZ];
  double         r_xx[OD_INTRA_NMODES][2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];
  double         scale[OD_INTRA_NMODES][2*B_SZ*2*B_SZ];
  double         beta[OD_INTRA_NMODES][B_SZ*B_SZ][2*B_SZ*2*B_SZ];
  double         freq[OD_INTRA_NMODES][NB_CONTEXTS][2];
  double         p0[OD_INTRA_NMODES];
  long long      n;
  double         satd_avg;
  double         total_bits;
  double         total_satd;
  double         total_count;
};


static int intra_xform_train_plane_start(void *_ctx,const char *_name,
 const th_info *_ti,int _pli,int _nxblocks,int _nyblocks){
  intra_xform_ctx *ctx;
  FILE            *map_file;
  char            *map_filename;
  FILE            *weights_file;
  char            *weights_filename;
  int              i;
  ctx=(intra_xform_ctx *)_ctx;
  ctx->map=(unsigned char *)malloc(_nxblocks*(size_t)_nyblocks);
  map_filename=get_map_filename(_name,_pli,_nxblocks,_nyblocks);
  map_file=fopen(map_filename,"rb");
  if(map_file==NULL){
    fprintf(stderr,"Error opening input file '%s'.\n",map_filename);
    return EXIT_FAILURE;
  }
  ctx->map_filename=map_filename;
  if(fread(ctx->map,_nxblocks*(size_t)_nyblocks,1,map_file)<1){
    fprintf(stderr,"Error reading from input file '%s'.\n",map_filename);
    return EXIT_FAILURE;
  }
  fclose(map_file);
  ctx->weights=(unsigned *)malloc(
   _nxblocks*(size_t)_nyblocks*sizeof(*ctx->weights));
  weights_filename=get_weights_filename(_name,_pli,_nxblocks,_nyblocks);
  weights_file=fopen(weights_filename,"rb");
  if(weights_file==NULL){
    fprintf(stderr,"Error opening input file '%s'.\n",weights_filename);
    return EXIT_FAILURE;
  }
  ctx->weights_filename=weights_filename;
  if(fread(ctx->weights,
   _nxblocks*(size_t)_nyblocks*sizeof(*ctx->weights),1,weights_file)<1){
    fprintf(stderr,"Error reading from input file '%s'.\n",weights_filename);
    return EXIT_FAILURE;
  }
  fclose(weights_file);
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  for(i=0;i<OD_INTRA_NMODES;i++)
    ctx->p0[i]=.1;
  return EXIT_SUCCESS;
}

#define APPLY_PREFILTER (1)
#define APPLY_DCT (1)

static od_coeff *xform_blocks(od_coeff _buf[3*B_SZ*3*B_SZ],
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
  return buf2;
}

static void intra_xform_train_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  intra_xform_ctx     *ctx;
  double               delta[2*B_SZ*2*B_SZ];
  double               dw;
  od_coeff             buf[3*B_SZ*3*B_SZ];
  od_coeff            *buf2;
  int                  mode;
  double               w;
  unsigned             wb;
  int                  j;
  int                  i;
  int                  k;
  int                  l;
  ctx=(intra_xform_ctx *)_ctx;
  wb=ctx->weights[_bj*ctx->nxblocks+_bi];
  if(wb<=0)return;
  buf2=xform_blocks(buf,_data,_stride);
  mode=ctx->map[_bj*ctx->nxblocks+_bi];
  w=ctx->r_w[mode];
  ctx->r_w[mode]+=wb;
  dw=wb/(w+wb);
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
          ctx->r_xx[mode][ci][cj]+=w*dw*delta[ci]*delta[cj];
        }
      }
    }
  }
  if (_bi>0 && _bj>0)
  {
    unsigned char *modes;
    int pos;
    int m;
    int width;
    modes=ctx->map;
    pos = _bj*ctx->nxblocks+_bi;
    width=ctx->nxblocks;
    for(m=0;m<OD_INTRA_NMODES;m++)
    {
      int c;
      c = GET_CONTEXT(modes,pos,m,width);
      ctx->freq[m][c][0]+=1;
      ctx->freq[m][c][1] += (mode==m);
    }
  }
}

static int intra_xform_train_plane_finish(void *_ctx){
  intra_xform_ctx *ctx;
  ctx=(intra_xform_ctx *)_ctx;
  free(ctx->weights_filename);
  free(ctx->weights);
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
    /*double    *r_x;*/
    r_xx_row  *r_xx;
    double    *scale;
    /*r_x=_ctx->r_x[mode];*/
    r_xx=_ctx->r_xx[mode];
    scale=_ctx->scale[mode];
    printf("  {\n");
    for(i=0;i<2*B_SZ*2*B_SZ;i++){
      scale[i]=sqrt(r_xx[i][i]);
      if(scale[i]<=0)scale[i]=1;
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
#if 1
        /*Include coefficients for the current block*/
        for(k=0;k<=i;k++){
          for(l=0;l<=j;l++){
            xi[nxi++]=2*B_SZ*(B_SZ+k)+B_SZ+l;
          }
        }
        nxi--;
#endif
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
  ctx->weights_filename=get_weights_filename(_name,_pli,_nxblocks,_nyblocks);
  ctx->map=(unsigned char *)malloc(_nxblocks*(size_t)_nyblocks);
  ctx->weights=(unsigned *)malloc(
   _nxblocks*(size_t)_nyblocks*sizeof(*ctx->weights));
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  return EXIT_SUCCESS;
}

static void intra_xform_update_block(void *_ctx,const unsigned char *_data,
 int _stride,int _bi,int _bj){
  intra_xform_ctx *ctx;
  od_coeff         buf[3*B_SZ*3*B_SZ];
  od_coeff        *buf2;
  double           best_satd;
  double           best_rlsatd;
  double           best_bits;
  double           next_best_satd;
  double           next_best_rlsatd;
  int              mode;
  int              best_mode;
  int              c0[OD_INTRA_NMODES]={0};
  double           bits=0;
  buf2=xform_blocks(buf,_data,_stride);
  ctx=(intra_xform_ctx *)_ctx;
  best_mode=0;
  best_satd=UINT_MAX;
  best_rlsatd=UINT_MAX;
  best_bits=0;
  next_best_satd=UINT_MAX;
  next_best_rlsatd=UINT_MAX;
  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double satd;
    double rlsatd;
    double diff;
    int      i;
    int      j;
    satd=0;
    for(i=0;i<B_SZ;i++){
      for(j=0;j<B_SZ;j++){
        const double *beta;
        double        p;
        int           k;
        int           l;
        beta=ctx->beta[mode][B_SZ*i+j];
        p=0;
        for(k=0;k<2*B_SZ;k++){
          for(l=0;l<2*B_SZ;l++){
            p+=beta[2*B_SZ*k+l]*buf2[3*B_SZ*k+l];
          }
        }
#ifdef INTRA_NO_RDO
        diff = fabs(buf2[3*B_SZ*(i+B_SZ)+j+B_SZ]-(od_coeff)floor(p+0.5));
        satd+=diff;
#else
        /* Simulates quantization with dead zone (without the annoying quantization effects) */
        diff = fabs(buf2[3*B_SZ*(i+B_SZ)+j+B_SZ]-(od_coeff)floor(p+0.5)) - 1;
        if (diff<0)diff=0;
        satd+=satd_weights[i*B_SZ+j]*diff;
#endif
        hist[i*B_SZ+j][4096+(buf2[3*B_SZ*(i+B_SZ)+j+B_SZ]-(od_coeff)floor(p+0.5))]++;
      }
    }
    rlsatd=satd;
    if (_bj>0){
      float p[OD_INTRA_NMODES];
      unsigned char *modes;
      int pos;
      int m;
      int width;
      float sum=0;
      float maxP=0;
      int maxM=0;
      modes=ctx->map;
      pos = _bj*ctx->nxblocks+_bi;
      width=ctx->nxblocks;
      for(m=0;m<OD_INTRA_NMODES;m++)
      {
        int c;
        if (_bi>0)
          c=GET_CONTEXT(modes,pos,m,width);
        else
          c=0;
        p[m] = ctx->freq[m][c][1]/(float)ctx->freq[m][c][0];
        /*p[m] = p[m]/(1.2-p[m]);*/
        p[m]+=1.e-5;
        if (p[m]<ctx->p0[m]) p[m]=ctx->p0[m];
        if (c==0 && _bi>0)
        {
          ctx->p0[m]*=(1-P0_COEF);
          c0[m]=1;
        }
        sum += p[m];
        if (p[m]>maxP)
        {
          maxP=p[m];
          maxM=m;
        }
      }
      /* Normalize all probabilities except the max */
      /*bits = -log(p[mode]/sum)/log(2);*/
      /* Normalize all probabilities except the max */
      bits = (mode==maxM) ? -log(p[mode])/log(2) : -log(p[mode]*(1-maxP)/(sum-maxP))/log(2);
#ifndef INTRA_NO_RDO
      satd += .5*bits;
#endif
      /* Bias towards DC mode */
      /*if (mode==0)satd-=.5;*/
    }
    if(satd<best_satd){
      next_best_satd=best_satd;
      next_best_rlsatd=best_rlsatd;
      best_satd=satd;
      best_rlsatd=rlsatd;
      best_mode=mode;
      best_bits=bits;
    }
    else if(satd<next_best_satd){next_best_satd=satd;next_best_rlsatd=rlsatd;}
  }
  if (c0[best_mode])
    ctx->p0[best_mode] += P0_COEF;
  /*fprintf(stderr,"%f\n", best_bits);*/
  ctx->total_bits += best_bits;
  ctx->total_satd+=best_satd;
  ctx->total_count += 1;
  ctx->satd_avg+=(best_satd-ctx->satd_avg)/++(ctx->n);
  ctx->map[_bj*ctx->nxblocks+_bi]=best_mode;
  ctx->weights[_bj*ctx->nxblocks+_bi]=floor((next_best_rlsatd-best_rlsatd)*1000.);
  if(next_best_rlsatd<=best_rlsatd)ctx->weights[_bj*ctx->nxblocks+_bi]=0;
  if(best_mode==0)ctx->weights[_bj*ctx->nxblocks+_bi]=1;
}

static int intra_xform_update_plane_finish(void *_ctx){
  intra_xform_ctx *ctx;
  FILE            *map_file;
  FILE            *weights_file;
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
  fclose(map_file);
  weights_file=fopen(ctx->weights_filename,"wb");
  if(weights_file==NULL){
    fprintf(stderr,"Error opening output file '%s'.\n",ctx->weights_filename);
    return EXIT_FAILURE;
  }
  if(fwrite(ctx->weights,
   ctx->nxblocks*(size_t)ctx->nyblocks*sizeof(*ctx->weights),1,
   weights_file)<1){
    fprintf(stderr,"Error writing to output file '%s'.\n",
     ctx->weights_filename);
    return EXIT_FAILURE;
  }
  fclose(weights_file);
  printf("Average SATD: %G\n",ctx->satd_avg);
  return intra_xform_train_plane_finish(_ctx);
}


int main(int _argc,const char **_argv){
  static intra_xform_ctx ctx;
  int                    ret;
  int                    i;
  for(i=0;i<OD_INTRA_NMODES;i++)
  {
    int j;
    for(j=0;j<NB_CONTEXTS;j++)
    {
      ctx.freq[i][j][0]=2;
      ctx.freq[i][j][1]=1;
    }
  }
  ctx.total_bits=0;
  ctx.total_count=0;
  ctx.total_satd=0;

  ret=apply_to_blocks(&ctx,intra_xform_train_plane_start,
   intra_xform_train_block,intra_xform_train_plane_finish,_argc,_argv);
  if(ret==EXIT_SUCCESS){
    update_intra_xforms(&ctx);
    ret=apply_to_blocks(&ctx,intra_xform_update_plane_start,
     intra_xform_update_block,intra_xform_update_plane_finish,_argc,_argv);
  }
  for(i=0;i<OD_INTRA_NMODES;i++)
  {
    int j;
    for(j=0;j<NB_CONTEXTS;j++)
      printf("%f ", ctx.freq[i][j][1]/(float)ctx.freq[i][j][0]);
    printf("\n");
  }
#if 0
  for(i=0;i<B_SZ*B_SZ;i++)
  {
    int j;
    printf("hist: ");
    for(j=0;j<8192;j++)
      printf("%lli ", hist[i][j]);
    printf("\n");
  }
#endif
  fprintf(stderr, "Average cost: %f bits/block, satd+cost: %f\n", ctx.total_bits/ctx.total_count, ctx.total_satd/ctx.total_count);
  return ret;
}
