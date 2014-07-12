/*Daala video codec
Copyright (c) 2012 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include "intra_fit_tools.h"
#include "svd.h"
#include "cholesky.h"
#include "../src/dct.h"
#include "../src/intra.h"

#define INTRA_NO_RDO (1)

int ExCount[3];
double Ex[3][B_SZ*B_SZ];
#define P0_COEF 0.002

#if !INTRA_NO_RDO
/* These weight are estimated from the entropy of the residual (in one older run) on subset1-y4m */
const float satd_weights[3][B_SZ*B_SZ] =
{{0.053046, 0.108601, 0.208930, 0.267489, 0.099610, 0.151836, 0.276093, 0.333044,
  0.180490, 0.256407, 0.420637, 0.485639, 0.224520, 0.295232, 0.464384, 0.519662,},
 {0.209979, 0.381558, 0.657913, 0.793790, 0.362177, 0.505781, 0.821064, 0.926135,
  0.603069, 0.779381, 1.078780, 1.151472, 0.725735, 0.864047, 1.135539, 1.175140,},
 {0.231557, 0.398671, 0.687661, 0.820069, 0.379319, 0.542037, 0.865300, 0.967672,
  0.643793, 0.824554, 1.124411, 1.191774, 0.760804, 0.909100, 1.181837, 1.223243,}
};

# if B_SZ==4
/* Less extreme weighting -- sqrt of the weights above */
const float satd_weights2[3][B_SZ*B_SZ] =
{{0.230317, 0.329547, 0.457088, 0.517193, 0.315611, 0.389662, 0.525445, 0.577099,
  0.424841, 0.506366, 0.648566, 0.696878, 0.473835, 0.543352, 0.681457, 0.720876,},
 {0.458235, 0.617704, 0.811118, 0.890949, 0.601811, 0.711183, 0.906126, 0.962359,
  0.776575, 0.882826, 1.038644, 1.073067, 0.851901, 0.929541, 1.065617, 1.084039,},
 {0.481204, 0.631404, 0.829253, 0.905577, 0.615889, 0.736232, 0.930215, 0.983703,
  0.802367, 0.908049, 1.060383, 1.091684, 0.872241, 0.953467, 1.087123, 1.106003,}
};
# else
#  error "No weights for B_SZ!=4 yet. #define INTRA_NO_RDO and try again."
# endif
#endif

#if 1
#define NB_CONTEXTS 8
#define GET_CONTEXT(modes,pos,m,width) (((modes)[(pos)-1]==(m))*4 + ((modes)[(pos)-(width)-1]==(m))*2 + ((modes)[(pos)-(width)]==(m))*1)
#else
#define NB_CONTEXTS 1000
#define GET_CONTEXT(modes,pos,m,width) (((modes)[(pos)-1])*100 + ((modes)[(pos)-(width)-1])*10 + ((modes)[(pos)-(width)]==(m))*1)
#endif

typedef struct intra_xform_ctx intra_xform_ctx;

struct intra_xform_ctx{
  char          *map_filename;
  unsigned char *map;
  char          *weights_filename;
  unsigned      *weights;
  int            nxblocks;
  int            nyblocks;
  int            pli;
  double         r_w[OD_INTRA_NMODES];
  double         r_x[OD_INTRA_NMODES][2*B_SZ*2*B_SZ];
  double         r_xx[OD_INTRA_NMODES][2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];
  double         scale[OD_INTRA_NMODES][2*B_SZ*2*B_SZ];
  double         beta[OD_INTRA_NMODES][B_SZ*B_SZ][2*B_SZ*2*B_SZ];
  double         freq[3][OD_INTRA_NMODES][NB_CONTEXTS][2];
  double         p0[OD_INTRA_NMODES];
  long long      n;
  double         satd_avg;
  double         total_bits;
  double         total_satd;
  double         total_count;
};

#define PRINT_BLOCKS (0)

static int intra_xform_train_plane_start(void *_ctx,const char *_name,
 const video_input_info *_info,int _pli,int _nxblocks,int _nyblocks){
  intra_xform_ctx *ctx;
  FILE            *map_file;
  char            *map_filename;
  FILE            *weights_file;
  char            *weights_filename;
  (void)_info;
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
#if PRINT_BLOCKS
  fprintf(stderr,"%i %i\n",_nxblocks,_nyblocks);
#endif
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  ctx->pli=_pli;
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
# if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
        (*OD_PRE_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(col,col);
# else
#  error "Need a prefilter implementation for this block size."
# endif
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
        (*OD_PRE_FILTER[B_SZ_LOG-OD_LOG_BSIZE0])(row,row);
      }
    }
  }
#endif
  buf2=_buf+3*B_SZ*(B_SZ>>1)+(B_SZ>>1);
#if APPLY_DCT
  for(by=0;by<2;by++){
    for(bx=0;bx<2;bx++){
# if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
      (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(buf2+B_SZ*(3*B_SZ*by+bx),3*B_SZ,
       buf2+B_SZ*(3*B_SZ*by+bx),3*B_SZ);
# else
#  error "Need an fDCT implementation for this block size."
# endif
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
  buf2=xform_blocks(buf,_data,_stride);
  mode=ctx->map[_bj*ctx->nxblocks+_bi];
#if PRINT_BLOCKS
  fprintf(stderr,"%i",mode);
  for (i=0;i<2*B_SZ;i++) {
    for (j=0;j<2*B_SZ;j++) {
      fprintf(stderr," %i",buf2[3*B_SZ*i+j]);
    }
  }
  fprintf(stderr,"\n");
#endif
  wb=ctx->weights[_bj*ctx->nxblocks+_bi];
  if(wb<=0)return;
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
      ctx->freq[ctx->pli][m][c][0]+=1;
      ctx->freq[ctx->pli][m][c][1] += (mode==m);
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
  int pli;
  /*Update the model for each coefficient in each mode.*/
  printf("/* This file is generated automatically by init_intra_xform */\n");

  printf("#include \"intra.h\"\n");
  printf("\n");
  printf("const double OD_INTRA_PRED_WEIGHTS_%ix%i"
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
        double  xty[2*B_SZ*2*B_SZ];
        double *beta;
        int     xii;
        int     xij;
        int     yi;
        nxi=3*B_SZ*B_SZ;
#if 0
        /*Include coefficients for the current block*/
        {
          int k;
          int l;
          for(k=0;k<=i;k++){
            for(l=0;l<=j;l++){
              xi[nxi++]=2*B_SZ*(B_SZ+k)+B_SZ+l;
            }
          }
          nxi--;
        }
#endif
        yi=2*B_SZ*(B_SZ+i)+B_SZ+j;
        for(xii=0;xii<nxi;xii++)xty[xii]=r_xx[xi[xii]][yi];
        beta=_ctx->beta[mode][B_SZ*i+j];
        memset(beta,0,2*B_SZ*2*B_SZ*sizeof(*beta));
#if defined(OD_USE_SVD)
        {
          double  xtx[2*2*B_SZ*2*B_SZ][2*B_SZ*2*B_SZ];
          double *xtxp[2*2*B_SZ*2*B_SZ];
          double  s[2*B_SZ*2*B_SZ];
          for(xii=0;xii<nxi;xii++){
            for(xij=0;xij<nxi;xij++){
              xtx[xii][xij]=r_xx[xi[xii]][xi[xij]];
            }
          }
          for(xii=0;xii<2*nxi;xii++)xtxp[xii]=xtx[xii];
          svd_pseudoinverse(xtxp,s,nxi,nxi);
          /*beta[yi]=r_x[yi];*/
          for(xii=0;xii<nxi;xii++){
            double beta_i;
            beta_i=0;
            for(xij=0;xij<nxi;xij++)beta_i+=xtx[xij][xii]*xty[xij];
            beta[xi[xii]]=beta_i*scale[yi]/scale[xi[xii]];
            /*beta[yi]-=beta_i*r_x[xi[xii]];*/
          }
        }
#else
        {
          double  xtx[UT_SZ(2*B_SZ*2*B_SZ,2*B_SZ*2*B_SZ)];
          double  tau[2*B_SZ*2*B_SZ];
          double  work[2*B_SZ*2*B_SZ];
          int     pivot[2*B_SZ*2*B_SZ];
          int     rank;
          for(xii=0;xii<nxi;xii++){
            for(xij=xii;xij<nxi;xij++){
              xtx[UT_IDX(xii,xij,nxi)]=r_xx[xi[xii]][xi[xij]];
            }
          }
          rank=cholesky(xtx,pivot,DBL_EPSILON,nxi);
          chdecomp(xtx,tau,rank,nxi);
          chsolve(xtx,pivot,tau,xty,xty,work,rank,nxi);
          for(xii=0;xii<nxi;xii++){
            beta[xi[xii]]=xty[xii]*scale[yi]/scale[xi[xii]];
            /*beta[yi]-=beta_i*r_x[xi[xii]];*/
          }
        }
#endif
        print_beta(mode,i,j,beta);
      }
      printf("    }%s\n",i<B_SZ-1?",":"");
    }
    printf("  }%s\n",mode<OD_INTRA_NMODES-1?",":"");
  }
  printf("};\n\n");

  printf("const unsigned char OD_INTRA_PRED_PROB_%dx%d[3][OD_INTRA_NMODES][OD_INTRA_NCONTEXTS]={\n",B_SZ,B_SZ);
  for(pli=0;pli<3;pli++)
  {
    int i;
    printf("{");
    for(i=0;i<OD_INTRA_NMODES;i++)
    {
      int j;
      printf("{");
      for(j=0;j<NB_CONTEXTS;j++)
        printf("%d, ", (int)floor(.5+256.*_ctx->freq[pli][i][j][1]/(float)_ctx->freq[pli][i][j][0]));
      printf("},\n");
    }
    printf("},\n");
  }
  printf("};\n\n");
}


static int intra_xform_update_plane_start(void *_ctx,const char *_name,
 const video_input_info *_info,int _pli,int _nxblocks,int _nyblocks){
  intra_xform_ctx *ctx;
  int i;
  (void)_info;
  ctx=(intra_xform_ctx *)_ctx;
  ctx->map_filename=get_map_filename(_name,_pli,_nxblocks,_nyblocks);
  ctx->weights_filename=get_weights_filename(_name,_pli,_nxblocks,_nyblocks);
  ctx->map=(unsigned char *)malloc(_nxblocks*(size_t)_nyblocks);
  ctx->weights=(unsigned *)malloc(
   _nxblocks*(size_t)_nyblocks*sizeof(*ctx->weights));
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  ctx->pli=_pli;
  for(i=0;i<OD_INTRA_NMODES;i++)
    ctx->p0[i]=ctx->freq[ctx->pli][i][0][1]/(float)ctx->freq[ctx->pli][i][0][0];
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
  double           best_error[B_SZ*B_SZ]={0};
  double           error[B_SZ*B_SZ];
  int              mode;
  int              best_mode;
  int              c0[OD_INTRA_NMODES]={0};
  double           bits=0;
  int              i;
  float sum=0;
  float sum2=0;
  unsigned char *modes;
  int pos;
  int m;
  int width;
  float p[OD_INTRA_NMODES];
  ogg_uint16_t cdf[OD_INTRA_NMODES];
/*If using this be sure to uncomment its assignment.*/
/*  ogg_uint32_t wsatd[OD_INTRA_NMODES];*/

  ctx=(intra_xform_ctx *)_ctx;
  modes=ctx->map;
  pos = _bj*ctx->nxblocks+_bi;
  width=ctx->nxblocks;

  buf2=xform_blocks(buf,_data,_stride);
  best_mode=0;
  best_satd=UINT_MAX;
  best_rlsatd=UINT_MAX;
  best_bits=0;
  next_best_satd=UINT_MAX;
  next_best_rlsatd=UINT_MAX;

  {
    int c;
    unsigned char probs[OD_INTRA_NMODES][OD_INTRA_NCONTEXTS];
    int left;
    int upleft;
    int up;

    left=(_bi==0)?0:modes[pos-1];
    up=(_bj==0)?0:modes[pos-width];
    upleft=(_bi==0||_bj==0)?0:modes[pos-width-1];
    for (m=0;m<OD_INTRA_NMODES;m++)
      for(c=0;c<OD_INTRA_NCONTEXTS;c++)
        probs[m][c] = 256.*ctx->freq[ctx->pli][m][c][1]/(float)ctx->freq[ctx->pli][m][c][0];
    od_intra_pred_cdf(cdf,probs,OD_INTRA_NMODES,left,upleft,up);
  }
  for(m=0;m<OD_INTRA_NMODES;m++)
  {
    int c;
    if (_bi>0 && _bj>0)
      c=GET_CONTEXT(modes,pos,m,width);
    else if (_bj>0)
#if 1
      c=(modes[pos-width]==m)+6*(m==0);
#else
    c=(modes[pos-width]);
#endif
    else if (_bi>0)
#if 1
      c=(modes[pos-1]==m)*4+3*(m==0);
#else
    c=(modes[pos-1])*100;
#endif
    else
      c=15*(m==0);
    p[m] = ctx->freq[ctx->pli][m][c][1]/(float)ctx->freq[ctx->pli][m][c][0];
    p[m]+=1.e-5;
    if (p[m]<ctx->p0[m]) p[m]=ctx->p0[m];
    if (c==0)
    {
      ctx->p0[m]*=(1-P0_COEF);
      c0[m]=1;
    }
    sum += p[m];
  }
  for(m=0;m<OD_INTRA_NMODES;m++)
  {
    p[m] *= (sum-p[m])/(1-p[m]);
    sum2+=p[m];
  }

  for(mode=0;mode<OD_INTRA_NMODES;mode++){
    double satd;
    double rlsatd;
    double diff;
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
        diff = fabs(buf2[3*B_SZ*(i+B_SZ)+j+B_SZ]-p) - 1;
        if (diff<0)diff=0;
        /*satd+=satd_weights2[i*B_SZ+j]*diff;*/
        satd+=satd_weights2[ctx->pli][i*B_SZ+j]*diff;
#endif
        error[i*B_SZ+j]=buf2[3*B_SZ*(i+B_SZ)+j+B_SZ]-p;
      }
    }
/*    wsatd[mode] = satd*64;*/
    rlsatd=satd;
    /* Normalize all probabilities except the max */
    /*bits = -log(p[mode]/sum)/log(2);*/
    /* Normalize all probabilities except the max */
    /*bits = (mode==maxM) ? -log(p[mode])/log(2) : -log(p[mode]*(1-maxP)/(sum-maxP))/log(2);*/

    bits = -log(p[mode]/sum2)/log(2);
    /*printf("{%f+l*%f= ", satd , bits);*/

#ifndef INTRA_NO_RDO
    satd += 1.1*bits;
#endif
    /*printf("%f} ", satd);*/
    /* Bias towards DC mode */
    /*if (mode==0)satd-=.5;*/
    if(satd<best_satd){
      next_best_satd=best_satd;
      next_best_rlsatd=best_rlsatd;
      best_satd=satd;
      best_rlsatd=rlsatd;
      best_mode=mode;
      best_bits=bits;
      for (i=0;i<B_SZ*B_SZ;i++)
        best_error[i]=error[i];
    }
    else if(satd<next_best_satd){next_best_satd=satd;next_best_rlsatd=rlsatd;}
  }
  /*printf("\n");*/
#if 0
  {
    int bmode;
    int left, up, upleft;

    left=(_bi==0)?0:modes[pos-1];
    up=(_bj==0)?0:modes[pos-width];
    upleft=(_bi==0||_bj==0)?0:modes[pos-width-1];
    bmode=od_intra_pred_search(cdf,wsatd,OD_INTRA_NMODES,64*1.1);
    od_intra_pred_update(p0,OD_INTRA_NMODES,bmode,left,upleft,up);
    /*if (bmode==best_mode)
      printf("+");
    else
      printf("-");
    printf("%d %d\n", bmode, best_mode);*/
  }
#endif
  for (i=0;i<B_SZ*B_SZ;i++)
    Ex[ctx->pli][i]+=fabs(best_error[i]);
  ExCount[ctx->pli]++;
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
  /*printf("Average SATD: %G\n",ctx->satd_avg);*/
  return intra_xform_train_plane_finish(_ctx);
}


int main(int _argc,const char **_argv){
  static intra_xform_ctx ctx;
  int                    ret;
  int                    i;
  int                    pli;
  for(pli=0;pli<3;pli++)
  {
    for(i=0;i<OD_INTRA_NMODES;i++)
    {
      int j;
      for(j=0;j<NB_CONTEXTS;j++)
      {
        ctx.freq[pli][i][j][0]=2;
        ctx.freq[pli][i][j][1]=1;
      }
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
#if 0
  for (pli=0;pli<3;pli++)
  {
    printf("Ex: ");
    for(i=0;i<B_SZ*B_SZ;i++)
    {
      printf("%f ", Ex[pli][i]/ExCount[pli]);
    }
    printf("\n");
  }
#endif
  fprintf(stderr, "Average cost: %f bits/block, satd+cost: %f\n", ctx.total_bits/ctx.total_count, ctx.total_satd/ctx.total_count);
  return ret;
}
