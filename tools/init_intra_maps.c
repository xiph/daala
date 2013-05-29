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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include "intra_fit_tools.h"
#include "../src/newdct.c"
#include "../src/internal.c"

/*For validation purposes only.
  Copied from libvpx.*/
#if B_SZ==4
#define B_DC_PRED (OD_INTRA_DC)
#define B_TM_PRED (OD_INTRA_TM)
#define B_VE_PRED (OD_INTRA_VE)
#define B_HE_PRED (OD_INTRA_HE)
#define B_LD_PRED (OD_INTRA_LD)
#define B_RD_PRED (OD_INTRA_RD)
#define B_VR_PRED (OD_INTRA_VR)
#define B_VL_PRED (OD_INTRA_VL)
#define B_HU_PRED (OD_INTRA_HU)
#define B_HD_PRED (OD_INTRA_HD)

void vp8_intra4x4_predict_c(const unsigned char *src, int src_stride,
                            int b_mode,
                            unsigned char *dst, int dst_stride)
{
    int i, r, c;
    const unsigned char *Above = src - src_stride;
    unsigned char Left[4];
    unsigned char top_left = Above[-1];
    Left[0] = src[-1];
    Left[1] = src[-1 + src_stride];
    Left[2] = src[-1 + 2 * src_stride];
    Left[3] = src[-1 + 3 * src_stride];
    switch (b_mode)
    {
    case B_DC_PRED:
    {
        int expected_dc = 0;
        for (i = 0; i < 4; i++)
        {
            expected_dc += Above[i];
            expected_dc += Left[i];
        }
        expected_dc = (expected_dc + 4) >> 3;
        for (r = 0; r < 4; r++)
        {
            for (c = 0; c < 4; c++)
            {
                dst[c] = expected_dc;
            }
            dst += dst_stride;
        }
    }
    break;
    case B_TM_PRED:
    {
        /* prediction similar to true_motion prediction */
        for (r = 0; r < 4; r++)
        {
            for (c = 0; c < 4; c++)
            {
                int pred = Above[c] - top_left + Left[r];
                if (pred < 0)
                    pred = 0;
                if (pred > 255)
                    pred = 255;
                dst[c] = pred;
            }
            dst += dst_stride;
        }
    }
    break;
    case B_VE_PRED:
    {
        unsigned int ap[4];
        ap[0] = (top_left  + 2 * Above[0] + Above[1] + 2) >> 2;
        ap[1] = (Above[0] + 2 * Above[1] + Above[2] + 2) >> 2;
        ap[2] = (Above[1] + 2 * Above[2] + Above[3] + 2) >> 2;
        ap[3] = (Above[2] + 2 * Above[3] + Above[4] + 2) >> 2;
        for (r = 0; r < 4; r++)
        {
            for (c = 0; c < 4; c++)
            {
                dst[c] = ap[c];
            }
            dst += dst_stride;
        }
    }
    break;
    case B_HE_PRED:
    {
        unsigned int lp[4];
        lp[0] = (top_left + 2 * Left[0] + Left[1] + 2) >> 2;
        lp[1] = (Left[0] + 2 * Left[1] + Left[2] + 2) >> 2;
        lp[2] = (Left[1] + 2 * Left[2] + Left[3] + 2) >> 2;
        lp[3] = (Left[2] + 2 * Left[3] + Left[3] + 2) >> 2;
        for (r = 0; r < 4; r++)
        {
            for (c = 0; c < 4; c++)
            {
                dst[c] = lp[r];
            }
            dst += dst_stride;
        }
    }
    break;
    case B_LD_PRED:
    {
        const unsigned char *ptr = Above;
        dst[0 * dst_stride + 0] = (ptr[0] + ptr[1] * 2 + ptr[2] + 2) >> 2;
        dst[0 * dst_stride + 1] =
            dst[1 * dst_stride + 0] = (ptr[1] + ptr[2] * 2 + ptr[3] + 2) >> 2;
        dst[0 * dst_stride + 2] =
            dst[1 * dst_stride + 1] =
                dst[2 * dst_stride + 0] = (ptr[2] + ptr[3] * 2 + ptr[4] + 2) >> 2;
        dst[0 * dst_stride + 3] =
            dst[1 * dst_stride + 2] =
                dst[2 * dst_stride + 1] =
                    dst[3 * dst_stride + 0] = (ptr[3] + ptr[4] * 2 + ptr[5] + 2) >> 2;
        dst[1 * dst_stride + 3] =
            dst[2 * dst_stride + 2] =
                dst[3 * dst_stride + 1] = (ptr[4] + ptr[5] * 2 + ptr[6] + 2) >> 2;
        dst[2 * dst_stride + 3] =
            dst[3 * dst_stride + 2] = (ptr[5] + ptr[6] * 2 + ptr[7] + 2) >> 2;
        dst[3 * dst_stride + 3] = (ptr[6] + ptr[7] * 2 + ptr[7] + 2) >> 2;
    }
    break;
    case B_RD_PRED:
    {
        unsigned char pp[9];
        pp[0] = Left[3];
        pp[1] = Left[2];
        pp[2] = Left[1];
        pp[3] = Left[0];
        pp[4] = top_left;
        pp[5] = Above[0];
        pp[6] = Above[1];
        pp[7] = Above[2];
        pp[8] = Above[3];
        dst[3 * dst_stride + 0] = (pp[0] + pp[1] * 2 + pp[2] + 2) >> 2;
        dst[3 * dst_stride + 1] =
            dst[2 * dst_stride + 0] = (pp[1] + pp[2] * 2 + pp[3] + 2) >> 2;
        dst[3 * dst_stride + 2] =
            dst[2 * dst_stride + 1] =
                dst[1 * dst_stride + 0] = (pp[2] + pp[3] * 2 + pp[4] + 2) >> 2;
        dst[3 * dst_stride + 3] =
            dst[2 * dst_stride + 2] =
                dst[1 * dst_stride + 1] =
                    dst[0 * dst_stride + 0] = (pp[3] + pp[4] * 2 + pp[5] + 2) >> 2;
        dst[2 * dst_stride + 3] =
            dst[1 * dst_stride + 2] =
                dst[0 * dst_stride + 1] = (pp[4] + pp[5] * 2 + pp[6] + 2) >> 2;
        dst[1 * dst_stride + 3] =
            dst[0 * dst_stride + 2] = (pp[5] + pp[6] * 2 + pp[7] + 2) >> 2;
        dst[0 * dst_stride + 3] = (pp[6] + pp[7] * 2 + pp[8] + 2) >> 2;
    }
    break;
    case B_VR_PRED:
    {
        unsigned char pp[9];
        pp[0] = Left[3];
        pp[1] = Left[2];
        pp[2] = Left[1];
        pp[3] = Left[0];
        pp[4] = top_left;
        pp[5] = Above[0];
        pp[6] = Above[1];
        pp[7] = Above[2];
        pp[8] = Above[3];
        dst[3 * dst_stride + 0] = (pp[1] + pp[2] * 2 + pp[3] + 2) >> 2;
        dst[2 * dst_stride + 0] = (pp[2] + pp[3] * 2 + pp[4] + 2) >> 2;
        dst[3 * dst_stride + 1] =
            dst[1 * dst_stride + 0] = (pp[3] + pp[4] * 2 + pp[5] + 2) >> 2;
        dst[2 * dst_stride + 1] =
            dst[0 * dst_stride + 0] = (pp[4] + pp[5] + 1) >> 1;
        dst[3 * dst_stride + 2] =
            dst[1 * dst_stride + 1] = (pp[4] + pp[5] * 2 + pp[6] + 2) >> 2;
        dst[2 * dst_stride + 2] =
            dst[0 * dst_stride + 1] = (pp[5] + pp[6] + 1) >> 1;
        dst[3 * dst_stride + 3] =
            dst[1 * dst_stride + 2] = (pp[5] + pp[6] * 2 + pp[7] + 2) >> 2;
        dst[2 * dst_stride + 3] =
            dst[0 * dst_stride + 2] = (pp[6] + pp[7] + 1) >> 1;
        dst[1 * dst_stride + 3] = (pp[6] + pp[7] * 2 + pp[8] + 2) >> 2;
        dst[0 * dst_stride + 3] = (pp[7] + pp[8] + 1) >> 1;
    }
    break;
    case B_VL_PRED:
    {
        const unsigned char *pp = Above;
        dst[0 * dst_stride + 0] = (pp[0] + pp[1] + 1) >> 1;
        dst[1 * dst_stride + 0] = (pp[0] + pp[1] * 2 + pp[2] + 2) >> 2;
        dst[2 * dst_stride + 0] =
            dst[0 * dst_stride + 1] = (pp[1] + pp[2] + 1) >> 1;
        dst[1 * dst_stride + 1] =
            dst[3 * dst_stride + 0] = (pp[1] + pp[2] * 2 + pp[3] + 2) >> 2;
        dst[2 * dst_stride + 1] =
            dst[0 * dst_stride + 2] = (pp[2] + pp[3] + 1) >> 1;
        dst[3 * dst_stride + 1] =
            dst[1 * dst_stride + 2] = (pp[2] + pp[3] * 2 + pp[4] + 2) >> 2;
        dst[0 * dst_stride + 3] =
            dst[2 * dst_stride + 2] = (pp[3] + pp[4] + 1) >> 1;
        dst[1 * dst_stride + 3] =
            dst[3 * dst_stride + 2] = (pp[3] + pp[4] * 2 + pp[5] + 2) >> 2;
        dst[2 * dst_stride + 3] = (pp[4] + pp[5] + 1) >> 1;
        dst[3 * dst_stride + 3] = (pp[4] + pp[5] * 2 + pp[6] + 2) >> 2;
    }
    break;
    case B_HD_PRED:
    {
        unsigned char pp[9];
        pp[0] = Left[3];
        pp[1] = Left[2];
        pp[2] = Left[1];
        pp[3] = Left[0];
        pp[4] = top_left;
        pp[5] = Above[0];
        pp[6] = Above[1];
        pp[7] = Above[2];
        pp[8] = Above[3];
        dst[3 * dst_stride + 0] = (pp[0] + pp[1] + 1) >> 1;
        dst[3 * dst_stride + 1] = (pp[0] + pp[1] * 2 + pp[2] + 2) >> 2;
        dst[2 * dst_stride + 0] =
            dst[3 * dst_stride + 2] = (pp[1] + pp[2] + 1) >> 1;
        dst[2 * dst_stride + 1] =
            dst[3 * dst_stride + 3] = (pp[1] + pp[2] * 2 + pp[3] + 2) >> 2;
        dst[2 * dst_stride + 2] =
            dst[1 * dst_stride + 0] = (pp[2] + pp[3] + 1) >> 1;
        dst[2 * dst_stride + 3] =
            dst[1 * dst_stride + 1] = (pp[2] + pp[3] * 2 + pp[4] + 2) >> 2;
        dst[1 * dst_stride + 2] =
            dst[0 * dst_stride + 0] = (pp[3] + pp[4] + 1) >> 1;
        dst[1 * dst_stride + 3] =
            dst[0 * dst_stride + 1] = (pp[3] + pp[4] * 2 + pp[5] + 2) >> 2;
        dst[0 * dst_stride + 2] = (pp[4] + pp[5] * 2 + pp[6] + 2) >> 2;
        dst[0 * dst_stride + 3] = (pp[5] + pp[6] * 2 + pp[7] + 2) >> 2;
    }
    break;
    case B_HU_PRED:
    {
        unsigned char *pp = Left;
        dst[0 * dst_stride + 0] = (pp[0] + pp[1] + 1) >> 1;
        dst[0 * dst_stride + 1] = (pp[0] + pp[1] * 2 + pp[2] + 2) >> 2;
        dst[0 * dst_stride + 2] =
            dst[1 * dst_stride + 0] = (pp[1] + pp[2] + 1) >> 1;
        dst[0 * dst_stride + 3] =
            dst[1 * dst_stride + 1] = (pp[1] + pp[2] * 2 + pp[3] + 2) >> 2;
        dst[1 * dst_stride + 2] =
            dst[2 * dst_stride + 0] = (pp[2] + pp[3] + 1) >> 1;
        dst[1 * dst_stride + 3] =
            dst[2 * dst_stride + 1] = (pp[2] + pp[3] * 2 + pp[3] + 2) >> 2;
        dst[2 * dst_stride + 2] =
            dst[2 * dst_stride + 3] =
                dst[3 * dst_stride + 0] =
                    dst[3 * dst_stride + 1] =
                        dst[3 * dst_stride + 2] =
                            dst[3 * dst_stride + 3] = pp[3];
    }
    break;
    }
}
#endif

#define USE_DCT_SATD (0)

#if !(USE_DCT_SATD&&B_SZ==4)
static void od_diff_hadamard(short _buf[B_SZ*B_SZ],
 const unsigned char *_src,int _src_stride,
 const unsigned char *_ref,int _ref_stride){
  int i;
  for(i=0;i<B_SZ;i++){
    int t[B_SZ];
    int r;
    int j;
    int k;
    int l;
    /*Hadamard stage 1:*/
    for(j=0;j<(B_SZ>>1);j++){
      t[j]=_src[j]-_ref[j]+_src[j+(B_SZ>>1)]-_ref[j+(B_SZ>>1)];
      t[j+(B_SZ>>1)]=_src[j]-_ref[j]-_src[j+(B_SZ>>1)]+_ref[j+(B_SZ>>1)];
    }
    /*Hadamard stages 2...N-1:*/
    for(l=B_SZ>>2;l>1;l>>=1){
      for(k=0;k<B_SZ;k+=2*l){
        for(j=0;j<l;j++){
          r=t[j+k];
          t[j+k]=r+t[j+k+l];
          t[j+k+l]=r-t[j+k+l];
        }
      }
    }
    /*Hadamard stage N:*/
    for(k=0;k<B_SZ;k+=2){
      _buf[k*B_SZ+i]=(short)(t[k]+t[k+1]);
      _buf[(k+1)*B_SZ+i]=(short)(t[k]-t[k+1]);
    }
    _src+=_src_stride;
    _ref+=_ref_stride;
  }
}

static unsigned od_hadamard_sad(const short _buf[B_SZ*B_SZ]){
  unsigned sad;
  int      i;
  sad=0;
  for(i=0;i<B_SZ;i++){
    int t[B_SZ];
    int j;
    int k;
    int l;
    int r;
    /*Hadamard stage 1:*/
    for(j=0;j<(B_SZ>>1);j++){
      t[j]=_buf[i*B_SZ+j]+_buf[i*B_SZ+j+(B_SZ>>1)];
      t[j+(B_SZ>>1)]=_buf[i*B_SZ+j]-_buf[i*B_SZ+j+(B_SZ>>1)];
    }
    /*Hadamard stages 2...N-1:*/
    for(l=B_SZ>>2;l>1;l>>=1){
      for(k=0;k<B_SZ;k+=2*l){
        for(j=0;j<l;j++){
          r=t[j+k];
          t[j+k]=r+t[j+k+l];
          t[j+k+l]=r-t[j+k+l];
        }
      }
    }
    /*Hadamard stage N:*/
    for(k=0;k<B_SZ;k+=2){
      sad+=abs(t[k]+t[k+1]);
      sad+=abs(t[k]-t[k+1]);
    }
  }
  return sad;
}

static unsigned od_satd(const unsigned char *_src,int _src_stride,
 const unsigned char *_ref,int _ref_stride){
  short buf[B_SZ*B_SZ];
  od_diff_hadamard(buf,_src,_src_stride,_ref,_ref_stride);
  return od_hadamard_sad(buf);
}
#else

/*Computes the SATD using the actual DCT, instead of the Hadamard transform.
  This is useful for comparing SATD numbers to the frequency domain techniques
   without having to worry about scaling.*/
static unsigned od_satd(const unsigned char *_src,int _src_stride,
 const unsigned char *_ref,int _ref_stride){
  od_coeff buf[B_SZ*B_SZ];
  unsigned satd;
  int      i;
  int      j;
  OD_ASSERT(B_SZ==4);
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++){
      buf[B_SZ*i+j]=_src[j]-_ref[j];
    }
    _src+=_src_stride;
    _ref+=_ref_stride;
  }
#if B_SZ_LOG>=OD_LOG_BSIZE0&&B_SZ_LOG<OD_LOG_BSIZE0+OD_NBSIZES
  (*OD_FDCT_2D[B_SZ_LOG-OD_LOG_BSIZE0])(buf,B_SZ,buf,B_SZ);
#else
# error "Need an fDCT implementation for this block size."
#endif
  satd=0;
  for(i=0;i<B_SZ;i++){
    for(j=0;j<B_SZ;j++)satd+=abs(buf[B_SZ*i+j]);
  }
  return satd;
}
#endif



typedef struct init_intra_maps_ctx init_intra_maps_ctx;



struct init_intra_maps_ctx{
  const char    *name;
  unsigned char *map;
  unsigned      *weights;
  int            pli;
  int            nxblocks;
  int            nyblocks;
  long long      n;
  double         satd_avg;
};



static int init_intra_plane_start(void *_ctx,const char *_name,
 const th_info *_ti,int _pli,int _nxblocks,int _nyblocks){
  init_intra_maps_ctx *ctx;
  (void)_ti;
  ctx=(init_intra_maps_ctx *)_ctx;
  ctx->name=_name;
  ctx->map=(unsigned char *)malloc(_nxblocks*(size_t)_nyblocks);
  ctx->weights=(unsigned *)malloc(
   _nxblocks*(size_t)_nyblocks*sizeof(*ctx->weights));
  ctx->pli=_pli;
  ctx->nxblocks=_nxblocks;
  ctx->nyblocks=_nyblocks;
  return EXIT_SUCCESS;
}

static void init_intra_block(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj){
  init_intra_maps_ctx *ctx;
  unsigned             best_satd;
  unsigned             next_best_satd;
  int                  mode;
  int                  best_mode;
  ctx=(init_intra_maps_ctx *)_ctx;
  best_mode=0;
  best_satd=UINT_MAX;
  next_best_satd=UINT_MAX;
  for(mode=0;mode<10;mode++){
    unsigned char        block[B_SZ*B_SZ];
    unsigned             satd;
    memset(block,0,B_SZ*B_SZ);
    vp8_intra_predict(block,B_SZ,_data,_stride,mode);
#if B_SZ==4
    {
      unsigned char block2[B_SZ*B_SZ];
      vp8_intra4x4_predict_c(_data,_stride,mode,block2,B_SZ);
      if(memcmp(block,block2,sizeof(block2))!=0){
        fprintf(stderr,"Intra prediction mismatch.\n");
        exit(EXIT_FAILURE);
      }
    }
#endif
    satd=od_satd(block,B_SZ,_data,_stride);
    if(satd<best_satd){
      next_best_satd=best_satd;
      best_satd=satd;
      best_mode=mode;
    }
    else if(satd<next_best_satd)next_best_satd=satd;
  }
  ctx->satd_avg+=(best_satd-ctx->satd_avg)/++(ctx->n);
  ctx->map[_bj*ctx->nxblocks+_bi]=(unsigned char)best_mode;
  ctx->weights[_bj*ctx->nxblocks+_bi]=next_best_satd-best_satd;
}

static int init_intra_plane_finish(void *_ctx){
  init_intra_maps_ctx *ctx;
  FILE                *map_file;
  char                *map_filename;
  FILE                *weights_file;
  char                *weights_filename;
  ctx=(init_intra_maps_ctx *)_ctx;
  map_filename=get_map_filename(ctx->name,
   ctx->pli,ctx->nxblocks,ctx->nyblocks);
  map_file=fopen(map_filename,"wb");
  if(map_file==NULL){
    fprintf(stderr,"Error opening output file '%s'.\n",map_filename);
    return EXIT_FAILURE;
  }
  if(fwrite(ctx->map,ctx->nxblocks*(size_t)ctx->nyblocks,1,map_file)<1){
    fprintf(stderr,"Error writing to output file '%s'.\n",map_filename);
    return EXIT_FAILURE;
  }
  fclose(map_file);
  free(map_filename);
  weights_filename=get_weights_filename(ctx->name,
   ctx->pli,ctx->nxblocks,ctx->nyblocks);
  weights_file=fopen(weights_filename,"wb");
  if(weights_file==NULL){
    fprintf(stderr,"Error opening output file '%s'.\n",weights_filename);
    return EXIT_FAILURE;
  }
  if(fwrite(ctx->weights,
   ctx->nxblocks*(size_t)ctx->nyblocks*sizeof(*ctx->weights),1,
   weights_file)<1){
    fprintf(stderr,"Error writing to output file '%s'.\n",weights_filename);
    return EXIT_FAILURE;
  }
  fclose(weights_file);
  free(weights_filename);
  free(ctx->weights);
  free(ctx->map);
  printf("Average SATD: %G\n",ctx->satd_avg);
  return EXIT_SUCCESS;
}

int main(int _argc,const char **_argv){
  init_intra_maps_ctx ctx;
  ctx.n=0;
  ctx.satd_avg=0;
  return apply_to_blocks(&ctx,init_intra_plane_start,init_intra_block,
   init_intra_plane_finish,_argc,_argv);
}
