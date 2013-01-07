/*
    Daala video codec
    Copyright (C) 2012 Daala project contributors

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

*/

#include "../generic_code.c"
#include "../pvq_encoder.c"
#include "../generic_encoder.c"
#include "../pvq_decoder.c"
#include "../generic_decoder.c"
#include "../entenc.c"
#include "../entdec.c"
#include "../entcode.c"
#include "../cdf_table.c"
#include "../internal.c"
#include <stdlib.h>
#include <math.h>

#define EC_BUF_SIZE (32)
#define MAX_VECTORS 10000
#define MAXN 256

int run_pvq(int *X,int len,int N){
  od_pvq_adapt_ctx adapt;
  int i, j;
  int num, den;
  od_ec_enc enc;
  od_ec_dec dec;
  unsigned char *buf;
  ogg_uint32_t   buf_sz;
  int *Ki;
  GenericEncoder model;
  int EK=65536;
  int bits_used;

  Ki = malloc(sizeof(*Ki)*len);
  adapt.mean_k_q8=163;
  adapt.mean_sum_ex_q8=64;
#if !defined(OD_DISABLE_PVQ_CODE1)
  adapt.pos_q4=30<<4;
#endif
  num = 100*4;
  den = 256*4;
  od_ec_enc_init(&enc, EC_BUF_SIZE);
  generic_model_init(&model);
  for(i=0;i<len;i++){
    int K=0;
    for (j=0;j<N;j++)
      K += abs(X[i*N+j]);
    Ki[i] = K;
    generic_encode(&enc, &model, K, &EK, 4);
    if (K<15)
      pvq_encode_delta(&enc,&X[i*N],N,K,&num,&den);
    else{
      pvq_encoder(&enc,&X[i*N],N,K,&adapt);
      if(adapt.k>=0){
        adapt.mean_k_q8+=(adapt.k<<8)-adapt.mean_k_q8>>OD_K_ADAPT_SPEED;
      }
      if(adapt.sum_ex_q8>=0){
        adapt.mean_sum_ex_q8+=
         adapt.sum_ex_q8-adapt.mean_sum_ex_q8>>OD_SUM_EX_ADAPT_SPEED;
      }
#if !defined(OD_DISABLE_PVQ_CODE1)
      if(adapt.pos>=0){
        adapt.mean_pos_q4=OD_MAXI(adapt.mean_pos_q4+
         ((adapt.pos<<4)-adapt.mean_pos_q4>>OD_POS_ADAPT_SPEED),N/8);
      }
#endif
    }
    /*if (i==0)
    {
      fprintf(stderr, "enc: K[%d]=%d, num=%d, den=%d, u=%d\n", i, Ki[i], num, den, u);
    }*/
    OD_ASSERT(!enc.error);
  }
  buf = od_ec_enc_done(&enc, &buf_sz);

  bits_used = od_ec_enc_tell(&enc);

  adapt.mean_k_q8=163;
  adapt.mean_sum_ex_q8=64;
#if !defined(OD_DISABLE_PVQ_CODE1)
  adapt.pos_q4=30<<4;
#endif
  num = 100*4;
  den = 256*4;
  od_ec_dec_init(&dec, buf, buf_sz);
  generic_model_init(&model);
  EK=65536;

  for (i=0;i<len;i++)
  {
    int y[MAXN];
    int K;
    K=generic_decode(&dec, &model, &EK, 4);
    if (K!=Ki[i]){
      fprintf(stderr, "mismatch for K of vector %d (N=%d)\n", i, N);
    }
    if (Ki[i]<15)
      pvq_decode_delta(&dec, y, N, Ki[i], &num, &den);
    else{
      pvq_decoder(&dec, y, N, Ki[i], &adapt);
      if(adapt.k>=0){
        adapt.mean_k_q8+=(adapt.k<<8)-adapt.mean_k_q8>>OD_K_ADAPT_SPEED;
      }
      if(adapt.sum_ex_q8>=0){
        adapt.mean_sum_ex_q8+=
         adapt.sum_ex_q8-adapt.mean_sum_ex_q8>>OD_SUM_EX_ADAPT_SPEED;
      }
#if !defined(OD_DISABLE_PVQ_CODE1)
      if(adapt.pos>=0){
        adapt.mean_pos_q4=OD_MAXI(adapt.mean_pos_q4+
         ((adapt.pos<<4)-adapt.mean_pos_q4>>OD_POS_ADAPT_SPEED),N/8);
      }
#endif
    }
    for (j=0;j<N;j++){
      if(y[j]!=X[i*N+j]){
        int k;
        fprintf(stderr, "mismatch for coef %d of vector %d (N=%d)\n", j, i, N);
        fprintf(stderr, "orig vector:\n");
        for(k=0;k<N;k++)
          fprintf(stderr,"%d ", X[i*N+k]);
        fprintf(stderr,"\n");
        fprintf(stderr, "decoded vector:\n");
        for(k=0;k<N;k++)
          fprintf(stderr,"%d ", y[k]);
        fprintf(stderr,"\n");
#if !defined(OD_DISABLE_PVQ_CODE1)
        fprintf(stderr, "K[%d]=%d, num=%d, den=%d, u=%d\n", i, Ki[i],
         adapt.mean_k_q8, adapt.mean_sum_ex_q8, adapt.mean_pos_q4);
#else
        fprintf(stderr, "K[%d]=%d, num=%d, den=%d\n", i, Ki[i],
         adapt.mean_k_q8, adapt.mean_sum_ex_q8);
#endif
        abort();
      }
    }
  }

  free(Ki);
  od_ec_enc_clear(&enc);
  return bits_used;
}

void test_pvq_sequence(int len,int N,float std)
{
  int i,j;
  int bits;
  int *X;
  X = malloc(sizeof(*X)*len*N);
  for(i=0;i<len;i++){
    for(j=0;j<N;j++){
      X[i*N+j]=-floor(.5+std*exp(-j*3./N)*log((rand()+1)/(float)RAND_MAX));
      if (rand()&1)
        X[i*N+j]=-X[i*N+j];
      /*printf("%d ", X[i][j]);*/
    }
    /*printf("\n");*/
  }

  bits = run_pvq(X,len,N);

  fprintf(stderr, "Coded %d dim, std=%1.1f with %f bits/sample (%1.4f bits/vector)\n",N,std,bits/(float)len/N,bits/(float)len);

  free(X);
}


int main(int argc, char **argv){
  if(argc==4){
    int *X;
    int i,j;
    int len,N;
    int bits;
    FILE *file;
    len=atoi(argv[1]);
    if(len<1){
      fprintf(stderr, "len must be at least 1\n");
      return 1;
    }
    N=atoi(argv[2]);
    if(N<2){
      fprintf(stderr, "N must be at least 2\n");
      return 1;
    }
    X=malloc(sizeof(*X)*len*N);
    if(X==NULL){
      fprintf(stderr, "cannot allocate buffer\n");
      return 1;
    }
    file=fopen(argv[3],"r");
    if(file==NULL){
      fprintf(stderr, "cannot open input file\n");
      free(X);
      return 1;
    }
    for(i=0;i<len;i++)
      for(j=0;j<N;j++)
        fscanf(file,"%d ",&X[i*N+j]);
    bits = run_pvq(X,len,N);
    fprintf(stderr, "Coded file with %f bits/sample (%f bits/vector)\n",bits/(float)len/N,bits/(float)len);
    fclose(file);
    free(X);
  } else {
    test_pvq_sequence(10000,128,.03);
    test_pvq_sequence(10000,128,.1);
    test_pvq_sequence(10000,128,.3);
    test_pvq_sequence(10000,128,1.);
    test_pvq_sequence(10000,128,3);
    test_pvq_sequence(10000,128,10);
    test_pvq_sequence(10000,128,30);
    test_pvq_sequence(10000,128,100);
    test_pvq_sequence(10000,128,300);
    test_pvq_sequence(10000,128,1000);
    test_pvq_sequence(10000,128,3000);
    test_pvq_sequence(10000,128,6000);

    test_pvq_sequence(10000,16,.03);
    test_pvq_sequence(10000,16,.1);
    test_pvq_sequence(10000,16,.3);
    test_pvq_sequence(10000,16,1.);
    test_pvq_sequence(10000,16,3);
    test_pvq_sequence(10000,16,10);
    test_pvq_sequence(10000,16,30);
    test_pvq_sequence(10000,16,100);
    test_pvq_sequence(10000,16,300);
    test_pvq_sequence(10000,16,1000);
    test_pvq_sequence(10000,16,3000);
    test_pvq_sequence(10000,16,10000);
    test_pvq_sequence(10000,16,30000);
  }
  return 0;
}
