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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "../generic_code.c"
#include "../pvq_encoder.c"
#include "../generic_encoder.c"
#include "../pvq_decoder.c"
#include "../generic_decoder.c"
#include "../entenc.c"
#include "../entdec.c"
#include "../entcode.c"
#include "../laplace_tables.c"
#include "../internal.c"
#include <stdlib.h>
#include <math.h>

#define EC_BUF_SIZE (32)
#define MAX_VECTORS 10000
#define MAXN 256


typedef struct od_pvq_adapt_ctx od_pvq_adapt_ctx;

/* FIXME: This is an old version of the pvq adaptation code,
   it should be rewritten using the 2D code (not yet sure what
   it should look like) */

/*TODO: These should be adapted based on target quality/bitrate.*/
# define OD_K_INIT_Q7 (2031)

# define OD_COUNT_INIT_Q7 (104)
# define OD_COUNT_EX_INIT_Q7 (128)

#  define OD_SUM_EX_INIT_Q7 (216)

#  define OD_K_ADAPT_SPEED (2)
#  define OD_SUM_EX_ADAPT_SPEED (2)

# define OD_DELTA_ADAPT_SPEED (1)

/*The scaling for the row contexts is different.*/
# define OD_K_ROW_INIT_Q8 (2*OD_K_INIT_Q7-(OD_K_INIT_Q7>>OD_K_ADAPT_SPEED))
# define OD_SUM_EX_ROW_INIT_Q8 \
 (2*OD_SUM_EX_INIT_Q7-(OD_SUM_EX_INIT_Q7>>OD_SUM_EX_ADAPT_SPEED))
# define OD_COUNT_ROW_INIT_Q8 \
 (2*OD_COUNT_INIT_Q7-(OD_COUNT_INIT_Q7>>OD_DELTA_ADAPT_SPEED))
# define OD_COUNT_EX_ROW_INIT_Q8 \
 (2*OD_COUNT_EX_INIT_Q7-(OD_COUNT_EX_INIT_Q7>>OD_DELTA_ADAPT_SPEED))
# define OD_POS_ROW_INIT_Q4 \
 (2*OD_POS_INIT_Q3-(OD_POS_INIT_Q3>>OD_POS_ADAPT_SPEED+1))

struct od_pvq_adapt_ctx{
  /*Mean value of K.*/
  int mean_k_q8;
  /*Mean value of the sum of (pulses left)/(positions left) for all
     positions.*/
  int mean_sum_ex_q8;
  /*For the delta version.*/
  int mean_count_q8;
  int mean_count_ex_q8;
  /*Actual value of K for this context.*/
  int k;
  /*Actual sum of (pulses left)/(positions left) for all positions for this
     context.*/
  int sum_ex_q8;
  int count_q8;
  int count_ex_q8;
};



int run_pvq(int *X,int len,int N){
  od_pvq_adapt_ctx pvq_adapt;
  od_adapt_ctx adapt;
  int i, j;
  od_ec_enc enc;
  od_ec_dec dec;
  unsigned char *buf;
  ogg_uint32_t   buf_sz;
  int *Ki;
  GenericEncoder model;
  int EK=65536;
  int bits_used;

  Ki = malloc(sizeof(*Ki)*len);
  pvq_adapt.mean_k_q8=163;
  pvq_adapt.mean_sum_ex_q8=64;
  pvq_adapt.mean_count_q8=100*4;
  pvq_adapt.mean_count_ex_q8=256*4;
  od_ec_enc_init(&enc, EC_BUF_SIZE);
  generic_model_init(&model);
  for(i=0;i<len;i++){
    int K=0;
    for (j=0;j<N;j++)
      K += abs(X[i*N+j]);
    Ki[i] = K;
    generic_encode(&enc, &model, K, &EK, 4);
    adapt.mean[OD_ADAPT_K_Q8] = pvq_adapt.mean_k_q8;
    adapt.mean[OD_ADAPT_SUM_EX_Q8] = pvq_adapt.mean_sum_ex_q8;
    adapt.mean[OD_ADAPT_COUNT_Q8] = pvq_adapt.mean_count_q8;
    adapt.mean[OD_ADAPT_COUNT_EX_Q8] = pvq_adapt.mean_count_ex_q8;
    pvq_encoder(&enc,&X[i*N],N,K,&adapt);
    pvq_adapt.k = adapt.curr[OD_ADAPT_K_Q8];
    pvq_adapt.sum_ex_q8 = adapt.curr[OD_ADAPT_SUM_EX_Q8];
    pvq_adapt.count_q8 = adapt.curr[OD_ADAPT_COUNT_Q8];
    pvq_adapt.count_ex_q8 = adapt.curr[OD_ADAPT_COUNT_EX_Q8];
    if(pvq_adapt.k>=0){
      pvq_adapt.mean_k_q8+=(pvq_adapt.k<<8)-pvq_adapt.mean_k_q8>>OD_K_ADAPT_SPEED;
      pvq_adapt.mean_sum_ex_q8+=
       pvq_adapt.sum_ex_q8-pvq_adapt.mean_sum_ex_q8>>OD_SUM_EX_ADAPT_SPEED;
    }
    if(pvq_adapt.count_q8>=0){
      pvq_adapt.mean_count_q8+=(pvq_adapt.count_q8<<8)-pvq_adapt.mean_count_q8>>OD_DELTA_ADAPT_SPEED;
      pvq_adapt.mean_count_ex_q8+=
       pvq_adapt.count_ex_q8-pvq_adapt.mean_count_ex_q8>>OD_DELTA_ADAPT_SPEED;
    }
    /*if (i==0)
    {
      fprintf(stderr, "enc: K[%d]=%d, num=%d, den=%d, u=%d\n", i, Ki[i], num, den, u);
    }*/
    OD_ASSERT(!enc.error);
  }
  buf = od_ec_enc_done(&enc, &buf_sz);

  bits_used = od_ec_enc_tell(&enc);

  pvq_adapt.mean_k_q8=163;
  pvq_adapt.mean_sum_ex_q8=64;
  pvq_adapt.mean_count_q8=100*4;
  pvq_adapt.mean_count_ex_q8=256*4;
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
    adapt.mean[OD_ADAPT_K_Q8] = pvq_adapt.mean_k_q8;
    adapt.mean[OD_ADAPT_SUM_EX_Q8] = pvq_adapt.mean_sum_ex_q8;
    adapt.mean[OD_ADAPT_COUNT_Q8] = pvq_adapt.mean_count_q8;
    adapt.mean[OD_ADAPT_COUNT_EX_Q8] = pvq_adapt.mean_count_ex_q8;
    pvq_decoder(&dec, y, N, Ki[i], &adapt);
    pvq_adapt.k = adapt.curr[OD_ADAPT_K_Q8];
    pvq_adapt.sum_ex_q8 = adapt.curr[OD_ADAPT_SUM_EX_Q8];
    pvq_adapt.count_q8 = adapt.curr[OD_ADAPT_COUNT_Q8];
    pvq_adapt.count_ex_q8 = adapt.curr[OD_ADAPT_COUNT_EX_Q8];
    if(pvq_adapt.k>=0){
      pvq_adapt.mean_k_q8+=(pvq_adapt.k<<8)-pvq_adapt.mean_k_q8>>OD_K_ADAPT_SPEED;
      pvq_adapt.mean_sum_ex_q8+=
       pvq_adapt.sum_ex_q8-pvq_adapt.mean_sum_ex_q8>>OD_SUM_EX_ADAPT_SPEED;
    }
    if(pvq_adapt.count_q8>=0){
      pvq_adapt.mean_count_q8+=(pvq_adapt.count_q8<<8)-pvq_adapt.mean_count_q8>>OD_DELTA_ADAPT_SPEED;
      pvq_adapt.mean_count_ex_q8+=
       pvq_adapt.count_ex_q8-pvq_adapt.mean_count_ex_q8>>OD_DELTA_ADAPT_SPEED;
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
        fprintf(stderr, "K[%d]=%d, num=%d, den=%d\n", i, Ki[i],
         pvq_adapt.mean_k_q8, pvq_adapt.mean_sum_ex_q8);
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

void test_pvq_basic(int N) {
  int i;
  int j;
  int bits;
  int *X;
  int len;
  len = 4*N + 1;
  X = malloc(sizeof(*X)*N*len);
  for (i = 0; i<N; i++) X[i] = 0;
  for (i = 0; i<N; i++) {
    for (j = 0; j < N; j++) {
      X[(i + 1)*N + j] = (j == i);
    }
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      X[(i + N + 1)*N + j] = -(j == i);
    }
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++){
      X[(i+ N*2 + 1)*N + j] = j==i ?(rand() & 1? -2 : 2): 0;
    }
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      X[(i + N*3 + 1)*N + j] = j==i ?(rand() & 1? -3 : 3): 0;
    }
  }
  bits = run_pvq(X, len, N);
  fprintf(stderr,
   "Coded %d dim, few pulses with %f bits/sample (%1.4f bits/vector)\n",
   N,bits/(float)(len)/N,bits/(float)(len));
  free(X);
}

void test_pvq_huge(int N) {
  int i;
  int j;
  int bits;
  int *X;
  int len;
  len = 4*N;
  X = malloc(sizeof(*X)*N*len);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      X[i*N + j] = (j == i)*65535;
    }
  }
  for (i = 0; i<N; i++) {
    for (j = 0; j<N; j++) {
      X[(i + N)*N + j] = (j == i)*-65536;
    }
  }
  for (i = 0; i<N; i++) {
    for (j = 0; j<N; j++) {
      X[(i + N*2)*N + j] = (j == i)*-65535;
    }
  }
  for (i = 0; i<N; i++) {
    for (j = 0; j<N; j++) {
      X[(i + N*3)*N + j] = (j == i)*65534;
    }
  }
  bits = run_pvq(X, len, N);
  fprintf(stderr,
   "Coded %d dim, 65536 pulses with %f bits/sample (%1.4f bits/vector)\n",
   N,bits/(float)(len)/N,bits/(float)(len));
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
        if(fscanf(file,"%d",&X[i*N+j])!=1)
          return 1;
    bits = run_pvq(X,len,N);
    fprintf(stderr, "Coded file with %f bits/sample (%f bits/vector)\n",bits/(float)len/N,bits/(float)len);
    fclose(file);
    free(X);
  } else {
    int i;
    for (i=2; i<18; i++) test_pvq_basic(i);
    test_pvq_basic(64);
    test_pvq_basic(128);
    test_pvq_basic(256);
    for (i=2; i<18; i++) test_pvq_huge(i);
    test_pvq_huge(64);
#ifdef FAILING_TESTS
    /*Currently failing due to signed integer overflow in expectation.*/
    test_pvq_huge(128);
    test_pvq_huge(256);
#endif
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
