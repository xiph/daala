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

#include "pvq_encoder.c"
#include "generic_encoder.c"
#include "pvq_decoder.c"
#include "generic_decoder.c"
#include "entenc.c"
#include "entdec.c"
#include "icdf_table.c"
#include "internal.c"
#include <stdlib.h>
#include <math.h>

#define EC_BUF_SIZE (1<<24)
#define MAX_VECTORS 10000
#define MAXN 256

void test_pvq_sequence(int len,int N,float std){
  int i, j;
  int num, den,u;
  ec_enc enc;
  ec_dec dec;
  unsigned char *buf;
  int *X;
  int *Ki;

  X = malloc(sizeof(*X)*len*N);
  Ki = malloc(sizeof(*Ki)*len);
  buf = malloc(sizeof(*buf)*EC_BUF_SIZE);
  for(i=0;i<len;i++){
    for(j=0;j<N;j++){
      X[i*N+j]=-floor(.5+std*log((rand()+1)/(float)RAND_MAX));
      if (rand()&1)
        X[i*N+j]=-X[i*N+j];
      /*printf("%d ", X[i][j]);*/
    }
    /*printf("\n");*/
  }
  num = 650*4;
  den = 256*4;
  u = 30<<4;
  ec_enc_init(&enc, buf, EC_BUF_SIZE);
  for(i=0;i<len;i++){
    int K=0;
    for (j=0;j<N;j++)
      K += abs(X[i*N+j]);
    Ki[i] = K;
    pvq_encoder(&enc,&X[i*N],N,K,&num,&den,&u);
    /*if (i==0)
    {
      fprintf(stderr, "enc: K[%d]=%d, num=%d, den=%d, u=%d\n", i, Ki[i], num, den, u);
    }*/
    od_assert(!ec_get_error(&enc));
    od_assert(ec_tell(&enc)<EC_BUF_SIZE<<3);
  }
  ec_enc_done(&enc);

  printf("tell()'s average for %d dim, std=%f is %f bits/sample\n", N, std, ec_tell(&enc)/(float)len/N);

  num = 650*4;
  den = 256*4;
  u = 30<<4;
  ec_dec_init(&dec, buf, EC_BUF_SIZE);

  for (i=0;i<len;i++)
  {
    int y[MAXN];
    pvq_decoder(&dec, y, N, Ki[i], &num, &den, &u);
    for (j=0;j<N;j++){
      if(y[j]!=X[i*N+j]){
        int k;
        fprintf(stderr, "mismatch for coef %d of vector %d (N=%d, std=%f)\n", j, i, N, std);
        fprintf(stderr, "orig vector:\n");
        for(k=0;k<N;k++)
          fprintf(stderr,"%d ", X[i*N+k]);
        fprintf(stderr,"\n");
        fprintf(stderr, "decoded vector:\n");
        for(k=0;k<N;k++)
          fprintf(stderr,"%d ", y[k]);
        fprintf(stderr,"\n");
        fprintf(stderr, "K[%d]=%d, num=%d, den=%d, u=%d\n", i, Ki[i], num, den, u);
        abort();
      }
    }
  }

  free(X);
  free(Ki);
  free(buf);
}

int main(void){
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

  return 0;
}
