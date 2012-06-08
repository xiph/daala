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

#include <stdlib.h>
#include <stdio.h>
#include "internal.h"
#include "entenc.h"
#include "entcode.h"
#include "pvq_code.h"

double bits_used=0;

extern const unsigned short icdf_table[][16];
extern const unsigned short expectA[];

void laplace_encode_special(ec_enc *enc, int pos, unsigned decay)
{
  unsigned decay2, decay4, decay8, decay16;
  unsigned short decay_icdf[2];
  if (decay>255)
    decay=255;
  if (decay==0)
    decay=1;
  decay2=decay*decay;
  decay4=decay2*decay2>>16;
  decay8=decay4*decay4>>16;
  decay16=decay8*decay8>>16;
  if (decay16<2)
    decay16=2;
  decay_icdf[0]=decay16>>1;
  decay_icdf[1]=0;
  while (pos>15)
  {
    ec_enc_icdf16(enc, 1, decay_icdf, 15);
    pos-=16;
    bits_used += -log2(decay_icdf[0]/32768.);
  }
  bits_used += -log2(1-decay_icdf[0]/32768.);
  ec_enc_icdf16(enc, 0, decay_icdf, 15);
#if 0
  ec_enc_bits(enc, pos, 4);
#else
  decay_icdf[0]=decay8>>1;
  ec_enc_icdf16(enc, (pos&0x8)!=0, decay_icdf, 15);
  decay_icdf[0]=decay4>>1;
  ec_enc_icdf16(enc, (pos&0x4)!=0, decay_icdf, 15);
  decay_icdf[0]=decay2>>1;
  ec_enc_icdf16(enc, (pos&0x2)!=0, decay_icdf, 15);
  decay_icdf[0]=decay<<7;
  ec_enc_icdf16(enc, (pos&0x1)!=0, decay_icdf, 15);
#endif
  bits_used += -log2(1-decay/256.)-(pos%16)*log2(decay/256.);

}


void laplace_encode(ec_enc *enc, int x, int Ex, int K)
{
  int j;
  int shift;
  int xs;
  unsigned short icdf[16];
  const unsigned short *icdf0, *icdf1;
  int sym;

  shift=od_ilog(Ex)-11;
  if(shift<0)
    shift=0;
  Ex=(Ex+(1<<shift>>1))>>shift;
  xs=(x+(1<<shift>>1))>>shift;
  icdf0=icdf_table[Ex>>4];
  icdf1=icdf_table[(Ex>>4)+1];
  for(j=0;j<16;j++)
    icdf[j]=((Ex&0xF)*icdf1[j]+(16-(Ex&0xF))*icdf0[j]+8)>>4;

  sym=xs;
  if (sym>15)
    sym=15;

  float tmp = 32768;
  if (K<15)
    tmp = 32768-icdf[K];
  bits_used += sym==0 ? -log2((32768-icdf[sym])/tmp) : -log2((icdf[sym-1]-icdf[sym])/tmp);
  if (K<15){
    /* Simple way of truncating the pdf when we have a bound */
    for (j=0;j<=K;j++)
      icdf[j]-=icdf[K];
  }
  ec_enc_icdf16(enc, sym, icdf, 15);

  if (shift)
  {
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(xs==0);
    ec_enc_bits(enc, x-(xs<<shift)+!special, shift-special);
    bits_used += shift-special;
  }

  if (xs>=15)
  {
    int a;
    unsigned decay;
    a=expectA[(Ex+8)>>4];
    decay = 256*exp(-256./a);

    laplace_encode_special(enc, xs-15, decay);

    bits_used += log2(1-exp(-256./a)) + (xs-15)*log2(exp(-256./a));
  }
}

static void pvq_encoder1(ec_enc *enc, const int *y,int N,int *u)
{
  int pos;
  int i;
  unsigned decay;
  int sign;

  pos=0;
  for(i=0;i<N;i++)
  {
    if (y[i]!=0)
    {
      pos=i;
      break;
    }
  }
  sign = y[pos]<0;

  decay = 256 - 4096 / *u; /* Approximates 256*exp(-16./ *u); */
  *u += pos - (*u>>4);
  if (*u<N/8)
    *u=N/8;

  laplace_encode_special(enc, pos, decay);

  ec_enc_bits(enc, sign, 1);
  bits_used += 1;
}

void pvq_encoder(ec_enc *enc, const int *y,int N,int K,int *num, int *den, int *u)
{
  int i;
  int sumEx;
  int Kn;
  sumEx=0;
  Kn=K;
  int expQ8;

  if (K==1)
  {
    pvq_encoder1(enc, y, N, u);
    return;
  }
  expQ8 = floor(.5+*num/(1+*den/256));

  for(i=0;i<N;i++){
    int Ex;
    int x;
    if (Kn==0)
      break;
    x = abs(y[i]);
    /* Expected value of x (round-to-nearest) */
#if 1
    Ex=(2*expQ8*Kn+(N-i))/(2*(N-i));
    if (Ex>Kn*256)
      Ex=Kn*256;
    sumEx += (2*256*Kn+(N-i))/(2*(N-i));
#else
    float r = 1-3.05/128;
    Ex = 256*K*(1-r)/(1-pow(r, N-i));
#endif
    /*printf("(%d %d %d 0x%x)\n", Kn, Ex, expQ8, enc->rng);*/
    /*fprintf(stderr, "%d %d %d\n", K, N-i, Ex);*/
    /* no need to encode the magnitude for the last bin */
    if (i!=N-1){
      laplace_encode(enc, x, Ex, Kn);
    }
    if (x!=0)
    {
      ec_enc_bits(enc, y[i]<0, 1);
      bits_used += 1;
    }
    Kn-=x;
  }
  if (K!=0)
  {
    *num += 256*K - (*num>>4);
    *den += sumEx - (*den>>4);
  }
}

#ifdef PVQ_MAIN

#if 1
#define SIZE 16
#define NB_VECTORS 60000
#else
#define SIZE 128
#define NB_VECTORS 213840
#endif

int main()
{
  int i, j;
  int num, den,u;
  ec_enc enc;
  ec_dec dec;
  int Ki[NB_VECTORS];
  unsigned char buf[1<<22];

  num = 650*4;
  den = 256*4;
  u = 30<<4;
  ec_enc_init(&enc, buf, 1<<22);
  for (i=0;i<NB_VECTORS;i++)
  {
    int K;
    int y[SIZE];
    K=0;
    /*if (i%400==0)
    {
      num=650*4;
      den=256*4;
    }*/
    for (j=0;j<SIZE;j++)
    {
      scanf("%d ", y+j);
      K += abs(y[j]);
    }
    Ki[i] = K;
    pvq_encoder(&enc, y, SIZE, K, &num, &den, &u);
  }
  ec_enc_done(&enc);
#if 1
  printf("DECODE\n");
  num = 650*4;
  den = 256*4;
  u = 30<<4;
  ec_dec_init(&dec, buf, 1<<22);
  for (i=0;i<NB_VECTORS;i++)
  {
    int y[SIZE];
    /*if (i%400==0)
    {
      num=650*4;
      den=256*4;
    }*/
    pvq_decoder(&dec, y, SIZE, Ki[i], &num, &den, &u);
    for (j=0;j<SIZE;j++)
      printf("%d ", y[j]);
    printf("\n");
  }
#endif

  printf("average bits = %f\n", bits_used/(float)NB_VECTORS);
  printf("tell()'s average = %f\n", ec_tell(&enc)/(float)NB_VECTORS);
  return 0;
}
#endif
