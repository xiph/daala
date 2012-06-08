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
#include "entdec.h"
#include "entcode.h"
#include "pvq_code.h"

extern const unsigned short icdf_table[][16];
extern const unsigned short expectA[];


int laplace_decode_special(ec_dec *dec, unsigned decay)
{
  unsigned decay2, decay4, decay8, decay16;
  int pos;
  unsigned short decay_icdf[2];
  if (decay>255)
    decay=255;
  decay2=decay*decay;
  decay4=decay2*decay2>>16;
  decay8=decay4*decay4>>16;
  decay16=decay8*decay8>>16;
  if (decay<2)
    decay=2;
  if (decay2<2)
    decay2=2;
  if (decay4<2)
    decay4=2;
  if (decay8<2)
    decay8=2;
  if (decay16<2)
    decay16=2;
  decay_icdf[0]=decay16>>1;
  decay_icdf[1]=0;
  pos = 0;
  while (ec_dec_icdf16(dec, decay_icdf, 15)==1)
  {
    pos+=16;
  }
#if 0
  pos += ec_dec_bits(dec, 4);
#else
  //printf("0x%x ", dec->rng);
  decay_icdf[0]=decay8>>1;
  pos += 8*ec_dec_icdf16(dec, decay_icdf, 15);
  decay_icdf[0]=decay4>>1;
  pos += 4*ec_dec_icdf16(dec, decay_icdf, 15);
  decay_icdf[0]=decay2>>1;
  pos += 2*ec_dec_icdf16(dec, decay_icdf, 15);
  decay_icdf[0]=decay<<7;
  pos += ec_dec_icdf16(dec, decay_icdf, 15);
  //printf("0x%x\n", dec->rng);
#endif
  return pos;
}


int laplace_decode(ec_dec *dec, int Ex, int K)
{
  int j;
  int shift;
  unsigned short icdf[16];
  const unsigned short *icdf0, *icdf1;
  int sym;
  int lsb=0;

  shift=od_ilog(Ex)-11;
  if(shift<0)
    shift=0;
  Ex=(Ex+(1<<shift>>1))>>shift;
  icdf0=icdf_table[Ex>>4];
  icdf1=icdf_table[(Ex>>4)+1];
  for(j=0;j<16;j++)
    icdf[j]=((Ex&0xF)*icdf1[j]+(16-(Ex&0xF))*icdf0[j]+8)>>4;

  if (K<15){
    /* Simple way of truncating the pdf when we have a bound */
    for (j=0;j<=K;j++)
      icdf[j]-=icdf[K];
  }
  sym = ec_dec_icdf16(dec, icdf, 15);
  if (shift)
  {
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(sym==0);
    lsb = ec_dec_bits(dec, shift-special)-!special;
  }

  if (sym==15)
  {
    int a;
    unsigned decay;
    a=expectA[(Ex+8)>>4];
    decay = 256*exp(-256./a);

    sym += laplace_decode_special(dec, decay);

  }

  return (sym<<shift)+lsb;
}


static void pvq_decoder1(ec_dec *dec, int *y,int N,int *u)
{
  int j;
  int pos;
  unsigned decay;

  for (j=0;j<N;j++)
    y[j]=0;
  decay = 256 - 4096 / *u; /* Approximates 256*exp(-16./ *u); */

  pos = laplace_decode_special(dec, decay);

  if (ec_dec_bits(dec, 1))
    y[pos]=-1;
  else
    y[pos]=1;
  *u += pos - (*u>>4);
  if (*u<N/8)
    *u=N/8;
  //bits_used += 1;
}

void pvq_decoder(ec_dec *dec, int *y,int N,int K,int *num, int *den, int *u)
{
  int i;
  int sumEx;
  int Kn;
  sumEx=0;
  Kn=K;
  int expQ8;

  if (K==1)
  {
    pvq_decoder1(dec, y, N, u);
    return;
  }
  expQ8 = floor(.5+*num/(1+*den/256));

  for(i=0;i<N;i++){
    int Ex;
    int x;
    if (Kn==0)
      break;
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
    /*printf("(%d %d %d 0x%x)\n", Kn, Ex, expQ8, dec->rng);*/
    /*fprintf(stderr, "%d %d %d\n", K, N-i, Ex);*/
    /* no need to encode the magnitude for the last bin */
    if (i!=N-1){
      x = laplace_decode(dec, Ex, Kn);
    } else {
      x = Kn;
    }
    if (x!=0)
    {
      if (ec_dec_bits(dec, 1))
        x=-x;
      //bits_used += 1;
    }
    y[i] = x;
    Kn-=abs(x);
  }
  if (K!=0)
  {
    *num += 256*K - (*num>>4);
    *den += sumEx - (*den>>4);
  }
  for(;i<N;i++)
    y[i] = 0;
}
