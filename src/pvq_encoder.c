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

void laplace_encode_special(ec_enc *enc, int pos, unsigned decay)
{
  unsigned decay2, decay4, decay8, decay16;
  unsigned short decay_icdf[2];
  unsigned char decay_icdf2[2];
  decay = OD_MINI(255,decay);
  decay2=decay*decay;
  decay4=decay2*decay2>>16;
  decay8=decay4*decay4>>16;
  decay16=decay8*decay8>>16;
  decay_icdf[0]=OD_MAXI(1,decay16>>1);
  decay_icdf[1]=0;
  while(pos>15){
    ec_enc_icdf16(enc,1,decay_icdf,15);
    pos-=16;
  }
  ec_enc_icdf16(enc,0,decay_icdf,15);
#if 0
  ec_enc_bits(enc, pos, 4);
#else
  decay_icdf2[1]=0;
  decay_icdf2[0]=OD_MAXI(1,decay8>>9);
  ec_enc_icdf_ft(enc, (pos&0x8)!=0,decay_icdf2,decay_icdf2[0]+128);
  decay_icdf2[0]=OD_MAXI(1,decay4>>9);
  ec_enc_icdf_ft(enc, (pos&0x4)!=0,decay_icdf2,decay_icdf2[0]+128);
  decay_icdf2[0]=OD_MAXI(1,decay2>>9);
  ec_enc_icdf_ft(enc, (pos&0x2)!=0,decay_icdf2,decay_icdf2[0]+128);
  decay_icdf2[0]=OD_MAXI(1,decay>>1);
  ec_enc_icdf_ft(enc, (pos&0x1)!=0,decay_icdf2,decay_icdf2[0]+128);
#endif

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

  if (K<15){
    /* Simple way of truncating the pdf when we have a bound */
    for (j=0;j<=K;j++)
      icdf[j]-=icdf[K];
  }
  ec_enc_icdf16(enc, sym, icdf, 15);

  if(shift){
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(xs==0);
    ec_enc_bits(enc,x-(xs<<shift)+(!special<<(shift-1)),shift-special);
  }

  if(xs>=15){
    unsigned decay;
    decay=decayE[(Ex+8)>>4];

    laplace_encode_special(enc,xs-15,decay);
  }
}

static void pvq_encoder1(ec_enc *enc, const int *y,int N,int *u)
{
  int pos;
  int i;
  unsigned decay;
  int sign;

  pos=0;
  for(i=0;i<N;i++){
    if(y[i]!=0){
      pos=i;
      break;
    }
  }
  sign=y[pos]<0;

  if(N>1){
    decay=256-4096/ *u; /* Approximates 256*exp(-16./ *u); */
    *u+=pos-(*u>>4);
    if (*u<N/8)
      *u=N/8;

    laplace_encode_special(enc,pos,decay);
  }

  ec_enc_bits(enc,sign,1);
}

void pvq_encoder(ec_enc *enc, const int *y,int N,int K,int *num, int *den, int *u)
{
  int i;
  int sumEx;
  int Kn;
  int expQ8;

  if(K==0)
    return;
  sumEx=0;
  Kn=K;

  if(K==1){
    pvq_encoder1(enc,y,N,u);
    return;
  }
  if (*num < 1<<23)
    expQ8=256**num/(1+*den);
  else
    expQ8=*num/(1+(*den>>8));

  for(i=0;i<N;i++){
    int Ex;
    int x;
    if(Kn==0)
      break;
    x=abs(y[i]);
    /* Expected value of x (round-to-nearest) */
    Ex=(2*expQ8*Kn+(N-i))/(2*(N-i));
    if(Ex>Kn*256)
      Ex=Kn*256;
    sumEx+=(2*256*Kn+(N-i))/(2*(N-i));

    /* no need to encode the magnitude for the last bin */
    if(i!=N-1){
      laplace_encode(enc,x,Ex,Kn);
    }
    if(x!=0){
      ec_enc_bits(enc,y[i]<0,1);
    }
    Kn-=x;
  }
  *num+=256*K-(*num>>4);
  *den+=sumEx-(*den>>4);
}
