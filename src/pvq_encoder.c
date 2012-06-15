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

/** Encodes the tail of a Laplace-distributed variable, i.e. it doesn't
 * do anything special for the zero case.
 *
 * @param [in,out] enc     range encoder
 * @param [in]     x       variable to encode (has to be positive)
 * @param [in]     decay   decay factor of the distribution in Q8 format, i.e. pdf ~= decay^x
 * @param [in]     max     maximum possible value of x (used to truncate the pdf)
 *
 */
void laplace_encode_special(ec_enc *enc,int x,unsigned decay,int max)
{
  unsigned decay2, decay4, decay8, decay16;
  unsigned short decay_icdf[2];
  unsigned char decay_icdf2[2];
  decay = OD_MINI(255,decay);
  /* powers of decay */
  decay2=decay*decay;
  decay4=decay2*decay2>>16;
  decay8=decay4*decay4>>16;
  decay16=decay8*decay8>>16;
  decay_icdf[0]=OD_MAXI(1,decay16>>1);
  decay_icdf[1]=0;
  if(max<0)
    max=0x7FFFFFFF;
  /* Encoding jumps of 16 with probability decay^16 */
  while(x>15 && max>=16){
    ec_enc_icdf16(enc,1,decay_icdf,15);
    x-=16;
    max-=16;
  }
  if(max>=16)
    ec_enc_icdf16(enc,0,decay_icdf,15);
#if 0
  ec_enc_bits(enc, x, 4);
#else
  decay_icdf2[1]=0;
  if(max>=8){
    /* p(x%16>8) = decay^8/(decay^8+1) */
    decay_icdf2[0]=OD_MAXI(1,decay8>>9);
    ec_enc_icdf_ft(enc, (x&0x8)!=0,decay_icdf2,decay_icdf2[0]+128);
    max-=8;
  }
  if(max>=4){
    /* p(x%8>4) = decay^4/(decay^4+1) */
    decay_icdf2[0]=OD_MAXI(1,decay4>>9);
    ec_enc_icdf_ft(enc, (x&0x4)!=0,decay_icdf2,decay_icdf2[0]+128);
    max-=4;
  }
  if(max>=2){
    /* p(x%4>2) = decay^2/(decay^2+1) */
    decay_icdf2[0]=OD_MAXI(1,decay2>>9);
    ec_enc_icdf_ft(enc, (x&0x2)!=0,decay_icdf2,decay_icdf2[0]+128);
    max-=2;
  }
  if(max>=1){
    /* p(x%2>1) = decay/(decay+1) */
    decay_icdf2[0]=OD_MAXI(1,decay>>1);
    ec_enc_icdf_ft(enc, (x&0x1)!=0,decay_icdf2,decay_icdf2[0]+128);
  }
#endif

}

/** Encodes a Laplace-distributed variable for use in PVQ
 *
 * @param [in,out] enc  range encoder
 * @param [in]     x    variable to encode (including sign)
 * @param [in]     ExQ8 expectation of the absolute value of x in Q8
 * @param [in]     K    maximum value of |x|
 */
static void laplace_encode(ec_enc *enc, int x, int ExQ8, int K)
{
  int j;
  int shift;
  int xs;
  unsigned short icdf[16];
  const unsigned short *icdf0, *icdf1;
  int sym;

  /* shift down x if expectation is too high */
  shift=od_ilog(ExQ8)-11;
  if(shift<0)
    shift=0;

  /* Apply the shift with rounding to Ex, K and xs */
  ExQ8=(ExQ8+(1<<shift>>1))>>shift;
  K=(K+(1<<shift>>1))>>shift;
  xs=(x+(1<<shift>>1))>>shift;

  /* Interpolate pre-computed icdfs based on Ex */
  icdf0=icdf_table[ExQ8>>4];
  icdf1=icdf_table[(ExQ8>>4)+1];
  for(j=0;j<16;j++)
    icdf[j]=((ExQ8&0xF)*icdf1[j]+(16-(ExQ8&0xF))*icdf0[j]+8)>>4;

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
    /* Because of the rounding, there's only half the number of possibilities for xs=0 */
    special=(xs==0);
    ec_enc_bits(enc,x-(xs<<shift)+(!special<<(shift-1)),shift-special);
  }

  if(xs>=15){
    unsigned decay;
    decay=decayE[(ExQ8+8)>>4];

    /* Handle the exponentially-decaying tail of the distribution */
    laplace_encode_special(enc,xs-15,decay,K);
  }
}

/** Encodes the position and sign of a single pulse in a vector of N.
 * The position is assumed to be Laplace-distributed.
 *
 * @param [in,out] enc range encoder
 * @param [in]     y   vector to encode
 * @param [in]     N   dimension of the vector
 * @param [in,out] u   mean position of the pulse (adapted)
 */
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
    /* Compute decay */
    decay=256-4096/ *u; /* Approximates 256*exp(-16./ *u); */
    /* Update mean position */
    *u+=pos-(*u>>4);
    if (*u<N/8)
      *u=N/8;

    laplace_encode_special(enc,pos,decay,N-1);
  }

  ec_enc_bits(enc,sign,1);
}

/** Encodes a vector of integers assumed to come from rounding a sequence of
 * Laplace-distributed real values in decreasing order of variance.
 *
 * @param [in,out] enc range encoder
 * @param [in]     y   vector to encode
 * @param [in]     N   dimension of the vector
 * @param [in]     K   sum of the absolute value of components of y
 * @param [in,out] num mean value of K (adapted)
 * @param [in,out] den mean value of remaining pulses/(N-i) (adapted)
 * @param [in,out] u   mean position of single-pulse sequences (adapted)
 */
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

  /* Special K==1 case to save CPU (should be roughly equivalent in terms of coding efficiency) */
  if(K==1){
    pvq_encoder1(enc,y,N,u);
    return;
  }
  /* Estimates the factor relating pulses_left and positions_left to E(|x|) */
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
    /* Expected value of x (round-to-nearest) is expQ8*pulses_left/positions_left */
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
  /* Adapting the estimates for expQ8 */
  *num+=256*K-(*num>>4);
  *den+=sumEx-(*den>>4);
}
