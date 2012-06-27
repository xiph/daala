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


/** Decodes the tail of a Laplace-distributed variable, i.e. it doesn't
 * do anything special for the zero case.
 *
 * @param [dec] range decoder
 * @param [decay] decay factor of the distribution, i.e. pdf ~= decay^x
 * @param [max] maximum possible value of x (used to truncate the pdf)
 *
 * @retval decoded variable x
 */
int laplace_decode_special(od_ec_dec *dec,unsigned decay,int max)
{
  int pos;
  unsigned decay2, decay4, decay8, decay16;
  unsigned decay_f;
  decay = OD_MINI(255,decay);
  /* powers of decay */
  decay2=decay*decay;
  decay4=decay2*decay2>>16;
  decay8=decay4*decay4>>16;
  decay16=decay8*decay8>>16;
  decay_f=32768U-OD_MAXI(1,decay16>>1);
  if(max<0)
    max=0x7FFFFFFF;
  pos=0;
  /* Decoding jumps of 16 with probability decay^16 */
  while(max>=16&&od_ec_decode_bool_q15(dec,decay_f)){
    pos+=16;
    max-=16;
  }
#if 0
  pos += od_ec_dec_bits(dec, 4);
#else
  if (max>=8){
    /* p(x%16>=8) = decay^8/(decay^8+1) */
    decay_f=OD_MAXI(1,decay8>>2);
    pos += 8*od_ec_decode_bool(dec,16384,16384+decay_f);
    max-=8;
  }
  if (max>=4){
    /* p(x%8>=4) = decay^4/(decay^4+1) */
    decay_f=OD_MAXI(1,decay4>>2);
    pos += 4*od_ec_decode_bool(dec,16384,16384+decay_f);
    max-=4;
  }
  if (max>=2){
    /* p(x%4>=2) = decay^2/(decay^2+1) */
    decay_f=OD_MAXI(1,decay2>>2);
    pos += 2*od_ec_decode_bool(dec,16384,16384+decay_f);
    max-=2;
  }
  if (max>=1){
    /* p(x%2>=1) = decay/(decay+1) */
    decay_f=OD_MAXI(1,decay<<6);
    pos += od_ec_decode_bool(dec,16384,16384+decay_f);
  }
#endif
  return pos;
}

/** Decodes a Laplace-distributed variable for use in PVQ
 *
 * @param [in,out] dec  range decoder
 * @param [in]     ExQ8 expectation of the absolute value of x
 * @param [in]     K    maximum value of |x|
 *
 * @retval decoded variable (including sign)
 */
static int laplace_decode(od_ec_dec *dec, int ExQ8, int K)
{
  int j;
  int shift;
  ogg_uint16_t cdf[16];
  const ogg_uint16_t *cdf0, *cdf1;
  int sym;
  int lsb=0;

  /* shift down x if expectation is too high */
  shift=od_ilog(ExQ8)-11;
  if(shift<0)
    shift=0;

  /* Apply the shift with rounding to Ex, K and xs */
  ExQ8=(ExQ8+(1<<shift>>1))>>shift;
  K=(K+(1<<shift>>1))>>shift;

  /* Interpolate pre-computed icdfs based on Ex */
  cdf0=cdf_table[ExQ8>>4];
  cdf1=cdf_table[(ExQ8>>4)+1];
  for(j=0;j<16;j++)
    cdf[j]=((ExQ8&0xF)*cdf1[j]+(16-(ExQ8&0xF))*cdf0[j]+8)>>4;

  if(K<15){
    /* Simple way of truncating the pdf when we have a bound */
    od_ec_decode_cdf_unscaled(dec,cdf,K+1);
  }
  sym=od_ec_decode_cdf_q15(dec,cdf,16);
  if(shift){
    int special;
    /* Because of the rounding, there's only half the number of possibilities for xs=0 */
    special=(sym==0);
    if (shift-special>0)
      lsb=od_ec_dec_bits(dec,shift-special);
    lsb-=(!special<<(shift-1));
  }

  if(sym==15){
    unsigned decay;
    decay=decayE[(ExQ8+8)>>4];

    /* Handle the exponentially-decaying tail of the distribution */
    sym+=laplace_decode_special(dec,decay,K);
  }

  return (sym<<shift)+lsb;
}


/** Decodes the position and sign of a single pulse in a vector of N.
 * The position is assumed to be Laplace-distributed.
 *
 * @param [in,out] dec range decoder
 * @param [in]     y   decoded vector
 * @param [in]     N   dimension of the vector
 * @param [in,out] u   mean position of the pulse (adapted)
 */
static void pvq_decoder1(od_ec_dec *dec, int *y,int N,int *u)
{
  int j;
  int pos;
  unsigned decay;

  for(j=0;j<N;j++)
    y[j]=0;

  if(N>1){
    /* Compute decay */
    decay=256-4096 / *u; /* Approximates 256*exp(-16./ *u); */

    pos=laplace_decode_special(dec,decay,N-1);

    /* Update mean position */
    *u+=pos-(*u>>4);
    if(*u<N/8)
      *u=N/8;
  } else {
    pos=0;
  }
  if(od_ec_dec_bits(dec,1))
    y[pos]=-1;
  else
    y[pos]=1;
}

/** Decodes a vector of integers assumed to come from rounding a sequence of
 * Laplace-distributed real values in decreasing order of variance.
 *
 * @param [in,out] dec range decoder
 * @param [in]     y   decoded vector
 * @param [in]     N   dimension of the vector
 * @param [in]     K   sum of the absolute value of components of y
 * @param [in,out] num mean value of K (adapted)
 * @param [in,out] den mean value of remaining pulses/(N-i) (adapted)
 * @param [in,out] u   mean position of single-pulse sequences (adapted)
 */
void pvq_decoder(od_ec_dec *dec, int *y,int N,int K,int *num, int *den, int *u)
{
  int i;
  int sumEx;
  int Kn;
  int expQ8;

  if(K==0){
    for(i=0;i<N;i++)
      y[i]=0;
    return;
  }
  sumEx=0;
  Kn=K;

  /* Special K==1 case to save CPU (should be roughly equivalent in terms of coding efficiency) */
  if(K==1){
    pvq_decoder1(dec,y,N,u);
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
    /* Expected value of x (round-to-nearest) is expQ8*pulses_left/positions_left */
    Ex=(2*expQ8*Kn+(N-i))/(2*(N-i));
    if(Ex>Kn*256)
      Ex=Kn*256;
    sumEx+=(2*256*Kn+(N-i))/(2*(N-i));
    /* no need to encode the magnitude for the last bin */
    if (i!=N-1){
      x=laplace_decode(dec, Ex, Kn);
    } else {
      x=Kn;
    }
    if(x!=0)
    {
      if(od_ec_dec_bits(dec, 1))
        x=-x;
    }
    y[i]=x;
    Kn-=abs(x);
  }
  /* Adapting the estimates for expQ8 */
  *num+= 256*K-(*num>>4);
  *den+= sumEx-(*den>>4);
  for(;i<N;i++)
    y[i]=0;
}

void pvq_decode_delta(od_ec_dec *dec, int *y,int N,int K,int *num, int *den)
{
  int i;
  int prev=0;
  int sumEx=0;
  int sumC=0;
  int coef = 256**num/ *den;
  int pos=0;
  int K0;
  int sign=0;
  for(i=0;i<N;i++)
    y[i]=0;
  K0=K;
  for(i=0;i<K0;i++){
    int count;

    count = laplace_decode(dec, coef*(N-prev)/K, N-prev-1);
    sumEx+=256*(N-prev);
    sumC+=count*K;
    pos += count;
    if (y[pos]==0)
      sign = od_ec_dec_bits(dec,1);
    y[pos]+=sign?-1:1;
    prev=pos;
    K--;
    if (K==0)
      break;
  }

  *num+=256*sumC-(*num>>6);
  *den+=sumEx-(*den>>6);
}
