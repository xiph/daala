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
void laplace_encode_special(od_ec_enc *enc,int x,unsigned decay,int max)
{
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
  /* Encoding jumps of 16 with probability decay^16 */
  while(x>15 && max>=16){
    od_ec_encode_bool_q15(enc,1,decay_f);
    x-=16;
    max-=16;
  }
  if(max>=16)
    od_ec_encode_bool_q15(enc,0,decay_f);
#if 0
  od_ec_enc_bits(enc, x, 4);
#else
  if(max>=8){
    /* p(x%16>8) = decay^8/(decay^8+1) */
    decay_f=OD_MAXI(1,decay8>>2);
    od_ec_encode_bool(enc, (x&0x8)!=0,16384,16384+decay_f);
    max-=8;
  }
  if(max>=4){
    /* p(x%8>4) = decay^4/(decay^4+1) */
    decay_f=OD_MAXI(1,decay4>>2);
    od_ec_encode_bool(enc, (x&0x4)!=0,16384,16384+decay_f);
    max-=4;
  }
  if(max>=2){
    /* p(x%4>2) = decay^2/(decay^2+1) */
    decay_f=OD_MAXI(1,decay2>>2);
    od_ec_encode_bool(enc, (x&0x2)!=0,16384,16384+decay_f);
    max-=2;
  }
  if(max>=1){
    /* p(x%2>1) = decay/(decay+1) */
    decay_f=OD_MAXI(1,decay<<6);
    od_ec_encode_bool(enc, (x&0x1)!=0,16384,16384+decay_f);
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
static void laplace_encode(od_ec_enc *enc, int x, int ExQ8, int K)
{
  int j;
  int shift;
  int xs;
  ogg_uint16_t cdf[16];
  const ogg_uint16_t *cdf0, *cdf1;
  int sym;

  /* shift down x if expectation is too high */
  shift=od_ilog(ExQ8)-11;
  if(shift<0)
    shift=0;

  /* Apply the shift with rounding to Ex, K and xs */
  ExQ8=(ExQ8+(1<<shift>>1))>>shift;
  K=(K+(1<<shift>>1))>>shift;
  xs=(x+(1<<shift>>1))>>shift;

  /* Interpolate pre-computed cdfs based on Ex */
  cdf0=cdf_table[ExQ8>>4];
  cdf1=cdf_table[(ExQ8>>4)+1];
  for(j=0;j<16;j++)
    cdf[j]=((ExQ8&0xF)*cdf1[j]+(16-(ExQ8&0xF))*cdf0[j]+8)>>4;

  sym=xs;
  if (sym>15)
    sym=15;

  if (K<15){
    /* Simple way of truncating the pdf when we have a bound */
    od_ec_encode_cdf_unscaled(enc, sym, cdf, K+1);
  } else {
    od_ec_encode_cdf_q15(enc, sym, cdf, 16);
  }
  if(shift){
    int special;
    /* Because of the rounding, there's only half the number of possibilities for xs=0 */
    special=(xs==0);
    if (shift-special>0)
      od_ec_enc_bits(enc,x-(xs<<shift)+(!special<<(shift-1)),shift-special);
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
static void pvq_encoder1(od_ec_enc *enc, const int *y,int N,int *u)
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

  od_ec_enc_bits(enc,sign,1);
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
void pvq_encoder(od_ec_enc *enc, const int *y,int N,int K,int *num, int *den, int *u)
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
      od_ec_enc_bits(enc,y[i]<0,1);
    }
    Kn-=x;
  }
  /* Adapting the estimates for expQ8 */
  *num+=256*K-(*num>>4);
  *den+=sumEx-(*den>>4);
}


#if 0

GenericEncoder model;
int init=0;
int ExQ16=65536;

float cwrs_table8[] = {
    0.0000000, 4.000000, 7.000000, 9.426265, 11.459432, 13.202124, 14.721099, 16.063058,
    17.262095, 18.343810, 19.327833, 20.229476, 21.060865, 21.831734, 22.549995, 23.222151,
    23.853602, 24.448876, 25.011801, 25.545641, 26.053199, 26.536899, 26.998848, 27.440890,
    27.864649, 28.271557, 28.662890, 29.039782, 29.403252, 29.754214, 30.093492, 30.421834,
    30.739917, 31.048359, 31.347722, 31.638523, 31.921234, 32.196291, 32.464096, 32.725019,
    32.979403, 33.227566, 33.469804, 33.706393, 33.937588, 34.163629, 34.384740, 34.601131,
    34.812998, 35.020526, 35.223889,};
float cwrs_table4[] = {
    0.0000000, 3.000000, 5.000000, 6.459432, 7.584963, 8.491853, 9.247928, 9.894818,
    10.459432, 10.960002, 11.409391, 11.816984, 12.189825, 12.533330, 12.851749, 13.148477,
    13.426265, 13.687376, 13.933691, 14.166791, 14.388017, 14.598518, 14.799282, 14.991168,
    15.174926, 15.351215, 15.520619, 15.683653, 15.840778, 15.992407, 16.138912, 16.280626,
    16.417853, 16.550867, 16.679920, 16.805240, 16.927037, 17.045504, 17.160817, 17.273140,
    17.382624, 17.489409, 17.593625, 17.695391, 17.794822, 17.892021, 17.987086, 18.080110,
    18.171177, 18.260368, 18.347760,};
float cwrs_table16[] = {
    0.0000000, 5.000000, 9.000000, 12.417853, 15.426265, 18.121048, 20.563673, 22.797199,
    24.853602, 26.757622, 28.528982, 30.183768, 31.735318, 33.194854, 34.571911, 35.874671,
    37.110202, 38.284657, 39.403416, 40.471212, 41.492224, 42.470163, 43.408332, 44.309685,
    45.176872, 46.012277, 46.818052, 47.596148, 48.348332, 49.076215, 49.781267, 50.464831,
    51.128140, 51.772326, 52.398430, 53.007415, 53.600167, 54.177511, 54.740209, 55.288969,
    55.824451, 56.347271, 56.858002, 57.357180, 57.845307, 58.322854, 58.790263, 59.247948,
    59.696300, 60.135686, 60.566454,};

void pvq_encoder2(od_ec_enc *enc, const int *y,int N0,int K,int *num, int *den, int *u)
{
  int i;
  int sumEx;
  int Kn;
  int expQ8;
  int N;
  int M=8;

  if(K==0)
    return;
  sumEx=0;
  Kn=K;

  /* Special K==1 case to save CPU (should be roughly equivalent in terms of coding efficiency) */
  if(K==1){
    pvq_encoder1(enc,y,N0,u);
    return;
  }
  N = N0/M;
  /* Estimates the factor relating pulses_left and positions_left to E(|x|) */
  if (*num < 1<<23)
    expQ8=256**num/(1+*den);
  else
    expQ8=*num/(1+(*den>>8));

  for(i=0;i<N;i++){
    int j;
    int Ex;
    int x;
    if(Kn==0)
      break;
    x=0;
    for(j=0;j<M;j++)
      x+=abs(y[M*i+j]);
    /* Expected value of x (round-to-nearest) is expQ8*pulses_left/positions_left */
    Ex=(2*expQ8*Kn+(N-i))/(2*(N-i));
    if(Ex>Kn*256)
      Ex=Kn*256;
    sumEx+=(2*256*Kn+(N-i))/(2*(N-i));

    /* no need to encode the magnitude for the last bin */
    if(i!=N-1){
      if (!init)
      {
        init=1;
        generic_model_init(&model);
      }
      ExQ16 = Ex*256;
      if (Ex>=256)
        generic_encode(enc, &model, x, &ExQ16, 4);
      else
        laplace_encode(enc,x,Ex,Kn);
      printf("%f\n", cwrs_table8[x]);
    }
    Kn-=x;
  }
  /* Adapting the estimates for expQ8 */
  *num+=256*K-(*num>>4);
  *den+=sumEx-(*den>>4);
}
#endif

#if 0
void pvq_encoder3(od_ec_enc *enc, const int *_y,int N,int K,int *num, int *den, int *u)
{
  int i;
  int count=0;
  int conseq=0;
  int prev=0;
  int y[1024];
  int sumEx=0;
  int sumC=0;
  for (i=0;i<N;i++)
    y[i]=abs(_y[i]);
  i=0;
  while(K!=0)
  {
    if (y[i]!=0)
    {
      //printf("%d %f %d %d %d\n", count, (N-prev)/(float)K, conseq, K, i);
      laplace_encode(enc, count, *num/(float)*den*256*(N-prev)/K, N-prev);
      sumEx += 256*(N-prev);
      sumC += count*K;
      if (!conseq)
        od_ec_enc_bits(enc, _y[i]<0, 1);
      y[i]--;
      K--;
      conseq=1;
      count=0;
      prev=i;
    } else {
      i++;
      count++;
      conseq=0;
    }
  }
  *num += 256*sumC-(*num>>6);
  *den += sumEx-(*den>>6);
  //printf("%f\n", *num/(float)*den);
}
#endif


void pvq_encode_delta(od_ec_enc *enc, const int *y,int N,int K,int *num,
 int *den)
{
  int i;
  int prev=0;
  int sumEx=0;
  int sumC=0;
  int coef = 256**num/ *den;
  for(i=0;i<N;i++){
    if(y[i]!=0){
      int j;
      int count;
      int mag;

      mag = abs(y[i]);
      count = i-prev;
      laplace_encode(enc, count, coef*(N-prev)/K, N-prev-1);
      sumEx+=256*(N-prev)+256*(mag-1)*(N-i);
      sumC+=count*K;
      od_ec_enc_bits(enc,y[i]<0,1);
      for(j=0;j<mag-1;j++)
        laplace_encode(enc,0,coef*(N-i)/(K-1-j),N-i-1);
      K-=mag;
      prev=i;
      if (K==0)
        break;
    }
  }

  *num+=256*sumC-(*num>>6);
  *den+=sumEx-(*den>>6);
}
