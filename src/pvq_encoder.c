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

#include <stdlib.h>
#include <stdio.h>
#include "internal.h"
#include "entenc.h"
#include "entcode.h"
#include "pvq_code.h"
#include "adapt.h"

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
#if 1
  int shift;
  int xs;
  int ms;
  int sym;
  const ogg_uint16_t *cdf;
  shift=0;
  /* We don't want a large decay value because that would require too many symbols.
     However, it's OK if the max is below 15. */
  while ( ((max>>shift)>=15||max==-1) && decay>235)
  {
    decay=(decay*decay+128)>>8;
    shift++;
  }
  decay = OD_MINI(decay, 254);
  decay = OD_MAXI(decay, 2);
  xs=x>>shift;
  ms=max>>shift;
  cdf = exp_cdf_table[(decay+1)>>1];
  /*printf("decay = %d\n", decay);*/
  do {
    sym = OD_MINI(xs, 15);
    /*{int i; printf("%d %d %d %d %d\n", x, xs, shift, sym, max);
    for (i=0;i<16;i++)
      printf("%d ", cdf[i]);
    printf("\n");}*/
    if (ms>0 && ms<15){
      /* Simple way of truncating the pdf when we have a bound */
      od_ec_encode_cdf_unscaled(enc, sym, cdf, ms+1);
    } else {
      od_ec_encode_cdf_q15(enc, sym, cdf, 16);
    }

    xs -= 15;
    ms -= 15;
  } while (sym>=15);
  if (shift)
    od_ec_enc_bits(enc, x&((1<<shift)-1), shift);
#else
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
    if(x&0x8)max-=8;
  }
  if(max>=4){
    /* p(x%8>4) = decay^4/(decay^4+1) */
    decay_f=OD_MAXI(1,decay4>>2);
    od_ec_encode_bool(enc, (x&0x4)!=0,16384,16384+decay_f);
    if(x&0x4)max-=4;
  }
  if(max>=2){
    /* p(x%4>2) = decay^2/(decay^2+1) */
    decay_f=OD_MAXI(1,decay2>>2);
    od_ec_encode_bool(enc, (x&0x2)!=0,16384,16384+decay_f);
    if (x&0x2)max-=2;
  }
  if(max>=1){
    /* p(x%2>1) = decay/(decay+1) */
    decay_f=OD_MAXI(1,decay<<6);
    od_ec_encode_bool(enc, (x&0x1)!=0,16384,16384+decay_f);
  }
#endif
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

/** Encodes a vector of integers assumed to come from rounding a sequence of
 * Laplace-distributed real values in decreasing order of variance.
 *
 * @param [in,out] enc range encoder
 * @param [in]     y   vector to encode
 * @param [in]     N   dimension of the vector
 * @param [in]     K   sum of the absolute value of components of y
 * @param [in,out] _adapt Adaptation context.
 */
void pvq_encoder(od_ec_enc *enc, const int *y,int N,int K,
 od_adapt_ctx *_adapt){
  int i;
  int sumEx;
  int Kn;
  int expQ8;
  int mean_k_q8;
  int mean_sum_ex_q8;
  _adapt->curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
  _adapt->curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
  if (K<=1)
  {
    pvq_encode_delta(enc, y, N, K, _adapt);
    return;
  }
  if(K==0){
    _adapt->curr[OD_ADAPT_K_Q8]=0;
    _adapt->curr[OD_ADAPT_SUM_EX_Q8]=0;
    return;
  }
  sumEx=0;
  Kn=K;
  /* Estimates the factor relating pulses_left and positions_left to E(|x|) */
  mean_k_q8=_adapt->mean[OD_ADAPT_K_Q8];
  mean_sum_ex_q8=_adapt->mean[OD_ADAPT_SUM_EX_Q8];
  if(mean_k_q8<1<<23)expQ8=256*mean_k_q8/(1+mean_sum_ex_q8);
  else expQ8=mean_k_q8/(1+(mean_sum_ex_q8>>8));
  for(i=0;i<N;i++){
    int Ex;
    int x;
    if(Kn==0)
      break;
    if (Kn<=1&&i!=N-1){
      pvq_encode_delta(enc, y+i, N-i, Kn, _adapt);
      i=N;
      break;
    }
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
  _adapt->curr[OD_ADAPT_K_Q8]=K-Kn;
  _adapt->curr[OD_ADAPT_SUM_EX_Q8]=sumEx;
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


void pvq_encode_delta(od_ec_enc *enc, const int *y,int N,int _K,
 od_adapt_ctx *_adapt)
{
  int i;
  int prev=0;
  int sumEx=0;
  int sumC=0;
  int first = 1;
  int K=_K;
  int coef = 256*_adapt->mean[OD_ADAPT_COUNT_Q8]/(1+_adapt->mean[OD_ADAPT_COUNT_EX_Q8]);
  coef=OD_MAXI(coef,1);
  for(i=0;i<N;i++){
    if(y[i]!=0){
      int j;
      int count;
      int mag;

      mag = abs(y[i]);
      count = i-prev;
      if(first)
      {
        int decay;
        int X = coef*(N-prev)/K;
        decay = OD_MINI(255,(int)((256*X/(X+256) + 8*X*X/(256*(N+1)*(N-1)*(N-1)))));
        /*Update mean position.*/
        laplace_encode_special(enc,count,decay,N-1);
        first = 0;
      } else {
        laplace_encode(enc, count, coef*(N-prev)/K, N-prev-1);
      }
      sumEx+=256*(N-prev);
      sumC+=count*K;
      od_ec_enc_bits(enc,y[i]<0,1);
      for(j=0;j<mag-1;j++)
      {
        laplace_encode(enc,0,coef*(N-i)/(K-1-j),N-i-1);
        sumEx += 256*(N-i);
      }
      K-=mag;
      prev=i;
      if (K==0)
        break;
    }
  }

  if (_K>0)
  {
    _adapt->curr[OD_ADAPT_COUNT_Q8]=256*sumC;
    _adapt->curr[OD_ADAPT_COUNT_EX_Q8]=sumEx;
  } else {
    _adapt->curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
    _adapt->curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
  }
  _adapt->curr[OD_ADAPT_K_Q8]=0;
  _adapt->curr[OD_ADAPT_SUM_EX_Q8]=0;
}
