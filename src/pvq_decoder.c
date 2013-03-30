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
  ms=max>>shift;
  cdf = exp_cdf_table[(decay+1)>>1];
  /*printf("decay = %d\n", decay);*/
  xs=0;
  do {
    sym = OD_MINI(xs, 15);
    /*{int i; printf("%d %d %d %d %d\n", x, xs, shift, sym, max);
    for (i=0;i<16;i++)
      printf("%d ", cdf[i]);
    printf("\n");}*/
    if (ms>0 && ms<15){
      /* Simple way of truncating the pdf when we have a bound */
      sym = od_ec_decode_cdf_unscaled(dec, cdf, ms+1);
    } else {
      sym = od_ec_decode_cdf_q15(dec, cdf, 16);
    }

    xs += sym;
    ms -= 15;
  } while (sym>=15);
  if (shift)
    pos = (xs<<shift) + od_ec_dec_bits(dec, shift);
  else
    pos = xs;

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
    if(pos&0x8)max-=8;
  }
  if (max>=4){
    /* p(x%8>=4) = decay^4/(decay^4+1) */
    decay_f=OD_MAXI(1,decay4>>2);
    pos += 4*od_ec_decode_bool(dec,16384,16384+decay_f);
    if(pos&0x4)max-=4;
  }
  if (max>=2){
    /* p(x%4>=2) = decay^2/(decay^2+1) */
    decay_f=OD_MAXI(1,decay2>>2);
    pos += 2*od_ec_decode_bool(dec,16384,16384+decay_f);
    if(pos&0x2)max-=2;
  }
  if (max>=1){
    /* p(x%2>=1) = decay/(decay+1) */
    decay_f=OD_MAXI(1,decay<<6);
    pos += od_ec_decode_bool(dec,16384,16384+decay_f);
  }
#endif
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
    sym=od_ec_decode_cdf_unscaled(dec,cdf,K+1);
  } else {
    sym=od_ec_decode_cdf_q15(dec,cdf,16);
  }
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

/** Decodes a vector of integers assumed to come from rounding a sequence of
 * Laplace-distributed real values in decreasing order of variance.
 *
 * @param [in,out] dec range decoder
 * @param [in]     y   decoded vector
 * @param [in]     N   dimension of the vector
 * @param [in]     K   sum of the absolute value of components of y
 * @param [in,out] _adapt Adaptation context.
 */
void pvq_decoder(od_ec_dec *dec, int *y,int N,int K,od_adapt_ctx *_adapt)
{
  int i;
  int sumEx;
  int Kn;
  int expQ8;
  int mean_k_q8;
  int mean_sum_ex_q8;
  _adapt->curr[OD_ADAPT_COUNT_Q8]=OD_ADAPT_NO_VALUE;
  _adapt->curr[OD_ADAPT_COUNT_EX_Q8]=OD_ADAPT_NO_VALUE;
  if(K<=1)
  {
    pvq_decode_delta(dec, y, N, K, _adapt);
    return;
  }
  if(K==0){
    _adapt->curr[OD_ADAPT_K_Q8]=0;
    _adapt->curr[OD_ADAPT_SUM_EX_Q8]=0;
    for(i=0;i<N;i++)
      y[i]=0;
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
    if(Kn<=1&&i!=N-1){
      pvq_decode_delta(dec, y+i, N-i, Kn, _adapt);
      i=N;
      break;
    }
    /* Expected value of x (round-to-nearest) is expQ8*pulses_left/positions_left */
    Ex=(2*expQ8*Kn+(N-i))/(2*(N-i));
    if(Ex>Kn*256)
      Ex=Kn*256;
    sumEx+=(2*256*Kn+(N-i))/(2*(N-i));
    /* no need to encode the magnitude for the last bin */
    if(i!=N-1){
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
  _adapt->curr[OD_ADAPT_K_Q8]=K-Kn;
  _adapt->curr[OD_ADAPT_SUM_EX_Q8]=sumEx;
  for(;i<N;i++)
    y[i]=0;
}

void pvq_decode_delta(od_ec_dec *dec, int *y,int N,int _K,
 od_adapt_ctx *_adapt)
{
  int i;
  int prev=0;
  int sumEx=0;
  int sumC=0;
  int coef = 256*_adapt->mean[OD_ADAPT_COUNT_Q8]/(1+_adapt->mean[OD_ADAPT_COUNT_EX_Q8]);
  int pos=0;
  int K0;
  int sign=0;
  int first = 1;
  int K=_K;

  for(i=0;i<N;i++)
    y[i]=0;
  K0=K;
  coef=OD_MAXI(coef,1);
  for(i=0;i<K0;i++){
    int count;

    if(first)
    {
      int decay;
      int X = coef*(N-prev)/K;
      decay = OD_MINI(255,(int)((256*X/(X+256) + 8*X*X/(256*(N+1)*(N-1)*(N-1)))));
      /*Update mean position.*/
      count = laplace_decode_special(dec,decay,N-1);
      first = 0;
    } else {
      count = laplace_decode(dec, coef*(N-prev)/K, N-prev-1);
    }
    sumEx+=256*(N-prev);
    sumC+=count*K;
    pos += count;
    if (y[pos]==0)
      sign = od_ec_dec_bits(dec,1);
    y[pos]+=sign?-1:1;
    prev=pos;
    K--;
    if(K==0)
      break;
  }

  if (_K>0)
  {
    _adapt->curr[OD_ADAPT_COUNT_Q8]=256*sumC;
    _adapt->curr[OD_ADAPT_COUNT_EX_Q8]=sumEx;
  } else {
    _adapt->curr[OD_ADAPT_COUNT_Q8]=-1;
    _adapt->curr[OD_ADAPT_COUNT_EX_Q8]=0;
  }
  _adapt->curr[OD_ADAPT_K_Q8]=0;
  _adapt->curr[OD_ADAPT_SUM_EX_Q8]=0;
}
