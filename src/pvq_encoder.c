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

double bits_used=0;

extern const unsigned short icdf_table[][16];
extern const unsigned short expectA[];


void laplace_encode(int x, int Ex)
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
  Ex>>=shift;
  xs=(x+(1<<shift>>1))>>shift;
  icdf0=icdf_table[Ex>>4];
  icdf1=icdf_table[(Ex>>4)+1];
  for(j=0;j<16;j++)
    icdf[j]=((Ex&0xF)*icdf1[j]+(16-(Ex&0xF))*icdf0[j]+8)>>4;

  sym=xs;
  if (sym>15)
    sym=15;

  /*ec_enc_icdf(enc, sym, icdf);*/
  bits_used += sym==0 ? -log2((32768-icdf[sym])/32768.) : -log2((icdf[sym-1]-icdf[sym])/32768.);

  if (shift)
  {
    int special;
    /* There's something special around zero after shift because of the rounding */
    special=(xs==0);
    /*ec_enc_bits(enc, x-(xs<<shift), shift-special);*/
    bits_used += shift-special;
  }

  if (xs>=15)
  {
    int a=expectA[(Ex+8)>>4];
    bits_used += log2(1-exp(-256./a)) + (xs-15)*log2(exp(-256./a));
  }
}

void pvq_encoder(const int *y,int N,int K,int expQ8)
{
  int i;

  for(i=0;i<N;i++){
    int Ex;
    int x;
    if (K==0)
      break;
    x = abs(y[i]);
    /* Expected value of x (round-to-nearest) */
    Ex=(2*expQ8*K+(N-i))/(2*(N-i));
    /*fprintf(stderr, "%d %d %d\n", K, N-i, Ex);*/
    laplace_encode(x, Ex);
    if (x!=0)
    {
      /*ec_enc_bits(enc, y[i]<0, 1);*/
      bits_used += 1;
    }
    K-=x;
  }
}

int main()
{
  int i, j;
  for (i=0;i<213840;i++)
  {
    int K;
    int y[128];
    K=0;
    for (j=0;j<128;j++)
    {
      scanf("%d ", y+j);
      K += abs(y[j]);
    }
    pvq_encoder(y, 128, K, 650);
  }
  printf("total bits = %f\n", bits_used);
  return 0;
}
