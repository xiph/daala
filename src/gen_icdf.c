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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define MAXN 32768

int main(int argc, char **argv)
{
  int i;
  int N;
  int shift;
  int j;
  float p[16];
  int pi[16];
  int sum;
  int aN[MAXN];

  N=atoi(argv[1]);
  shift=atoi(argv[2]);
  printf("/* This file is auto-generated using gen_icdf.c */\n\n");
  printf("const unsigned short icdf_table[%d][16] = {\n", N+1);
  printf("{15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0},\n");
  aN[0]=1;
  for (i=1;i<=N;i++){
    float Ex;
    float gamma;
    float a;
    int icdf;
    float maxp;
    int maxj;
    Ex=(float)i/(1<<shift);
    gamma = (sqrt(1+4*Ex*Ex)-1)/(2*Ex);
    a=-.5/log(gamma);
    aN[i] = floor(.5+256*exp(-1./a));
    /*printf("%f %f ", Ex, a);*/
    p[0]=1-exp(-.5/a);
    for(j=1;j<15;j++)
      p[j]=exp(-.5/a)*(exp(-((float)j-1.)/a)-exp(-(float)j/a));
    p[15]=exp(-.5/a)*exp(-14./a);

    sum=0;
    maxp=0;
    maxj=0;
    for(j=0;j<16;j++)
    {
      if (p[j]>maxp)
      {
        maxp=p[j];
        maxj=j;
      }
      pi[j] = floor(.5+32768*p[j]);
      if (pi[j]==0)
        pi[j]=1;
      sum += pi[j];
    }
    pi[maxj] += 32768-sum;


    icdf = 32768;
    printf("{");
    for(j=0;j<16;j++)
    {
      icdf -= pi[j];
      printf("%d,", icdf);
      /*printf("%d ", icdf);*/
    }
    printf("},\n");
    /*printf("\n");*/
  }
  printf("};\n");
  printf("\n\n");
  printf("const unsigned char decayE[%d] = {\n", N+1);
  for(i=0;i<=N;i++)
    printf("%d,\n", aN[i]);
  printf("};\n");
  return 0;
}
