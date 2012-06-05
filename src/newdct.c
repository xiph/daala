/*Daala video codec
  Copyright (C) 2002-2012 Daala project contributors

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
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.*/
#include "filter.h"

/*This should translate directly to 3 or 4 instructions for a constant _b:
#define OD_UNBIASED_RSHIFT(_a,_b) ((_a)+(((1<<(_b))-1)&-((_a)<0))>>(_b))*/
/*This version relies on a smart compiler:*/
#define OD_UNBIASED_RSHIFT(_a,_b) ((_a)/(1<<(_b)))

#if 1
# define OD_DCT_RSHIFT(_a,_b) OD_DIV_POW2_RE(_a,_b)
#elif 1
# define OD_DCT_RSHIFT(_a,_b) OD_UNBIASED_RSHIFT(_a,_b)
#else
# define OD_DCT_RSHIFT(_a,_b) ((_a)>>(_b))
#endif

#if 0
void od_bin_fdct4(
od_coeff _y[],const od_coeff _x[]){
  /*9 adds, 4 shifts, 3 "muls".*/
  int t0;
  int t1;
  int t2;
  int t3;
  /*+1/-1 butterflies and initial permutation:*/
  t0=_x[0]+_x[3];
  t2=_x[1]+_x[2];
  t3=t2-(_x[2]<<1);
  t1=t0-(_x[3]<<1);
  /*69/64~=2\sqrt{2}*sin(\frac{\pi}{8})~=1.0823922002923939687994464107328*/
  t3=(t3*69>>6)+(t3>0);
  /*21/16~=\sqrt{2}*cos(\frac{\pi}{8})~=1.3065629648763765278566431734272*/
  t1=(t1*21>>4)+(t1>0);
  /*+ Embedded 2-point type-II DCT.*/
  t0+=t2;
  t2=t0-(t2<<1);
  /*+ Embedded 2-point type-II DCT.*/
  t1+=OD_DCT_RSHIFT(t3,1);
  t3=t1-t3;
  /*91/64~=\sqrt{2}~=1.4142135623730950488016887242097*/
  t3=(t3*91>>6)+(t3>0);
  /*Clean-up.*/
  t3-=t1;
  _y[0]=(od_coeff)t0;
  _y[1]=(od_coeff)t1;
  _y[2]=(od_coeff)t2;
  _y[3]=(od_coeff)t3;
}

void od_bin_idct4(od_coeff _x[],const od_coeff _y[]){
  int t0;
  int t1;
  int t2;
  int t3;
  t3=_y[1]-(_y[3]+_y[1]<<6)/91;
  t1=(_y[1]-(OD_DCT_RSHIFT(t3,1))<<4)/21;
  t2=OD_DCT_RSHIFT(_y[0]-_y[2],1);
  t0=_y[0]-t2;
  t1=OD_DCT_RSHIFT(t0-t1,1);
  t3=OD_DCT_RSHIFT(t2-(t3<<6)/69,1);
  _x[3]=(od_coeff)t1;
  _x[2]=(od_coeff)t3;
  _x[1]=(od_coeff)(t2-t3);
  _x[0]=(od_coeff)(t0-t1);
}
#else

void od_bin_fdct4(od_coeff _y[],const od_coeff _x[]){
  /*9 adds, 2 shifts, 3 "muls".*/
  int t0;
  int t1;
  int t2;
  int t2h;
  int t3;
  /*+1/-1 butterflies and initial permutation:*/
  t3=_x[0]-_x[3];
  t2=_x[1]+_x[2];
  t2h=OD_DCT_RSHIFT(t2,1);
  t1=t2h-_x[2];
  t0=_x[0]-OD_DCT_RSHIFT(t3,1);
  /*+ Embedded 2-point type-II DCT.*/
  t0+=t2h;
  t2=t0-t2;
  /*+ Embedded 2-point type-IV DST.*/
  /*45/64~=4*sin(\frac{\pi}{8})-2*tan(\frac{\pi}{8})~=
     0.70230660471416898931046248770220*/
  t3-=t1*45+32>>6;
  /*21/32~=\sqrt{1/2}*cos(\frac{\pi}{8}))~=0.65328148243818826392832158671359*/
  t1+=t3*21+16>>5;
  /*71/64~=4*sin(\frac{\pi}{8})-tan(\frac{\pi}{8})~=
     1.1165201670872640381121512119119*/
  t3-=t1*71+32>>6;
  _y[0]=(od_coeff)t0;
  _y[1]=(od_coeff)t1;
  _y[2]=(od_coeff)t2;
  _y[3]=(od_coeff)t3;
}

void od_bin_idct4(od_coeff _x[],const od_coeff _y[]){
  int t0;
  int t1;
  int t2;
  int t2h;
  int t3;
  t3=_y[3]+(_y[1]*71+32>>6);
  t1=_y[1]-(t3*21+16>>5);
  t3+=t1*45+32>>6;
  t2=_y[0]-_y[2];
  t2h=OD_DCT_RSHIFT(t2,1);
  t0=_y[0]-t2h+OD_DCT_RSHIFT(t3,1);
  t1=t2h-t1;
  _x[0]=(od_coeff)t0;
  _x[1]=(od_coeff)(t2-t1);
  _x[2]=(od_coeff)t1;
  _x[3]=(od_coeff)(t0-t3);
}
#endif

void od_bin_fdct8(od_coeff _y[],const od_coeff _x[]){
  /*29 adds, 11 shifts, 11 "muls".*/
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t4h;
  int t5;
  int t6;
  int t7;
  /*+1/-1 butterflies and initial permutation:*/
  t0=_x[0]+_x[7];
  t4=_x[1]+_x[6];
  t2=_x[2]+_x[5];
  t6=_x[3]+_x[4];
  t3=_x[4]-OD_DCT_RSHIFT(t6,1);
  t7=OD_DCT_RSHIFT(t2,1)-_x[5];
  t5=OD_DCT_RSHIFT(t4,1)-_x[6];
  t1=OD_DCT_RSHIFT(t0,1)-_x[7];
  /*25/16~=8*sin(\frac{\pi}{16})~=1.5607225761290261427862789478162*/
  t3=(t3*25>>4)+(t3>0);
  /*71/64~=2*sin(\frac{3\pi}{16})~=1.1111404660392044494856616278971*/
  t7=(t7*71>>6)+(t7>0);
  /*53/32~=2*cos(\frac{3\pi}{16})~=1.6629392246050904741575767552358*/
  t5=(t5*53>>5)+(t5>0);
  /*63/32~=2*cos(\frac{\pi}{16})~=1.9615705608064608982523644722685*/
  t1=(t1*63>>5)+(t1>0);
  /*+ Embedded 4-point type-II DCT.*/
  t6=t0-t6;
  t4+=t2;
  t4h=OD_DCT_RSHIFT(t4,1);
  t2=t4h-t2;
  t0-=OD_DCT_RSHIFT(t6,1);
  /*|-+ Embedded 2-point type-II DCT.*/
  t0+=t4h;
  t4=t0-t4;
  /*|-+ Embedded 2-point type-IV DST.*/
  /*45/64~=4*sin(\frac{\pi}{8})-2*tan(\frac{\pi}{8})~=
     0.70230660471416898931046248770220*/
  t6-=t2*45+32>>6;
  /*21/32~=\sqrt{1/2}*cos(\frac{\pi}{8}))~=0.65328148243818826392832158671359*/
  t2+=t6*21+16>>5;
  /*71/64~=4*sin(\frac{\pi}{8})-tan(\frac{\pi}{8})~=
     1.1165201670872640381121512119119*/
  t6-=t2*71+32>>6;
  /*+ Embedded 4-point type-II DCT.*/
  t1-=OD_DCT_RSHIFT(t3,2);
  t5+=t7;
  t7=OD_DCT_RSHIFT(t5,1)-t7;
  t3+=t1<<1;
  /*69/64~=2*\sqrt{2}*sin(\frac{\pi}{8})~=1.0823922002923939687994464107328*/
  t7=(t7*69>>6)+(t7>0);
  /*21/16~=\sqrt{2}*cos(\frac{\pi}{8})~=1.3065629648763765278566431734272*/
  t3=(t3*21>>4)+(t3>0);
  /*|-+ Embedded 2-point type-II DCT.*/
  t5=t1-t5;
  t1-=OD_DCT_RSHIFT(t5,1);
  /*91/64~=\sqrt{2}~=1.4142135623730950488016887242097*/
  t1=(t1*91>>6)+(t1>0);
  /*|-+ Embedded 2-point type-II DCT.*/
  t7=OD_DCT_RSHIFT(t3,1)-t7;
  t3-=t7;
  /*91/64~=\sqrt{2}~=1.4142135623730950488016887242097*/
  t7=(t7*91>>6)+(t7>0);
  /*Clean-up.*/
  t3-=t1;
  t7-=t5;
  t5-=t3;
  t7-=t1;
  _y[0]=(od_coeff)t0;
  _y[1]=(od_coeff)t1;
  _y[2]=(od_coeff)t2;
  _y[3]=(od_coeff)t3;
  _y[4]=(od_coeff)t4;
  _y[5]=(od_coeff)t5;
  _y[6]=(od_coeff)t6;
  _y[7]=(od_coeff)t7;
}

void od_bin_idct8(od_coeff _y[],const od_coeff _x[]){
  int t0;
  int t1;
  int t2;
  int t3;
  int t4;
  int t4h;
  int t5;
  int t6;
  int t7;
  t7=_x[7]+_x[1];
  t5=_x[5]+_x[3];
  t7=(t7+t5<<6)/91;
  t3=_x[3]+_x[1]+t7;
  t1=(_x[1]<<6)/91+OD_DCT_RSHIFT(t5,1);
  t5=t1-t5;
  t7=OD_DCT_RSHIFT(t5,1)-(OD_DCT_RSHIFT(t3,1)-t7<<6)/69;
  t3=(t3<<4)/21-(t1<<1);
  t5=(t5-t7<<5)/53;
  t1=(t1+OD_DCT_RSHIFT(t3,2)<<5)/63;
  t6=_x[6]+(_x[2]*71+32>>6);
  t2=_x[2]-(t6*21+16>>5);
  t6+=t2*45+32>>6;
  t4=_x[0]-_x[4];
  t4h=OD_DCT_RSHIFT(t4,1);
  t0=_x[0]-t4h+OD_DCT_RSHIFT(t6,1);
  t6=t0-t6;
  t2=t4h-t2;
  t4-=t2;
  t3=(t3<<4)/25+OD_DCT_RSHIFT(t6,1);
  t7=OD_DCT_RSHIFT(t2,1)-(t7<<6)/71;
  t5=OD_DCT_RSHIFT(t4,1)-t5;
  t1=OD_DCT_RSHIFT(t0,1)-t1;
  _y[0]=t0-t1;
  _y[1]=t4-t5;
  _y[2]=t2-t7;
  _y[3]=t6-t3;
  _y[4]=t3;
  _y[5]=t7;
  _y[6]=t5;
  _y[7]=t1;
}

#if 0
/*Test code.*/
#include <stdio.h>
#include <math.h>
#include <string.h>


/*The auto-correlation coefficent. 0.95 is a common value.*/
#define INPUT_AUTOCORR (0.95)
#define INPUT_AUTOCORR_2 (INPUT_AUTOCORR*INPUT_AUTOCORR)
#define INPUT_AUTOCORR_4 (INPUT_AUTOCORR_2*INPUT_AUTOCORR_2)
#define INPUT_AUTOCORR_8 (INPUT_AUTOCORR_4*INPUT_AUTOCORR_4)

/*An autocorrelation table.
  A common model for natural-image input is an AR-0 process with an
   autocorrelation coefficient of 0.95.
  This table contains various powers of 0.95, so that
   AUTOCORR[i-j+15]==(0.95)**abs(i-j), for i,j in [0..15].
  This makes it easy to compute element i,j of the covariance matrix.*/
static const double AUTOCORR[31]={
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_2,
  INPUT_AUTOCORR,
  1,
  INPUT_AUTOCORR,
  INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_2*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2,
  INPUT_AUTOCORR_8*INPUT_AUTOCORR_4*INPUT_AUTOCORR_2*INPUT_AUTOCORR
};



/*State information for the IEEE 1180 pseudo-random number generator.*/
static int ieee1180_rand_x;

/*Initializes the IEEE 1180 random number generator.*/
static void ieee1180_srand(int _seed){
  ieee1180_rand_x=_seed;
}

/*Computes a random number between -l and h, inclusive, accoring to the
   specification in IEEE Std 1180-1990, "IEEE Standard Specifications for the
   Implementations of 8x8 Inverse Discrete Cosine Transform."*/
static int ieee1180_rand(int _l,int _h){
  double x;
  int    i;
  ieee1180_rand_x=(ieee1180_rand_x*1103515245)+12345;
  i=ieee1180_rand_x&0x7FFFFFFE;
  x=i/(double)0x7FFFFFFF*(_l+_h+1);
  return (int)x-_l;
}

static char *ieee1180_meets(double _val,double _limit){
   return fabs(_val)<=_limit?"meets":"FAILS";
}

/*The number of different input ranges.*/
#define IEEE1180_NRANGES (3)

/*The number of blocks of data to generate.*/
#define IEEE1180_NBLOCKS (10000)

static const int IEEE1180_L[IEEE1180_NRANGES]={256,5,300};
static const int IEEE1180_H[IEEE1180_NRANGES]={255,5,300};


/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT4_ISCALE[4]={
  1,1,1,1
};

/*The true forward 4-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT4_BASIS[4][4]={
  {
     0.5,                                 0.5,
     0.5,                                 0.5
  },
  {
     0.65328148243818826392832158671359,  0.27059805007309849219986160268319,
    -0.27059805007309849219986160268319, -0.65328148243818826392832158671359
  },
  {
     0.5,                                -0.5,
    -0.5,                                 0.5
  },
  {
     0.27059805007309849219986160268319, -0.65328148243818826392832158671359,
     0.65328148243818826392832158671359, -0.27059805007309849219986160268319
  },
};

void idct4(double _x[],const double _y[]){
  double t[8];
  int    i;
  int    j;
  for(j=0;j<4;j++){
    t[j]=0;
    for(i=0;i<4;i++)t[j]+=DCT4_BASIS[i][j]*_y[i];
  }
  for(j=0;j<4;j++)_x[j]=t[j];
}

void fdct4(double _x[],const double _y[]){
  double t[8];
  int    i;
  int    j;
  for(j=0;j<4;j++){
    t[j]=0;
    for(i=0;i<4;i++)t[j]+=DCT4_BASIS[j][i]*_y[i];
  }
  for(j=0;j<4;j++)_x[j]=t[j];
}

static void ieee1180_print_results4(long _sumerrs[4][4],long _sumsqerrs[4][4],
 int _maxerr[4][4],int _l,int _h,int _sign){
  double max;
  double total;
  int    m;
  int    i;
  int    j;
  printf("IEEE1180-1990 test results:\n");
  printf("Input range: [%i,%i]\n",-_l,_h);
  printf("Sign: %i\n",_sign);
  printf("Iterations: %i\n\n",IEEE1180_NBLOCKS);
  printf("Peak absolute values of errors:\n");
  for(i=0,m=0;i<4;i++){
    for(j=0;j<4;j++){
      if(_maxerr[i][j]>m)m=_maxerr[i][j];
      printf("%4i",_maxerr[i][j]);
    }
    printf("\n");
  }
  printf("Worst peak error = %i (%s spec limit 1)\n\n",m,
   ieee1180_meets((double)m,1.0));
  printf("Mean square errors:\n");
  max=total=0;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      double err;
      err=_sumsqerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f",err);
      total+=_sumsqerrs[i][j];
      if(max<err)max=err;
    }
    printf("\n");
  }
  printf("Worst pmse = %.6f (%s spec limit 0.06)\n",max,
   ieee1180_meets(max,0.06));
  total/=4*4*(double)IEEE1180_NBLOCKS;
  printf("Overall mse = %.6f (%s spec limit 0.02)\n\n",total,
   ieee1180_meets(total,0.02));
  printf("Mean errors:\n");
  max=total=0;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      double err;
      err=_sumerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f",err);
      total+=_sumerrs[i][j];
      if(err<0)err=-err;
      if(max<err)max=err;
    }
    printf("\n");
  }
  total/=4*4*(double)IEEE1180_NBLOCKS;
  printf("Worst mean error = %.6f (%s spec limit 0.015)\n",max,
   ieee1180_meets(max,0.015));
  printf("Overall mean error = %.6f (%s spec limit 0.0015)\n\n",total,
   ieee1180_meets(total,0.0015));
}

static void ieee1180_test_block4(long _sumerrs[4][4],long _sumsqerrs[4][4],
 int _maxerr[4][4],int _l,int _h,int _sign){
  od_coeff block[4][4];
  od_coeff refcoefs[4][4];
  od_coeff refout[4][4];
  od_coeff testout[4][4];
  double   floatcoefs[4][4];
  int      maxerr;
  int      i;
  int      j;
  for(i=0;i<4;i++)for(j=0;j<4;j++)block[i][j]=ieee1180_rand(_l,_h)*_sign;
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  for(i=0;i<4;i++)od_bin_fdct4(refcoefs[i],block[i]);
  for(j=0;j<4;j++){
    od_coeff x[4];
    for(i=0;i<4;i++)x[i]=refcoefs[i][j];
    od_bin_fdct4(x,x);
    for(i=0;i<4;i++)refcoefs[i][j]=x[i];
  }
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for(i=0;i<4;i++)for(j=0;j<4;j++){
    /*Modification of IEEE1180: inputs to reference iDCT are scaled to match
       the scaling factors introduced by the forward integer transform.*/
    floatcoefs[i][j]=refcoefs[i][j]/(DCT4_ISCALE[i]*DCT4_ISCALE[j]);
  }
  for(i=0;i<4;i++)idct4(floatcoefs[i],floatcoefs[i]);
  for(j=0;j<4;j++){
    double x[4];
    for(i=0;i<4;i++)x[i]=floatcoefs[i][j];
    idct4(x,x);
    for(i=0;i<4;i++)floatcoefs[i][j]=x[i];
  }
  for(i=0;i<4;i++)for(j=0;j<4;j++){
    refout[i][j]=(od_coeff)(floatcoefs[i][j]+0.5);
    if(refout[i][j]>255)refout[i][j]=255;
    else if(refout[i][j]<-256)refout[i][j]=-256;
  }
  for(j=0;j<4;j++){
    od_coeff x[4];
    for(i=0;i<4;i++)x[i]=refcoefs[i][j];
    od_bin_idct4(x,x);
    for(i=0;i<4;i++)testout[i][j]=x[i];
  }
  for(i=0;i<4;i++){
    od_bin_idct4(testout[i],testout[i]);
    for(j=0;j<4;j++){
      if(testout[i][j]>255)testout[i][j]=255;
      else if(testout[i][j]<-256)testout[i][j]=-256;
    }
  }
  for(i=0;i<4;i++)for(j=0;j<4;j++){
    int err;
    err=testout[i][j]-refout[i][j];
    _sumerrs[i][j]+=err;
    _sumsqerrs[i][j]+=err*err;
    if(err<0)err=-err;
    if(_maxerr[i][j]<err)_maxerr[i][j]=err;
  }
  for(i=0,maxerr=0;i<4;i++)for(j=0;j<4;j++){
    int err;
    err=testout[i][j]-refout[i][j];
    if(err<0)err=-err;
    if(err>maxerr)maxerr=err;
  }
  /*if(maxerr>1){
    int u;
    int v;
    printf("Excessive peak error: %i\n",maxerr);
    printf("Input:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",block[u][v]);
      }
      printf("\n");
    }
    printf("Forward transform coefficients:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",refcoefs[u][v]);
      }
      printf("\n");
    }
    printf("Reference inverse:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",refout[u][v]);
      }
      printf("\n");
    }
    printf("Integerized inverse:\n");
    for(u=0;u<4;u++){
      for(v=0;v<4;v++){
        printf("%5i",testout[u][v]);
      }
      printf("\n");
    }
  }*/
}

static void ieee1180_test4(void){
  long sumerrs[4][4];
  long sumsqerrs[4][4];
  int  maxerr[4][4];
  int  i;
  int  j;
  ieee1180_srand(1);
  for(i=0;i<IEEE1180_NRANGES;i++){
    memset(sumerrs,0,sizeof(sumerrs));
    memset(sumsqerrs,0,sizeof(sumsqerrs));
    memset(maxerr,0,sizeof(maxerr));
    for(j=0;j<IEEE1180_NBLOCKS;j++){
      ieee1180_test_block4(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
       IEEE1180_H[i],1);
    }
    ieee1180_print_results4(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
     IEEE1180_H[i],1);
  }
  ieee1180_srand(1);
  for(i=0;i<IEEE1180_NRANGES;i++){
    memset(sumerrs,0,sizeof(sumerrs));
    memset(sumsqerrs,0,sizeof(sumsqerrs));
    memset(maxerr,0,sizeof(maxerr));
    for(j=0;j<IEEE1180_NBLOCKS;j++){
      ieee1180_test_block4(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
       IEEE1180_H[i],-1);
    }
    ieee1180_print_results4(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
     IEEE1180_H[i],-1);
  }
}



static void print_basis4(double _basis[4][4]){
  int i;
  int j;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
       printf("%8.5lf%c",_basis[i][j],j==3?'\n':' ');
    }
  }
}

static void compute_fbasis4(double _basis[4][4]){
  int i;
  int j;
  for(i=0;i<4;i++){
    od_coeff x[4];
    for(j=0;j<4;j++)x[j]=(i==j)<<8;
    od_bin_fdct4(x,x);
    for(j=0;j<4;j++)_basis[j][i]=x[j]/256.0;
  }
}

static void compute_ftrue_basis4(double _basis[4][4]){
  int i;
  int j;
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
      _basis[j][i]=sqrt(2.0/4)*cos((i+0.5)*j*M_PI/4)*DCT4_ISCALE[j];
      if(j==0)_basis[j][i]*=M_SQRT1_2;
    }
  }
}

static double compute_mse4(double _basis[4][4],double _tbasis[4][4]){
  double e[4][4];
  double ret;
  int    i;
  int    j;
  int    k;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      e[i][j]=0;
      for(k=0;k<4;k++){
        e[i][j]+=(_basis[i][k]-_tbasis[i][k])*AUTOCORR[k-j+15];
      }
    }
  }
  ret=0;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      ret+=e[i][j]*(_basis[i][j]-_tbasis[i][j]);
    }
  }
  return ret/4;
}

static void check_bias4(){
  double rtacc[4];
  double facc[4];
  double q8acc[4];
  double q7acc[4];
  int    i;
  int    j;
  for(j=0;j<4;j++)q7acc[j]=q8acc[j]=facc[j]=rtacc[j]=0;
  ieee1180_srand(1);
  for(i=0;i<10000000;i++){
    od_coeff x[4];
    od_coeff x2[4];
    od_coeff y[4];
    od_coeff y2[4];
    for(j=0;j<4;j++)x[j]=ieee1180_rand(255,255);
    od_bin_fdct4(y,x);
    for(j=0;j<4;j++)facc[j]+=y[j];
    od_bin_idct4(x2,y);
    for(j=0;j<4;j++)if(x[j]!=x2[j]){
      printf("Mismatch:\n");
      printf("in:    ");
      for(j=0;j<4;j++)printf(" %i",x[j]);
      printf("\nxform: ");
      for(j=0;j<4;j++)printf(" %i",y[j]);
      printf("\nout:   ");
      for(j=0;j<4;j++)printf(" %i",x2[j]);
      printf("\n\n");
      break;
    }
    for(j=0;j<4;j++)y2[j]=y[j]+ieee1180_rand(1,1);
    od_bin_idct4(x2,y2);
    for(j=0;j<4;j++)rtacc[j]+=x2[j]-x[j];
    for(j=0;j<4;j++)y2[j]=y[j]/8<<3;
    od_bin_idct4(x2,y2);
    for(j=0;j<4;j++)q8acc[j]+=x2[j]-x[j];
    for(j=0;j<4;j++)y2[j]=y[j]/7*7;
    od_bin_idct4(x2,y2);
    for(j=0;j<4;j++)q7acc[j]+=x2[j]-x[j];
  }
  printf("1-D Forward Bias:\n");
  for(j=0;j<4;j++)printf("% -18.15G%s",facc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
  printf("1-D Round-Trip Bias:\n");
  for(j=0;j<4;j++)printf("% -18.15G%s",rtacc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
  printf("1-D Q=8 Bias:\n");
  for(j=0;j<4;j++)printf("% -18.15G%s",q8acc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
  printf("1-D Q=7 Bias:\n");
  for(j=0;j<4;j++)printf("% -18.15G%s",q7acc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
}

static void check4(void){
  od_coeff min[4];
  od_coeff max[4];
  double   basis[4][4];
  double   tbasis[4][4];
  int      i;
  int      j;
  for(j=0;j<4;j++)min[j]=max[j]=0;
  for(i=0;i<1<<4;i++){
    od_coeff x[4];
    od_coeff y[4];
    od_coeff x2[4];
    for(j=0;j<4;j++)x[j]=(i>>j&1)?255:-256;
    od_bin_fdct4(y,x);
    od_bin_idct4(x2,y);
    for(j=0;j<4;j++){
      if(y[j]<min[j])min[j]=y[j];
      else if(y[j]>max[j])max[j]=y[j];
    }
    for(j=0;j<4;j++)if(x[j]!=x2[j]){
      printf("Mismatch:\n");
      printf("in:    ");
      for(j=0;j<4;j++)printf(" %i",x[j]);
      printf("\nxform: ");
      for(j=0;j<4;j++)printf(" %i",y[j]);
      printf("\nout:   ");
      for(j=0;j<4;j++)printf(" %i",x2[j]);
      printf("\n\n");
      break;
    }
  }
  printf("Min:");
  for(j=0;j<4;j++)printf(" %5i",min[j]);
  printf("\nMax:");
  for(j=0;j<4;j++)printf(" %5i",max[j]);
  printf("\nod_bin_fdct4 basis:\n");
  compute_fbasis4(basis);
  print_basis4(basis);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis4(tbasis);
  print_basis4(tbasis);
  printf("MSE: %.32lg\n\n",compute_mse4(basis,tbasis));
  ieee1180_test4();
  check_bias4();
}


/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT8_ISCALE[8]={
  M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2,M_SQRT2
};

/*The true forward 8-point type-II DCT basis, to 32-digit (100 bit) precision.
  The inverse is merely the transpose.*/
static const double DCT8_BASIS[8][8]={
  {
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242,  0.35355339059327376220042218105242
  },
  {
     0.49039264020161522456309111806712,  0.41573480615127261853939418880895,
     0.27778511650980111237141540697427,  0.097545161008064133924142434238511,
    -0.097545161008064133924142434238511,-0.27778511650980111237141540697427,
    -0.41573480615127261853939418880895, -0.49039264020161522456309111806712
  },
  {
     0.46193976625574337806409159469839,  0.19134171618254488586422999201520,
    -0.19134171618254488586422999201520, -0.46193976625574337806409159469839,
    -0.46193976625574337806409159469839, -0.19134171618254488586422999201520,
     0.19134171618254488586422999201520,  0.46193976625574337806409159469839
  },
  {
     0.41573480615127261853939418880895, -0.097545161008064133924142434238511,
    -0.49039264020161522456309111806712, -0.27778511650980111237141540697427,
     0.27778511650980111237141540697427,  0.49039264020161522456309111806712,
     0.097545161008064133924142434238511,-0.41573480615127261853939418880895
  },
  {
     0.35355339059327376220042218105242, -0.35355339059327376220042218105242,
    -0.35355339059327376220042218105242,  0.35355339059327376220042218105242,
     0.35355339059327376220042218105242, -0.35355339059327376220042218105242,
    -0.35355339059327376220042218105242,  0.35355339059327376220042218105242
  },
  {
     0.27778511650980111237141540697427, -0.49039264020161522456309111806712,
     0.097545161008064133924142434238511, 0.41573480615127261853939418880895,
    -0.41573480615127261853939418880895, -0.097545161008064133924142434238511,
     0.49039264020161522456309111806712, -0.27778511650980111237141540697427
  },
  {
     0.19134171618254488586422999201520, -0.46193976625574337806409159469839,
     0.46193976625574337806409159469839, -0.19134171618254488586422999201520,
    -0.19134171618254488586422999201520,  0.46193976625574337806409159469839,
    -0.46193976625574337806409159469839,  0.19134171618254488586422999201520
  },
  {
     0.097545161008064133924142434238511,-0.27778511650980111237141540697427,
     0.41573480615127261853939418880895, -0.49039264020161522456309111806712,
     0.49039264020161522456309111806712, -0.41573480615127261853939418880895,
     0.27778511650980111237141540697427, -0.097545161008064133924142434238511
  }
};

void idct8(double _x[],const double _y[]){
  double t[8];
  int    i;
  int    j;
  for(j=0;j<8;j++){
    t[j]=0;
    for(i=0;i<8;i++)t[j]+=DCT8_BASIS[i][j]*_y[i];
  }
  for(j=0;j<8;j++)_x[j]=t[j];
}

static void ieee1180_print_results8(long _sumerrs[8][8],long _sumsqerrs[8][8],
 int _maxerr[8][8],int _l,int _h,int _sign){
  double max;
  double total;
  int    m;
  int    i;
  int    j;
  printf("IEEE1180-1990 test results:\n");
  printf("Input range: [%i,%i]\n",-_l,_h);
  printf("Sign: %i\n",_sign);
  printf("Iterations: %i\n\n",IEEE1180_NBLOCKS);
  printf("Peak absolute values of errors:\n");
  for(i=0,m=0;i<8;i++){
    for(j=0;j<8;j++){
      if(_maxerr[i][j]>m)m=_maxerr[i][j];
      printf("%4i",_maxerr[i][j]);
    }
    printf("\n");
  }
  printf("Worst peak error = %i (%s spec limit 1)\n\n",m,
   ieee1180_meets((double)m,1.0));
  printf("Mean square errors:\n");
  max=total=0;
  for(i=0;i<8;i++){
    for(j=0;j<8;j++){
      double err;
      err=_sumsqerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f",err);
      total+=_sumsqerrs[i][j];
      if(max<err)max=err;
    }
    printf("\n");
  }
  printf("Worst pmse = %.6f (%s spec limit 0.06)\n",max,
   ieee1180_meets(max,0.06));
  total/=8*8*(double)IEEE1180_NBLOCKS;
  printf("Overall mse = %.6f (%s spec limit 0.02)\n\n",total,
   ieee1180_meets(total,0.02));
  printf("Mean errors:\n");
  max=total=0;
  for(i=0;i<8;i++){
    for(j=0;j<8;j++){
      double err;
      err=_sumerrs[i][j]/(double)IEEE1180_NBLOCKS;
      printf(" %8.4f",err);
      total+=_sumerrs[i][j];
      if(err<0)err=-err;
      if(max<err)max=err;
    }
    printf("\n");
  }
  total/=8*8*(double)IEEE1180_NBLOCKS;
  printf("Worst mean error = %.6f (%s spec limit 0.015)\n",max,
   ieee1180_meets(max,0.015));
  printf("Overall mean error = %.6f (%s spec limit 0.0015)\n\n",total,
   ieee1180_meets(total,0.0015));
}

static void ieee1180_test_block8(long _sumerrs[8][8],long _sumsqerrs[8][8],
 int _maxerr[8][8],int _l,int _h,int _sign){
  od_coeff block[8][8];
  od_coeff refcoefs[8][8];
  od_coeff refout[8][8];
  od_coeff testout[8][8];
  double   floatcoefs[8][8];
  int      maxerr;
  int      i;
  int      j;
  for(i=0;i<8;i++)for(j=0;j<8;j++)block[i][j]=ieee1180_rand(_l,_h)*_sign;
  /*Modification of IEEE1180: use our integerized DCT, not a true DCT.*/
  for(i=0;i<8;i++)od_bin_fdct8(refcoefs[i],block[i]);
  for(j=0;j<8;j++){
    od_coeff x[8];
    for(i=0;i<8;i++)x[i]=refcoefs[i][j];
    od_bin_fdct8(x,x);
    for(i=0;i<8;i++)refcoefs[i][j]=x[i];
  }
  /*Modification of IEEE1180: no rounding or range clipping (coefficients
     are always in range with our integerized DCT).*/
  for(i=0;i<8;i++)for(j=0;j<8;j++){
    /*Modification of IEEE1180: inputs to reference iDCT are scaled to match
       the scaling factors introduced by the forward integer transform.*/
    floatcoefs[i][j]=refcoefs[i][j]/(DCT8_ISCALE[i]*DCT8_ISCALE[j]);
  }
  for(i=0;i<8;i++)idct8(floatcoefs[i],floatcoefs[i]);
  for(j=0;j<8;j++){
    double x[8];
    for(i=0;i<8;i++)x[i]=floatcoefs[i][j];
    idct8(x,x);
    for(i=0;i<8;i++)floatcoefs[i][j]=x[i];
  }
  for(i=0;i<8;i++)for(j=0;j<8;j++){
    refout[i][j]=(od_coeff)(floatcoefs[i][j]+0.5);
    if(refout[i][j]>255)refout[i][j]=255;
    else if(refout[i][j]<-256)refout[i][j]=-256;
  }
  for(j=0;j<8;j++){
    od_coeff x[8];
    for(i=0;i<8;i++)x[i]=refcoefs[i][j];
    od_bin_idct8(x,x);
    for(i=0;i<8;i++)testout[i][j]=x[i];
  }
  for(i=0;i<8;i++){
    od_bin_idct8(testout[i],testout[i]);
    for(j=0;j<8;j++){
      if(testout[i][j]>255)testout[i][j]=255;
      else if(testout[i][j]<-256)testout[i][j]=-256;
    }
  }
  for(i=0;i<8;i++)for(j=0;j<8;j++){
    int err;
    err=testout[i][j]-refout[i][j];
    _sumerrs[i][j]+=err;
    _sumsqerrs[i][j]+=err*err;
    if(err<0)err=-err;
    if(_maxerr[i][j]<err)_maxerr[i][j]=err;
  }
  for(i=0,maxerr=0;i<8;i++)for(j=0;j<8;j++){
    int err;
    err=testout[i][j]-refout[i][j];
    if(err<0)err=-err;
    if(err>maxerr)maxerr=err;
  }
  /*if(maxerr>1){
    int u;
    int v;
    printf("Excessive peak error: %i\n",maxerr);
    printf("Input:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",block[u][v]);
      }
      printf("\n");
    }
    printf("Forward transform coefficients:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",refcoefs[u][v]);
      }
      printf("\n");
    }
    printf("Reference inverse:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",refout[u][v]);
      }
      printf("\n");
    }
    printf("Integerized inverse:\n");
    for(u=0;u<8;u++){
      for(v=0;v<8;v++){
        printf("%5i",testout[u][v]);
      }
      printf("\n");
    }
  }*/
}

static void ieee1180_test8(void){
  long sumerrs[8][8];
  long sumsqerrs[8][8];
  int  maxerr[8][8];
  int  i;
  int  j;
  ieee1180_srand(1);
  for(i=0;i<IEEE1180_NRANGES;i++){
    memset(sumerrs,0,sizeof(sumerrs));
    memset(sumsqerrs,0,sizeof(sumsqerrs));
    memset(maxerr,0,sizeof(maxerr));
    for(j=0;j<IEEE1180_NBLOCKS;j++){
      ieee1180_test_block8(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
       IEEE1180_H[i],1);
    }
    ieee1180_print_results8(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
     IEEE1180_H[i],1);
  }
  ieee1180_srand(1);
  for(i=0;i<IEEE1180_NRANGES;i++){
    memset(sumerrs,0,sizeof(sumerrs));
    memset(sumsqerrs,0,sizeof(sumsqerrs));
    memset(maxerr,0,sizeof(maxerr));
    for(j=0;j<IEEE1180_NBLOCKS;j++){
      ieee1180_test_block8(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
       IEEE1180_H[i],-1);
    }
    ieee1180_print_results8(sumerrs,sumsqerrs,maxerr,IEEE1180_L[i],
     IEEE1180_H[i],-1);
  }
}

static void print_basis8(double _basis[8][8]){
  int i;
  int j;
  for(i=0;i<8;i++){
    for(j=0;j<8;j++){
      printf("%8.5lf%c",_basis[i][j],j==8-1?'\n':' ');
    }
  }
}

static void compute_fbasis8(double _basis[8][8]){
  int i;
  int j;
  for(i=0;i<8;i++){
    od_coeff x[8];
    for(j=0;j<8;j++)x[j]=(i==j)*256;
    od_bin_fdct8(x,x);
    for(j=0;j<8;j++)_basis[j][i]=x[j]/256.0;
  }
}

static void compute_ibasis8(double _basis[8][8]){
  int i;
  int j;
  for(i=0;i<8;i++){
    od_coeff x[8];
    for(j=0;j<8;j++)x[j]=(i==j)*256;
    od_bin_idct8(x,x);
    for(j=0;j<8;j++)_basis[j][i]=x[j]/256.0;
  }
}

static void compute_ftrue_basis8(double _basis[8][8]){
  int i;
  int j;
  for(j=0;j<8;j++){
    for(i=0;i<8;i++){
      _basis[j][i]=sqrt(2.0/8)*cos((i+0.5)*j*M_PI/8)*DCT8_ISCALE[j];
      if(j==0)_basis[j][i]*=M_SQRT1_2;
    }
  }
}

static void compute_itrue_basis8(double _basis[8][8]){
  int i;
  int j;
  for(i=0;i<8;i++){
    double x[8];
    for(j=0;j<8;j++)x[j]=(i==j);
    idct8(x,x);
    for(j=0;j<8;j++)_basis[j][i]=x[j]/DCT8_ISCALE[i];
  }
}

static double compute_mse8(double _basis[8][8],double _tbasis[8][8]){
  double e[8][8];
  double ret;
  int    i;
  int    j;
  int    k;
  for(i=0;i<8;i++){
    for(j=0;j<8;j++){
      e[i][j]=0;
      for(k=0;k<8;k++){
        e[i][j]+=(_basis[i][k]-_tbasis[i][k])*AUTOCORR[k-j+15];
      }
    }
  }
  ret=0;
  for(i=0;i<8;i++){
    for(j=0;j<8;j++){
      ret+=e[i][j]*(_basis[i][j]-_tbasis[i][j]);
    }
  }
  return ret/8;
}

static void check_bias8(){
  double rtacc[8];
  double facc[8];
  double q8acc[8];
  double q7acc[8];
  int    i;
  int    j;
  for(j=0;j<8;j++)q7acc[j]=q8acc[j]=facc[j]=rtacc[j]=0;
  ieee1180_srand(1);
  for(i=0;i<10000000;i++){
    od_coeff x[8];
    od_coeff x2[8];
    od_coeff y[8];
    od_coeff y2[8];
    for(j=0;j<8;j++)x[j]=ieee1180_rand(255,255);
    od_bin_fdct8(y,x);
    for(j=0;j<8;j++)facc[j]+=y[j];
    od_bin_idct8(x2,y);
    for(j=0;j<8;j++)if(x[j]!=x2[j]){
      printf("Mismatch:\n");
      printf("in:    ");
      for(j=0;j<8;j++)printf(" %i",x[j]);
      printf("\nxform: ");
      for(j=0;j<8;j++)printf(" %i",y[j]);
      printf("\nout:   ");
      for(j=0;j<8;j++)printf(" %i",x2[j]);
      printf("\n\n");
      break;
    }
    for(j=0;j<8;j++)y2[j]=y[j]+ieee1180_rand(1,1);
    od_bin_idct8(x2,y2);
    for(j=0;j<8;j++)rtacc[j]+=x2[j]-x[j];
    for(j=0;j<8;j++)y2[j]=y[j]/8<<3;
    od_bin_idct8(x2,y2);
    for(j=0;j<8;j++)q8acc[j]+=x2[j]-x[j];
    for(j=0;j<8;j++)y2[j]=y[j]/7*7;
    od_bin_idct8(x2,y2);
    for(j=0;j<8;j++)q7acc[j]+=x2[j]-x[j];
  }
  printf("1-D Forward Bias:\n");
  for(j=0;j<8;j++)printf("% -18.15G%s",facc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
  printf("1-D Round-Trip Bias:\n");
  for(j=0;j<8;j++)printf("% -18.15G%s",rtacc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
  printf("1-D Q=8 Bias:\n");
  for(j=0;j<8;j++)printf("% -18.15G%s",q8acc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
  printf("1-D Q=7 Bias:\n");
  for(j=0;j<8;j++)printf("% -18.15G%s",q7acc[j]/i,(j&3)==3?"\n":"  ");
  printf("\n");
}

#if 0
static void bin_fxform_2d8(od_coeff _x[8*2][8*2]){
  od_coeff y[8*2];
  int      u;
  int      v;
  /*Perform pre-filtering.*/
  for(u=0;u<8*2;u++){
    od_pre_filter8(_x[u],_x[u]);
    od_pre_filter8(_x[u]+8,_x[u]+8);
  }
  for(v=0;v<8*2;v++){
    for(u=0;u<8*2;u++)y[u]=_x[u][v];
    od_pre_filter8(y,y);
    od_pre_filter8(y+8,y+8);
    for(u=0;u<8*2;u++)_x[u][v]=y[u];
  }
  /*Perform DCT.*/
  for(u=8/2;u<8*3/2;u++)od_bin_fdct8(_x[u]+8/2,_x[u]+8/2);
  for(v=8/2;v<8*3/2;v++){
    for(u=8/2;u<8*3/2;u++)y[u]=_x[u][v];
    od_bin_fdct8(y+8/2,y+8/2);
    for(u=8/2;u<8*3/2;u++)_x[u][v]=y[u];
  }
}

static void dynamic_range8(void){
  double   basis2[8][8][8*2][8*2];
  od_coeff min2[8][8];
  od_coeff max2[8][8];
  int      i;
  int      j;
  int      u;
  int      v;
  for(i=0;i<8*2;i++){
    for(j=0;j<8*2;j++){
      od_coeff x[8*2][8*2];
      /*Generate impulse.*/
      for(u=0;u<8*2;u++){
        for(v=0;v<8*2;v++){
          x[u][v]=(u==i&&v==j)*256;
        }
      }
      bin_fxform_2d8(x);
      /*Retrieve basis elements.*/
      for(u=0;u<8;u++){
        for(v=0;v<8;v++){
          basis2[u][v][i][j]=x[u+8/2][v+8/2]/256.0;
        }
      }
    }
  }
  for(u=0;u<8;u++){
    for(v=0;v<8;v++){
      od_coeff x[8*2][8*2];
      for(i=0;i<8*2;i++){
        for(j=0;j<8*2;j++){
          x[i][j]=basis2[u][v][i][j]<0?-255:255;
        }
      }
      bin_fxform_2d8(x);
      max2[u][v]=x[u+8/2][v+8/2];
      for(i=0;i<8*2;i++){
        for(j=0;j<8*2;j++){
          x[i][j]=basis2[u][v][i][j]>0?-255:255;
        }
      }
      bin_fxform_2d8(x);
      min2[u][v]=x[u+8/2][v+8/2];
    }
  }
  printf("2-D ranges:\n");
  for(u=0;u<8;u++){
    printf("Min %2i:",u);
    for(v=0;v<8;v++)printf(" %6i",min2[u][v]);
    printf("\nMax %2i:",u);
    for(v=0;v<8;v++)printf(" %6i",max2[u][v]);
    printf("\n");
  }
}
#endif

static void check8(void){
  od_coeff min[8];
  od_coeff max[8];
  double   basis[8][8];
  double   tbasis[8][8];
  int      i;
  int      j;
  /*dynamic_range8();*/
  for(j=0;j<8;j++)min[j]=max[j]=0;
  for(i=0;i<1<<8;i++){
    od_coeff x[8];
    od_coeff y[8];
    od_coeff x2[8];
    for(j=0;j<8;j++)x[j]=(i>>j&1)?255:-256;
    od_bin_fdct8(y,x);
    od_bin_idct8(x2,y);
    for(j=0;j<8;j++){
      if(y[j]<min[j])min[j]=y[j];
      else if(y[j]>max[j])max[j]=y[j];
    }
    for(j=0;j<8;j++)if(x[j]!=x2[j]){
      printf("Mismatch:\n");
      printf("in:    ");
      for(j=0;j<8;j++)printf(" %i",x[j]);
      printf("\nxform: ");
      for(j=0;j<8;j++)printf(" %i",y[j]);
      printf("\nout:   ");
      for(j=0;j<8;j++)printf(" %i",x2[j]);
      printf("\n\n");
      break;
    }
  }
  printf("Min:");
  for(j=0;j<8;j++)printf(" %5i",min[j]);
  printf("\nMax:");
  for(j=0;j<8;j++)printf(" %5i",max[j]);
  printf("\nod_bin_idct8 basis:\n");
  compute_ibasis8(basis);
  print_basis8(basis);
  printf("Scaled type-II iDCT basis:\n");
  compute_itrue_basis8(tbasis);
  print_basis8(tbasis);
  printf("\nod_bin_fdct8 basis:\n");
  compute_fbasis8(basis);
  print_basis8(basis);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis8(tbasis);
  print_basis8(tbasis);
  printf("MSE: %.32lg\n\n",compute_mse8(basis,tbasis));
  ieee1180_test8();
  check_bias8();
}


int main(void){
  check4();
  check8();
  return 0;
}
#endif
