/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "od_defs.h"
#include "od_filter.h"

int NE_FILTER_PARAMS4[3];
int NE_FILTER_PARAMS8[10];
int NE_FILTER_PARAMS16[22];

void ne_filter_params_init(){
  int i;
  for(i=0;i<3;i++){
    NE_FILTER_PARAMS4[i]=OD_FILTER_PARAMS4[i];
  }
  for(i=0;i<10;i++){
    NE_FILTER_PARAMS8[i]=OD_FILTER_PARAMS8[i];
  }
  for(i=0;i<22;i++){
    NE_FILTER_PARAMS16[i]=OD_FILTER_PARAMS16[i];
  }
}

const od_filter_func NE_PRE_FILTER[OD_NBSIZES]={
  od_pre_filter4,
  od_pre_filter8,
  od_pre_filter16
};

const od_filter_func NE_POST_FILTER[OD_NBSIZES]={
  od_post_filter4,
  od_post_filter8,
  od_post_filter16
};

static void ne_pre_filter4_double(double _y[4],const double _x[4],
 const int _f[3]){
  double t[4];
  t[3]= _x[0] - _x[3];
  t[2]= _x[1] - _x[2];
  t[1]= _x[1] - (t[2]/2);
  t[0]= _x[0] - (t[3]/2);
  t[3] += t[2]*_f[0]/(1 << FILTER_BITS);
  t[2] = t[3]*_f[1]/(1 << FILTER_BITS) - t[2];
  t[3] += t[2]*_f[2]/(1 << FILTER_BITS);
  t[0] += t[2]/2;
  _y[0] = t[0];
  t[1] += t[3]/2;
  _y[1] = t[1];
  _y[2] = (t[1] - t[3]);
  _y[3] = (t[0] - t[2]);
}

static void ne_post_filter4_double(double _x[4],const double _y[4],
 const int _f[3]){
  double t[4];
  t[2] = _y[0] - _y[3];
  t[3] = _y[1] - _y[2];
  t[1] = _y[1] - (t[3]/2);
  t[0] = _y[0] - (t[2]/2);
  t[3] -= t[2]*_f[2]/(1 << FILTER_BITS);
  t[2] = (t[3]*_f[1]/(1 << FILTER_BITS)) - t[2];
  t[3] -= t[2]*_f[0]/(1 << FILTER_BITS);
  t[0] += t[3]/2;
  _x[0] = t[0];
  t[1] += t[2]/2;
  _x[1] = t[1];
  _x[2] = (t[1] - t[2]);
  _x[3] = (t[0] - t[3]);
}

static void ne_pre_filter8_double(double _y[8],const double _x[8],
 const int _f[10]){
  double t[8];
  t[7]=_x[0]-_x[7];
  t[6]=_x[1]-_x[6];
  t[5]=_x[2]-_x[5];
  t[4]=_x[3]-_x[4];
  t[3]=_x[3]-(t[4]/2);
  t[2]=_x[2]-(t[5]/2);
  t[1]=_x[1]-(t[6]/2);
  t[0]=_x[0]-(t[7]/2);
  t[4]=t[4]*_f[0]/(1<<FILTER_BITS);
  t[5]=t[5]*_f[1]/(1<<FILTER_BITS);
  t[6]=t[6]*_f[2]/(1<<FILTER_BITS);
  t[7]=t[7]*_f[3]/(1<<FILTER_BITS);
#if USE_TYPE3
  t[7]+=t[6]*_f[6]/(1<<FILTER_BITS);
  t[6]+=t[7]*_f[9]/(1<<FILTER_BITS);
  t[6]+=t[5]*_f[5]/(1<<FILTER_BITS);
  t[5]+=t[6]*_f[8]/(1<<FILTER_BITS);
  t[5]+=t[4]*_f[4]/(1<<FILTER_BITS);
  t[4]+=t[5]*_f[7]/(1<<FILTER_BITS);
#else
  t[5]+=t[4]*_f[4]/(1<<FILTER_BITS);
  t[6]+=t[5]*_f[5]/(1<<FILTER_BITS);
  t[7]+=t[6]*_f[6]/(1<<FILTER_BITS);
  t[6]+=t[7]*_f[9]/(1<<FILTER_BITS);
  t[5]+=t[6]*_f[8]/(1<<FILTER_BITS);
  t[4]+=t[5]*_f[7]/(1<<FILTER_BITS);
#endif
  t[0]+=t[7]/2;
  _y[0]=t[0];
  t[1]+=t[6]/2;
  _y[1]=t[1];
  t[2]+=t[5]/2;
  _y[2]=t[2];
  t[3]+=t[4]/2;
  _y[3]=t[3];
  _y[4]=(t[3]-t[4]);
  _y[5]=(t[2]-t[5]);
  _y[6]=(t[1]-t[6]);
  _y[7]=(t[0]-t[7]);
}

static void ne_post_filter8_double(double _x[8],const double _y[8],
 const int _f[10]){
  double t[8];
  t[7]=_y[0]-_y[7];
  t[6]=_y[1]-_y[6];
  t[5]=_y[2]-_y[5];
  t[4]=_y[3]-_y[4];
  t[3]=_y[3]-(t[4]/2);
  t[2]=_y[2]-(t[5]/2);
  t[1]=_y[1]-(t[6]/2);
  t[0]=_y[0]-(t[7]/2);
#if USE_TYPE3
  t[4]-=t[5]*_f[7]/(1<<FILTER_BITS);
  t[5]-=t[4]*_f[4]/(1<<FILTER_BITS);
  t[5]-=t[6]*_f[8]/(1<<FILTER_BITS);
  t[6]-=t[5]*_f[5]/(1<<FILTER_BITS);
  t[6]-=t[7]*_f[9]/(1<<FILTER_BITS);
  t[7]-=t[6]*_f[6]/(1<<FILTER_BITS);
#else
  t[4]-=t[5]*_f[7]/(1<<FILTER_BITS);
  t[5]-=t[6]*_f[8]/(1<<FILTER_BITS);
  t[6]-=t[7]*_f[9]/(1<<FILTER_BITS);
  t[7]-=t[6]*_f[6]/(1<<FILTER_BITS);
  t[6]-=t[5]*_f[5]/(1<<FILTER_BITS);
  t[5]-=t[4]*_f[4]/(1<<FILTER_BITS);
#endif
  t[7]=t[7]*(1<<FILTER_BITS)/_f[3];
  t[6]=t[6]*(1<<FILTER_BITS)/_f[2];
  t[5]=t[5]*(1<<FILTER_BITS)/_f[1];
  t[4]=t[4]*(1<<FILTER_BITS)/_f[0];
  t[0]+=t[7]/2;
  _x[0]=t[0];
  t[1]+=t[6]/2;
  _x[1]=t[1];
  t[2]+=t[5]/2;
  _x[2]=t[2];
  t[3]+=t[4]/2;
  _x[3]=t[3];
  _x[4]=(t[3]-t[4]);
  _x[5]=(t[2]-t[5]);
  _x[6]=(t[1]-t[6]);
  _x[7]=(t[0]-t[7]);
}

static void ne_pre_filter16_double(double _y[16],const double _x[16],
 const int _f[22]){
  double t[16];
  t[15]=_x[0]-_x[15];
  t[14]=_x[1]-_x[14];
  t[13]=_x[2]-_x[13];
  t[12]=_x[3]-_x[12];
  t[11]=_x[4]-_x[11];
  t[10]=_x[5]-_x[10];
  t[9]=_x[6]-_x[9];
  t[8]=_x[7]-_x[8];
  t[7]=_x[7]-(t[8]/2);
  t[6]=_x[6]-(t[9]/2);
  t[5]=_x[5]-(t[10]/2);
  t[4]=_x[4]-(t[11]/2);
  t[3]=_x[3]-(t[12]/2);
  t[2]=_x[2]-(t[13]/2);
  t[1]=_x[1]-(t[14]/2);
  t[0]=_x[0]-(t[15]/2);
  t[8]=t[8]*_f[0]/(1<<FILTER_BITS);
  t[9]=t[9]*_f[1]/(1<<FILTER_BITS);
  t[10]=t[10]*_f[2]/(1<<FILTER_BITS);
  t[11]=t[11]*_f[3]/(1<<FILTER_BITS);
  t[12]=t[12]*_f[4]/(1<<FILTER_BITS);
  t[13]=t[13]*_f[5]/(1<<FILTER_BITS);
  t[14]=t[14]*_f[6]/(1<<FILTER_BITS);
  t[15]=t[15]*_f[7]/(1<<FILTER_BITS);
#if USE_TYPE3
  t[15]+=t[14]*_f[14]/(1<<FILTER_BITS);
  t[14]+=t[15]*_f[21]/(1<<FILTER_BITS);
  t[14]+=t[13]*_f[13]/(1<<FILTER_BITS);
  t[13]+=t[14]*_f[20]/(1<<FILTER_BITS);
  t[13]+=t[12]*_f[12]/(1<<FILTER_BITS);
  t[12]+=t[13]*_f[19]/(1<<FILTER_BITS);
  t[12]+=t[11]*_f[11]/(1<<FILTER_BITS);
  t[11]+=t[12]*_f[18]/(1<<FILTER_BITS);
  t[11]+=t[10]*_f[10]/(1<<FILTER_BITS);
  t[10]+=t[11]*_f[17]/(1<<FILTER_BITS);
  t[10]+=t[9]*_f[9]/(1<<FILTER_BITS);
  t[9]+=t[10]*_f[16]/(1<<FILTER_BITS);
  t[9]+=t[8]*_f[8]/(1<<FILTER_BITS);
  t[8]+=t[9]*_f[15]/(1<<FILTER_BITS);
#else
  t[9]+=t[8]*_f[8]/(1<<FILTER_BITS);
  t[10]+=t[9]*_f[9]/(1<<FILTER_BITS);
  t[11]+=t[10]*_f[10]/(1<<FILTER_BITS);
  t[12]+=t[11]*_f[11]/(1<<FILTER_BITS);
  t[13]+=t[12]*_f[12]/(1<<FILTER_BITS);
  t[14]+=t[13]*_f[13]/(1<<FILTER_BITS);
  t[15]+=t[14]*_f[14]/(1<<FILTER_BITS);
  t[14]+=t[15]*_f[21]/(1<<FILTER_BITS);
  t[13]+=t[14]*_f[20]/(1<<FILTER_BITS);
  t[12]+=t[13]*_f[19]/(1<<FILTER_BITS);
  t[11]+=t[12]*_f[18]/(1<<FILTER_BITS);
  t[10]+=t[11]*_f[17]/(1<<FILTER_BITS);
  t[9]+=t[10]*_f[16]/(1<<FILTER_BITS);
  t[8]+=t[9]*_f[15]/(1<<FILTER_BITS);
#endif
  t[0]+=t[15]/2;
  _y[0]=t[0];
  t[1]+=t[14]/2;
  _y[1]=t[1];
  t[2]+=t[13]/2;
  _y[2]=t[2];
  t[3]+=t[12]/2;
  _y[3]=t[3];
  t[4]+=t[11]/2;
  _y[4]=t[4];
  t[5]+=t[10]/2;
  _y[5]=t[5];
  t[6]+=t[9]/2;
  _y[6]=t[6];
  t[7]+=t[8]/2;
  _y[7]=t[7];
  _y[8]=(t[7]-t[8]);
  _y[9]=(t[6]-t[9]);
  _y[10]=(t[5]-t[10]);
  _y[11]=(t[4]-t[11]);
  _y[12]=(t[3]-t[12]);
  _y[13]=(t[2]-t[13]);
  _y[14]=(t[1]-t[14]);
  _y[15]=(t[0]-t[15]);
}

static void ne_post_filter16_double(double _x[16],const double _y[16],
 const int _f[22]){
  double t[16];
  t[15]=_y[0]-_y[15];
  t[14]=_y[1]-_y[14];
  t[13]=_y[2]-_y[13];
  t[12]=_y[3]-_y[12];
  t[11]=_y[4]-_y[11];
  t[10]=_y[5]-_y[10];
  t[9]=_y[6]-_y[9];
  t[8]=_y[7]-_y[8];
  t[7]=_y[7]-(t[8]/2);
  t[6]=_y[6]-(t[9]/2);
  t[5]=_y[5]-(t[10]/2);
  t[4]=_y[4]-(t[11]/2);
  t[3]=_y[3]-(t[12]/2);
  t[2]=_y[2]-(t[13]/2);
  t[1]=_y[1]-(t[14]/2);
  t[0]=_y[0]-(t[15]/2);
#if USE_TYPE3
  t[8]-=t[9]*_f[15]/(1<<FILTER_BITS);
  t[9]-=t[8]*_f[8]/(1<<FILTER_BITS);
  t[9]-=t[10]*_f[16]/(1<<FILTER_BITS);
  t[10]-=t[9]*_f[9]/(1<<FILTER_BITS);
  t[10]-=t[11]*_f[17]/(1<<FILTER_BITS);
  t[11]-=t[10]*_f[10]/(1<<FILTER_BITS);
  t[11]-=t[12]*_f[18]/(1<<FILTER_BITS);
  t[12]-=t[11]*_f[11]/(1<<FILTER_BITS);
  t[12]-=t[13]*_f[19]/(1<<FILTER_BITS);
  t[13]-=t[12]*_f[12]/(1<<FILTER_BITS);
  t[13]-=t[14]*_f[20]/(1<<FILTER_BITS);
  t[14]-=t[13]*_f[13]/(1<<FILTER_BITS);
  t[14]-=t[15]*_f[21]/(1<<FILTER_BITS);
  t[15]-=t[14]*_f[14]/(1<<FILTER_BITS);
#else
  t[8]-=t[9]*_f[15]/(1<<FILTER_BITS);
  t[9]-=t[10]*_f[16]/(1<<FILTER_BITS);
  t[10]-=t[11]*_f[17]/(1<<FILTER_BITS);
  t[11]-=t[12]*_f[18]/(1<<FILTER_BITS);
  t[12]-=t[13]*_f[19]/(1<<FILTER_BITS);
  t[13]-=t[14]*_f[20]/(1<<FILTER_BITS);
  t[14]-=t[15]*_f[21]/(1<<FILTER_BITS);
  t[15]-=t[14]*_f[14]/(1<<FILTER_BITS);
  t[14]-=t[13]*_f[13]/(1<<FILTER_BITS);
  t[13]-=t[12]*_f[12]/(1<<FILTER_BITS);
  t[12]-=t[11]*_f[11]/(1<<FILTER_BITS);
  t[11]-=t[10]*_f[10]/(1<<FILTER_BITS);
  t[10]-=t[9]*_f[9]/(1<<FILTER_BITS);
  t[9]-=t[8]*_f[8]/(1<<FILTER_BITS);
#endif
  t[15]=t[15]*(1<<FILTER_BITS)/_f[7];
  t[14]=t[14]*(1<<FILTER_BITS)/_f[6];
  t[13]=t[13]*(1<<FILTER_BITS)/_f[5];
  t[12]=t[12]*(1<<FILTER_BITS)/_f[4];
  t[11]=t[11]*(1<<FILTER_BITS)/_f[3];
  t[10]=t[10]*(1<<FILTER_BITS)/_f[2];
  t[9]=t[9]*(1<<FILTER_BITS)/_f[1];
  t[8]=t[8]*(1<<FILTER_BITS)/_f[0];
  t[0]+=t[15]/2;
  _x[0]=t[0];
  t[1]+=t[14]/2;
  _x[1]=t[1];
  t[2]+=t[13]/2;
  _x[2]=t[2];
  t[3]+=t[12]/2;
  _x[3]=t[3];
  t[4]+=t[11]/2;
  _x[4]=t[4];
  t[5]+=t[10]/2;
  _x[5]=t[5];
  t[6]+=t[9]/2;
  _x[6]=t[6];
  t[7]+=t[8]/2;
  _x[7]=t[7];
  _x[8]=(t[7]-t[8]);
  _x[9]=(t[6]-t[9]);
  _x[10]=(t[5]-t[10]);
  _x[11]=(t[4]-t[11]);
  _x[12]=(t[3]-t[12]);
  _x[13]=(t[2]-t[13]);
  _x[14]=(t[1]-t[14]);
  _x[15]=(t[0]-t[15]);
}

const ne_filter_func_double NE_PRE_FILTER_DOUBLE[OD_NBSIZES]={
  ne_pre_filter4_double,
  ne_pre_filter8_double,
  ne_pre_filter16_double
};

const ne_filter_func_double NE_POST_FILTER_DOUBLE[OD_NBSIZES]={
  ne_post_filter4_double,
  ne_post_filter8_double,
  ne_post_filter16_double
};
