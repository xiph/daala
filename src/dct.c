/*Daala video codec
Copyright (c) 2002-2010 Daala project contributors.  All rights reserved.

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

#include "dct.h"

void od_bin_fdct4(od_coeff _y[],const od_coeff _x[]){
  /*8 adds, 5 shifts.*/
  int t[4];
  t[0]=_x[0]+_x[3];
  t[2]=_x[1]+_x[2];
  t[3]=(t[2]+1>>1)-_x[2];
  t[1]=(t[0]>>1)-_x[3];
  /*Embedded 2-point type-II DCT.*/
  t[0]+=t[2];
  _y[0]=(od_coeff)t[0];
  /*This is scaled by 1/2.*/
  _y[2]=(od_coeff)((t[0]+1>>1)-t[2]);
  /*Embedded 2-point type-IV DCT.*/
  /*13/32 ~= \tan(\frac{\pi}{8}) = \sqrt{2}-1 ~=
    0.41421356237309504880168872420970*/
  t[3]=(t[1]*13+16>>5)-t[3];
  _y[3]=(od_coeff)t[3];
  /*11/32 ~= \sin(\frac{\pi}{8})\cos(\frac{\pi}{8}) = \sqrt(\frac{1}{8}) ~=
    0.3535533905927376220042218105242 */
  _y[1]=(od_coeff)(t[1]-(t[3]*11+16>>5));
}

void od_bin_fdct8(od_coeff _y[],const od_coeff _x[]){
  /*29 adds, 15 shifts.*/
  int t[8];
  /*+1/-1 Butterflies and initial permutation:*/
  t[1]=_x[0]-_x[7];
  t[5]=_x[1]-_x[6];
  t[3]=_x[2]-_x[5];
  t[2]=_x[3]+_x[4];
  /*These are scaled by 1/2.*/
  t[7]=_x[3]-(t[2]+1>>1);
  t[6]=_x[2]-(t[3]+1>>1);
  /*Well, except this one.*/
  t[4]=_x[1]+_x[6];
  t[0]=_x[0]-(t[1]+1>>1);
  /*Embedded 4-point type-II DCT:*/
  /*These are scaled by 1/2.*/
  t[0]+=t[2]+1>>1;
  t[6]=(t[4]+1>>1)-t[6];
  /*+1/-1 Butterflies (no additional permutation).*/
  t[4]-=t[6];
  t[2]=t[0]-t[2];
  /*Embedded 2-point type-II DCT.*/
  t[4]=t[0]-t[4];
  _y[4]=(od_coeff)t[4];
  /*This is scaled by 1/2.*/
  _y[0]=(od_coeff)(t[0]-(t[4]+1>>1));
  /*Embedded 2-point type-IV DCT.*/
  /*13/32 ~= \tan(\frac{\pi}{8}) = \sqrt{2}-1 ~=
    0.41421356237309504880168872420970*/
  t[6]=(t[2]*13+16>>5)-t[6];
  _y[6]=(od_coeff)t[6];
  /*11/32 ~= \sin(\frac{\pi}{8})\cos(\frac{\pi}{8}) = \sqrt(\frac{1}{8}) ~=
    0.3535533905927376220042218105242 */
  _y[2]=(od_coeff)(t[2]-(t[6]*11+16>>5));
  /*Embedded 4-point type-IV DST:*/
  /*39/64 ~= 2\frac{1-\cos(\frac{3\pi}{16})}{\sin(\frac{3\pi}{16})} =
    2\frac{2-\sqrt{2+2\sqrt{1+\sqrt{\frac{1}{2}}}-\sqrt{2+\sqrt{2}}}}
     {\sqrt{2-\sqrt{2+\sqrt{2}}}+\sqrt{(2+\sqrt{2})(2-\sqrt{2+\sqrt{2}})}} ~=
    2*0.30334668360734239167588394694130*/
  t[1]-=t[7]*39+32>>6;
  /*9/32 ~= \frac{1}{2}\sin(\frac{3\pi}{16}) =
    \frac{1}{4}(\sqrt{2-\sqrt{2+\sqrt{2}}}+
     \sqrt{(2+\sqrt{2})(2-\sqrt{2+\sqrt{2}})}) ~=
    0.5*0.55557023301970222474283081394853*/
  t[7]+=t[1]*9+16>>5;
  t[1]-=t[7]*39+32>>6;
  /*3/32 ~= \frac{1-\cos(\frac{\pi}{16})}{\sin(\frac{\pi}{16})} =
    \frac{2-\sqrt{2+\sqrt{2+sqrt{2}}}}{\sqrt{2-\sqrt{2+\sqrt{2}}}} ~=
    0.098491403357164253077197521291327*/
  t[5]-=t[3]*3+16>>5;
  /*3/16 ~= \sin(\frac{\pi}{16}) =
    \frac{1}{2}\sqrt{2-\sqrt{2+\sqrt{2}}} ~=
    0.19509032201612826784727476747702*/
  t[3]+=t[5]*3+8>>4;
  t[5]-=t[3]*3+16>>5;
  t[7]+=t[5]>>1;
  _y[5]=(od_coeff)(t[7]-t[5]);
  t[1]+=t[3];
  _y[3]=(od_coeff)((t[1]+1>>1)-t[3]);
  t[7]=(t[1]>>1)-t[7];
  _y[1]=(od_coeff)(t[1]-t[7]);
  _y[7]=(od_coeff)t[7];
}

void od_bin_fdct16(od_coeff _y[],const od_coeff _x[]){
  int t[16];
  int c00;
  int c04;
  int c08;
  int c12;
  t[11]=_x[0]-_x[15];
  t[5]=_x[1]-_x[14];
  t[7]=_x[2]-_x[13];
  t[1]=_x[3]-_x[12];
  t[15]=_x[4]-_x[11];
  t[9]=_x[5]-_x[10];
  t[3]=_x[6]-_x[9];
  t[13]=_x[7]-_x[8];
  t[2]=_x[7]-(t[13]+1>>1);
  t[10]=_x[6]-(t[3]+1>>1);
  t[6]=_x[5]-(t[9]+1>>1);
  t[14]=_x[4]-(t[15]+1>>1);
  t[4]=_x[3]-(t[1]+1>>1);
  t[12]=_x[2]-(t[7]+1>>1);
  t[8]=_x[1]-(t[5]+1>>1);
  t[0]=_x[0]-(t[11]+1>>1);
  /*Embedded 8-point type-II DCT:*/
  /*+1/-1 Butterflies (no additional permutation).*/
  c00=t[0];
  t[0]+=t[2];
  c08=t[8];
  t[8]+=t[10];
  c12=t[12];
  t[12]+=t[6];
  c04=t[4];
  t[4]+=t[14];
  t[14]=c04-t[14];
  t[10]=c08-t[10];
  t[6]=c12-t[6];
  t[2]=c00-t[2];
  /*Embedded 4-point type-II DCT:*/
  /*+1/-1 Butterflies (no additional permutation).*/
  c00=t[0];
  t[0]+=t[4];
  t[12]=t[8]-t[12];
  t[8]-=t[12]+1>>1;
  t[4]=c00-t[4];
  /*Embedded 2-point type-II DCT.*/
  t[8]=(t[0]+1>>1)-t[8];
  _y[0]=(od_coeff)(t[0]-t[8]);
  _y[8]=(od_coeff)t[8];
  /*Embedded 2-point type-IV DCT.*/
  /*13/32 ~= \tan(\frac{\pi}{8}) = \sqrt{2}-1 ~=
    0.41421356237309504880168872420970*/
  t[12]=(t[4]*13+16>>5)-t[12];
  _y[12]=(od_coeff)t[12];
  /*11/32 ~= \sin(\frac{\pi}{8})\cos(\frac{\pi}{8}) = \sqrt(\frac{1}{8}) ~=
    0.3535533905927376220042218105242 */
  _y[4]=(od_coeff)(t[4]-(t[12]*11+16>>5));
  /*Embedded 4-point type-IV DCT:*/
  /*19/64 ~= \frac{1-\cos(\frac{3\pi}{16})}{\sin(\frac{3\pi}{16})} =
    \frac{2-\sqrt{2+2\sqrt{1+\sqrt{\frac{1}{2}}}-\sqrt{2+\sqrt{2}}}}
     {\sqrt{2-\sqrt{2+\sqrt{2}}}+\sqrt{(2+\sqrt{2})(2-\sqrt{2+\sqrt{2}})}} ~=
    0.30334668360734239167588394694130*/
  t[2]-=t[14]*19+32>>6;
  /*9/16 ~= \sin(\frac{3\pi}{16}) =
    \frac{1}{2}(\sqrt{2-\sqrt{2+\sqrt{2}}}+
     \sqrt{(2+\sqrt{2})(2-\sqrt{2+\sqrt{2}})}) ~=
    0.55557023301970222474283081394853*/
  t[14]+=t[2]*9+8>>4;
  t[2]-=t[14]*19+32>>6;
  /*3/32 ~= \frac{1-\cos(\frac{\pi}{16})}{\sin(\frac{\pi}{16})} =
    \frac{2-\sqrt{2+\sqrt{2+sqrt{2}}}}{\sqrt{2-\sqrt{2+\sqrt{2}}}} ~=
    0.098491403357164253077197521291327*/
  t[10]-=t[6]*3+16>>5;
  /*3/16 ~= \sin(\frac{\pi}{16}) =
    \frac{1}{2}\sqrt{2-\sqrt{2+\sqrt{2}}} ~=
    0.19509032201612826784727476747702*/
  t[6]+=t[10]*3+8>>4;
  t[10]-=t[6]*3+16>>5;
  t[10]=t[14]-t[10];
  _y[10]=(od_coeff)t[10];
  t[14]-=t[10]>>1;
  t[6]=t[2]-t[6];
  _y[6]=(od_coeff)t[6];
  t[2]-=t[6]>>1;
  _y[2]=(od_coeff)(t[2]+t[14]);
  _y[14]=(od_coeff)(t[2]-t[14]);
  /*Embedded 8-point type-IV DCT:*/
  /*23/64 ~= \frac{1-\cos(\frac{7\pi}{32})}{\sin(\frac{7\pi}{32})} =
    \frac{2+(1+\sqrt{2+\sqrt{2}}(1+\sqrt{2+\sqrt{2+\sqrt{2}}}))
     \sqrt{2+\sqrt{2+\sqrt{2+\sqrt{2}}}}}
    {(1+\sqrt{2+\sqrt{2}}(1+\sqrt{2+\sqrt{2+\sqrt{2}}}))
     \sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}}} ~=
    0.35780572131452410467248774377447*/
  t[11]-=t[13]*23+32>>6;
  /*41/64 ~= \sin(\frac{7\pi}{32}) =
    \frac{1}{2}(1+\sqrt{2+\sqrt{2}}(1+\sqrt{2+\sqrt{2+\sqrt{2}}}))
    \sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}} ~=
    0.63439328416364549821517161322549*/
  t[13]+=t[11]*41+32>>6;
  t[11]-=t[13]*23+32>>6;
  /*19/32 ~= \frac{1-\cos(\frac{11\pi}{32})}{\sin(\frac{11\pi}{32})} =
    \frac{2+(1+\sqrt{2}+\sqrt{2+\sqrt{2}}-
      (1+\sqrt{2})\sqrt{2+\sqrt{2+\sqrt{2}}})
     \sqrt{2+\sqrt{2+\sqrt{2+\sqrt{2}}}}}
    {(1+\sqrt{2}+\sqrt{2+\sqrt{2}}+(1+\sqrt{2})\sqrt{2+\sqrt{2+\sqrt{2}}})
     sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}}} ~=
    0.59937693368192376627138986901440*/
  t[3]-=t[5]*19+16>>5;
  /*7/8 ~= \sin(\frac{11\pi}{32}) =
    \frac{1}{2}(1+\sqrt{2}+\sqrt{2+\sqrt{2}}+
     (1+\sqrt{2})\sqrt{2+\sqrt{2+\sqrt{2}}})
    \sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}} ~=
    0.88192126434835502971275686366039*/
  t[5]+=t[3]*7+4>>3;
  t[3]-=t[5]*19+16>>5;
  /*9/64 ~= \frac{1-\cos(\frac{3\pi}{32})}{\sin(\frac{3\pi}{32}} =
    \frac{2+(1-\sqrt{2+\sqrt{2+\sqrt{2}}})\sqrt{2+\sqrt{2+\sqrt{2+\sqrt{2}}}}}
    {(1+\sqrt{2+\sqrt{2+\sqrt{2}}})\sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}}} ~=
    0.14833598753834742875367651148691*/
  t[7]-=t[9]*9+32>>6;
  /*9/32 ~= \sin(\frac{3\pi}{32}) =
    \frac{1}{2}(1+\sqrt{2+\sqrt{2+\sqrt{2}}})
    \sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}} ~=
    0.29028467725446236763619237581740*/
  t[9]+=t[7]*19+32>>6;
  t[7]-=t[9]*9+32>>6;
  /*29/32 ~= \frac{1-\cos(\frac{15\pi}{32})}{\sin(\frac{15\pi}{32})} =
    \frac{2-\sqrt{2-\sqrt{2+\sqrt{2+\sqrt{2}}}}}
    {\sqrt{2+\sqrt{2+\sqrt{2+\sqrt{2}}}}} ~=
    0.90634716901914715794614271726891*/
  t[15]-=t[1]*29+16>>5;
  /*1 ~= \sin(\frac{15\pi}{11}) =
    \frac{1}{2}\sqrt{2+\sqrt{2+\sqrt{2+\sqrt{2}}}} ~=
    0.99518472667219688624483695310948*/
  t[1]+=t[15];
  t[15]-=t[1]*29+16>>5;
  t[15]=t[13]-t[15];
  t[13]-=t[15]>>1;
  t[3]+=t[9];
  t[9]=(t[3]>>1)-t[9];
  t[1]+=t[11];
  t[11]=(t[1]>>1)-t[11];
  t[5]=t[7]-t[5];
  t[7]-=t[5]+1>>1;
  t[13]+=t[5]>>1;
  t[5]=t[13]-t[5];
  t[11]=(t[3]>>1)-t[11];
  t[3]-=t[11];
  t[13]=t[13]+(t[3]*13+16>>5);
  _y[13]=(od_coeff)t[13];
  _y[3]=(od_coeff)((t[13]*11+16>>5)-t[3]);
  t[9]+=t[15]+1>>1;
  _y[9]=(od_coeff)t[9];
  t[15]=t[9]-t[15];
  t[7]=(t[1]>>1)-t[7];
  t[1]-=t[7];
  t[15]+=t[1];
  _y[15]=(od_coeff)t[15];
  _y[1]=(od_coeff)(t[1]-(t[15]+1>>1));
  _y[7]=(od_coeff)t[7];
  t[5]=t[5]+(t[11]*13+16>>5);
  _y[5]=(od_coeff)t[5];
  _y[11]=(od_coeff)(t[11]-(t[5]*11+16>>5));
}

void od_bin_idct4(od_coeff _x[],const od_coeff _y[]){
  int t[4];
  t[1]=(_y[3]*11+16>>5)+_y[1];
  t[3]=(t[1]*13+16>>5)-_y[3];
  t[2]=(_y[0]+1>>1)-_y[2];
  t[0]=_y[0]-t[2];
  _x[3]=(t[0]>>1)-t[1];
  _x[2]=(t[2]+1>>1)-t[3];
  _x[1]=t[2]-_x[2];
  _x[0]=t[0]-_x[3];
}

void od_bin_idct8(od_coeff _x[],const od_coeff _y[]){
  int t[8];
  t[1]=_y[1]+_y[7];
  t[7]=(t[1]>>1)-_y[7];
  t[3]=(t[1]+1>>1)-_y[3];
  t[1]-=t[3];
  t[5]=t[7]-_y[5];
  t[7]-=t[5]>>1;
  t[5]+=t[3]*3+16>>5;
  t[3]-=t[5]*3+8>>4;
  t[5]+=t[3]*3+16>>5;
  t[1]+=t[7]*39+32>>6;
  t[7]-=t[1]*9+16>>5;
  t[1]+=t[7]*39+32>>6;
  t[2]=_y[2]+(_y[6]*11+16>>5);
  t[6]=(t[2]*13+16>>5)-_y[6];
  t[0]=_y[0]+(_y[4]+1>>1);
  t[4]=t[0]-_y[4];
  t[2]=t[0]-t[2];
  t[4]+=t[6];
  t[6]=(t[4]+1>>1)-t[6];
  t[0]-=t[2]+1>>1;
  t[0]+=t[1]+1>>1;
  _x[0]=(od_coeff)t[0];
  t[4]=t[4]+t[5]>>1;
  _x[1]=(od_coeff)t[4];
  t[6]+=t[3]+1>>1;
  _x[2]=(od_coeff)t[6];
  t[7]+=t[2]+1>>1;
  _x[3]=(od_coeff)t[7];
  _x[4]=t[2]-t[7];
  _x[5]=-t[3]+t[6];
  _x[6]=-t[5]+t[4];
  _x[7]=-t[1]+t[0];
}

void od_bin_idct16(od_coeff _x[],const od_coeff _y[]){
  int t[16];
  t[11]=_y[11]+(_y[5]*11+16>>5);
  t[5]=_y[5]-(t[11]*13+16>>5);
  t[1]=_y[1]+(_y[15]+1>>1);
  t[15]=_y[15]-t[1];
  t[1]+=_y[7];
  t[7]=(t[1]>>1)-_y[7];
  t[15]=_y[9]-t[15];
  t[9]=_y[9]-(t[15]+1>>1);
  t[3]=(_y[13]*11+16>>5)-_y[3];
  t[13]=_y[13]-(t[3]*13+16>>5);
  t[3]+=t[11];
  t[11]=(t[3]>>1)-t[11];
  t[5]=t[13]-t[5];
  t[13]-=t[5]>>1;
  t[7]+=t[5]+1>>1;
  t[5]=t[7]-t[5];
  t[11]=(t[1]>>1)-t[11];
  t[1]-=t[11];
  t[9]=(t[3]>>1)-t[9];
  t[3]-=t[9];
  t[13]+=t[15]>>1;
  t[15]=t[13]-t[15];
  t[15]+=t[1]*29+16>>5;
  t[1]-=t[15];
  t[15]+=t[1]*29+16>>5;
  t[7]+=t[9]*9+32>>6;
  t[9]-=t[7]*19+32>>6;
  t[7]+=t[9]*9+32>>6;
  t[3]+=t[5]*19+16>>5;
  t[5]-=t[3]*7+4>>3;
  t[3]+=t[5]*19+16>>5;
  t[11]+=t[13]*23+32>>6;
  t[13]-=t[11]*41+32>>6;
  t[11]+=t[13]*23+32>>6;
  t[2]=_y[2]+_y[14]>>1;
  t[14]=t[2]-_y[14];
  t[2]+=_y[6]>>1;
  t[6]=t[2]-_y[6];
  t[14]+=_y[10]>>1;
  t[10]=t[14]-_y[10];
  t[10]+=t[6]*3+16>>5;
  t[6]-=t[10]*3+8>>4;
  t[10]+=t[6]*3+16>>5;
  t[2]+=t[14]*19+32>>6;
  t[14]-=t[2]*9+8>>4;
  t[2]+=t[14]*19+32>>6;
  t[4]=_y[4]+(_y[12]*11+16>>5);
  t[12]=(t[4]*13+16>>5)-_y[12];
  t[0]=_y[0]+_y[8];
  t[8]=(t[0]+1>>1)-_y[8];
  t[0]=t[0]+t[4]>>1;
  t[8]+=t[12]+1>>1;
  t[12]=t[8]-t[12];
  t[4]=t[0]-t[4];
  t[0]=t[0]+t[2]>>1;
  t[8]=t[8]+t[10]>>1;
  t[12]=t[12]+t[6]>>1;
  t[4]=t[4]+t[14]>>1;
  t[14]=t[4]-t[14];
  t[6]=t[12]-t[6];
  t[10]=t[8]-t[10];
  t[2]=t[0]-t[2];
  t[0]=t[0]+(t[11]+1>>1);
  _x[0]=(od_coeff)t[0];
  t[8]=t[8]+(t[5]+1>>1);
  _x[1]=(od_coeff)t[8];
  t[12]=t[12]+(t[7]+1>>1);
  _x[2]=(od_coeff)t[12];
  t[4]=t[4]+(t[1]+1>>1);
  _x[3]=(od_coeff)t[4];
  t[14]=t[14]+(t[15]+1>>1);
  _x[4]=(od_coeff)t[14];
  t[6]=t[6]+(t[9]+1>>1);
  _x[5]=(od_coeff)t[6];
  t[10]=t[10]+(t[3]+1>>1);
  _x[6]=(od_coeff)t[10];
  t[2]=t[2]+(t[13]+1>>1);
  _x[7]=(od_coeff)t[2];
  _x[8]=t[2]-t[13];
  _x[9]=t[10]-t[3];
  _x[10]=t[6]-t[9];
  _x[11]=t[14]-t[15];
  _x[12]=t[4]-t[1];
  _x[13]=t[12]-t[7];
  _x[14]=t[8]-t[5];
  _x[15]=t[0]-t[11];
}

#if 0
/*Test code.*/
#include <stdio.h>
#include <math.h>
#include <string.h>

/*The (2-D) scaling factors that make a true DCT approximation out of the
   integer transform.*/
static const double DCT8x8_FSCALE[8][8]={
  {
    2,1,M_SQRT2/SIN_3PI_8,M_SQRT2,
    1,M_SQRT2,M_SQRT2*SIN_3PI_8,1
  },
  {
    1,0.5,M_SQRT1_2/SIN_3PI_8,M_SQRT1_2,
    0.5,M_SQRT1_2,M_SQRT1_2*SIN_3PI_8,0.5
  },
  {
    M_SQRT2/SIN_3PI_8,M_SQRT1_2/SIN_3PI_8,
    1/(SIN_3PI_8*SIN_3PI_8),1/SIN_3PI_8,
    M_SQRT1_2/SIN_3PI_8,1/SIN_3PI_8,
    1,M_SQRT1_2/SIN_3PI_8
  },
  {
    M_SQRT2,M_SQRT1_2,1/SIN_3PI_8,1,M_SQRT1_2,1,SIN_3PI_8,M_SQRT1_2
  },
  {
    1,0.5,M_SQRT1_2/SIN_3PI_8,M_SQRT1_2,
    0.5,M_SQRT1_2,M_SQRT1_2*SIN_3PI_8,0.5
  },
  {
    M_SQRT2,M_SQRT1_2,1/SIN_3PI_8,1,M_SQRT1_2,1,SIN_3PI_8,M_SQRT1_2
  },
  {
    M_SQRT2*SIN_3PI_8,M_SQRT1_2*SIN_3PI_8,1,SIN_3PI_8,
    M_SQRT1_2*SIN_3PI_8,SIN_3PI_8,SIN_3PI_8*SIN_3PI_8,M_SQRT1_2*SIN_3PI_8
  },
  {
    1,0.5,M_SQRT1_2/SIN_3PI_8,M_SQRT1_2,
    0.5,M_SQRT1_2,M_SQRT1_2*SIN_3PI_8,0.5
  }
};

/*The (2-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT8x8_ISCALE[8][8]={
  {
    0.5,1,M_SQRT1_2*SIN_3PI_8,M_SQRT1_2,1,M_SQRT1_2,M_SQRT1_2/SIN_3PI_8,1
  },
  {
    1,2,M_SQRT2*SIN_3PI_8,M_SQRT2,2,M_SQRT2,M_SQRT2/SIN_3PI_8,2
  },
  {
    M_SQRT1_2*SIN_3PI_8,M_SQRT2*SIN_3PI_8,SIN_3PI_8*SIN_3PI_8,SIN_3PI_8,
    M_SQRT1_2*SIN_3PI_8,SIN_3PI_8,1,M_SQRT2*SIN_3PI_8
  },
  {
    M_SQRT1_2,M_SQRT2,SIN_3PI_8,1,M_SQRT2,1,1/SIN_3PI_8,M_SQRT2
  },
  {
    1,2,M_SQRT2*SIN_3PI_8,M_SQRT2,2,M_SQRT2,M_SQRT2/SIN_3PI_8,2
  },
  {
    M_SQRT1_2,M_SQRT2,SIN_3PI_8,1,M_SQRT2,1,1/SIN_3PI_8,M_SQRT2
  },
  {
    M_SQRT1_2/SIN_3PI_8,M_SQRT2/SIN_3PI_8,1,1/SIN_3PI_8,
    M_SQRT1_2/SIN_3PI_8,1/SIN_3PI_8,1/(SIN_3PI_8*SIN_3PI_8),M_SQRT2/SIN_3PI_8
  },
  {
    1,2,M_SQRT2*SIN_3PI_8,M_SQRT2,2,M_SQRT2,M_SQRT2/SIN_3PI_8,2
  }
};



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



#if 0
/*A floating-point implementation of an inverse 8-point DCT.
  The scaling factors of the input are the same as in the integerized version.
  x: The destination vector (of size 8).
     This may be the same as the source.
  y: The source vector (of size 8).*/
void idct8(double _x[],const double _y[]){
  double t[8];
  double a0;
  double a1;
  double c01;
  double c05;
  t[7]=0.5*(_y[1]-_y[7]);
  t[1]=0.5*(_y[1]+_y[7]);
  t[3]=0.5*(t[1]-_y[3]);
  t[1]-=t[3];
  t[5]=0.5*(t[7]-_y[5]);
  t[7]-=t[5];
  /*These rotations are done with Loeffler's original factorization, not the
     lifting method, as it has less accumulated error.*/
  a0=(t[3]+t[5])*0.98078528040323044912618223613424/*cos(pi/16)*/;
  c05=t[5];
  t[5]=t[3]*-0.78569495838710218127789736765722/*sin(pi/16)-cos(pi/16)*/+a0;
  t[3]=c05*-1.1758756024193587169744671046113/*-sin(pi/16)-cos(pi/16)*/+a0;
  a1=(t[1]+t[7])*0.83146961230254523707878837761791/*cos(3pi/16)*/;
  c01=t[1];
  t[1]=t[7]*-0.27589937928294301233595756366937/*sin(3pi/16)-cos(3pi/16)*/+a1;
  t[7]=c01*-1.3870398453221474618216191915664/*-sin(3pi/16)-cos(3pi/16)*/+a1;
  /*This rotation is still done using the scaled lifting method used in the
     integer transform (with more accurate multipliers, of course), to keep
     the scaling factors the same.*/
  t[2]=_y[2]+_y[6]*0.3535533905927376220042218105242/*sin(pi/8)cos(pi/8)*/;
  t[6]=t[2]*0.41421356237309504880168872420970/*tan(pi/8)*/-_y[6];
  t[4]=0.5*_y[0]-_y[4];
  t[0]=_y[0]-t[4];
  t[6]=0.5*t[4]-t[6];
  t[2]=0.5*t[0]-t[2];
  t[4]-=t[6];
  t[0]-=t[2];
  _x[0]=t[1]+0.5*t[0];
  _x[1]=t[5]+0.5*t[4];
  _x[2]=t[3]+0.5*t[6];
  _x[3]=t[7]+0.5*t[2];
  _x[4]=t[2]-_x[3];
  _x[5]=t[6]-_x[2];
  _x[6]=t[4]-_x[1];
  _x[7]=t[0]-_x[0];
}
#else

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
#endif

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
    floatcoefs[i][j]=refcoefs[i][j]*DCT8x8_FSCALE[i][j];
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
    for(j=0;j<8;j++)_basis[j][i]=x[j]*DCT8_FSCALE[i];
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

static void check8(void){
  od_coeff min[8];
  od_coeff max[8];
  double   basis[8][8];
  double   tbasis[8][8];
  int      i;
  int      j;
  dynamic_range8();
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
}



static void print_basis16(double _basis[16][16]){
  int i;
  int j;
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      printf("%8.5lf%c",_basis[i][j],j==16-1?'\n':' ');
    }
  }
}

static void compute_fbasis16(double _basis[16][16]){
  int i;
  int j;
  for(i=0;i<16;i++){
    od_coeff x[16];
    for(j=0;j<16;j++)x[j]=(i==j)*256;
    od_bin_fdct16(x,x);
    for(j=0;j<16;j++)_basis[j][i]=x[j]/256.0;
  }
}

static void compute_ftrue_basis16(double _basis[16][16]){
  int i;
  int j;
  for(j=0;j<16;j++){
    for(i=0;i<16;i++){
      _basis[j][i]=sqrt(2.0/16)*cos((i+0.5)*j*M_PI/16)*DCT16_ISCALE[j];
      if(j==0)_basis[j][i]*=M_SQRT1_2;
    }
  }
}

static double compute_mse16(double _basis[16][16],double _tbasis[16][16]){
  double e[16][16];
  double ret;
  int    i;
  int    j;
  int    k;
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      e[i][j]=0;
      for(k=0;k<16;k++){
        e[i][j]+=(_basis[i][k]-_tbasis[i][k])*AUTOCORR[k-j+15];
      }
    }
  }
  ret=0;
  for(i=0;i<16;i++){
    for(j=0;j<16;j++){
      ret+=e[i][j]*(_basis[i][j]-_tbasis[i][j]);
    }
  }
  return ret/16;
}

/*
static void check16(void){
  od_coeff min[16];
  od_coeff max[16];
  double   basis[16][16];
  double   tbasis[16][16];
  double   lbasis[16][16];
  int      i;
  int      j;
  for(j=0;j<16;j++)min[j]=max[j]=0;
  for(i=0;i<1<<16;i++){
    od_coeff x[16];
    od_coeff y[16];
    od_coeff x2[16];
    for(j=0;j<16;j++)x[j]=(i>>j&1)?255:-256;
    od_bin_fdct16(y,x);
    od_bin_idct16(x2,y);
    for(j=0;j<16;j++){
      if(y[j]<min[j])min[j]=y[j];
      else if(y[j]>max[j])max[j]=y[j];
    }
    for(j=0;j<16;j++)if(x[j]!=x2[j]){
      printf("Mismatch:\n");
      printf("in:    ");
      for(j=0;j<16;j++)printf(" %i",x[j]);
      printf("\nxform: ");
      for(j=0;j<16;j++)printf(" %i",y[j]);
      printf("\nout:   ");
      for(j=0;j<16;j++)printf(" %i",x2[j]);
      printf("\n\n");
      break;
    }
  }
  printf("Min:");
  for(j=0;j<16;j++)printf(" %5i",min[j]);
  printf("\nMax:");
  for(j=0;j<16;j++)printf(" %5i",max[j]);
  printf("\nod_bin_fdct16 basis:\n");
  compute_fbasis16(basis);
  print_basis16(basis);
  printf("Scaled type-II DCT basis:\n");
  compute_ftrue_basis16(tbasis);
  print_basis16(tbasis);
  printf("MSE: %.32lg\n\n",compute_mse6(basis,tbasis));
}
*/

static void bin_fxform_2d16(od_coeff _x[16*2][16*2]){
   od_coeff y[16*2];
   int      u;
   int      v;
   /*Perform pre-filtering.*/
   for(u=0;u<16*2;u++){
      od_pre_filter16(_x[u],_x[u]);
      od_pre_filter16(_x[u]+16,_x[u]+16);
   }
   for(v=0;v<16*2;v++){
     for(u=0;u<16*2;u++)y[u]=_x[u][v];
     od_pre_filter16(y,y);
     od_pre_filter16(y+16,y+16);
     for(u=0;u<16*2;u++)_x[u][v]=y[u];
   }
   /*Perform DCT.*/
   for(u=16/2;u<16*3/2;u++)od_bin_fdct16(_x[u]+16/2,_x[u]+16/2);
   for(v=16/2;v<16*3/2;v++){
     for(u=16/2;u<16*3/2;u++)y[u]=_x[u][v];
     od_bin_fdct16(y+16/2,y+16/2);
     for(u=16/2;u<16*3/2;u++)_x[u][v]=y[u];
   }
}



static void check16(void){
  od_coeff min2[16][16];
  od_coeff max2[16][16];
  double   basis2[16][16][16*2][16*2];
  int      i;
  int      j;
  int      u;
  int      v;
  for(i=0;i<16*2;i++){
     for(j=0;j<16*2;j++){
        od_coeff x[16*2][16*2];
        /*Generate impulse.*/
        for(u=0;u<16*2;u++){
           for(v=0;v<16*2;v++){
              x[u][v]=(u==i&&v==j)*256;
           }
        }
        bin_fxform_2d16(x);
        /*Retrieve basis elements.*/
        for(u=0;u<16;u++){
           for(v=0;v<16;v++){
              basis2[u][v][i][j]=x[u+16/2][v+16/2]/256.0;
           }
        }
     }
  }
  for(u=0;u<16;u++){
     for(v=0;v<16;v++){
        od_coeff x[16*2][16*2];
        for(i=0;i<16*2;i++){
           for(j=0;j<16*2;j++){
              x[i][j]=basis2[u][v][i][j]<0?-255:255;
           }
        }
        bin_fxform_2d16(x);
        max2[u][v]=x[u+16/2][v+16/2];
        for(i=0;i<16*2;i++){
           for(j=0;j<16*2;j++){
              x[i][j]=basis2[u][v][i][j]>0?-255:255;
           }
        }
        bin_fxform_2d16(x);
        min2[u][v]=x[u+16/2][v+16/2];
     }
  }
  printf("2-D ranges:\n");
  for(u=0;u<16;u++){
     printf("Min %2i:",u);
     for(v=0;v<16;v++)printf(" %6i",min2[u][v]);
     printf("\nMax %2i:",u);
     for(v=0;v<16;v++)printf(" %6i",max2[u][v]);
     printf("\n");
  }
  {
    od_coeff min[16];
    od_coeff max[16];
    double   basis[16][16];
    double   tbasis[16][16];
    for(i=0;i<1<<16;i++){
      od_coeff x[16];
      od_coeff y[16];
      od_coeff x2[16];
      for(j=0;j<16;j++)x[j]=(i>>j&1)?255:-256;
      od_bin_fdct16(y,x);
      od_bin_idct16(x2,y);
      for(j=0;j<16;j++){
        if(y[j]<min[j])min[j]=y[j];
        else if(y[j]>max[j])max[j]=y[j];
      }
      for(j=0;j<16;j++)if(x[j]!=x2[j]){
        printf("Mismatch:\n");
        printf("in:    ");
        for(j=0;j<16;j++)printf(" %i",x[j]);
        printf("\nxform: ");
        for(j=0;j<16;j++)printf(" %i",y[j]);
        printf("\nout:   ");
        for(j=0;j<16;j++)printf(" %i",x2[j]);
        printf("\n\n");
        break;
      }
    }
    printf("Min:");
    for(j=0;j<16;j++)printf(" %5i",min[j]);
    printf("\nMax:");
    for(j=0;j<16;j++)printf(" %5i",max[j]);
    printf("\nod_bin_fdct16 basis:\n");
    compute_fbasis16(basis);
    print_basis16(basis);
    printf("Scaled type-II DCT basis:\n");
    compute_ftrue_basis16(tbasis);
    print_basis16(tbasis);
    printf("MSE: %.32lg\n\n",compute_mse16(basis,tbasis));
  }
}


int main(void){
  check4();
  check8();
  check16();
  return 0;
}
#endif
