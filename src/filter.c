/*
    Daala video codec
    Copyright (C) 2010 Timothy B. Terriberry and Daala project contributors

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


#include "filter.h"



/*Pre-/post-filter pairs of various sizes.
  For a FIR, PR, LP filter bank, the pre-filter must have the structure:

   0.5*{{I,J},{J,-I}}*{{U,0},{0,V}},{{I,J},{J,-I}}

  To create a pre-/post-filter for an M-sample block transform, I, J, U and V
   are blocks of size floor(M/2).
  I is the identity, J is the reversal matrix, and U and V are arbitrary
   invertible matrices.
  If U and V are orthonormal, then the pre-/post-filters are orthonormal.
  Otherwise they are merely biorthogonal.
  For filters with larger inputs, multiple pre-/post-filter stages can be used.
  However, one level should be sufficient to eliminate blocking artifacts, and
   additional levels would expand the influence of ringing artifacts.

  All the filters here are an even narrower example of the above structure.
  U is taken to be the identity, as it does not help much with energy
   compaction.
  V is chosen so that the filter pair is (1,2) regular, assuming a block filter
   with 0 DC leakage, such as the DCT.
  This is particularly important for the way that motion-compensated prediction
   is done.
  A 2-regular synthesis filter can reproduce a linear ramp from just the DC
   coefficient, which matches the linearly interpolated offsets that were
   subtracted out of the motion-compensation phase.

  In order to yield a fast implementation, the V filter is chosen to be of the
   form:
    x0 -s0----------...----+---- y0
            p0 |           | u0
    x1 -s1-----+----...--+------ y1
              p1 |       | u1
    x2 -s2-------+--...+-------- y2
                p2 |   | u2
    x3 -s3---------+...--------- y3
                     .
                     .
                     .
  Here, s(i) is a scaling factor, such that s(i) >= 1, to preserve perfect
   reconstruction given an integer implementation.
  p(i) and u(i) are arbitrary, so long as p(i)u(i) != -1, to keep the transform
   invertible.
  In order to make V (1,2) regular, we have the conditions:
    s0+M*u0==M
    (2*i+1)*s(i)+M*u(i)+M*(1-u(i-1))*p(i-1)==M, i in [1..K-2]
    (M-1)*s(K-1)+M*(1-u(K-2))*p(K-2)==M
   where K=floor(M/2).
  These impose additonal constraints on u(i) and p(i), derived from the
   constraints above, such as u(0) <= (M-1)/M.
  It is desirable to have u(i), p(i) and 1/s(i) be dyadic rationals, the latter
   to provide for a fast inverse transform.
  However, as can be seen by the constraints, it is very easy to make u(i) and
   p(i) be dyadic, or 1/s(i) be dyadic, but solutions with all of them dyadic
   are very sparse, and require at least s(0) to be a power of 2.
  Such solutions do not have a good coding gain, and so we settle for having
   s(i) be dyadic instead of 1/s(i).

  As an aside, the form
    x0 -s0---------------+---- y0
                    p0 | | u0
    x1 -s1-----------+-+------ y1
                p1 | | u1
    x2 -s2-------+-+---------- y2
            p2 | | u2
    x3 -s3-----+-------------- y3
                 .
                 .
                 .
   yields slightly higher coding gains, but the conditions for (1,2) regularity
    s0+M*u0==M
    (2*i+1)*s(i)+M*u(i)+(2*i-1)*s(i-1)*p(i-1)==M, i in [1..K-2]
    (M-1)*s(K-1)+(M-3)*s(K-2)*p(K-2)==M
   make dyadic rational approximations more sparse, such that the coding gains
   of the approximations are actually lower.

  Our solutions are found by numerically solving for the u(i) and p(i) that
   yield the optimal coding gain, given the above restrictions.
  Unfortunately, finding dyadic rational approximations of the optimal solution
   is non-trivial, since dyadic rational u(i)'s and p(i)'s that yield good,
   dyadic rational s(i)'s are sparse.
  The neighborhood of the optimal point is searched for such a solution by
   brute force.
  This search continues in an expanding area around the optimal point, until
   the boundaries of this region all have a smaller coding gain than the best
   acceptable solution found so far.
  The maximum denominator for all coefficients was allowed to be 64.*/

#define TEST

#if defined(TEST)
extern int minv[32];
extern int maxv[32];
#define _Check(_val,_idx) \
   if((_val)<minv[(_idx)])minv[(_idx)]=(_val); \
   if((_val)>maxv[(_idx)])maxv[(_idx)]=(_val);
#else
#define _Check(_val,_idx)
#endif



void od_pre_filter4(ogg_int16_t _y[4],const ogg_int16_t _x[4]){
   /*Optimal coding gain without dyadic rational restrictions: 8.60382 dB.
     S={1.4688242104187466,1.4228174370096593},
     0={0.6327939473953134,-0.18276680703263349}
     8-bit coding gain: 8.59848 dB
     S={47/32,357/256}, 0={81/128,-1/8}
     7-bit coding gain: same as 6-bit
     6-bit coding gain: 8.59689 dB
     S={23/16,93/64}, 0={41/64,-1/4}
     5-bit coding gain: 8.55956 dB
     S={13/8,47/32}, 0={19/32,-1/4}
     4-bit coding gain: 8.53007 dB
     S={5/4,23/16}, 0={11/16,-1/4}*/
   int t[4];
   /*+1/-1 butterflies (required for FIR, PR, LP).*/
   t[3]=_x[0]-_x[3];
   t[2]=_x[1]-_x[2];
   t[1]=_x[1]-(t[2]>>1);
   t[0]=_x[0]-(t[3]>>1);
   /*U filter (arbitrary invertible, omitted).*/
   /*V filter (arbitrary invertible).*/
   /*Scaling factors: the biorthogonal part.*/
   /*Note: t[i]+=t[i]>>15&1 is equivalent to: if(t[i]>0)t[i]++
     This step ensures that the scaling is trivially invertible on the decoder's
      side, with perfect reconstruction.*/
   _Check(t[2],0);
   t[2]=t[2]*5>>2;      /*s0= 1.47093*/
   t[2]+=-t[2]>>15&1;
   _Check(t[3],1);
   t[3]=t[3]*23>>4;     /*s1= 1.42128*/
   t[3]+=-t[3]>>15&1;
   /*Rotation:*/
   _Check(t[2],2);
   t[3]-=t[2]+2>>2;     /*p0=-0.179364*/
   _Check(t[3],3);
   t[2]+=t[3]*11+8>>4;  /*u0= 0.632267*/
   /*More +1/-1 butterflies (required for FIR, PR, LP).*/
   t[0]+=t[3]>>1;
   _y[0]=(ogg_int16_t)t[0];
   t[1]+=t[2]>>1;
   _y[1]=(ogg_int16_t)t[1];
   _y[2]=(ogg_int16_t)(t[1]-t[2]);
   _y[3]=(ogg_int16_t)(t[0]-t[3]);
}

void od_post_filter4(ogg_int16_t _x[4],const ogg_int16_t _y[4]){
   int t[4];
   t[3]=_y[0]-_y[3];
   t[2]=_y[1]-_y[2];
   t[1]=_y[1]-(t[2]>>1);
   t[0]=_y[0]-(t[3]>>1);
   t[2]-=t[3]*11+8>>4;
   t[3]+=t[2]+2>>2;
   t[3]=(t[3]<<4)/23;
   t[2]=(t[2]<<2)/5;
   t[0]+=t[3]>>1;
   _x[0]=(ogg_int16_t)t[0];
   t[1]+=t[2]>>1;
   _x[1]=(ogg_int16_t)t[1];
   _x[2]=(ogg_int16_t)(t[1]-t[2]);
   _x[3]=(ogg_int16_t)(t[0]-t[3]);
}



void od_pre_filter8(ogg_int16_t _y[8],const ogg_int16_t _x[8]){
   int t[8];
   /*+1/-1 butterflies (required for FIR, PR, LP).*/
   t[7]=_x[0]-_x[7];
   t[6]=_x[1]-_x[6];
   t[5]=_x[2]-_x[5];
   t[4]=_x[3]-_x[4];
   t[3]=_x[3]-(t[4]>>1);
   t[2]=_x[2]-(t[5]>>1);
   t[1]=_x[1]-(t[6]>>1);
   t[0]=_x[0]-(t[7]>>1);
   /*U filter (arbitrary invertible, omitted).*/
   /*V filter (arbitrary invertible, can be optimized for speed or accuracy).*/
#if 0
   /*Non-regularity constrained type-III coding gain: ? dB
      (optimal without dyadic rational restrictions: 9.6115 dB).*/
   t[4]=t[4]*23>>5;     /* 1.40*/
   t[5]=t[5]*7>>3;      /* 1.12*/
   t[6]=t[6]*7>>3;      /* 1.14*/
   t[7]=t[7]*27>>5;     /* 1.19*/
   t[7]-=t[6]*7+32>>6;  /*-0.11*/
   t[6]+=t[7]*11+16>>5; /* 0.34*/
   t[6]-=t[5]*21+32>>6; /*-0.33*/
   t[5]+=t[6]*9+8>>4;   /* 0.56*/
   t[5]-=t[4]*13+16>>5; /*-0.40*/
   t[4]+=t[5]*25+16>>5; /* 0.78*/
#elif 0
   /*(1,2) regular type-III coding gain: 9.55331 dN
      (optimal without dyadic rational restrictions: 9.57059),
     S={3/2,5/4,37/32,319/256}, 0={13/16,-1/2}, 1={5/8,-3/8}, 2={29/64,-1/8}*/
   /*Scaling factors: the biorthogonal part.*/
   /*Note: t[i]+=t[i]>>15&1 is equivalent to: if(t[i]>0)t[i]++
     This step ensures that the scaling is trivially invertible on the decoder's
      side, with perfect reconstruction.*/
   t[4]=t[4]*3>>1;      /* 1.39901*/
   t[4]+=-t[4]>>15&1;
   t[5]=t[5]*5>>2;      /* 1.21533*/
   t[5]+=-t[5]>>15&1;
   t[6]=t[6]*37>>5;     /* 1.21348*/
   t[6]+=-t[6]>>15&1;
   t[7]=t[7]*319>>8;    /* 1.23047*/
   t[7]+=-t[7]>>15&1;
   /*Rotations:*/
   t[7]-=t[6]+4>>3;     /*-0.10108*/
   t[6]+=t[7]*29+32>>6; /* 0.38514*/
   t[6]-=t[5]*3+4>>3;   /*-0.31500*/
   t[5]+=t[6]*5+4>>3;   /* 0.62298*/
   t[5]-=t[4]+1>>1;     /*-0.45021*/
   t[4]+=t[5]*13+8>>4;  /* 0.82512*/
#else
   /*Optimal (1,2) regular coding gain with aribtrary precision: 9.56126 dB,
     6-bit (1,2) regular type-IV coding gain: 9.56051 dB
     S={11/8,75/64,19/16,19/16}, 0={53/64,-3/8}, 1={5/8,-5/16}, 2={3/8,-1/16}
     4-bit (1,2) regular type-IV coding gain: 9.55846 dB
     S={3/2,19/16,19/16,19/16}, 0={13/16,-3/8}, 1={5/8,-5/16}, 2={3/8,-1/16}
     */
   /*Scaling factors: the biorthogonal part.*/
   /*Note: t[i]+=t[i]>>15&1; is equivalent to: if(t[i]>0)t[i]++;
     This step ensures that the scaling is trivially invertible on the
      decoder's side, with perfect reconstruction.*/
   t[4]=t[4]*3>>2;      /* 1.40617*/
   t[4]+=-t[4]>>15&1;
   t[5]=t[5]*19>>4;     /* 1.21864*/
   t[5]+=-t[5]>>15&1;
   t[6]=t[6]*19>>4;     /* 1.19258*/
   t[6]+=-t[6]>>15&1;
   t[7]=t[7]*19>>4;     /* 1.21660*/
   t[7]+=-t[7]>>15&1;
   /*Rotations:*/
   t[5]-=t[4]*3+4>>3;   /*-0.45208*/
   t[6]-=t[5]*5+8>>4;   /*-0.30440*/
   t[7]-=t[6]+8>>4;     /*-0.10258*/
   t[6]+=t[7]*3+4>>3;   /* 0.37100*/
   t[5]+=t[6]*5+4>>3;   /* 0.61773*/
   t[4]+=t[5]*13+8>>4;  /* 0.82423*/
#endif
   /*More +1/-1 butterflies (required for FIR, PR, LP).*/
   t[0]+=t[7]>>1;
   _y[0]=(ogg_int16_t)t[0];
   t[1]+=t[6]>>1;
   _y[1]=(ogg_int16_t)t[1];
   t[2]+=t[5]>>1;
   _y[2]=(ogg_int16_t)t[2];
   t[3]+=t[4]>>1;
   _y[3]=(ogg_int16_t)t[3];
   _y[4]=(ogg_int16_t)(t[3]-t[4]);
   _y[5]=(ogg_int16_t)(t[2]-t[5]);
   _y[6]=(ogg_int16_t)(t[1]-t[6]);
   _y[7]=(ogg_int16_t)(t[0]-t[7]);
}



void od_post_filter8(ogg_int16_t _x[8],const ogg_int16_t _y[8]){
   int t[8];
   t[7]=_y[0]-_y[7];
   t[6]=_y[1]-_y[6];
   t[5]=_y[2]-_y[5];
   t[4]=_y[3]-_y[4];
   t[3]=_y[3]-(t[4]>>1);
   t[2]=_y[2]-(t[5]>>1);
   t[1]=_y[1]-(t[6]>>1);
   t[0]=_y[0]-(t[7]>>1);
   /*Optimal (1,2) regular coding gain with aribtrary precision: 9.56126 dB,
     6-bit (1,2) regular type-IV coding gain: 9.56051 dB
     S={11/8,75/64,19/16,19/16}, 0={53/64,-3/8}, 1={5/8,-5/16}, 2={3/8,-1/16}
     4-bit (1,2) regular type-IV coding gain: 9.55846 dB
     S={3/2,19/16,19/16,19/16}, 0={13/16,-3/8}, 1={5/8,-5/16}, 2={3/8,-1/16}*/
   t[4]-=t[5]*13+8>>4;  /* 0.82423*/
   t[5]-=t[6]*5+4>>3;   /* 0.61773*/
   t[6]-=t[7]*3+4>>3;   /* 0.37100*/
   t[7]+=t[6]+8>>4;     /*-0.10258*/
   t[6]+=t[5]*5+8>>4;   /*-0.30440*/
   t[5]+=t[4]*3+4>>3;   /*-0.45208*/
   t[7]=(t[7]<<4)/19;   /* 1.21660*/
   t[6]=(t[6]<<4)/19;   /* 1.19258*/
   t[5]=(t[5]<<4)/19;   /* 1.21864*/
   t[4]=(t[4]<<2)/3;    /* 1.40617*/
   t[0]+=t[7]>>1;
   _x[0]=(ogg_int16_t)t[0];
   t[1]+=t[6]>>1;
   _x[1]=(ogg_int16_t)t[1];
   t[2]+=t[5]>>1;
   _x[2]=(ogg_int16_t)t[2];
   t[3]+=t[4]>>1;
   _x[3]=(ogg_int16_t)t[3];
   _x[4]=(ogg_int16_t)(t[3]-t[4]);
   _x[5]=(ogg_int16_t)(t[2]-t[5]);
   _x[6]=(ogg_int16_t)(t[1]-t[6]);
   _x[7]=(ogg_int16_t)(t[0]-t[7]);
}



void od_pre_filter16(ogg_int16_t _y[16],const ogg_int16_t _x[16]){
   int t[16];
   /*
"Found: "{9.788267439364699, {3/2, 39/32, 77/64, 39/32, 19/16, 87/64, 35/32,`
    1}, {29/32, -7/16, 13/16, -27/64, 45/64, -3/8, 37/64, -1/4, 7/16, -5/64,`
    7/64, 1/8, 0, 1/16}}
Found: {9.791121179618907, {3/2, 39/32, 35/32, 37/32, 77/64, 87/64, 35/32,`
    1}, {29/32, -7/16, 13/16, -13/32, 47/64, -3/8, 19/32, -9/32, 7/16, -5/64,`
    7/64, 1/8, 0, 1/16}}
   */
   /*+1/-1 butterflies (required for FIR, PR, LP).*/
   t[15]=_x[0]-_x[15];
   t[14]=_x[1]-_x[14];
   t[13]=_x[2]-_x[13];
   t[12]=_x[3]-_x[12];
   t[11]=_x[4]-_x[11];
   t[10]=_x[5]-_x[10];
   t[9]=_x[6]-_x[9];
   t[8]=_x[7]-_x[8];
   t[7]=_x[7]-(t[8]>>1);
   t[6]=_x[6]-(t[9]>>1);
   t[5]=_x[5]-(t[10]>>1);
   t[4]=_x[4]-(t[11]>>1);
   t[3]=_x[3]-(t[12]>>1);
   t[2]=_x[2]-(t[13]>>1);
   t[1]=_x[1]-(t[14]>>1);
   t[0]=_x[0]-(t[15]>>1);
   /*U filter (arbitrary invertible, omitted).*/
   /*V filter (arbitrary invertible, can be optimized for speed or accuracy).*/
   /*Scaling factors: the biorthogonal part.*/
   t[8]=(t[8]<<5)/23;     /*1.41*/
   t[8]+=-t[8]>>15&1;
   /*if(t[8]>0)t[8]++;*/
   t[9]=(t[9]<<5)/29;     /*1.09*/
   t[9]+=-t[9]>>15&1;
   /*if(t[9]>0)t[9]++;*/
   t[10]=(t[10]<<4)/15;   /*1.06*/
   t[10]+=-t[10]>>15&1;
   /*if(t[10]>0)t[10]++;*/
   t[11]=(t[11]<<4)/15;   /*1.06*/
   t[11]+=-t[11]>>15&1;
   /*if(t[11]>0)t[11]++;*/
   t[12]=(t[12]<<4)/15;   /*1.06*/
   t[12]+=-t[12]>>15&1;
   /*if(t[12]>0)t[12]++;*/
   t[13]=(t[13]<<4)/15;   /*1.06*/
   t[13]+=-t[13]>>15&1;
   /*if(t[13]>0)t[13]++;*/
   t[14]=(t[14]<<4)/15;   /*1.08*/
   t[14]+=-t[14]>>15&1;
   /*if(t[14]>0)t[14]++;*/
   t[15]=(t[15]<<5)/29;   /*1.11*/
   t[15]+=-t[15]>>15&1;
   /*if(t[15]>0)t[15]++;*/
   /*Rotations:*/
   t[15]-=t[14]*3+16>>5;  /*-0.09*/
   t[14]+=t[15]*13+32>>6; /*0.21*/
   t[14]-=t[13]*15+32>>6; /*-0.24*/
   t[13]+=t[14]*3+4>>3;   /*0.37*/
   t[13]-=t[12]*23+32>>6; /*-0.36*/
   t[12]+=t[13]+1>>1;     /*0.50*/
   t[12]-=t[11]*15+16>>5; /*-0.46*/
   t[11]+=t[12]*19+16>>5; /*0.60*/
   t[11]-=t[10]*17+16>>5; /*-0.53*/
   t[10]+=t[11]*43+32>>6; /*0.67*/
   t[10]-=t[9]*9+8>>4;    /*-0.56*/
   t[9]+=t[10]*47+32>>6;  /*0.74*/
   t[9]-=t[8]*15+16>>5;   /*-0.47*/
   t[8]+=t[9]*55+32>>6;   /*0.86*/
   /*More +1/-1 butterflies (required for FIR, PR, LP).*/
   t[0]+=t[15]>>1;
   _y[0]=(ogg_int16_t)t[0];
   t[1]+=t[14]>>1;
   _y[1]=(ogg_int16_t)t[1];
   t[2]+=t[13]>>1;
   _y[2]=(ogg_int16_t)t[2];
   t[3]+=t[12]>>1;
   _y[3]=(ogg_int16_t)t[3];
   t[4]+=t[11]>>1;
   _y[4]=(ogg_int16_t)t[4];
   t[5]+=t[10]>>1;
   _y[5]=(ogg_int16_t)t[5];
   t[6]+=t[9]>>1;
   _y[6]=(ogg_int16_t)t[6];
   t[7]+=t[8]>>1;
   _y[7]=(ogg_int16_t)t[7];
   _y[8]=(ogg_int16_t)(t[7]-t[8]);
   _y[9]=(ogg_int16_t)(t[6]-t[9]);
   _y[10]=(ogg_int16_t)(t[5]-t[10]);
   _y[11]=(ogg_int16_t)(t[4]-t[11]);
   _y[12]=(ogg_int16_t)(t[3]-t[12]);
   _y[13]=(ogg_int16_t)(t[2]-t[13]);
   _y[14]=(ogg_int16_t)(t[1]-t[14]);
   _y[15]=(ogg_int16_t)(t[0]-t[15]);
}

void od_post_filter16(ogg_int16_t _x[16],const ogg_int16_t _y[16]){
   int t[16];
   t[15]=_y[0]-_y[15];
   t[14]=_y[1]-_y[14];
   t[13]=_y[2]-_y[13];
   t[12]=_y[3]-_y[12];
   t[11]=_y[4]-_y[11];
   t[10]=_y[5]-_y[10];
   t[9]=_y[6]-_y[9];
   t[8]=_y[7]-_y[8];
   t[7]=_y[7]-(t[8]>>1);
   t[6]=_y[6]-(t[9]>>1);
   t[5]=_y[5]-(t[10]>>1);
   t[4]=_y[4]-(t[11]>>1);
   t[3]=_y[3]-(t[12]>>1);
   t[2]=_y[2]-(t[13]>>1);
   t[1]=_y[1]-(t[14]>>1);
   t[0]=_y[0]-(t[15]>>1);
   t[8]-=t[9]*55+32>>6;
   t[9]+=t[8]*15+16>>5;
   t[9]-=t[10]*47+32>>6;
   t[10]+=t[9]*9+8>>4;
   t[10]-=t[11]*43+32>>6;
   t[11]+=t[10]*17+16>>5;
   t[11]-=t[12]*19+16>>5;
   t[12]+=t[11]*15+16>>5;
   t[12]-=t[13]+1>>1;
   t[13]+=t[12]*23+32>>6;
   t[13]-=t[14]*3+4>>3;
   t[14]+=t[13]*15+32>>6;
   t[14]-=t[15]*13+32>>6;
   t[15]+=t[14]*3+16>>5;
   t[15]=t[15]*29>>5;
   t[14]=t[14]*15>>4;
   t[13]=t[13]*15>>4;
   t[12]=t[12]*15>>4;
   t[11]=t[11]*15>>4;
   t[10]=t[10]*15>>4;
   t[9]=t[9]*29>>5;
   t[8]=t[8]*23>>5;
   t[0]+=t[15]>>1;
   _x[0]=(ogg_int16_t)t[0];
   t[1]+=t[14]>>1;
   _x[1]=(ogg_int16_t)t[1];
   t[2]+=t[13]>>1;
   _x[2]=(ogg_int16_t)t[2];
   t[3]+=t[12]>>1;
   _x[3]=(ogg_int16_t)t[3];
   t[4]+=t[11]>>1;
   _x[4]=(ogg_int16_t)t[4];
   t[5]+=t[10]>>1;
   _x[5]=(ogg_int16_t)t[5];
   t[6]+=t[9]>>1;
   _x[6]=(ogg_int16_t)t[6];
   t[7]+=t[8]>>1;
   _x[7]=(ogg_int16_t)t[7];
   _x[8]=(ogg_int16_t)(t[7]-t[8]);
   _x[9]=(ogg_int16_t)(t[6]-t[9]);
   _x[10]=(ogg_int16_t)(t[5]-t[10]);
   _x[11]=(ogg_int16_t)(t[4]-t[11]);
   _x[12]=(ogg_int16_t)(t[3]-t[12]);
   _x[13]=(ogg_int16_t)(t[2]-t[13]);
   _x[14]=(ogg_int16_t)(t[1]-t[14]);
   _x[15]=(ogg_int16_t)(t[0]-t[15]);
}


#if defined(TEST)
#include <stdio.h>

int minv[32];
int maxv[32];

int main(void){
   ogg_int32_t min[16];
   ogg_int32_t max[16];
   int         mini[16];
   int         maxi[16];
   int         i;
   int         j;
   for(j=0;j<16;j++)min[j]=max[j]=mini[j]=maxi[j]=0;
   for(i=0;i<256;i++){
      ogg_int16_t x[16];
      for(j=0;j<16;j++)x[j]=i>>j&1?1325:-1325;
      od_pre_filter8(x,x);
      /*for(j=0;j<16;j++){
         if(x[j]<min[j]){
            min[j]=x[j];
            mini[j]=i;
         }
         else if(x[j]>max[j]){
            max[j]=x[j];
            maxi[j]=i;
         }
      }*/
   }
   for(j=0;j<18;j++){
      printf("Min: %5i  ",minv[j]);
      printf("Max: %5i\n",maxv[j]);
   }
   return 0;
}
#endif
