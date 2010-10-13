/*
    Daala video codec
    Copyright (C) 2006-2010 Daala project contributors

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


#include <stdio.h>
#include <stddef.h>
#include "mc.h"



/*A table of indices used to set up the rotated versions of each vector
   interpolation formula.*/
static const int MIDXS[][4]={
  {0,1,2,3},/*0*/
  {0,0,2,3},/*1*/
  {1,1,2,3},/*2*/
  {0,0,3,3},/*3*/
  {1,1,2,2},/*4*/
  {0,1,0,0},/*5*/
  {2,1,2,2},/*6*/
  {0,1,1,1},/*7*/
};

/*Set up the finite differences needed to interpolate a motion vector
   component.
  _dmv:         Returns the motion vector deltas.
  _dmv[0]:      The initial value.
  _dmv[1]:      The initial amount to increment per unit change in i.
  _dmv[2]:      The amount to increment per unit change in j.
  _dmv[3]:      The amount to increment _dxdi by per unit change in j.
  _mvx:         The component value of the 4 motion vectors.
  _m:           The index of the motion vector to use for each corner in the
                 base orientation.
  _r:           The amount to rotate (clockwise) the formulas by (0...3).
  _log_xblk_sz: The log base 2 of the horizontal block dimension.
  _log_yblk_sz: The log base 2 of the vertical block dimension.*/
void od_mc_setup_mvc(ogg_int32_t _dmv[4],const ogg_int32_t _mvs[4],
 const int _m[4],int _r,int _log_xblk_sz,int _log_yblk_sz){
  _dmv[0]=_mvs[_m[0-_r&3]+_r&3];
  _dmv[1]=_mvs[_m[1-_r&3]+_r&3]-_dmv[0]>>_log_xblk_sz;
  _dmv[2]=_mvs[_m[3-_r&3]+_r&3]-_dmv[0]>>_log_yblk_sz;
  _dmv[3]=_mvs[_m[0-_r&3]+_r&3]+_mvs[_m[2-_r&3]+_r&3]-
   _mvs[_m[1-_r&3]+_r&3]-_mvs[_m[3-_r&3]+_r&3]>>_log_xblk_sz+_log_yblk_sz;
  /*Advance the vector to the (0.5,0.5) position.*/
  _dmv[0]+=_dmv[2]>>1;
  _dmv[1]+=_dmv[3]>>1;
  _dmv[0]+=_dmv[1]>>1;
  _dmv[2]+=_dmv[3]>>1;
}

/*Form the prediction given by one interpolated motion vector.
  _dst:         The destination buffer (xstride must be 1).
  _dystride:    The byte offset between destination pixel rows.
  _src:         The source buffer (xstride must be 1).
  _systride:    The byte offset between source pixel rows.
  _mvx:         The X component of the motion vectors.
  _mvy:         The Y component of the motion vectors.
  _m:           The index of the motion vector to use for each corner in the
                 base orientation.
  _r:           The amount to rotate (clockwise) the formulas by (0...3).
  _log_xblk_sz: The log base 2 of the horizontal block dimension.
  _log_yblk_sz: The log base 2 of the vertical block dimension.*/
void od_mc_predict1imv8_c(unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],const int _m[4],int _r,int _log_xblk_sz,
 int _log_yblk_sz){
  ogg_int32_t x;
  ogg_int32_t y;
  ogg_int32_t dmvx[4];
  ogg_int32_t dmvy[4];
  int         xblk_sz;
  int         yblk_sz;
  int         i;
  int         j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  od_mc_setup_mvc(dmvx,_mvx,_m,_r,_log_xblk_sz,_log_yblk_sz);
  od_mc_setup_mvc(dmvy,_mvy,_m,_r,_log_xblk_sz,_log_yblk_sz);
  if(_dystride==xblk_sz&&
   !dmvx[1]&&!dmvy[1]&&!dmvx[2]&&!dmvy[2]&&!dmvx[3]&&!dmvy[3]){
    od_mc_predict1fmv8_c(_dst,_src,_systride,dmvx[0],dmvy[0],
     _log_xblk_sz,_log_yblk_sz);
  }
  dmvx[1]+=1<<17;
  dmvy[2]+=1<<17;
  for(j=0;j<yblk_sz;j++){
    x=dmvx[0];
    y=dmvy[0];
    for(i=0;i<xblk_sz;i++){
      const unsigned char *p;
      ogg_int32_t          xf;
      ogg_int32_t          yf;
      unsigned             a;
      unsigned             b;
      /*printf("<%16.12f,%16.12f>%s",(x-(i<<17))/(double)0x20000,
       (y-(j<<17))/(double)0x20000,i+1<xblk_sz?"::":"\n");*/
      xf=x&0xFFFF;
      yf=y&0xFFFF;
      p=_src+(x>>16)+(y>>16)*_systride;
      a=(unsigned)(p[0]+((p[1]-p[0])*xf>>16));
      b=(unsigned)((p+_systride)[0]+
       (((p+_systride)[1]-(p+_systride)[0])*xf>>16));
      _dst[i]=(unsigned char)(a+((b-a)*yf>>16));
      x+=dmvx[1];
      y+=dmvy[1];
    }
    dmvx[0]+=dmvx[2];
    dmvy[0]+=dmvy[2];
    dmvx[1]+=dmvx[3];
    dmvy[1]+=dmvy[3];
    _dst+=_dystride;
  }
  /*_dst-=_dystride*yblk_sz;
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++)printf("%2X ",*(_dst+i+j*_dystride));
    printf("\n");
  }*/
}

static void od_mc_predict1imv8(od_state *_state,unsigned char *_dst,
 int _dystride,const unsigned char *_src,int _systride,
 const ogg_int32_t _mvx[4],const ogg_int32_t _mvy[4],const int _m[4],int _r,
 int _log_xblk_sz,int _log_yblk_sz){
  (*_state->opt_vtbl.mc_predict1imv8)(_dst,_dystride,_src,_systride,
   _mvx,_mvy,_m,_r,_log_xblk_sz,_log_yblk_sz);
}

/*Form the prediction given by one fixed motion vector.
  _dst:         The destination buffer (xstride must be 1).
  _src:         The source buffer (xstride must be 1).
  _systride:    The byte offset between source pixel rows.
  _mvx:         The X component of the motion vector.
  _mvy:         The Y component of the motion vector.
  _log_xblk_sz: The log base 2 of the horizontal block dimension.
  _log_yblk_sz: The log base 2 of the vertical block dimension.*/
void od_mc_predict1fmv8_c(unsigned char *_dst,const unsigned char *_src,
 int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz){
  ogg_uint32_t mvxf;
  ogg_uint32_t mvyf;
  int         xblk_sz;
  int         yblk_sz;
  int         i;
  int         j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  _src+=(_mvx>>16)+(_mvy>>16)*_systride;
  mvxf=(ogg_uint32_t)(_mvx&0xFFFF);
  mvyf=(ogg_uint32_t)(_mvy&0xFFFF);
  if(mvxf!=0){
    if(mvyf!=0){
      for(j=0;j<yblk_sz;j++){
        for(i=0;i<xblk_sz;i++){
          ogg_uint32_t a;
          ogg_uint32_t b;
          /*printf("<%16.12f,%16.12f>%s",_mvx/(double)0x40000,
           _mvy/(double)0x40000,i+1<xblk_sz?"::":"\n");*/
          a=((ogg_uint32_t)_src[i<<1]<<16)+(_src[i<<1|1]-_src[i<<1])*mvxf>>16;
          b=((ogg_uint32_t)(_src+_systride)[i<<1]<<16)+
           ((_src+_systride)[i<<1|1]-(_src+_systride)[i<<1])*mvxf>>16;
          _dst[j*xblk_sz+i]=(unsigned char)((a<<16)+(b-a)*mvyf>>16);
        }
        _src+=_systride<<1;
      }
    }
    else{
      for(j=0;j<yblk_sz;j++){
        for(i=0;i<xblk_sz;i++){
          /*printf("<%16.12f,%16.12f>%s",_mvx/(double)0x40000,
           _mvy/(double)0x40000,i+1<xblk_sz?"::":"\n");*/
          _dst[j*xblk_sz+i]=(unsigned char)(
           ((ogg_uint32_t)_src[i<<1]<<16)+(_src[i<<1|1]-_src[i<<1])*mvxf>>16);
        }
        _src+=_systride<<1;
      }
    }
  }
  else{
    if(mvyf!=0){
      for(j=0;j<yblk_sz;j++){
        for(i=0;i<xblk_sz;i++){
          /*printf("<%16.12f,%16.12f>%s",_mvx/(double)0x40000,
           _mvy/(double)0x40000,i+1<xblk_sz?"::":"\n");*/
          _dst[j*xblk_sz+i]=(unsigned char)(((ogg_uint32_t)_src[i<<1]<<16)+
           ((_src+_systride)[i<<1]-_src[i<<1])*mvyf>>16);
        }
        _src+=_systride<<1;
      }
    }
    else{
      for(j=0;j<yblk_sz;j++){
        for(i=0;i<xblk_sz;i++){
          /*printf("<%16.12f,%16.12f>%s",_mvx/(double)0x40000,
           _mvy/(double)0x40000,i+1<xblk_sz?"::":"\n");*/
          _dst[j*xblk_sz+i]=_src[i<<1];
        }
        _src+=_systride<<1;
      }
    }
  }
  /*_dst-=xblk_sz*yblk_sz;
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++)printf("%2X ",*(_dst+i+j*xblk_sz));
    printf("\n");
  }*/
}

static void od_mc_predict1fmv8(od_state *_state,unsigned char *_dst,
 const unsigned char *_src,int _systride,ogg_int32_t _mvx,ogg_int32_t _mvy,
 int _log_xblk_sz,int _log_yblk_sz){
  (*_state->opt_vtbl.mc_predict1fmv8)(_dst,_src,_systride,_mvx,_mvy,
   _log_xblk_sz,_log_yblk_sz);
}

/*Perform normal bilinear blending.*/
void od_mc_blend_full8_c(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  unsigned a;
  unsigned b;
  int      log_blk_sz2;
  int      xblk_sz;
  int      yblk_sz;
  int      round;
  int      i;
  int      j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  log_blk_sz2=_log_xblk_sz+_log_yblk_sz;
  round=1<<log_blk_sz2-1;
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      a=_src[0][j*xblk_sz+i];
      b=_src[3][j*xblk_sz+i];
      a=(a<<_log_xblk_sz)+(_src[1][j*xblk_sz+i]-a)*i;
      b=(b<<_log_xblk_sz)+(_src[2][j*xblk_sz+i]-b)*i;
      _dst[i]=(unsigned char)((a<<_log_yblk_sz)+(b-a)*j+round>>log_blk_sz2);
    }
    _dst+=_dystride;
  }
}

/*Perform normal bilinear blending.*/
static void od_mc_blend_full8(od_state *_state,unsigned char *_dst,
 int _dystride,const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  (*_state->opt_vtbl.mc_blend_full8)(_dst,_dystride,_src,
   _log_xblk_sz,_log_yblk_sz);
}

#if 0
/*Perform multiresolution bilinear blending.*/
static void od_mc_blend_multi8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  const unsigned char *p;
  unsigned char       *dst0;
  unsigned char       *dst;
  ptrdiff_t            o;
  ptrdiff_t            o0;
  int                  ll[4];
  int                  lh;
  int                  hl;
  int                  hh;
  int                  a;
  int                  b;
  int                  c;
  int                  d;
  int                  log_blk_sz2;
  int                  xblk_sz;
  int                  yblk_sz;
  int                  xblk_sz_2;
  int                  yblk_sz_2;
  int                  i;
  int                  j;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  log_blk_sz2=_log_xblk_sz+_log_yblk_sz;
  o0=0;
  dst0=_dst;
  /*Perform multiresolution blending.*/
  xblk_sz_2=xblk_sz>>1;
  yblk_sz_2=yblk_sz>>1;
  for(j=1;j<yblk_sz_2;j+=2){
    o=o0;
    dst=dst0;
    /*Upper-left quadrant.*/
    for(i=1;i<xblk_sz_2;i+=2){
      p=_src[0]+o;
      /*Forward Haar wavelet.*/
      ll[0]=p[0]+p[1];
      lh=p[0]-p[1];
      hl=(p+xblk_sz)[0]+(p+xblk_sz)[1];
      hh=(p+xblk_sz)[0]-(p+xblk_sz)[1];
      c=ll[0]-hl;
      ll[0]+=hl;
      hl=c;
      /*No need to finish the transform; we'd just invert it later.*/
      p=_src[1]+o;
      ll[1]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[2]+o;
      ll[2]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[3]+o;
      ll[3]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      /*LL blending.*/
      a=(ll[0]<<_log_xblk_sz)+(ll[1]-ll[0])*i;
      b=(ll[3]<<_log_xblk_sz)+(ll[2]-ll[3])*i;
      a=(int)(((ogg_int32_t)a<<_log_yblk_sz)+
       (ogg_int32_t)(b-a)*j>>log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c=a-hl+1>>1;
      a=a+hl+1>>1;
      d=c-hh+1>>1;
      c=c+hh+1>>1;
      b=a-lh+1>>1;
      a=a+lh+1>>1;
      dst[0]=OD_CLAMP255(a);
      dst[1]=OD_CLAMP255(b);
      (dst+_dystride)[0]=OD_CLAMP255(c);
      (dst+_dystride)[1]=OD_CLAMP255(d);
      o+=2;
      dst+=2;
    }
    /*Upper-right quadrant.*/
    for(;i<xblk_sz;i+=2){
      p=_src[0]+o;
      ll[0]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[1]+o;
      /*Forward Haar wavelet.*/
      ll[1]=p[0]+p[1];
      lh=p[0]-p[1];
      hl=(p+xblk_sz)[0]+(p+xblk_sz)[1];
      hh=(p+xblk_sz)[0]-(p+xblk_sz)[1];
      c=ll[1]-hl;
      ll[1]+=hl;
      hl=c;
      /*No need to finish the transform; we'd just invert it later.*/
      p=_src[2]+o;
      ll[2]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[3]+o;
      ll[3]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      /*LL blending.*/
      a=(ll[0]<<_log_xblk_sz)+(ll[1]-ll[0])*i;
      b=(ll[3]<<_log_xblk_sz)+(ll[2]-ll[3])*i;
      a=(int)(((ogg_int32_t)a<<_log_yblk_sz)+
       (ogg_int32_t)(b-a)*j>>log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c=a-hl+1>>1;
      a=a+hl+1>>1;
      d=c-hh+1>>1;
      c=c+hh+1>>1;
      b=a-lh+1>>1;
      a=a+lh+1>>1;
      dst[0]=OD_CLAMP255(a);
      dst[1]=OD_CLAMP255(b);
      (dst+_dystride)[0]=OD_CLAMP255(c);
      (dst+_dystride)[1]=OD_CLAMP255(d);
      o+=2;
      dst+=2;
    }
    o0+=xblk_sz<<1;
    dst0+=_dystride<<1;
  }
  for(;j<yblk_sz;j+=2){
    o=o0;
    dst=dst0;
    /*Lower-left quadrant.*/
    for(i=1;i<xblk_sz_2;i+=2){
      p=_src[0]+o;
      ll[0]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[1]+o;
      ll[1]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[2]+o;
      ll[2]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[3]+o;
      /*Forward Haar wavelet.*/
      ll[3]=p[0]+p[1];
      lh=p[0]-p[1];
      hl=(p+xblk_sz)[0]+(p+xblk_sz)[1];
      hh=(p+xblk_sz)[0]-(p+xblk_sz)[1];
      c=ll[3]-hl;
      ll[3]+=hl;
      hl=c;
      /*No need to finish the transform; we'd just invert it later.*/
      /*LL blending.*/
      a=(ll[0]<<_log_xblk_sz)+(ll[1]-ll[0])*i;
      b=(ll[3]<<_log_xblk_sz)+(ll[2]-ll[3])*i;
      a=(int)(((ogg_int32_t)a<<_log_yblk_sz)+
       (ogg_int32_t)(b-a)*j>>log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c=a-hl+1>>1;
      a=a+hl+1>>1;
      d=c-hh+1>>1;
      c=c+hh+1>>1;
      b=a-lh+1>>1;
      a=a+lh+1>>1;
      dst[0]=OD_CLAMP255(a);
      dst[1]=OD_CLAMP255(b);
      (dst+_dystride)[0]=OD_CLAMP255(c);
      (dst+_dystride)[1]=OD_CLAMP255(d);
      o+=2;
      dst+=2;
    }
    /*Lower-right quadrant.*/
    for(;i<xblk_sz;i+=2){
      p=_src[0]+o;
      ll[0]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[1]+o;
      ll[1]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[2]+o;
      /*Forward Haar wavelet.*/
      ll[2]=p[0]+p[1];
      lh=p[0]-p[1];
      hl=(p+xblk_sz)[0]+(p+xblk_sz)[1];
      hh=(p+xblk_sz)[0]-(p+xblk_sz)[1];
      c=ll[2]-hl;
      ll[2]+=hl;
      hl=c;
      /*No need to finish the transform; we'd just invert it later.*/
      p=_src[3]+o;
      ll[3]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      /*LL blending.*/
      a=(ll[0]<<_log_xblk_sz)+(ll[1]-ll[0])*i;
      b=(ll[3]<<_log_xblk_sz)+(ll[2]-ll[3])*i;
      a=(int)(((ogg_int32_t)a<<_log_yblk_sz)+
       (ogg_int32_t)(b-a)*j>>log_blk_sz2);
      /*Inverse Haar wavelet.*/
      c=a-hl+1>>1;
      a=a+hl+1>>1;
      d=c-hh+1>>1;
      c=c+hh+1>>1;
      b=a-lh+1>>1;
      a=a+lh+1>>1;
      dst[0]=OD_CLAMP255(a);
      dst[1]=OD_CLAMP255(b);
      (dst+_dystride)[0]=OD_CLAMP255(c);
      (dst+_dystride)[1]=OD_CLAMP255(d);
      o+=2;
      dst+=2;
    }
    o0+=xblk_sz<<1;
    dst0+=_dystride<<1;
  }
}
#else

/*Perform multiresolution bilinear blending.*/
static void od_mc_blend_multi8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _log_xblk_sz,int _log_yblk_sz){
  unsigned char        src_ll[4][8][8];
  int                  dst_ll[8][8];
  const unsigned char *p;
  unsigned char       *dst;
  ptrdiff_t            o;
  int                  a;
  int                  b;
  int                  c;
  int                  d;
  int                  e;
  int                  f;
  int                  g;
  int                  h;
  int                  log_blk_sz2;
  int                  xblk_sz;
  int                  yblk_sz;
  int                  xblk_sz_2;
  int                  yblk_sz_2;
  int                  xblk_sz_4;
  int                  yblk_sz_4;
  int                  round;
  int                  i;
  int                  i2;
  int                  j;
  int                  j2;
  int                  k;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  /*Perform multiresolution blending.*/
  xblk_sz_2=xblk_sz>>1;
  yblk_sz_2=yblk_sz>>1;
  log_blk_sz2=_log_xblk_sz+_log_yblk_sz;
  round=1<<log_blk_sz2-1;
  /*Compute the low-pass band for each src block.*/
  for(k=0;k<4;k++){
    unsigned lh[4][8];
    p=_src[k];
    src_ll[k][0][0]=p[0];
    for(i=1;i<xblk_sz_2;i++){
      i2=i<<1;
      src_ll[k][0][i]=(unsigned char)(p[i2-1]+2*p[i2]+p[i2+1]+2>>2);
    }
    p+=xblk_sz;
    lh[1][0]=p[0]<<2;
    for(i=1;i<xblk_sz_2;i++){
      i2=i<<1;
      lh[1][i]=p[i2-1]+2*p[i2]+p[i2+1];
    }
    p+=xblk_sz;
    for(j=1;j<yblk_sz_2;j++){
      j2=j<<1&3;
      lh[j2][0]=p[0]<<2;
      for(i=1;i<xblk_sz_2;i++){
        i2=i<<1;
        lh[j2][i]=p[i2-1]+2*p[i2]+p[i2+1];
      }
      p+=xblk_sz;
      lh[j2+1][0]=p[0]<<2;
      for(i=1;i<xblk_sz_2;i++){
        i2=i<<1;
        lh[j2+1][i]=p[i2-1]+2*p[i2]+p[i2+1];
      }
      p+=xblk_sz;
      for(i=0;i<xblk_sz_2;i++){
        src_ll[k][j][i]=(unsigned char)
         (lh[j2-1&3][i]+2*lh[j2][i]+lh[j2+1][i]+8>>4);
      }
    }
  }
  /*Blend the low-pass bands.*/
  for(j=0;j<xblk_sz_2;j++){
    for(i=0;i<xblk_sz_2;i++){
      a=(src_ll[0][j][i]<<_log_xblk_sz-1)+(src_ll[1][j][i]-src_ll[0][j][i])*i;
      b=(src_ll[3][j][i]<<_log_xblk_sz-1)+(src_ll[2][j][i]-src_ll[3][j][i])*i;
      dst_ll[j][i]=(a<<_log_yblk_sz-1)+(b-a)*j;
    }
  }
  /*Perform the high-pass filtering for each quadrant.*/
  xblk_sz_4=xblk_sz>>2;
  yblk_sz_4=yblk_sz>>2;
  o=0;
  dst=_dst;
  for(j=0;j<yblk_sz_4;j++){
    /*Upper-left quadrant.*/
    for(i=0;i<xblk_sz_4;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[0][j][i]<<log_blk_sz2;
      f=src_ll[0][j][i]+src_ll[0][j][i+1]<<log_blk_sz2-1;
      g=src_ll[0][j][i]+src_ll[0][j+1][i]<<log_blk_sz2-1;
      h=src_ll[0][j][i]+src_ll[0][j][i+1]+
       src_ll[0][j+1][i]+src_ll[0][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(
       ((_src[0]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(
       ((_src[0]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[0]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[0]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
    }
    /*Upper-right quadrant.*/
    for(;i<xblk_sz_2-1;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[1][j][i]<<log_blk_sz2;
      f=src_ll[1][j][i]+src_ll[1][j][i+1]<<log_blk_sz2-1;
      g=src_ll[1][j][i]+src_ll[1][j+1][i]<<log_blk_sz2-1;
      h=src_ll[1][j][i]+src_ll[1][j][i+1]+
       src_ll[1][j+1][i]+src_ll[1][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(
       ((_src[1]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(
       ((_src[1]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[1]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[1]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
    }
    /*Upper-right quadrant, last column.*/
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=3*dst_ll[j][i]-dst_ll[j][i-1]<<1;
    c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j+1][i])-(dst_ll[j][i-1]+dst_ll[j+1][i-1]);
    e=src_ll[1][j][i]<<log_blk_sz2;
    f=3*src_ll[1][j][i]-src_ll[1][j][i-1]<<log_blk_sz2-1;
    g=src_ll[1][j][i]+src_ll[1][j+1][i]<<log_blk_sz2-1;
    h=3*(src_ll[1][j][i]+src_ll[1][j+1][i])-
     (src_ll[1][j][i-1]+src_ll[1][j+1][i-1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(
     ((_src[1]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(
     ((_src[1]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[1]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[1]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
    o+=xblk_sz<<1;
    dst+=_dystride<<1;
  }
  for(;j<yblk_sz_2-1;j++){
    /*Lower-left quadrant.*/
    for(i=0;i<xblk_sz_4;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[3][j][i]<<log_blk_sz2;
      f=src_ll[3][j][i]+src_ll[3][j][i+1]<<log_blk_sz2-1;
      g=src_ll[3][j][i]+src_ll[3][j+1][i]<<log_blk_sz2-1;
      h=src_ll[3][j][i]+src_ll[3][j][i+1]+
       src_ll[3][j+1][i]+src_ll[3][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(
       ((_src[3]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(
       ((_src[3]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[3]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[3]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
    }
    /*Lower-right quadrant.*/
    for(;i<xblk_sz_2-1;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[2][j][i]<<log_blk_sz2;
      f=src_ll[2][j][i]+src_ll[2][j][i+1]<<log_blk_sz2-1;
      g=src_ll[2][j][i]+src_ll[2][j+1][i]<<log_blk_sz2-1;
      h=src_ll[2][j][i]+src_ll[2][j][i+1]+
       src_ll[2][j+1][i]+src_ll[2][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(
       ((_src[2]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(
       ((_src[2]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[2]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[2]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
    }
    /*Lower-right quadrant, last column.*/
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=3*dst_ll[j][i]-dst_ll[j][i-1]<<1;
    c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j+1][i])-(dst_ll[j][i-1]+dst_ll[j+1][i-1]);
    e=src_ll[2][j][i]<<log_blk_sz2;
    f=3*src_ll[2][j][i]-src_ll[2][j][i-1]<<log_blk_sz2-1;
    g=src_ll[2][j][i]+src_ll[2][j+1][i]<<log_blk_sz2-1;
    h=3*(src_ll[2][j][i]+src_ll[2][j+1][i])-
     (src_ll[2][j][i-1]+src_ll[2][j+1][i-1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(
     ((_src[2]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(
     ((_src[2]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
    o+=xblk_sz<<1;
    dst+=_dystride<<1;
  }
  /*Lower-left quadrant, last row.*/
  for(i=0;i<xblk_sz_4;i++){
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
    c=3*dst_ll[j][i]-dst_ll[j-1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j][i+1])-(dst_ll[j-1][i]+dst_ll[j-1][i+1]);
    e=src_ll[3][j][i]<<log_blk_sz2;
    f=src_ll[3][j][i]+src_ll[3][j][i+1]<<log_blk_sz2-1;
    g=3*src_ll[3][j][i]-src_ll[3][j-1][i]<<log_blk_sz2-1;
    h=3*(src_ll[3][j][i]+src_ll[3][j][i+1])-
     (src_ll[3][j-1][i]+src_ll[3][j-1][i+1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(
     ((_src[3]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(
     ((_src[3]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[3]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[3]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
  }
  /*Lower-right quadrant, last row.*/
  for(;i<xblk_sz_2-1;i++){
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
    c=3*dst_ll[j][i]-dst_ll[j-1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j][i+1])-(dst_ll[j-1][i]+dst_ll[j-1][i+1]);
    e=src_ll[2][j][i]<<log_blk_sz2;
    f=src_ll[2][j][i]+src_ll[2][j][i+1]<<log_blk_sz2-1;
    g=3*src_ll[2][j][i]-src_ll[2][j-1][i]<<log_blk_sz2-1;
    h=3*(src_ll[2][j][i]+src_ll[2][j][i+1])-
     (src_ll[2][j-1][i]+src_ll[2][j-1][i+1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(
     ((_src[2]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(
     ((_src[2]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
  }
  /*Lower-right quadrant, last row and column.*/
  i2=i<<1;
  a=dst_ll[j][i]<<2;
  b=3*dst_ll[j][i]-dst_ll[j][i-1]<<1;
  c=3*dst_ll[j][i]-dst_ll[j-1][i]<<1;
  d=9*dst_ll[j][i]-3*(dst_ll[j-1][i]+dst_ll[j][i-1])+dst_ll[j-1][i-1];
  e=src_ll[2][j][i]<<log_blk_sz2;
  f=3*src_ll[2][j][i]-src_ll[2][j][i-1]<<log_blk_sz2-1;
  g=3*src_ll[2][j][i]-src_ll[2][j-1][i]<<log_blk_sz2-1;
  h=9*src_ll[2][j][i]-3*(src_ll[2][j][i-1]+src_ll[2][j-1][i])+
   src_ll[2][j-1][i-1]<<log_blk_sz2-2;
  dst[i2]=OD_CLAMP255(
   ((_src[2]+o)[i2]<<log_blk_sz2)+a-e+round>>log_blk_sz2);
  dst[i2+1]=OD_CLAMP255(
   ((_src[2]+o)[i2+1]<<log_blk_sz2)+b-f+round>>log_blk_sz2);
  (dst+_dystride)[i2]=OD_CLAMP255(
   ((_src[2]+o+xblk_sz)[i2]<<log_blk_sz2)+c-g+round>>log_blk_sz2);
  (dst+_dystride)[i2+1]=OD_CLAMP255(
   ((_src[2]+o+xblk_sz)[i2+1]<<log_blk_sz2)+d-h+round>>log_blk_sz2);
}
#endif

static void od_mc_setup_s_split(int _s0[4],int _dsdi[4],int _dsdj[4],
 int _ddsdidj[4],int _c,int _s,int _log_xblk_sz,int _log_yblk_sz){
  int log_blk_sz2;
  int k;
  log_blk_sz2=_log_xblk_sz+_log_yblk_sz;
  _s0[0]=2<<log_blk_sz2;
  _s0[1]=_s0[2]=_s0[3]=0;
  _dsdi[0]=-2<<_log_xblk_sz;
  _dsdi[1]=2<<_log_xblk_sz;
  _dsdi[2]=_dsdi[3]=0;
  _dsdj[0]=-2<<_log_yblk_sz;
  _dsdj[1]=_dsdj[2]=0;
  _dsdj[3]=2<<_log_yblk_sz;
  _ddsdidj[0]=_ddsdidj[2]=2;
  _ddsdidj[1]=_ddsdidj[3]=-2;
  if(!(_s&1)){
    k=_c+1&3;
    _s0[k]>>=1;
    _s0[_c]+=_s0[k];
    _dsdi[k]>>=1;
    _dsdi[_c]+=_dsdi[k];
    _dsdj[k]>>=1;
    _dsdj[_c]+=_dsdj[k];
    _ddsdidj[k]>>=1;
    _ddsdidj[_c]+=_ddsdidj[k];
  }
  if(!(_s&2)){
    k=_c+3&3;
    _s0[k]>>=1;
    _s0[_c]+=_s0[k];
    _dsdi[k]>>=1;
    _dsdi[_c]+=_dsdi[k];
    _dsdj[k]>>=1;
    _dsdj[_c]+=_dsdj[k];
    _ddsdidj[k]>>=1;
    _ddsdidj[_c]+=_ddsdidj[k];
  }
  /*Advance the weights to the (0.5,0.5) position.
    LOOP VECTORIZES.
  for(k=0;k<4;k++){
    _s0[k]+=_dsdj[k]>>1;
    _dsdi[k]+=_ddsdidj[k]>>1;
    _s0[k]+=_dsdi[k]>>1;
    _dsdj[k]+=_ddsdidj[k]>>1;
  }*/
}

/*Perform normal blending with bilinear weights modified for unsplit edges.*/
void od_mc_blend_full_split8_c(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz){
  int a;
  int b;
  int c;
  int d;
  int s[4];
  int s0[4];
  int dsdi[4];
  int dsdj[4];
  int ddsdidj[4];
  int xblk_sz;
  int yblk_sz;
  int log_blk_sz2p1;
  int round;
  int i;
  int j;
  int k;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  /*The block is too small; perform normal blending.*/
  log_blk_sz2p1=_log_xblk_sz+_log_yblk_sz+1;
  round=1<<log_blk_sz2p1-1;
  od_mc_setup_s_split(s0,dsdi,dsdj,ddsdidj,_c,_s,_log_xblk_sz,_log_yblk_sz);
  /*LOOP VECTORIZES.*/
  for(k=0;k<4;k++)s[k]=s0[k];
  for(j=0;j<yblk_sz;j++){
    for(i=0;i<xblk_sz;i++){
      a=_src[0][j*xblk_sz+i];
      b=(_src[1][j*xblk_sz+i]-a)*s[1];
      c=(_src[2][j*xblk_sz+i]-a)*s[2];
      d=(_src[3][j*xblk_sz+i]-a)*s[3];
      _dst[i]=(unsigned char)((a<<log_blk_sz2p1)+b+c+d+round>>log_blk_sz2p1);
      /*LOOP VECTORIZES.*/
      for(k=0;k<4;k++)s[k]+=dsdi[k];
    }
    _dst+=_dystride;
    /*LOOP VECTORIZES.*/
    for(k=0;k<4;k++){
      s0[k]+=dsdj[k];
      s[k]=s0[k];
      dsdi[k]+=ddsdidj[k];
    }
  }
}

/*Perform normal blending with bilinear weights modified for unsplit edges.*/
void od_mc_blend_full_split8(od_state *_state,unsigned char *_dst,
 int _dystride,const unsigned char *_src[4],int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz){
  (*_state->opt_vtbl.mc_blend_full_split8)(_dst,_dystride,_src,_c,_s,
   _log_xblk_sz,_log_yblk_sz);
}

#if 0
/*There are other ways to implement multiresolution blending for the modified
   bilinear weights, but using a table lookup to select which predictor to
   draw the high-frequency coefficients from moves all the complexity into the
   tables, and leaves the code dead simple.*/

/*The MV from which to use the high-frequency coefficients for a 2x2 LL band.*/
static const unsigned char OD_MC_SIDXS_22[3][4][4]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,
      0,2
    },{
      /*Corner: 1; split: none*/
      1,1,
      3,1
    },{
      /*Corner: 2; split: none*/
      0,2,
      2,2
    },{
      /*Corner: 3; split: none*/
      3,1,
      3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,1,
      0,2
    },{
      /*Corner: 1; split: 2*/
      1,1,
      3,2
    },{
      /*Corner: 2; split: 3*/
      0,2,
      3,2
    },{
      /*Corner: 3; split: 0*/
      1,1,
      3,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,
      3,2
    },{
      /*Corner: 1; split: 0*/
      0,1,
      3,1
    },{
      /*Corner: 2; split: 1*/
      0,1,
      2,2
    },{
      /*Corner: 3; split: 2*/
      3,1,
      3,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 2x4 LL band.*/
static const unsigned char OD_MC_SIDXS_24[3][4][8]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,0,0,
      0,0,2,2
    },{
      /*Corner: 1; split: none*/
      1,1,1,1,
      3,3,1,1
    },{
      /*Corner: 2; split: none*/
      0,0,2,2,
      2,2,2,2
    },{
      /*Corner: 3; split: none*/
      3,3,1,1,
      3,3,3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,0,1,1,
      0,0,2,2
    },{
      /*Corner: 1; split: 2*/
      1,1,1,1,
      3,3,2,2
    },{
      /*Corner: 2; split: 3*/
      0,0,2,2,
      3,3,2,2
    },{
      /*Corner: 3; split: 0*/
      1,1,1,1,
      3,3,2,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,0,0,
      3,3,2,2
    },{
      /*Corner: 1; split: 0*/
      0,0,1,1,
      3,3,1,1
    },{
      /*Corner: 2; split: 1*/
      0,0,1,1,
      2,2,2,2
    },{
      /*Corner: 3; split: 2*/
      3,3,1,1,
      3,3,2,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 2x8 LL band.*/
static const unsigned char OD_MC_SIDXS_28[3][4][16]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,0,0,0,0,0,0,
      0,0,0,0,2,2,2,2
    },{
      /*Corner: 1; split: none*/
      1,1,1,1,1,1,1,1,
      3,3,3,3,1,1,1,1
    },{
      /*Corner: 2; split: none*/
      0,0,0,0,2,2,2,2,
      2,2,2,2,2,2,2,2
    },{
      /*Corner: 3; split: none*/
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,3,3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,0,0,0,1,1,1,1,
      0,0,0,0,2,2,2,2
    },{
      /*Corner: 1; split: 2*/
      1,1,1,1,1,1,1,1,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 2; split: 3*/
      0,0,0,0,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 3; split: 0*/
      1,1,1,1,1,1,1,1,
      3,3,3,3,2,2,2,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,0,0,0,0,0,0,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 1; split: 0*/
      0,0,0,0,1,1,1,1,
      3,3,3,3,1,1,1,1
    },{
      /*Corner: 2; split: 1*/
      0,0,0,0,1,1,1,1,
      2,2,2,2,2,2,2,2
    },{
      /*Corner: 3; split: 2*/
      3,3,3,3,1,1,1,1,
      3,3,3,3,2,2,2,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 4x2 LL band.*/
static const unsigned char OD_MC_SIDXS_42[3][4][8]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,
      0,0,
      0,2,
      0,2
    },{
      /*Corner: 1; split: none*/
      1,1,
      1,1,
      3,1,
      3,1
    },{
      /*Corner: 2; split: none*/
      0,2,
      0,2,
      2,2,
      2,2
    },{
      /*Corner: 3; split: none*/
      3,1,
      3,1,
      3,3,
      3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,1,
      0,1,
      0,2,
      0,2
    },{
      /*Corner: 1; split: 2*/
      1,1,
      1,1,
      3,2,
      3,2
    },{
      /*Corner: 2; split: 3*/
      0,2,
      0,2,
      3,2,
      3,2
    },{
      /*Corner: 3; split: 0*/
      1,1,
      1,1,
      3,2,
      3,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,
      0,0,
      3,2,
      3,2
    },{
      /*Corner: 1; split: 0*/
      0,1,
      0,1,
      3,1,
      3,1
    },{
      /*Corner: 2; split: 1*/
      0,1,
      0,1,
      2,2,
      2,2
    },{
      /*Corner: 3; split: 2*/
      3,1,
      3,1,
      3,2,
      3,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 4x4 LL band.*/
static const unsigned char OD_MC_SIDXS_44[3][4][16]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,0,0,
      0,0,0,0,
      0,0,2,2,
      0,0,2,2
    },{
      /*Corner: 1; split: none*/
      1,1,1,1,
      1,1,1,1,
      3,3,1,1,
      3,3,1,1
    },{
      /*Corner: 2; split: none*/
      0,0,2,2,
      0,0,2,2,
      2,2,2,2,
      2,2,2,2
    },{
      /*Corner: 3; split: none*/
      3,3,1,1,
      3,3,1,1,
      3,3,3,3,
      3,3,3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,0,1,1,
      0,0,1,1,
      0,0,2,2,
      0,0,2,2
    },{
      /*Corner: 1; split: 2*/
      1,1,1,1,
      1,1,1,1,
      3,3,2,2,
      3,3,2,2
    },{
      /*Corner: 2; split: 3*/
      0,0,2,2,
      0,0,2,2,
      3,3,2,2,
      3,3,2,2
    },{
      /*Corner: 3; split: 0*/
      1,1,1,1,
      1,1,1,1,
      3,3,2,2,
      3,3,2,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,0,0,
      0,0,0,0,
      3,3,2,2,
      3,3,2,2
    },{
      /*Corner: 1; split: 0*/
      0,0,1,1,
      0,0,1,1,
      3,3,1,1,
      3,3,1,1
    },{
      /*Corner: 2; split: 1*/
      0,0,1,1,
      0,0,1,1,
      2,2,2,2,
      2,2,2,2
    },{
      /*Corner: 3; split: 2*/
      3,3,1,1,
      3,3,1,1,
      3,3,2,2,
      3,3,2,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for a 4x8 LL band.*/
static const unsigned char OD_MC_SIDXS_48[3][4][32]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,2,
      0,0,0,0,0,2,2,2,
      0,0,0,2,2,2,2,2
    },{
      /*Corner: 1; split: none*/
      1,1,1,1,1,1,1,1,
      3,1,1,1,1,1,1,1,
      3,3,3,1,1,1,1,1,
      3,3,3,3,3,1,1,1
    },{
      /*Corner: 2; split: none*/
      0,0,0,0,0,2,2,2,
      0,0,0,2,2,2,2,2,
      0,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2
    },{
      /*Corner: 3; split: none*/
      3,3,3,1,1,1,1,1,
      3,3,3,3,3,1,1,1,
      3,3,3,3,3,3,3,1,
      3,3,3,3,3,3,3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,0,0,0,1,1,1,1,
      0,0,0,0,0,1,1,1,
      0,0,0,0,2,2,2,2,
      0,0,0,2,2,2,2,2
    },{
      /*Corner: 1; split: 2*/
      1,1,1,1,1,1,1,1,
      3,1,1,1,1,1,1,1,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 2; split: 3*/
      0,0,0,0,0,2,2,2,
      0,0,0,0,2,2,2,2,
      3,3,3,2,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 3; split: 0*/
      1,1,1,1,1,1,1,1,
      3,1,1,1,1,1,1,1,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 1; split: 0*/
      0,0,0,0,1,1,1,1,
      0,0,0,1,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,1,1,1
    },{
      /*Corner: 2; split: 1*/
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2
    },{
      /*Corner: 3; split: 2*/
      3,3,3,1,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,2,2,2,
      3,3,3,3,2,2,2,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for an 8x2 LL band.*/
static const unsigned char OD_MC_SIDXS_82[3][4][16]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,
      0,0,
      0,0,
      0,0,
      0,2,
      0,2,
      0,2,
      0,2
    },{
      /*Corner: 1; split: none*/
      1,1,
      1,1,
      1,1,
      1,1,
      3,1,
      3,1,
      3,1,
      3,1
    },{
      /*Corner: 2; split: none*/
      0,2,
      0,2,
      0,2,
      0,2,
      2,2,
      2,2,
      2,2,
      2,2
    },{
      /*Corner: 3; split: none*/
      3,1,
      3,1,
      3,1,
      3,1,
      3,3,
      3,3,
      3,3,
      3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,1,
      0,1,
      0,1,
      0,1,
      0,2,
      0,2,
      0,2,
      0,2
    },{
      /*Corner: 1; split: 2*/
      1,1,
      1,1,
      1,1,
      1,1,
      3,2,
      3,2,
      3,2,
      3,2
    },{
      /*Corner: 2; split: 3*/
      0,2,
      0,2,
      0,2,
      0,2,
      3,2,
      3,2,
      3,2,
      3,2
    },{
      /*Corner: 3; split: 0*/
      1,1,
      1,1,
      1,1,
      1,1,
      3,2,
      3,2,
      3,2,
      3,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,
      0,0,
      0,0,
      0,0,
      3,2,
      3,2,
      3,2,
      3,2
    },{
      /*Corner: 1; split: 0*/
      0,1,
      0,1,
      0,1,
      0,1,
      3,1,
      3,1,
      3,1,
      3,1
    },{
      /*Corner: 2; split: 1*/
      0,1,
      0,1,
      0,1,
      0,1,
      0,1,
      2,2,
      2,2,
      2,2
    },{
      /*Corner: 3; split: 2*/
      3,1,
      3,1,
      3,1,
      3,1,
      3,2,
      3,2,
      3,2,
      3,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for an 8x4 LL band.*/
static const unsigned char OD_MC_SIDXS_84[3][4][32]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,2,
      0,0,0,2,
      0,0,2,2,
      0,0,2,2,
      0,2,2,2
    },{
      /*Corner: 1; split: none*/
      1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      3,1,1,1,
      3,1,1,1,
      3,3,1,1,
      3,3,1,1,
      3,3,3,1
    },{
      /*Corner: 2; split: none*/
      0,0,0,2,
      0,0,2,2,
      0,0,2,2,
      0,2,2,2,
      0,2,2,2,
      2,2,2,2,
      2,2,2,2,
      2,2,2,2
    },{
      /*Corner: 3; split: none*/
      3,1,1,1,
      3,3,1,1,
      3,3,1,1,
      3,3,3,1,
      3,3,3,1,
      3,3,3,3,
      3,3,3,3,
      3,3,3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,0,1,1,
      0,0,1,1,
      0,0,1,1,
      0,0,1,1,
      0,0,2,2,
      0,0,2,2,
      0,0,2,2,
      0,2,2,2
    },{
      /*Corner: 1; split: 2*/
      1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      3,1,1,1,
      3,3,1,2,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2
    },{
      /*Corner: 2; split: 3*/
      0,0,0,2,
      0,0,2,2,
      0,0,2,2,
      0,0,2,2,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2
    },{
      /*Corner: 3; split: 0*/
      1,1,1,1,
      1,1,1,1,
      1,1,1,1,
      3,1,1,1,
      3,3,1,2,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,0,0,
      0,0,0,0,
      0,0,0,0,
      0,0,0,2,
      3,0,2,2,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2
    },{
      /*Corner: 1; split: 0*/
      0,0,1,1,
      0,0,1,1,
      0,0,1,1,
      0,0,1,1,
      3,3,1,1,
      3,3,1,1,
      3,3,1,1,
      3,3,3,1
    },{
      /*Corner: 2; split: 1*/
      0,0,1,1,
      0,0,1,1,
      0,0,1,1,
      0,0,2,1,
      0,2,2,2,
      2,2,2,2,
      2,2,2,2,
      2,2,2,2
    },{
      /*Corner: 3; split: 2*/
      3,1,1,1,
      3,3,1,1,
      3,3,1,1,
      3,3,1,1,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2,
      3,3,2,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients for an 8x8 LL band.*/
static const unsigned char OD_MC_SIDXS_88[3][4][64]={
  {
    {
      /*Corner: 0; split: none*/
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,2,
      0,0,0,0,0,2,2,2,
      0,0,0,0,2,2,2,2,
      0,0,0,0,2,2,2,2,
      0,0,0,2,2,2,2,2
    },{
      /*Corner: 1; split: none*/
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      3,1,1,1,1,1,1,1,
      3,3,3,1,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,1,1,1
    },{
      /*Corner: 2; split: none*/
      0,0,0,0,0,2,2,2,
      0,0,0,0,2,2,2,2,
      0,0,0,0,2,2,2,2,
      0,0,0,2,2,2,2,2,
      0,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2
    },{
      /*Corner: 3; split: none*/
      3,3,3,1,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,1,1,1,
      3,3,3,3,3,3,3,1,
      3,3,3,3,3,3,3,3,
      3,3,3,3,3,3,3,3,
      3,3,3,3,3,3,3,3
    }
  },{
    {
      /*Corner: 0; split: 1*/
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,0,0,0,0,1,1,1,
      0,0,0,0,2,2,2,2,
      0,0,0,0,2,2,2,2,
      0,0,0,2,2,2,2,2,
      0,0,0,2,2,2,2,2
    },{
      /*Corner: 1; split: 2*/
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      3,3,1,1,1,1,1,1,
      3,3,3,3,1,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 2; split: 3*/
      0,0,0,0,0,2,2,2,
      0,0,0,0,0,2,2,2,
      0,0,0,0,2,2,2,2,
      0,0,0,0,2,2,2,2,
      3,3,3,2,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 3; split: 0*/
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      1,1,1,1,1,1,1,1,
      3,3,1,1,1,1,1,1,
      3,3,3,3,1,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    }
  },{
    {
      /*Corner: 0; split: 3*/
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,2,2,
      3,3,3,0,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    },{
      /*Corner: 1; split: 0*/
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,0,0,1,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,1,1,1,
      3,3,3,3,3,1,1,1
    },{
      /*Corner: 2; split: 1*/
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,0,0,0,1,1,1,1,
      0,0,0,0,2,1,1,1,
      0,0,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2,
      2,2,2,2,2,2,2,2
    },{
      /*Corner: 3; split: 2*/
      3,3,3,1,1,1,1,1,
      3,3,3,1,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,1,1,1,1,
      3,3,3,3,3,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2,
      3,3,3,3,2,2,2,2
    }
  }
};

/*The MV from which to use the high-frequency coefficients, indexed by:
   [log_yblk_sz-2][log_xblk_sz-2][!!s3<<1|!!s1][c][j<<log_xblk_sz-1|i].*/
static const unsigned char *OD_MC_SIDXS[3][3][3][4]={
  {
    {
      {
        OD_MC_SIDXS_22[0][0],OD_MC_SIDXS_22[0][1],
        OD_MC_SIDXS_22[0][2],OD_MC_SIDXS_22[0][3]
      },{
        OD_MC_SIDXS_22[1][0],OD_MC_SIDXS_22[1][1],
        OD_MC_SIDXS_22[1][2],OD_MC_SIDXS_22[1][3]
      },{
        OD_MC_SIDXS_22[2][0],OD_MC_SIDXS_22[2][1],
        OD_MC_SIDXS_22[2][2],OD_MC_SIDXS_22[2][3]
      }
    },{
      {
        OD_MC_SIDXS_24[0][0],OD_MC_SIDXS_24[0][1],
        OD_MC_SIDXS_24[0][2],OD_MC_SIDXS_24[0][3]
      },{
        OD_MC_SIDXS_24[1][0],OD_MC_SIDXS_24[1][1],
        OD_MC_SIDXS_24[1][2],OD_MC_SIDXS_24[1][3]
      },{
        OD_MC_SIDXS_24[2][0],OD_MC_SIDXS_24[2][1],
        OD_MC_SIDXS_24[2][2],OD_MC_SIDXS_24[2][3]
      }
    },{
      {
        OD_MC_SIDXS_28[0][0],OD_MC_SIDXS_28[0][1],
        OD_MC_SIDXS_28[0][2],OD_MC_SIDXS_28[0][3]
      },{
        OD_MC_SIDXS_28[1][0],OD_MC_SIDXS_28[1][1],
        OD_MC_SIDXS_28[1][2],OD_MC_SIDXS_28[1][3]
      },{
        OD_MC_SIDXS_28[2][0],OD_MC_SIDXS_28[2][1],
        OD_MC_SIDXS_28[2][2],OD_MC_SIDXS_28[2][3]
      }
    }
  },{
    {
      {
        OD_MC_SIDXS_42[0][0],OD_MC_SIDXS_42[0][1],
        OD_MC_SIDXS_42[0][2],OD_MC_SIDXS_42[0][3]
      },{
        OD_MC_SIDXS_42[1][0],OD_MC_SIDXS_42[1][1],
        OD_MC_SIDXS_42[1][2],OD_MC_SIDXS_42[1][3]
      },{
        OD_MC_SIDXS_42[2][0],OD_MC_SIDXS_42[2][1],
        OD_MC_SIDXS_42[2][2],OD_MC_SIDXS_42[2][3]
      }
    },{
      {
        OD_MC_SIDXS_44[0][0],OD_MC_SIDXS_44[0][1],
        OD_MC_SIDXS_44[0][2],OD_MC_SIDXS_44[0][3]
      },{
        OD_MC_SIDXS_44[1][0],OD_MC_SIDXS_44[1][1],
        OD_MC_SIDXS_44[1][2],OD_MC_SIDXS_44[1][3]
      },{
        OD_MC_SIDXS_44[2][0],OD_MC_SIDXS_44[2][1],
        OD_MC_SIDXS_44[2][2],OD_MC_SIDXS_44[2][3]
      }
    },{
      {
        OD_MC_SIDXS_48[0][0],OD_MC_SIDXS_48[0][1],
        OD_MC_SIDXS_48[0][2],OD_MC_SIDXS_48[0][3]
      },{
        OD_MC_SIDXS_48[1][0],OD_MC_SIDXS_48[1][1],
        OD_MC_SIDXS_48[1][2],OD_MC_SIDXS_48[1][3]
      },{
        OD_MC_SIDXS_48[2][0],OD_MC_SIDXS_48[2][1],
        OD_MC_SIDXS_48[2][2],OD_MC_SIDXS_48[2][3]
      }
    }
  },{
    {
      {
        OD_MC_SIDXS_82[0][0],OD_MC_SIDXS_82[0][1],
        OD_MC_SIDXS_82[0][2],OD_MC_SIDXS_82[0][3]
      },{
        OD_MC_SIDXS_82[1][0],OD_MC_SIDXS_82[1][1],
        OD_MC_SIDXS_82[1][2],OD_MC_SIDXS_82[1][3]
      },{
        OD_MC_SIDXS_82[2][0],OD_MC_SIDXS_82[2][1],
        OD_MC_SIDXS_82[2][2],OD_MC_SIDXS_82[2][3]
      }
    },{
      {
        OD_MC_SIDXS_84[0][0],OD_MC_SIDXS_84[0][1],
        OD_MC_SIDXS_84[0][2],OD_MC_SIDXS_84[0][3]
      },{
        OD_MC_SIDXS_84[1][0],OD_MC_SIDXS_84[1][1],
        OD_MC_SIDXS_84[1][2],OD_MC_SIDXS_84[1][3]
      },{
        OD_MC_SIDXS_84[2][0],OD_MC_SIDXS_84[2][1],
        OD_MC_SIDXS_84[2][2],OD_MC_SIDXS_84[2][3]
      }
    },{
      {
        OD_MC_SIDXS_88[0][0],OD_MC_SIDXS_88[0][1],
        OD_MC_SIDXS_88[0][2],OD_MC_SIDXS_88[0][3]
      },{
        OD_MC_SIDXS_88[1][0],OD_MC_SIDXS_88[1][1],
        OD_MC_SIDXS_88[1][2],OD_MC_SIDXS_88[1][3]
      },{
        OD_MC_SIDXS_88[2][0],OD_MC_SIDXS_88[2][1],
        OD_MC_SIDXS_88[2][2],OD_MC_SIDXS_88[2][3]
      }
    }
  }
};

/*Perform multiresolution blending with bilinear weights modified for unsplit
   edges.*/
static void od_mc_blend_multi_split8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz){
  const unsigned char *p;
  const unsigned char *sidx0;
  const unsigned char *sidx;
  unsigned char       *dst0;
  unsigned char       *dst;
  ptrdiff_t            o;
  ptrdiff_t            o0;
  int                  ll[4];
  int                  lh;
  int                  hl;
  int                  hh;
  int                  a;
  int                  b;
  int                  c;
  int                  d;
  int                  s[4];
  int                  s0[4];
  int                  dsdi[4];
  int                  dsdj[4];
  int                  ddsdidj[4];
  int                  xblk_sz;
  int                  yblk_sz;
  int                  log_blk_sz2p1;
  int                  i;
  int                  j;
  int                  k;
  /*Perform multiresolution blending.*/
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  o0=0;
  dst0=_dst;
  log_blk_sz2p1=_log_xblk_sz+_log_yblk_sz+1;
  od_mc_setup_s_split(s0,dsdi,dsdj,ddsdidj,_c,_s,
   _log_xblk_sz-1,_log_yblk_sz-1);
  sidx0=OD_MC_SIDXS[_log_yblk_sz-2][_log_xblk_sz-2][_s][_c];
  for(k=0;k<4;k++)s[k]=s0[k];
  for(j=1;j<yblk_sz;j+=2){
    o=o0;
    dst=dst0;
    sidx=sidx0;
    /*Upper-left quadrant.*/
    for(i=1;i<xblk_sz;i+=2){
      k=*sidx++;
      p=_src[k]+o;
      /*Forward Haar wavelet.*/
      ll[k]=p[0]+p[1];
      lh=p[0]-p[1];
      hl=(p+xblk_sz)[0]+(p+xblk_sz)[1];
      hh=(p+xblk_sz)[0]-(p+xblk_sz)[1];
      c=ll[k]-hl;
      ll[k]+=hl;
      hl=c;
      /*No need to finish the transform; we'd just invert it later.*/
      p=_src[k+1&3]+o;
      ll[k+1&3]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[k+2&3]+o;
      ll[k+2&3]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      p=_src[k+3&3]+o;
      ll[k+3&3]=p[0]+p[1]+(p+xblk_sz)[0]+(p+xblk_sz)[1];
      /*LL blending.*/
      a=(int)((ogg_int32_t)ll[0]*s[0]+
       (ogg_int32_t)ll[1]*s[1]+(ogg_int32_t)ll[2]*s[2]+
       (ogg_int32_t)ll[3]*s[3]>>log_blk_sz2p1);
      /*Inverse Haar wavelet.*/
      c=a-hl+1>>1;
      a=a+hl+1>>1;
      d=c-hh+1>>1;
      c=c+hh+1>>1;
      b=a-lh+1>>1;
      a=a+lh+1>>1;
      dst[0]=OD_CLAMP255(a);
      dst[1]=OD_CLAMP255(b);
      (dst+_dystride)[0]=OD_CLAMP255(c);
      (dst+_dystride)[1]=OD_CLAMP255(d);
      o+=2;
      dst+=2;
      for(k=0;k<4;k++)s[k]+=dsdi[k];
    }
    o0+=xblk_sz<<1;
    dst0+=_dystride<<1;
    sidx0+=xblk_sz>>1;
    for(k=0;k<4;k++){
      s0[k]+=dsdj[k];
      s[k]=s0[k];
      dsdi[k]+=ddsdidj[k];
    }
  }
}
#else
/*Sets up a second set of image pointers based on the given split state to
   properly shift weight from one image to another.*/
static void od_mc_setup_split_ptrs(const unsigned char *_drc[4],
 const unsigned char *_src[4],int _c,int _s){
  int j;
  int k;
  _drc[_c]=_src[_c];
  j=_c+(_s&1)&3;
  k=_c+1&3;
  _drc[k]=_src[j];
  j=_c+(_s&2)+((_s&2)>>1)&3;
  k=_c+3&3;
  _drc[k]=_src[j];
  k=_c^2;
  _drc[k]=_src[k];
}

/*Perform multiresolution bilinear blending.*/
static void od_mc_blend_multi_split8(unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz){
  unsigned char        src_ll[4][8][8];
  int                  dst_ll[8][8];
  const unsigned char *drc[4];
  const unsigned char *p;
  const unsigned char *q;
  unsigned char       *dst;
  ptrdiff_t            o;
  int                  a;
  int                  b;
  int                  c;
  int                  d;
  int                  e;
  int                  f;
  int                  g;
  int                  h;
  int                  log_blk_sz2;
  int                  xblk_sz;
  int                  yblk_sz;
  int                  xblk_sz_2;
  int                  yblk_sz_2;
  int                  xblk_sz_4;
  int                  yblk_sz_4;
  int                  round;
  int                  i;
  int                  i2;
  int                  j;
  int                  j2;
  int                  k;
  xblk_sz=1<<_log_xblk_sz;
  yblk_sz=1<<_log_yblk_sz;
  od_mc_setup_split_ptrs(drc,_src,_c,_s);
  /*Perform multiresolution blending.*/
  xblk_sz_2=xblk_sz>>1;
  yblk_sz_2=yblk_sz>>1;
  log_blk_sz2=_log_xblk_sz+_log_yblk_sz;
  round=1<<log_blk_sz2-1;
  /*Compute the low-pass band for each src block.*/
  for(k=0;k<4;k++){
    unsigned lh[4][8];
    p=_src[k];
    q=drc[k];
    src_ll[k][0][0]=p[0]+q[0]+1>>1;
    for(i=1;i<xblk_sz_2;i++){
      i2=i<<1;
      src_ll[k][0][i]=(unsigned char)
       (p[i2-1]+q[i2-1]+2*(p[i2]+q[i2])+p[i2+1]+q[i2+1]+4>>3);
    }
    p+=xblk_sz;
    q+=xblk_sz;
    lh[1][0]=p[0]+q[0]<<2;
    for(i=1;i<xblk_sz_2;i++){
      i2=i<<1;
      lh[1][i]=p[i2-1]+q[i2-1]+2*(p[i2]+q[i2])+p[i2+1]+q[i2+1];
    }
    p+=xblk_sz;
    q+=xblk_sz;
    for(j=1;j<yblk_sz_2;j++){
      j2=j<<1&3;
      lh[j2][0]=p[0]+q[0]<<2;
      for(i=1;i<xblk_sz_2;i++){
        i2=i<<1;
        lh[j2][i]=p[i2-1]+q[i2-1]+2*(p[i2]+q[i2])+p[i2+1]+q[i2+1];
      }
      p+=xblk_sz;
      q+=xblk_sz;
      lh[j2+1][0]=p[0]<<2;
      for(i=1;i<xblk_sz_2;i++){
        i2=i<<1;
        lh[j2+1][i]=p[i2-1]+q[i2-1]+2*(p[i2]+q[i2])+p[i2+1]+q[i2+1];
      }
      p+=xblk_sz;
      q+=xblk_sz;
      for(i=0;i<xblk_sz_2;i++){
        src_ll[k][j][i]=(unsigned char)
         (lh[j2-1&3][i]+2*lh[j2][i]+lh[j2+1][i]+16>>5);
      }
    }
  }
  /*Blend the low-pass bands.*/
  for(j=0;j<xblk_sz_2;j++){
    for(i=0;i<xblk_sz_2;i++){
      a=(src_ll[0][j][i]<<_log_xblk_sz-1)+(src_ll[1][j][i]-src_ll[0][j][i])*i;
      b=(src_ll[3][j][i]<<_log_xblk_sz-1)+(src_ll[2][j][i]-src_ll[3][j][i])*i;
      dst_ll[j][i]=(a<<_log_yblk_sz-1)+(b-a)*j;
    }
  }
  /*Perform the high-pass filtering for each quadrant.*/
  xblk_sz_4=xblk_sz>>2;
  yblk_sz_4=yblk_sz>>2;
  o=0;
  dst=_dst;
  for(j=0;j<yblk_sz_4;j++){
    /*Upper-left quadrant.*/
    for(i=0;i<xblk_sz_4;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[0][j][i]<<log_blk_sz2;
      f=src_ll[0][j][i]+src_ll[0][j][i+1]<<log_blk_sz2-1;
      g=src_ll[0][j][i]+src_ll[0][j+1][i]<<log_blk_sz2-1;
      h=src_ll[0][j][i]+src_ll[0][j][i+1]+
       src_ll[0][j+1][i]+src_ll[0][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(((_src[0]+o)[i2]+(drc[0]+o)[i2]<<log_blk_sz2-1)+
       a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(((_src[0]+o)[i2+1]+(drc[0]+o)[i2+1]<<log_blk_sz2-1)+
       b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[0]+o+xblk_sz)[i2]+(drc[0]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
       c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[0]+o+xblk_sz)[i2+1]+(drc[0]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
       d-h+round>>log_blk_sz2);
    }
    /*Upper-right quadrant.*/
    for(;i<xblk_sz_2-1;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[1][j][i]<<log_blk_sz2;
      f=src_ll[1][j][i]+src_ll[1][j][i+1]<<log_blk_sz2-1;
      g=src_ll[1][j][i]+src_ll[1][j+1][i]<<log_blk_sz2-1;
      h=src_ll[1][j][i]+src_ll[1][j][i+1]+
       src_ll[1][j+1][i]+src_ll[1][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(((_src[1]+o)[i2]+(drc[1]+o)[i2]<<log_blk_sz2-1)+
       a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(((_src[1]+o)[i2+1]+(drc[1]+o)[i2+1]<<log_blk_sz2-1)+
       b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[1]+o+xblk_sz)[i2]+(drc[1]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
       c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[1]+o+xblk_sz)[i2+1]+(drc[1]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
       d-h+round>>log_blk_sz2);
    }
    /*Upper-right quadrant, last column.*/
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=3*dst_ll[j][i]-dst_ll[j][i-1]<<1;
    c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j+1][i])-(dst_ll[j][i-1]+dst_ll[j+1][i-1]);
    e=src_ll[1][j][i]<<log_blk_sz2;
    f=3*src_ll[1][j][i]-src_ll[1][j][i-1]<<log_blk_sz2-1;
    g=src_ll[1][j][i]+src_ll[1][j+1][i]<<log_blk_sz2-1;
    h=3*(src_ll[1][j][i]+src_ll[1][j+1][i])-
     (src_ll[1][j][i-1]+src_ll[1][j+1][i-1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(((_src[1]+o)[i2]+(drc[1]+o)[i2]<<log_blk_sz2-1)+
     a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(((_src[1]+o)[i2+1]+(drc[1]+o)[i2+1]<<log_blk_sz2-1)+
     b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[1]+o+xblk_sz)[i2]+(drc[1]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
     c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[1]+o+xblk_sz)[i2+1]+(drc[1]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
     d-h+round>>log_blk_sz2);
    o+=xblk_sz<<1;
    dst+=_dystride<<1;
  }
  for(;j<yblk_sz_2-1;j++){
    /*Lower-left quadrant.*/
    for(i=0;i<xblk_sz_4;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[3][j][i]<<log_blk_sz2;
      f=src_ll[3][j][i]+src_ll[3][j][i+1]<<log_blk_sz2-1;
      g=src_ll[3][j][i]+src_ll[3][j+1][i]<<log_blk_sz2-1;
      h=src_ll[3][j][i]+src_ll[3][j][i+1]+
       src_ll[3][j+1][i]+src_ll[3][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(((_src[3]+o)[i2]+(drc[3]+o)[i2]<<log_blk_sz2-1)+
       a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(((_src[3]+o)[i2+1]+(drc[3]+o)[i2+1]<<log_blk_sz2-1)+
       b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[3]+o+xblk_sz)[i2]+(drc[3]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
       c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[3]+o+xblk_sz)[i2+1]+(drc[3]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
       d-h+round>>log_blk_sz2);
    }
    /*Lower-right quadrant.*/
    for(;i<xblk_sz_2-1;i++){
      i2=i<<1;
      a=dst_ll[j][i]<<2;
      b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
      c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
      d=dst_ll[j][i]+dst_ll[j][i+1]+dst_ll[j+1][i]+dst_ll[j+1][i+1];
      e=src_ll[2][j][i]<<log_blk_sz2;
      f=src_ll[2][j][i]+src_ll[2][j][i+1]<<log_blk_sz2-1;
      g=src_ll[2][j][i]+src_ll[2][j+1][i]<<log_blk_sz2-1;
      h=src_ll[2][j][i]+src_ll[2][j][i+1]+
       src_ll[2][j+1][i]+src_ll[2][j+1][i+1]<<log_blk_sz2-2;
      dst[i2]=OD_CLAMP255(((_src[2]+o)[i2]+(drc[2]+o)[i2]<<log_blk_sz2-1)+
       a-e+round>>log_blk_sz2);
      dst[i2+1]=OD_CLAMP255(((_src[2]+o)[i2+1]+(drc[2]+o)[i2+1]<<log_blk_sz2-1)+
       b-f+round>>log_blk_sz2);
      (dst+_dystride)[i2]=OD_CLAMP255(
       ((_src[2]+o+xblk_sz)[i2]+(drc[2]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
       c-g+round>>log_blk_sz2);
      (dst+_dystride)[i2+1]=OD_CLAMP255(
       ((_src[2]+o+xblk_sz)[i2+1]+(drc[2]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
       d-h+round>>log_blk_sz2);
    }
    /*Lower-right quadrant, last column.*/
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=3*dst_ll[j][i]-dst_ll[j][i-1]<<1;
    c=dst_ll[j][i]+dst_ll[j+1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j+1][i])-(dst_ll[j][i-1]+dst_ll[j+1][i-1]);
    e=src_ll[2][j][i]<<log_blk_sz2;
    f=3*src_ll[2][j][i]-src_ll[2][j][i-1]<<log_blk_sz2-1;
    g=src_ll[2][j][i]+src_ll[2][j+1][i]<<log_blk_sz2-1;
    h=3*(src_ll[2][j][i]+src_ll[2][j+1][i])-
     (src_ll[2][j][i-1]+src_ll[2][j+1][i-1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(((_src[2]+o)[i2]+(drc[2]+o)[i2+1]<<log_blk_sz2-1)+
     a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(((_src[2]+o)[i2+1]+(drc[2]+o)[i2]<<log_blk_sz2-1)+
     b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2]+(drc[2]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
     c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2+1]+(drc[2]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
     d-h+round>>log_blk_sz2);
    o+=xblk_sz<<1;
    dst+=_dystride<<1;
  }
  /*Lower-left quadrant, last row.*/
  for(i=0;i<xblk_sz_4;i++){
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
    c=3*dst_ll[j][i]-dst_ll[j-1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j][i+1])-(dst_ll[j-1][i]+dst_ll[j-1][i+1]);
    e=src_ll[3][j][i]<<log_blk_sz2;
    f=src_ll[3][j][i]+src_ll[3][j][i+1]<<log_blk_sz2-1;
    g=3*src_ll[3][j][i]-src_ll[3][j-1][i]<<log_blk_sz2-1;
    h=3*(src_ll[3][j][i]+src_ll[3][j][i+1])-
     (src_ll[3][j-1][i]+src_ll[3][j-1][i+1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(((_src[3]+o)[i2]+(drc[3]+o)[i2]<<log_blk_sz2-1)+
     a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(((_src[3]+o)[i2+1]+(drc[3]+o)[i2+1]<<log_blk_sz2-1)+
     b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[3]+o+xblk_sz)[i2]+(drc[3]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
     c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[3]+o+xblk_sz)[i2+1]+(drc[3]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
     d-h+round>>log_blk_sz2);
  }
  /*Lower-right quadrant, last row.*/
  for(;i<xblk_sz_2-1;i++){
    i2=i<<1;
    a=dst_ll[j][i]<<2;
    b=dst_ll[j][i]+dst_ll[j][i+1]<<1;
    c=3*dst_ll[j][i]-dst_ll[j-1][i]<<1;
    d=3*(dst_ll[j][i]+dst_ll[j][i+1])-(dst_ll[j-1][i]+dst_ll[j-1][i+1]);
    e=src_ll[2][j][i]<<log_blk_sz2;
    f=src_ll[2][j][i]+src_ll[2][j][i+1]<<log_blk_sz2-1;
    g=3*src_ll[2][j][i]-src_ll[2][j-1][i]<<log_blk_sz2-1;
    h=3*(src_ll[2][j][i]+src_ll[2][j][i+1])-
     (src_ll[2][j-1][i]+src_ll[2][j-1][i+1])<<log_blk_sz2-2;
    dst[i2]=OD_CLAMP255(((_src[2]+o)[i2]+(drc[2]+o)[i2]<<log_blk_sz2-1)+
     a-e+round>>log_blk_sz2);
    dst[i2+1]=OD_CLAMP255(((_src[2]+o)[i2+1]+(drc[2]+o)[i2+1]<<log_blk_sz2-1)+
     b-f+round>>log_blk_sz2);
    (dst+_dystride)[i2]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2]+(drc[2]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
     c-g+round>>log_blk_sz2);
    (dst+_dystride)[i2+1]=OD_CLAMP255(
     ((_src[2]+o+xblk_sz)[i2+1]+(drc[2]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
     d-h+round>>log_blk_sz2);
  }
  /*Lower-right quadrant, last row and column.*/
  i2=i<<1;
  a=dst_ll[j][i]<<2;
  b=3*dst_ll[j][i]-dst_ll[j][i-1]<<1;
  c=3*dst_ll[j][i]-dst_ll[j-1][i]<<1;
  d=9*dst_ll[j][i]-3*(dst_ll[j-1][i]+dst_ll[j][i-1])+dst_ll[j-1][i-1];
  e=src_ll[2][j][i]<<log_blk_sz2;
  f=3*src_ll[2][j][i]-src_ll[2][j][i-1]<<log_blk_sz2-1;
  g=3*src_ll[2][j][i]-src_ll[2][j-1][i]<<log_blk_sz2-1;
  h=9*src_ll[2][j][i]-3*(src_ll[2][j][i-1]+src_ll[2][j-1][i])+
   src_ll[2][j-1][i-1]<<log_blk_sz2-2;
  dst[i2]=OD_CLAMP255(((_src[2]+o)[i2]+(drc[2]+o)[i2]<<log_blk_sz2-1)+
   a-e+round>>log_blk_sz2);
  dst[i2+1]=OD_CLAMP255(((_src[2]+o)[i2+1]+(drc[2]+o)[i2+1]<<log_blk_sz2-1)+
   b-f+round>>log_blk_sz2);
  (dst+_dystride)[i2]=OD_CLAMP255(
   ((_src[2]+o+xblk_sz)[i2]+(drc[2]+o+xblk_sz)[i2]<<log_blk_sz2-1)+
   c-g+round>>log_blk_sz2);
  (dst+_dystride)[i2+1]=OD_CLAMP255(
   ((_src[2]+o+xblk_sz)[i2+1]+(drc[2]+o+xblk_sz)[i2+1]<<log_blk_sz2-1)+
   d-h+round>>log_blk_sz2);
}
#endif

static void od_mc_blend8(od_state *_state,unsigned char *_dst,int _dystride,
 const unsigned char *_src[4],int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz){
  if(0&&_log_xblk_sz>1&&_log_yblk_sz>1){
    /*Perform multiresolution blending.*/
    if(_s==3)od_mc_blend_multi8(_dst,_dystride,_src,_log_xblk_sz,_log_yblk_sz);
    else{
      od_mc_blend_multi_split8(_dst,_dystride,_src,
       _c,_s,_log_xblk_sz,_log_yblk_sz);
    }
  }
  else{
    /*The block is too small; perform normal blending.*/
    if(_s==3){
      od_mc_blend_full8(_state,_dst,_dystride,_src,
       _log_xblk_sz,_log_yblk_sz);
    }
    else{
      od_mc_blend_full_split8(_state,_dst,_dystride,_src,_c,_s,
       _log_xblk_sz,_log_yblk_sz);
    }
  }
}

void od_mc_predict8(od_state *_state,unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],int _interp_type,int _c,int _s,
 int _log_xblk_sz,int _log_yblk_sz){
  const unsigned char *pred[4];
  unsigned char        __attribute__((aligned(16))) buf[4][16*16];
  int                  r;
  r=0;
  switch(_interp_type){
    case OD_MC_INTERP_VVVV:{
      od_mc_predict1imv8(_state,_dst,_dystride,_src,_systride,
       _mvx,_mvy,MIDXS[0]/*0,1,2,3*/,0,_log_xblk_sz,_log_yblk_sz);
    }break;
    case OD_MC_INTERP_VVVB:r++;
    case OD_MC_INTERP_VVBV:r++;
    case OD_MC_INTERP_VBVV:r++;
    case OD_MC_INTERP_BVVV:{
      od_mc_predict1imv8(_state,buf[0],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[1]/*0,0,2,3*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[0+r&3]=buf[0];
      od_mc_predict1imv8(_state,buf[1],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[2]/*1,1,2,3*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[1+r&3]=buf[1];
      od_mc_predict1imv8(_state,buf[2],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[0]/*0,1,2,3*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[2+r&3]=buf[2];
      pred[3+r&3]=buf[2];
      od_mc_blend8(_state,_dst,_dystride,pred,
       _c,_s|((_c-r&3)!=0)|((_c-r&3)!=1)<<1,_log_xblk_sz,_log_yblk_sz);
    }break;
    case OD_MC_INTERP_VBVB:r++;
    case OD_MC_INTERP_BVBV:{
      od_mc_predict1imv8(_state,buf[0],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[3]/*0,0,3,3*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[3+r&3]=buf[0];
      pred[0+r&3]=buf[0];
      od_mc_predict1imv8(_state,buf[1],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[4]/*1,1,2,2*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[1+r&3]=buf[1];
      pred[2+r&3]=buf[1];
      od_mc_blend8(_state,_dst,_dystride,pred,
       _c,_s|(_c-r&1)|!(_c-r&1)<<1,_log_xblk_sz,_log_yblk_sz);
    }break;
    case OD_MC_INTERP_VBBV:r++;
    case OD_MC_INTERP_BBVV:r++;
    case OD_MC_INTERP_BVVB:r++;
    case OD_MC_INTERP_VVBB:{
      od_mc_predict1imv8(_state,buf[0],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[5]/*0,1,0,0*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[0+r&3]=buf[0];
      od_mc_predict1imv8(_state,buf[1],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[0]/*0,1,2,3*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[1+r&3]=buf[1];
      od_mc_predict1imv8(_state,buf[2],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[6]/*2,1,2,2*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[2+r&3]=buf[2];
      od_mc_predict1fmv8(_state,buf[3],_src,_systride,
       _mvx[3+r&3],_mvy[3+r&3],_log_xblk_sz,_log_yblk_sz);
      pred[3+r&3]=buf[3];
      od_mc_blend8(_state,_dst,_dystride,pred,
       _c,_s|(_c+2-r&2)>>1|(_c+1-r&2),_log_xblk_sz,_log_yblk_sz);
    }break;
    case OD_MC_INTERP_BBBV:r++;
    case OD_MC_INTERP_BBVB:r++;
    case OD_MC_INTERP_BVBB:r++;
    case OD_MC_INTERP_VBBB:{
      od_mc_predict1imv8(_state,buf[0],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[5]/*0,1,0,0*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[0+r&3]=buf[0];
      od_mc_predict1imv8(_state,buf[1],1<<_log_xblk_sz,_src,_systride,
       _mvx,_mvy,MIDXS[7]/*0,1,1,1*/,r,_log_xblk_sz,_log_yblk_sz);
      pred[1+r&3]=buf[1];
      od_mc_predict1fmv8(_state,buf[2],_src,_systride,
       _mvx[2+r&3],_mvy[2+r&3],_log_xblk_sz,_log_yblk_sz);
      pred[2+r&3]=buf[2];
      if(_mvx[3+r&3]==_mvx[2+r&3]&&_mvy[3+r&3]==_mvy[2+r&3]){
        pred[3+r&3]=pred[2+r&3];
      }
      else{
        od_mc_predict1fmv8(_state,buf[3],_src,_systride,
         _mvx[3+r&3],_mvy[3+r&3],_log_xblk_sz,_log_yblk_sz);
        pred[3+r&3]=buf[3];
      }
      od_mc_blend8(_state,_dst,_dystride,pred,
       _c,_s|((_c-r&3)==0)|((_c-r&3)==1)<<1,_log_xblk_sz,_log_yblk_sz);
    }break;
    case OD_MC_INTERP_BBBB:{
      od_mc_predict1fmv8(_state,buf[0],_src,_systride,
       _mvx[0],_mvy[0],_log_xblk_sz,_log_yblk_sz);
      pred[0]=buf[0];
      if(_mvx[1]==_mvx[0]&&_mvy[1]==_mvy[1])pred[1]=pred[0];
      else{
        od_mc_predict1fmv8(_state,buf[1],_src,_systride,
         _mvx[1],_mvy[1],_log_xblk_sz,_log_yblk_sz);
        pred[1]=buf[1];
      }
      if(_mvx[2]==_mvx[0]&&_mvy[2]==_mvy[0])pred[2]=pred[0];
      else if(_mvx[2]==_mvx[1]&&_mvy[2]==_mvy[1])pred[2]=pred[1];
      else{
        od_mc_predict1fmv8(_state,buf[2],_src,_systride,
         _mvx[2],_mvy[2],_log_xblk_sz,_log_yblk_sz);
        pred[2]=buf[2];
      }
      if(_mvx[3]==_mvx[0]&&_mvy[3]==_mvy[0])pred[3]=pred[0];
      else if(_mvx[3]==_mvx[1]&&_mvy[3]==_mvy[1])pred[3]=pred[1];
      else if(_mvx[3]==_mvx[2]&&_mvy[3]==_mvy[2])pred[3]=pred[2];
      else{
        od_mc_predict1fmv8(_state,buf[3],_src,_systride,
         _mvx[3],_mvy[3],_log_xblk_sz,_log_yblk_sz);
        pred[3]=buf[3];
      }
      od_mc_blend8(_state,_dst,_dystride,pred,_c,_s,_log_xblk_sz,_log_yblk_sz);
    }
  }
}

#if 0
#include <stdio.h>

static unsigned char mask[4][4]={
  {0, 8, 2,10},
  {12,4,14, 6},
  {3, 11,1, 9},
  {15,7,13, 5}
};

static unsigned char img[16*7][16*7];

static ogg_int32_t mvs[4][4][2];
static ogg_int32_t mvs2[4][4][2];

static int edge_types[4][4]={
  {0x0,0x1,0x2,0x8},
  {0x4,0x6,0xA,0xC},
  {0x5,0x7,0xE,0xD},
  {0x3,0xF,0xB,0x9}
};

static int edge_types2[8][8]={
  {0x0,0x4,0x1,0x5,0x0,0x2,0xA,0x8},
  {0x2,0x9,0x0,0x1,0x2,0xA,0x8,0x0},
  {0x6,0x8,0x2,0xA,0xC,0x6,0xC,0x4},
  {0x5,0x4,0x6,0xE,0x9,0x3,0xF,0xD},
  {0x1,0x5,0x3,0xF,0xA,0xE,0xB,0x9},
  {0x6,0xD,0x4,0x7,0xE,0xF,0xE,0xC},
  {0x5,0x7,0xB,0xF,0xD,0x7,0xB,0xD},
  {0x3,0xB,0xC,0x7,0x9,0x3,0x8,0x1}
};

static void fill_mvs(int _log_blk_sz){
  int i;
  int j;
  for(j=0;j<4;j++){
    for(i=0;i<4;i++){
      mvs[j][i][0]=(mask[j+0&3][i+0&3]-8)<<_log_blk_sz+15^
       mask[j+2&3][i+2&3]<<_log_blk_sz+12;
      mvs[j][i][1]=(mask[j+1&3][i+1&3]-8)<<_log_blk_sz+15^
       mask[j+3&3][i+3&3]<<_log_blk_sz+12;
      mvs2[j][i][0]=(mask[j+3&3][i+3&3]-8)<<_log_blk_sz+15^
       mask[j+1&3][i+1&3]<<_log_blk_sz+12;
      mvs2[j][i][1]=(mask[j+2&3][i+2&3]-8)<<_log_blk_sz+15^
       mask[j+0&3][i+0&3]<<_log_blk_sz+12;
    }
  }
}

static void fill_img(int _log_blk_sz){
  int i;
  int j;
  for(j=0;j<7;j++){
    for(i=0;i<7;i++){
      int c;
      int x;
      int y;
      c=mask[j&3][i&3];
      c=c<<4|15;
      for(y=j<<_log_blk_sz;y<(j+1)<<_log_blk_sz;y++){
        for(x=i<<_log_blk_sz;x<(i+1)<<_log_blk_sz;x++){
          img[y][x]=c;
        }
      }
    }
  }
  for(j=0;j<6<<_log_blk_sz;j++){
    for(i=0;i<6<<_log_blk_sz;i++){
      printf("%2X%c",img[j][i],i+1<6<<_log_blk_sz?' ':'\n');
    }
  }
}

int main(void){
  int log_blk_sz;
  for(log_blk_sz=2;log_blk_sz<=4;log_blk_sz++){
    int blk_sz;
    int i;
    int j;
    blk_sz=1<<log_blk_sz;
    fill_img(log_blk_sz);
    fill_mvs(log_blk_sz);
    for(j=0;j<4;j++){
      for(i=0;i<4;i++){
        ogg_int32_t   mvx[4];
        ogg_int32_t   mvy[4];
        unsigned char dst[17][17];
        unsigned char dst2[4][9][9];
        unsigned      mismatch[4][9][9];
        int           etype;
        int           x;
        int           y;
        int           c;
        mvx[0]=mvs[j][i][0];
        mvy[0]=mvs[j][i][1];
        mvx[1]=mvs[j][i+1&3][0];
        mvy[1]=mvs[j][i+1&3][1];
        mvx[2]=mvs[j+1&3][i+1&3][0];
        mvy[2]=mvs[j+1&3][i+1&3][1];
        mvx[3]=mvs[j+1&3][i][0];
        mvy[3]=mvs[j+1&3][i][1];
        etype=edge_types[j][i];
        printf("Block (%i,%i): size %i, interpolation type: %c%c%c%c (0x%X)\n",
         i,j,1<<log_blk_sz,
         etype&1?'V':'B',etype&2?'V':'B',
         etype&4?'V':'B',etype&8?'V':'B',etype);
        printf("<%8.4lf,%8.4lf> <%8.4lf,%8.4lf>\n",
         mvx[0]/(double)0x40000,mvy[0]/(double)0x40000,
         mvx[1]/(double)0x40000,mvy[1]/(double)0x40000);
        printf("<%8.4lf,%8.4lf> <%8.4lf,%8.4lf>\n",
         mvx[3]/(double)0x40000,mvy[3]/(double)0x40000,
         mvx[2]/(double)0x40000,mvy[2]/(double)0x40000);
        od_mc_predict8(dst[0],sizeof(dst[0]),
         img[j+1<<log_blk_sz]+(i+1<<log_blk_sz),sizeof(img[0]),mvx,mvy,
         etype,0,0,0,log_blk_sz,log_blk_sz);
        for(y=0;y<blk_sz;y++){
          for(x=0;x<blk_sz;x++){
            printf("%2X%c",dst[y][x],x+1<blk_sz?' ':'\n');
          }
        }
        printf("\n");
        for(c=0;c<4;c++){
          int s1;
          int s3;
          mvx[0]=mvs[j][i][0];
          mvy[0]=mvs[j][i][1];
          mvx[1]=mvs[j][i+1&3][0];
          mvy[1]=mvs[j][i+1&3][1];
          mvx[2]=mvs[j+1&3][i+1&3][0];
          mvy[2]=mvs[j+1&3][i+1&3][1];
          mvx[3]=mvs[j+1&3][i][0];
          mvy[3]=mvs[j+1&3][i][1];
          mvx[c+2&3]=mvs2[j][i][0];
          mvy[c+2&3]=mvs2[j][i][1];
          etype=edge_types2[j<<1|(c>>1)][i<<1|((c+1&3)>>1)];
          if(1||!(c&1)){
            etype&=~(1<<(c+1&3));
            etype|=(etype|etype<<4)>>3&1<<(c+1&3);
            s1=0;
          }
          else s1=1;
          if(s1||(etype>>c&1)){
            mvx[c+1&3]=mvx[c]+mvx[c+1&3]>>1;
            mvy[c+1&3]=mvy[c]+mvy[c+1&3]>>1;
          }
          if(1||(c&1)){
            etype&=~(1<<(c+2&3));
            etype|=(etype|etype<<4)>>1&1<<(c+2&3);
            s3=0;
          }
          else s3=1;
          if(s3||(etype<<(-c&3)&8)){
            mvx[c+3&3]=mvx[c]+mvx[c+3&3]>>1;
            mvy[c+3&3]=mvy[c]+mvy[c+3&3]>>1;
          }
          printf("Block (%i.%i,%i.%i): size %i, "
           "interpolation type: %c%c%c%c (0x%X)\n",
           i,((c+1&3)>>1)*5,j,(c>>1)*5,1<<log_blk_sz-1,
           etype&1?'V':'B',etype&2?'V':'B',
           etype&4?'V':'B',etype&8?'V':'B',etype);
          printf("<%9.5lf,%9.5lf> <%9.5lf,%9.5lf>\n",
           mvx[0]/(double)0x40000,mvy[0]/(double)0x40000,
           mvx[1]/(double)0x40000,mvy[1]/(double)0x40000);
          printf("<%9.5lf,%9.5lf> <%9.5lf,%9.5lf>\n",
           mvx[3]/(double)0x40000,mvy[3]/(double)0x40000,
           mvx[2]/(double)0x40000,mvy[2]/(double)0x40000);
          od_mc_predict8(dst2[c][0],sizeof(dst2[c][0]),
           img[(j+1<<1|(c>>1))<<log_blk_sz-1]+
           ((i+1<<1|((c+1&3)>>1))<<log_blk_sz-1),
           sizeof(img[0]),mvx,mvy,etype,c,s1,s3,log_blk_sz-1,log_blk_sz-1);
          memset(mismatch[c][0],0,sizeof(mismatch[c]));
          switch(c){
            case 0:{
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst[0][x])mismatch[c][0][x]++;
              }
              for(y=1;y<blk_sz>>1;y++){
                if(dst2[c][y][0]!=dst[y][0])mismatch[c][y][0]++;
              }
            }break;
            case 1:{
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst[0][x+(blk_sz>>1)])mismatch[c][0][x]++;
              }
              for(y=1;y<blk_sz>>1;y++){
                if(dst2[c][y][blk_sz>>1]!=dst[y][blk_sz]){
                  mismatch[c][y][blk_sz>>1]++;
                }
              }
              for(y=0;y<blk_sz>>1;y++){
                if(dst2[c][y][0]!=dst2[0][y][blk_sz>>1])mismatch[c][y][0]++;
              }
            }break;
            case 2:{
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst2[1][blk_sz>>1][x])mismatch[c][0][x]++;
              }
              for(y=0;y<blk_sz>>1;y++){
                if(dst2[c][y][blk_sz>>1]!=dst[y+(blk_sz>>1)][blk_sz]){
                  mismatch[c][y][blk_sz>>1]++;
                }
              }
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][blk_sz>>1][x]!=dst[blk_sz][x+(blk_sz>>1)]){
                  mismatch[c][blk_sz>>1][x]++;
                }
              }
            }break;
            case 3:{
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst2[0][blk_sz>>1][x])mismatch[c][0][x]++;
              }
              for(y=1;y<blk_sz>>1;y++){
                if(dst2[c][y][blk_sz>>1]!=dst2[2][y][0]){
                  mismatch[c][y][blk_sz>>1]++;
                }
              }
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][blk_sz>>1][x]!=dst[blk_sz][x]){
                  mismatch[c][blk_sz>>1][x]++;
                }
              }
              for(y=0;y<blk_sz>>1;y++){
                if(dst2[c][y][0]!=dst[y+(blk_sz>>1)][0])mismatch[c][y][0]++;
              }
            }break;
          }
          for(y=0;y<blk_sz>>1;y++){
            for(x=0;x<blk_sz>>1;x++){
              printf("%c%2X",mismatch[c][y][x]?'!':' ',dst2[c][y][x]);
            }
            printf("\n");
          }
          printf("\n");
        }
      }
    }
  }
  return 0;
}

#endif
