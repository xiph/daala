#include <stddef.h>
#include "internal.h"

/*Interpolation along an edge by interpolating vectors.*/
#define OD_MC_INTERP_VECTOR (0)
/*Interpolation along an edge by blending samples from different vectors.*/
#define OD_MC_INTERP_BLEND  (1)

#define _MCInterpType(_a,_b,_c,_d) ((_a)|((_b)<<1)|((_c)<<2)|((_d)<<3))

#define OD_MC_INTERP_VVVV (0)
#define OD_MC_INTERP_BVVV (1)
#define OD_MC_INTERP_VBVV (2)
#define OD_MC_INTERP_BBVV (3)
#define OD_MC_INTERP_VVBV (4)
#define OD_MC_INTERP_BVBV (5)
#define OD_MC_INTERP_VBBV (6)
#define OD_MC_INTERP_BBBV (7)
#define OD_MC_INTERP_VVVB (8)
#define OD_MC_INTERP_BVVB (9)
#define OD_MC_INTERP_VBVB (10)
#define OD_MC_INTERP_BBVB (11)
#define OD_MC_INTERP_VVBB (12)
#define OD_MC_INTERP_BVBB (13)
#define OD_MC_INTERP_VBBB (14)
#define OD_MC_INTERP_BBBB (15)

/*
  Indexing of vetices and edges:
  0--0--1
  |     |
  3     1
  |     |
  3--2--2

  Interpolation formulas:
  i = horizontal position in the destination block, 0...blk_sz-1
  j = vertical position  in the destination block, 0...blk_sz-1
  blk_sz = length of an edge in the destination block in pixels
  blk_sz2 = blk_sz*blk_sz
  w[k]   = bilinear VECTOR weight of vertex k:
  w[0]   = (blk_sz-i)*(blk_sz-j)/blk_sz2
  w[1]   = i*(blk_sz-j)/blk_sz2
  w[2]   = i*j/blk_sz2
  w[3]   = (blk_sz-i)*j/blk_sz2
  u[k]   = bilinear weight of edge k:
  u[0]   = (blk_sz-j)/blk_sz
  u[1]   = i/blk_sz
  u[2]   = j/blk_sz
  u[3]   = (blk_sz-i)/blk_sz
  s[k]   = bilinear weight of vertex k:
  m[k]   = motion vector k
  s[k,c] = BLEND weight of vertex k in a split block whose outside corner is c:
  s[c+0,c] = w[c]+w[c+1]/2+w[c+3]/2
  s[c+1,c] = w[c+1]/2
  s[c+2,c] = w[c+2]
  s[c+3,c] = w[c+3]/2
  src[m] = source image value of pixel (i,j) offset by the given motion
   vector, bilinearly interpolated from the values at half-pel locations.
   First, horizontal interpolation is done, and the result is truncated to
   an integer, followed by vertical interpolation and truncation to an
   integer.

  Blending and interpolation functions are defined by the type of edges:
  V: vector interpolation edge.
     Along this edge, the motion vector is linearly interpolated from the
      value at one vertex to the other.
  B: blend edge
     Along this edge, the predicted value is linearly interpolated between
      that of the source image offset by the motion vector at one vertex
      to the value of the source image offset by the motion vector at the
      other vertex.
  This allows for C^0 switching between overlapped block motion compensation
   and smooth, per-pixel motion vector fields.
  The various interpolation formulas below are chosen so that:
   1) They have the proper form along each edge,
   2) They still match each other when a block is split to add a new motion
       vector in the middle (subject to the modifications to the blending
       weights listed above).
   2) All the vectors can be computed with simple bilinear interpolation,
   3) All the resulting source image pixels can be blended with simple
       bilinear interpolation.

  There are 6 similar classes of edge configurations:

  VVVV:
     src[m[0]*w[0]+m[1]*w[1]+m[2]*w[2]+m[3]*w[3]]
  BVVV (and rotations thereof):
     src[u[0]*m[0]+w[2]*m[2]+w[3]*m[3]]*w[0]+
      src[u[0]*m[1]+w[2]*m[2]+w[3]*m[3]]*w[1]+
      src[m[0]*w[0]+m[1]*w[1]+m[2]*w[2]+m[3]*w[3]]*u[2]
  BVBV (and rotation thereof):
     src[u[0]*m[0]+u[3]*m[2]]*u[2]+
      src[u[0]*m[1]+u[3]*m[3]]*u[1]
  VVBB (and rotations thereof):
     src[(w[0]+u[2])*m[0]+w[1]*m[1]]*w[0]+
      src[w[0]*m[0]+w[1]*m[1]+w[2]*m[2]+w[3]*m[3]]*w[1]+
      src[w[1]*m[1]+(w[2]+u[3])*m[2]]*w[2]+
      src[m[3]]*w[3]
  VBBB (and rotations thereof):
     src[(u[3]+w[2])*m[0]+w[1]*m[1]]*w[0]+
      src[w[0]*m[0]+(w[1]+u[2])*m[1]]*w[1]+
      src[m[2]]*w[2]+src[m[3]]*w[3]
  BBBB:
     src[m[0]]*w[0]+src[m[1]]*w[1]+src[m[2]]*w[2]+src[m[3]]*w[3]

  The remaining 10 edge configurations can be obtained by rotations of these
   formulas.

  Lots of indexing and finite differences where tiny errors may be made here.
  This is also a good candidate for SIMD optimization.
  Each case can be optimized separately.
  The first and last are expected to be more common than the others.
*/

/*A table of indices used to set up the rotated versions of each vector
   interpolation formula (see doc/mc-interpolation.nb for details).*/
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

/*These are for finite differences for the edge and vertex weights.*/
static const int SIDXS[][4]={
  { 1, 0, 0, 0},/*0*//*w0[k]*/
  {-1, 1, 0, 0},/*1*//*dwdi[k]*/
  {-1, 0, 0, 1},/*2*//*dwdj[k]*/
  { 1,-1, 1,-1},/*3*//*ddwdidj[k]*/
  { 0, 0, 1, 1},/*4*//*u0[k+3&3]*/
  { 1, 0,-1, 0},/*5*//*dudi[k+3&3*/
  { 0, 1, 0,-1},/*6*//*dudj[k+3&3*/
};

/*Set up the finite differences needed to interpolate a motion vector
   component.
  _x0:         The initial value.
  _dxdi:       The initial amount to increment per unit change in i.
  _dxdj:       The amount to increment per unit change in j.
  _ddxdidj:    The amount to increment _dxdi by per unit change in j.
  _mvx:        The component value of the 4 motion vectors.
  _m:          The index of the motion vector to use for each corner in the
                base orientation.
  _r:          The amount to rotate (clockwise) the formulas by (0...3).
  _log_blk_sz: The log base 2 of the block size.*/
static void od_mc_setup_mvc(ogg_int32_t *_x0,ogg_int32_t *_dxdi,
 ogg_int32_t *_dxdj,ogg_int32_t *_ddxdidj,const ogg_int32_t _mvx[4],
 const int _m[4],int _r,int _log_blk_sz){
  int log_blk_sz2;
  int k;
  log_blk_sz2=_log_blk_sz<<1;
  *_x0=_mvx[_m[0-_r&3]+_r&3];
  *_dxdi=_mvx[_m[1-_r&3]+_r&3]-*_x0>>_log_blk_sz;
  *_dxdj=_mvx[_m[3-_r&3]+_r&3]-*_x0>>_log_blk_sz;
  *_ddxdidj=_mvx[_m[0-_r&3]+_r&3]+_mvx[_m[2-_r&3]+_r&3]-
   _mvx[_m[1-_r&3]+_r&3]-_mvx[_m[3-_r&3]+_r&3]>>log_blk_sz2;
}

static void od_mc_setup_w(int _w0[4],int _dwdi[4],int _dwdj[4],
 int _ddwdidj[4],int _r,int _log_blk_sz){
  int log_blk_sz2;
  log_blk_sz2=_log_blk_sz<<1;
  _w0[0-_r&3]=1<<log_blk_sz2;
  _w0[1-_r&3]=_w0[2-_r&3]=_w0[3-_r&3]=0;
  _dwdi[0-_r&3]=-1<<_log_blk_sz;
  _dwdi[1-_r&3]=1<<_log_blk_sz;
  _dwdi[2-_r&3]=_dwdi[3-_r&3]=0;
  _dwdj[0-_r&3]=-1<<_log_blk_sz;
  _dwdj[1-_r&3]=_dwdj[2-_r&3]=0;
  _dwdj[3-_r&3]=1<<_log_blk_sz;
  _ddwdidj[0-_r&3]=_ddwdidj[2-_r&3]=1;
  _ddwdidj[1-_r&3]=_ddwdidj[3-_r&3]=-1;
}

static void od_mc_predict8(unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],int _interp_type,int _log_blk_sz){
  int log_blk_sz2;
  int blk_sz;
  int r;
  blk_sz=1<<_log_blk_sz;
  log_blk_sz2=_log_blk_sz<<1;
  r=0;
  switch(_interp_type){
    case OD_MC_INTERP_VVVV:{
      ogg_int32_t x0;
      ogg_int32_t y0;
      ogg_int32_t x;
      ogg_int32_t y;
      ogg_int32_t dxdi;
      ogg_int32_t dydi;
      ogg_int32_t dxdj;
      ogg_int32_t dydj;
      ogg_int32_t ddxdidj;
      ogg_int32_t ddydidj;
      int         i;
      int         j;
      od_mc_setup_mvc(&x0,&dxdi,&dxdj,&ddxdidj,_mvx,
       MIDXS[0]/*0,1,2,3*/,0,_log_blk_sz);
      od_mc_setup_mvc(&y0,&dydi,&dydj,&ddydidj,_mvy,
       MIDXS[0]/*0,1,2,3,*/,0,_log_blk_sz);
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        x=x0;
        y=y0;
        dst=_dst;
        for(i=0;i<=blk_sz;i++){
          const unsigned char *p;
          ogg_int32_t          xf;
          ogg_int32_t          yf;
          printf("<%8.4lf,%8.4lf>%s",x/(double)0x40000,y/(double)0x40000,
           i<blk_sz?"::":"\n");
          xf=x&0x3FFFF;
          yf=y&0x3FFFF;
          p=_src+(x>>18)+i+((y>>18)+j)*_systride;
          dst[0]=(unsigned char)((p[0]*(0x40000-xf)+p[1]*xf>>18)*
           (0x40000-yf)+(p[_systride]*(0x40000-xf)+
           p[_systride+1]*xf>>18)*yf>>18);
          x+=dxdi;
          y+=dydi;
          dst++;
        }
        x0+=dxdj;
        y0+=dydj;
        dxdi+=ddxdidj;
        dydi+=ddydidj;
        _dst+=_dystride;
      }
    }break;
    case OD_MC_INTERP_VVVB:r++;
    case OD_MC_INTERP_VVBV:r++;
    case OD_MC_INTERP_VBVV:r++;
    case OD_MC_INTERP_BVVV:{
      ogg_int32_t x0[3];
      ogg_int32_t y0[3];
      ogg_int32_t x[3];
      ogg_int32_t y[3];
      ogg_int32_t dxdi[3];
      ogg_int32_t dydi[3];
      ogg_int32_t dxdj[3];
      ogg_int32_t dydj[3];
      ogg_int32_t ddxdidj[3];
      ogg_int32_t ddydidj[3];
      int         w0[4];
      int         w[4];
      int         dwdi[4];
      int         i;
      int         j;
      int         k;
      od_mc_setup_mvc(x0+0,dxdi+0,dxdj+0,ddxdidj+0,_mvx,
       MIDXS[1]/*0,0,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+0,dydi+0,dydj+0,ddydidj+0,_mvy,
       MIDXS[1]/*0,0,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(x0+1,dxdi+1,dxdj+1,ddxdidj+1,_mvx,
       MIDXS[2]/*1,1,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+1,dydi+1,dydj+1,ddydidj+1,_mvy,
       MIDXS[2]/*1,1,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(x0+2,dxdi+2,dxdj+2,ddxdidj+2,_mvx,
       MIDXS[0]/*0,1,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+2,dydi+2,dydj+2,ddydidj+2,_mvy,
       MIDXS[0]/*0,1,2,3*/,r,_log_blk_sz);
      w0[0-r&3]=1<<log_blk_sz2;
      w0[1-r&3]=w0[2-r&3]=w0[3-r&3]=0;
      dwdi[0-r&3]=-blk_sz;
      dwdi[1-r&3]=blk_sz;
      dwdi[2-r&3]=dwdi[3-r&3]=0;
      for(k=0;k<3;k++){
        x[k]=x0[k];
        y[k]=y0[k];
      }
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        dst=_dst;
        for(k=0;k<4;k++)w[k]=w0[k];
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t c[3];
          for(k=0;k<3;k++){
            const unsigned char *p;
            ogg_int32_t          xf;
            ogg_int32_t          yf;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            p=_src+(x[k]>>18)+i+((y[k]>>18)+j)*_systride;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
          }
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           (x[0]-dxdi[0])/(double)0x40000,
           (y[0]-dydi[0])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           (x[1]-dxdi[1])/(double)0x40000,
           (y[1]-dydi[1])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[2]+w[3],
           (x[2]-dxdi[2])/(double)0x40000,
           (y[2]-dydi[2])/(double)0x40000,i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*(w[2]+w[3])>>log_blk_sz2);
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        for(k=0;k<3;k++){
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
          dxdi[k]+=ddxdidj[k];
          dydi[k]+=ddydidj[k];
          x[k]=x0[k];
          y[k]=y0[k];
        }
        w0[0-r&3]-=blk_sz;
        w0[3-r&3]+=blk_sz;
        dwdi[0-r&3]++;
        dwdi[1-r&3]--;
        dwdi[2-r&3]++;
        dwdi[3-r&3]--;
        _dst+=_dystride;
      }
    }break;
    case OD_MC_INTERP_VBVB:r++;
    case OD_MC_INTERP_BVBV:{
      ogg_int32_t x0[2];
      ogg_int32_t y0[2];
      ogg_int32_t x[2];
      ogg_int32_t y[2];
      ogg_int32_t dxdi[2];
      ogg_int32_t dydi[2];
      ogg_int32_t dxdj[2];
      ogg_int32_t dydj[2];
      int         i[2];
      int         k;
      od_mc_setup_mvc(x0+0,dxdi+0,dxdj+0,i,_mvx,
       MIDXS[3]/*0,0,3,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+0,dydi+0,dydj+0,i,_mvy,
       MIDXS[3]/*0,0,3,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(x0+1,dxdi+1,dxdj+1,i,_mvx,
       MIDXS[4]/*1,1,2,2*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+1,dydi+1,dydj+1,i,_mvy,
       MIDXS[4]/*1,1,2,2*/,r,_log_blk_sz);
      for(i[1]=0;i[1]<=blk_sz;i[1]++){
        unsigned char *dst;
        for(k=0;k<2;k++){
          x[k]=x0[k];
          y[k]=y0[k];
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
        }
        dst=_dst;
        for(i[0]=0;i[0]<=blk_sz;i[0]++){
          ogg_int32_t c[3];
          for(k=0;k<2;k++){
            const unsigned char *p;
            ogg_int32_t          xf;
            ogg_int32_t          yf;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            p=_src+(x[k]>>18)+i[0]+((y[k]>>18)+i[1])*_systride;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
          }
          printf("%2X<%8.4lf,%8.4lf>",blk_sz-i[r],
           (x[0]-dxdi[0])/(double)0x40000,
           (y[0]-dydi[0])/(double)0x40000);
          printf("%2X<%8.4lf,%8.4lf>%s",i[r],
           (x[1]-dxdi[1])/(double)0x40000,
           (y[1]-dydi[1])/(double)0x40000,i[0]<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(c[0]*(blk_sz-i[r])+c[1]*i[r]>>_log_blk_sz);
          dst++;
        }
        _dst+=_dystride;
      }
    }break;
    case OD_MC_INTERP_VBBV:r++;
    case OD_MC_INTERP_BBVV:r++;
    case OD_MC_INTERP_BVVB:r++;
    case OD_MC_INTERP_VVBB:{
      const unsigned char *mvp;
      ogg_int32_t          mvxf;
      ogg_int32_t          mvyf;
      ogg_int32_t          x0[3];
      ogg_int32_t          y0[3];
      ogg_int32_t          x[3];
      ogg_int32_t          y[3];
      ogg_int32_t          dxdi[3];
      ogg_int32_t          dydi[3];
      ogg_int32_t          dxdj[3];
      ogg_int32_t          dydj[3];
      ogg_int32_t          ddxdidj[3];
      ogg_int32_t          ddydidj[3];
      ptrdiff_t            o0;
      int                  w0[4];
      int                  w[4];
      int                  dwdi[4];
      int                  i;
      int                  j;
      int                  k;
      od_mc_setup_mvc(x0+0,dxdi+0,dxdj+0,ddxdidj+0,_mvx,
       MIDXS[5]/*0,1,0,0*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+0,dydi+0,dydj+0,ddydidj+0,_mvy,
       MIDXS[5]/*0,1,0,0*/,r,_log_blk_sz);
      od_mc_setup_mvc(x0+1,dxdi+1,dxdj+1,ddxdidj+1,_mvx,
       MIDXS[0]/*0,1,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+1,dydi+1,dydj+1,ddydidj+1,_mvy,
       MIDXS[0]/*0,1,2,3*/,r,_log_blk_sz);
      od_mc_setup_mvc(x0+2,dxdi+2,dxdj+2,ddxdidj+2,_mvx,
       MIDXS[6]/*2,1,2,2*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+2,dydi+2,dydj+2,ddydidj+2,_mvy,
       MIDXS[6]/*2,1,2,2*/,r,_log_blk_sz);
      mvp=_src+(_mvx[3+r&3]>>18)+(_mvy[3+r&3]>>18)*_systride;
      mvxf=_mvx[3+r&3]&0x3FFFF;
      mvyf=_mvy[3+r&3]&0x3FFFF;
      o0=0;
      w0[0-r&3]=1<<log_blk_sz2;
      w0[1-r&3]=w0[2-r&3]=w0[3-r&3]=0;
      dwdi[0-r&3]=-blk_sz;
      dwdi[1-r&3]=blk_sz;
      dwdi[2-r&3]=dwdi[3-r&3]=0;
      for(k=0;k<3;k++){
        x[k]=x0[k];
        y[k]=y0[k];
      }
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        ptrdiff_t      o;
        o=o0;
        dst=_dst;
        for(k=0;k<4;k++)w[k]=w0[k];
        for(i=0;i<=blk_sz;i++){
          const unsigned char *p;
          ogg_int32_t          c[4];
          ogg_int32_t          xf;
          ogg_int32_t          yf;
          for(k=0;k<3;k++){
            p=_src+o+(x[k]>>18)+(y[k]>>18)*_systride;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
          }
          p=mvp+o;
          c[3]=(p[0]*(0x40000-mvxf)+p[1]*mvxf>>18)*(0x40000-mvyf)+
           (p[_systride]*(0x40000-mvxf)+p[_systride+1]*mvxf>>18)*mvyf>>18;
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           (x[0]-dxdi[0])/(double)0x40000,
           (y[0]-dydi[0])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           (x[1]-dxdi[1])/(double)0x40000,
           (y[1]-dydi[1])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[2],
           (x[2]-dxdi[2])/(double)0x40000,
           (y[2]-dydi[2])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[3],
           _mvx[3+r&3]/(double)0x40000,
           _mvy[3+r&3]/(double)0x40000,i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*w[2]+c[3]*w[3]>>log_blk_sz2);
          o++;
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        for(k=0;k<3;k++){
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
          x[k]=x0[k];
          y[k]=y0[k];
          dxdi[k]+=ddxdidj[k];
          dydi[k]+=ddydidj[k];
        }
        o0+=_systride;
        _dst+=_dystride;
        w0[0-r&3]-=blk_sz;
        w0[3-r&3]+=blk_sz;
        dwdi[0-r&3]++;
        dwdi[1-r&3]--;
        dwdi[2-r&3]++;
        dwdi[3-r&3]--;
      }
    }break;
    case OD_MC_INTERP_BBBV:r++;
    case OD_MC_INTERP_BBVB:r++;
    case OD_MC_INTERP_BVBB:r++;
    case OD_MC_INTERP_VBBB:{
      const unsigned char *mvp[2];
      ogg_int32_t          mvxf[2];
      ogg_int32_t          mvyf[2];
      ogg_int32_t          x0[2];
      ogg_int32_t          y0[2];
      ogg_int32_t          x[2];
      ogg_int32_t          y[2];
      ogg_int32_t          dxdi[2];
      ogg_int32_t          dydi[2];
      ogg_int32_t          dxdj[2];
      ogg_int32_t          dydj[2];
      ogg_int32_t          ddxdidj;
      ogg_int32_t          ddydidj;
      ptrdiff_t            o0;
      int                  w0[4];
      int                  w[4];
      int                  dwdi[4];
      int                  i;
      int                  j;
      int                  k;
      od_mc_setup_mvc(x0+0,dxdi+0,dxdj+0,&ddxdidj,_mvx,
       MIDXS[5]/*0,1,0,0*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+0,dydi+0,dydj+0,&ddydidj,_mvy,
       MIDXS[5]/*0,1,0,0*/,r,_log_blk_sz);
      od_mc_setup_mvc(x0+1,dxdi+1,dxdj+1,&ddxdidj,_mvx,
       MIDXS[7]/*0,1,1,1*/,r,_log_blk_sz);
      od_mc_setup_mvc(y0+1,dydi+1,dydj+1,&ddydidj,_mvy,
       MIDXS[7]/*0,1,1,1*/,r,_log_blk_sz);
      for(k=0;k<2;k++){
        mvp[k]=_src+(_mvx[2+k+r&3]>>18)+(_mvy[2+k+r&3]>>18)*_systride;
        mvxf[k]=_mvx[2+k+r&3]&0x3FFFF;
        mvyf[k]=_mvy[2+k+r&3]&0x3FFFF;
      }
      o0=0;
      w0[0-r&3]=1<<log_blk_sz2;
      w0[1-r&3]=w0[2-r&3]=w0[3-r&3]=0;
      dwdi[0-r&3]=-blk_sz;
      dwdi[1-r&3]=blk_sz;
      dwdi[2-r&3]=dwdi[3-r&3]=0;
      for(k=0;k<2;k++){
        x[k]=x0[k];
        y[k]=y0[k];
      }
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        ptrdiff_t      o;
        o=o0;
        dst=_dst;
        for(k=0;k<4;k++)w[k]=w0[k];
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t c[4];
          for(k=0;k<2;k++){
            const unsigned char *p;
            ogg_int32_t          xf;
            ogg_int32_t          yf;
            p=_src+o+(x[k]>>18)+(y[k]>>18)*_systride;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            p=mvp[k]+o;
            c[k+2]=(p[0]*(0x40000-mvxf[k])+p[1]*mvxf[k]>>18)*
             (0x40000-mvyf[k])+(p[_systride]*(0x40000-mvxf[k])+
             p[_systride+1]*mvxf[k]>>18)*mvyf[k]>>18;
          }
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           (x[0]-dxdi[0])/(double)0x40000,(y[0]-dydi[0])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           (x[1]-dxdi[1])/(double)0x40000,(y[1]-dydi[1])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[2],
           _mvx[2+r&3]/(double)0x40000,_mvy[2+r&3]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[3],
           _mvx[3+r&3]/(double)0x40000,_mvy[3+r&3]/(double)0x40000,
           i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*w[2]+c[3]*w[3]>>log_blk_sz2);
          o++;
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        for(k=0;k<2;k++){
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
          x[k]=x0[k];
          y[k]=y0[k];
          dxdi[k]+=ddxdidj;
          dydi[k]+=ddydidj;
        }
        o0+=_systride;
        _dst+=_dystride;
        w0[0-r&3]-=blk_sz;
        w0[3-r&3]+=blk_sz;
        dwdi[0-r&3]++;
        dwdi[1-r&3]--;
        dwdi[2-r&3]++;
        dwdi[3-r&3]--;
      }
    }break;
    case OD_MC_INTERP_BBBB:{
      const unsigned char *mvp[4];
      ogg_int32_t          mvxf[4];
      ogg_int32_t          mvyf[4];
      ptrdiff_t            o0;
      int                  w0[4];
      int                  w[4];
      int                  dwdi[4];
      int                  i;
      int                  j;
      int                  k;
      for(k=0;k<4;k++){
        mvp[k]=_src+(_mvx[k]>>18)+(_mvy[k]>>18)*_systride;
        mvxf[k]=_mvx[k]&0x3FFFF;
        mvyf[k]=_mvy[k]&0x3FFFF;
      }
      o0=0;
      w0[0]=1<<log_blk_sz2;
      w0[1]=w0[2]=w0[3]=0;
      dwdi[0]=-blk_sz;
      dwdi[1]=blk_sz;
      dwdi[2]=dwdi[3]=0;
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        ptrdiff_t      o;
        o=o0;
        dst=_dst;
        for(k=0;k<4;k++)w[k]=w0[k];
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t c[4];
          for(k=0;k<4;k++){
            const unsigned char *p;
            p=mvp[k]+o;
            c[k]=(p[0]*(0x40000-mvxf[k])+p[1]*mvxf[k]>>18)*(0x40000-mvyf[k])+
             (p[_systride]*(0x40000-mvxf[k])+
             p[_systride+1]*mvxf[k]>>18)*mvyf[k]>>18;
          }
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           _mvx[0]/(double)0x40000,_mvy[0]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           _mvx[1]/(double)0x40000,_mvy[1]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[2],
           _mvx[2]/(double)0x40000,_mvy[2]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[3],
           _mvx[3]/(double)0x40000,_mvy[3]/(double)0x40000,
           i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*w[2]+c[3]*w[3]>>log_blk_sz2);
          o++;
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        o0+=_systride;
        _dst+=_dystride;
        w0[0]-=blk_sz;
        w0[3]+=blk_sz;
        dwdi[0]++;
        dwdi[1]--;
        dwdi[2]++;
        dwdi[3]--;
      }
    }
  }
}

static void od_mc_setup_mvc_split(ogg_int32_t *_x0,ogg_int32_t *_dxdi,
 ogg_int32_t *_dxdj,ogg_int32_t *_ddxdidj,const ogg_int32_t _mvx[4],
 const int _m[4],int _r,int _c,int _s1,int _s3,int _log_blk_sz){
  ogg_int32_t dx;
  int         log_blk_sz2;
  int         k;
  od_mc_setup_mvc(_x0,_dxdi,_dxdj,_ddxdidj,_mvx,_m,_r,_log_blk_sz);
  log_blk_sz2=_log_blk_sz<<1;
  if(_s1){
    k=_c+1&3;
    dx=_mvx[_m[_c-_r&3]+_r&3]-_mvx[_m[k-_r&3]+_r&3]>>1;
    *_x0+=dx*SIDXS[0][k];
    *_dxdi+=dx*SIDXS[1][k]>>_log_blk_sz;
    *_dxdj+=dx*SIDXS[2][k]>>_log_blk_sz;
    *_ddxdidj+=dx*SIDXS[3][k]>>log_blk_sz2;
  }
  if(_s3){
    k=_c+3&3;
    dx=_mvx[_m[_c-_r&3]+_r&3]-_mvx[_m[k-_r&3]+_r&3]>>1;
    *_x0+=dx*SIDXS[0][k];
    *_dxdi+=dx*SIDXS[1][k]>>_log_blk_sz;
    *_dxdj+=dx*SIDXS[2][k]>>_log_blk_sz;
    *_ddxdidj+=dx*SIDXS[3][k]>>log_blk_sz2;
  }
}

static void od_mc_setup_w_split(int _w0[4],int _dwdi[4],int _dwdj[4],
 int _ddwdidj[4],int _r,int _c,int _s1,int _s3,int _log_blk_sz){
  int log_blk_sz2;
  int k;
  log_blk_sz2=_log_blk_sz<<1;
  _w0[0-_r&3]=2<<log_blk_sz2;
  _w0[1-_r&3]=_w0[2-_r&3]=_w0[3-_r&3]=0;
  _dwdi[0-_r&3]=-2<<_log_blk_sz;
  _dwdi[1-_r&3]=2<<_log_blk_sz;
  _dwdi[2-_r&3]=_dwdi[3-_r&3]=0;
  _dwdj[0-_r&3]=-2<<_log_blk_sz;
  _dwdj[1-_r&3]=_dwdj[2-_r&3]=0;
  _dwdj[3-_r&3]=2<<_log_blk_sz;
  _ddwdidj[0-_r&3]=_ddwdidj[2-_r&3]=2;
  _ddwdidj[1-_r&3]=_ddwdidj[3-_r&3]=-2;
  _c=_c-_r&3;
  if(_s1){
    k=_c+1&3;
    _w0[k]>>=1;
    _w0[_c]+=_w0[k];
    _dwdi[k]>>=1;
    _dwdi[_c]+=_dwdi[k];
    _dwdj[k]>>=1;
    _dwdj[_c]+=_dwdj[k];
    _ddwdidj[k]>>=1;
    _ddwdidj[_c]+=_ddwdidj[k];
  }
  if(_s3){
    k=_c+3&3;
    _w0[k]>>=1;
    _w0[_c]+=_w0[k];
    _dwdi[k]>>=1;
    _dwdi[_c]+=_dwdi[k];
    _dwdj[k]>>=1;
    _dwdj[_c]+=_dwdj[k];
    _ddwdidj[k]>>=1;
    _ddwdidj[_c]+=_ddwdidj[k];
  }
}

static void od_mc_predict8_split(unsigned char *_dst,int _dystride,
 const unsigned char *_src,int _systride,const ogg_int32_t _mvx[4],
 const ogg_int32_t _mvy[4],int _interp_type,int _c,int _s1,int _s3,
 int _log_blk_sz){
  int log_blk_sz2;
  int blk_sz;
  int r;
  blk_sz=1<<_log_blk_sz;
  log_blk_sz2=_log_blk_sz<<1;
  r=0;
  switch(_interp_type){
    case OD_MC_INTERP_VVVV:{
      ogg_int32_t x0;
      ogg_int32_t y0;
      ogg_int32_t x;
      ogg_int32_t y;
      ogg_int32_t dxdi;
      ogg_int32_t dydi;
      ogg_int32_t dxdj;
      ogg_int32_t dydj;
      ogg_int32_t ddxdidj;
      ogg_int32_t ddydidj;
      int         i;
      int         j;
      od_mc_setup_mvc_split(&x0,&dxdi,&dxdj,&ddxdidj,_mvx,
       MIDXS[0]/*0,1,2,3*/,0,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(&y0,&dydi,&dydj,&ddydidj,_mvy,
       MIDXS[0]/*0,1,2,3*/,0,_c,0,0,_log_blk_sz);
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        x=x0;
        y=y0;
        dst=_dst;
        for(i=0;i<=blk_sz;i++){
          const unsigned char *p;
          ogg_int32_t          xf;
          ogg_int32_t          yf;
          printf("<%8.4lf,%8.4lf>%s",x/(double)0x40000,y/(double)0x40000,
           i<blk_sz?"::":"\n");
          xf=x&0x3FFFF;
          yf=y&0x3FFFF;
          p=_src+(x>>18)+i+((y>>18)+j)*_systride;
          dst[0]=(unsigned char)((p[0]*(0x40000-xf)+p[1]*xf>>18)*
           (0x40000-yf)+(p[_systride]*(0x40000-xf)+
           p[_systride+1]*xf>>18)*yf>>18);
          x+=dxdi;
          y+=dydi;
          dst++;
        }
        x0+=dxdj;
        y0+=dydj;
        dxdi+=ddxdidj;
        dydi+=ddydidj;
        _dst+=_dystride;
      }
    }break;
    case OD_MC_INTERP_VVVB:r++;
    case OD_MC_INTERP_VVBV:r++;
    case OD_MC_INTERP_VBVV:r++;
    case OD_MC_INTERP_BVVV:{
      ogg_int32_t x0[3];
      ogg_int32_t y0[3];
      ogg_int32_t x[3];
      ogg_int32_t y[3];
      ogg_int32_t dxdi[3];
      ogg_int32_t dydi[3];
      ogg_int32_t dxdj[3];
      ogg_int32_t dydj[3];
      ogg_int32_t ddxdidj[3];
      ogg_int32_t ddydidj[3];
      int         w0[4];
      int         w[4];
      int         dwdi[4];
      int         dwdj[4];
      int         ddwdidj[4];
      int         i;
      int         j;
      int         k;
      od_mc_setup_mvc_split(x0+0,dxdi+0,dxdj+0,ddxdidj+0,_mvx,
       MIDXS[1]/*0,0,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+0,dydi+0,dydj+0,ddydidj+0,_mvy,
       MIDXS[1]/*0,0,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(x0+1,dxdi+1,dxdj+1,ddxdidj+1,_mvx,
       MIDXS[2]/*1,1,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+1,dydi+1,dydj+1,ddydidj+1,_mvy,
       MIDXS[2]/*1,1,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(x0+2,dxdi+2,dxdj+2,ddxdidj+2,_mvx,
       MIDXS[0]/*0,1,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+2,dydi+2,dydj+2,ddydidj+2,_mvy,
       MIDXS[0]/*0,1,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_w_split(w0,dwdi,dwdj,ddwdidj,r,
       _c,_s1&&(_c+1-r&3)==1,_s3&&(_c+3-r&3)==0,_log_blk_sz);
      w0[2]+=w0[3];
      dwdi[2]+=dwdi[3];
      dwdj[2]+=dwdj[3];
      ddwdidj[2]+=ddwdidj[3];
      for(k=0;k<3;k++){
        x[k]=x0[k];
        y[k]=y0[k];
        w[k]=w0[k];
      }
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        dst=_dst;
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t c[3];
          for(k=0;k<3;k++){
            const unsigned char *p;
            ogg_int32_t          xf;
            ogg_int32_t          yf;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            p=_src+(x[k]>>18)+i+((y[k]>>18)+j)*_systride;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
          }
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           (x[0]-dxdi[0])/(double)0x40000,
           (y[0]-dydi[0])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           (x[1]-dxdi[1])/(double)0x40000,
           (y[1]-dydi[1])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[2],
           (x[2]-dxdi[2])/(double)0x40000,
           (y[2]-dydi[2])/(double)0x40000,i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*w[2]>>log_blk_sz2+1);
          dst++;
          for(k=0;k<3;k++)w[k]+=dwdi[k];
        }
        for(k=0;k<3;k++){
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
          w0[k]+=dwdj[k];
          dxdi[k]+=ddxdidj[k];
          dydi[k]+=ddydidj[k];
          dwdi[k]+=ddwdidj[k];
          x[k]=x0[k];
          y[k]=y0[k];
          w[k]=w0[k];
        }
        _dst+=_dystride;
      }
    }break;
    case OD_MC_INTERP_VBVB:r++;
    case OD_MC_INTERP_BVBV:{
      ogg_int32_t x0[2];
      ogg_int32_t y0[2];
      ogg_int32_t x[2];
      ogg_int32_t y[2];
      ogg_int32_t dxdi[2];
      ogg_int32_t dydi[2];
      ogg_int32_t dxdj[2];
      ogg_int32_t dydj[2];
      int         u0;
      int         u;
      int         dudi;
      int         dudj;
      int         ddudidj;
      int         i;
      int         j;
      int         k;
      od_mc_setup_mvc_split(x0+0,dxdi+0,dxdj+0,&i,_mvx,
       MIDXS[3]/*0,0,3,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+0,dydi+0,dydj+0,&i,_mvy,
       MIDXS[3]/*0,0,3,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(x0+1,dxdi+1,dxdj+1,&i,_mvx,
       MIDXS[4]/*1,1,2,2*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+1,dydi+1,dydj+1,&i,_mvy,
       MIDXS[4]/*1,1,2,2*/,r,_c,0,0,_log_blk_sz);
      u0=SIDXS[4][r]/*0,0,1,1*/<<log_blk_sz2+1;
      dudi=SIDXS[5][r]/*1,0,-1,0*/<<_log_blk_sz+1;
      dudj=SIDXS[6][r]/*0,1,0,-1*/<<_log_blk_sz+1;
      ddudidj=0;
      if(_s3&&(_c-r&1))k=_c+3&3;
      else if(_s1&&!(_c-r&1))k=_c+1&3;
      else k=-1;
      if(k>=0){
        if(k+1-r&2){
          u0-=SIDXS[0][k]<<log_blk_sz2;
          dudi-=SIDXS[1][k]/*-1,1,0,0*/<<_log_blk_sz;
          dudj-=SIDXS[2][k]/*-1,0,0,1*/<<_log_blk_sz;
          ddudidj-=SIDXS[3][k]/*1,-1,1,-1*/;
        }
        else{
          u0+=SIDXS[0][k]<<log_blk_sz2;
          dudi+=SIDXS[1][k]/*-1,1,0,0*/<<_log_blk_sz;
          dudj+=SIDXS[2][k]/*-1,0,0,1*/<<_log_blk_sz;
          ddudidj+=SIDXS[3][k]/*1,-1,1,-1*/;
        }
      }
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        for(k=0;k<2;k++){
          x[k]=x0[k];
          y[k]=y0[k];
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
        }
        dst=_dst;
        u=u0;
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t c[3];
          for(k=0;k<2;k++){
            const unsigned char *p;
            ogg_int32_t          xf;
            ogg_int32_t          yf;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            p=_src+(x[k]>>18)+i+((y[k]>>18)+j)*_systride;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
          }
          printf("%2X<%8.4lf,%8.4lf>",(1<<log_blk_sz2+1)-u,
           (x[0]-dxdi[0])/(double)0x40000,
           (y[0]-dydi[0])/(double)0x40000);
          printf("%2X<%8.4lf,%8.4lf>%s",u,
           (x[1]-dxdi[1])/(double)0x40000,
           (y[1]-dydi[1])/(double)0x40000,i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*((1<<log_blk_sz2+1)-u)+c[1]*u>>log_blk_sz2+1);
          dst++;
          u+=dudi;
        }
        _dst+=_dystride;
        u0+=dudj;
        dudi+=ddudidj;
      }
    }break;
    case OD_MC_INTERP_VBBV:r++;
    case OD_MC_INTERP_BBVV:r++;
    case OD_MC_INTERP_BVVB:r++;
    case OD_MC_INTERP_VVBB:{
      const unsigned char *mvp;
      ogg_int32_t          mvxf;
      ogg_int32_t          mvyf;
      ogg_int32_t          x0[3];
      ogg_int32_t          y0[3];
      ogg_int32_t          x[3];
      ogg_int32_t          y[3];
      ogg_int32_t          dxdi[3];
      ogg_int32_t          dydi[3];
      ogg_int32_t          dxdj[3];
      ogg_int32_t          dydj[3];
      ogg_int32_t          ddxdidj[3];
      ogg_int32_t          ddydidj[3];
      ptrdiff_t            o0;
      int                  w0[4];
      int                  w[4];
      int                  dwdi[4];
      int                  dwdj[4];
      int                  ddwdidj[4];
      int                  i;
      int                  j;
      int                  k;
      od_mc_setup_mvc_split(x0+0,dxdi+0,dxdj+0,ddxdidj+0,_mvx,
       MIDXS[5]/*0,1,0,0*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+0,dydi+0,dydj+0,ddydidj+0,_mvy,
       MIDXS[5]/*0,1,0,0*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(x0+1,dxdi+1,dxdj+1,ddxdidj+1,_mvx,
       MIDXS[0]/*0,1,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+1,dydi+1,dydj+1,ddydidj+1,_mvy,
       MIDXS[0]/*0,1,2,3*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(x0+2,dxdi+2,dxdj+2,ddxdidj+2,_mvx,
       MIDXS[6]/*2,1,2,2*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+2,dydi+2,dydj+2,ddydidj+2,_mvy,
       MIDXS[6]/*2,1,2,2*/,r,_c,0,0,_log_blk_sz);
      mvp=_src+(_mvx[3+r&3]>>18)+(_mvy[3+r&3]>>18)*_systride;
      mvxf=_mvx[3+r&3]&0x3FFFF;
      mvyf=_mvy[3+r&3]&0x3FFFF;
      o0=0;
      od_mc_setup_w_split(w0,dwdi,dwdj,ddwdidj,r,
       _c,_s1&&(_c-r&2),_s3&&(_c+3-r&2),_log_blk_sz);
      for(k=0;k<3;k++){
        x[k]=x0[k];
        y[k]=y0[k];
      }
      for(k=0;k<4;k++)w[k]=w0[k];
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        ptrdiff_t      o;
        o=o0;
        dst=_dst;
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t    c[4];
          const unsigned char *p;
          ogg_int32_t          xf;
          ogg_int32_t          yf;
          for(k=0;k<3;k++){
            p=_src+o+(x[k]>>18)+(y[k]>>18)*_systride;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
          }
          p=mvp+o;
          c[3]=(p[0]*(0x40000-mvxf)+p[1]*mvxf>>18)*(0x40000-mvyf)+
           (p[_systride]*(0x40000-mvxf)+p[_systride+1]*mvxf>>18)*mvyf>>18;
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           (x[0]-dxdi[0])/(double)0x40000,
           (y[0]-dydi[0])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           (x[1]-dxdi[1])/(double)0x40000,
           (y[1]-dydi[1])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[2],
           (x[2]-dxdi[2])/(double)0x40000,
           (y[2]-dydi[2])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[3],
           _mvx[3+r&3]/(double)0x40000,
           _mvy[3+r&3]/(double)0x40000,i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*w[2]+c[3]*w[3]>>log_blk_sz2+1);
          o++;
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        for(k=0;k<3;k++){
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
          dxdi[k]+=ddxdidj[k];
          dydi[k]+=ddydidj[k];
          x[k]=x0[k];
          y[k]=y0[k];
        }
        o0+=_systride;
        _dst+=_dystride;
        for(k=0;k<4;k++){
          w0[k]+=dwdj[k];
          dwdi[k]+=ddwdidj[k];
          w[k]=w0[k];
        }
      }
    }break;
    case OD_MC_INTERP_BBBV:r++;
    case OD_MC_INTERP_BBVB:r++;
    case OD_MC_INTERP_BVBB:r++;
    case OD_MC_INTERP_VBBB:{
      const unsigned char *mvp[2];
      ogg_int32_t          mvxf[2];
      ogg_int32_t          mvyf[2];
      ogg_int32_t          x0[2];
      ogg_int32_t          y0[2];
      ogg_int32_t          x[2];
      ogg_int32_t          y[2];
      ogg_int32_t          dxdi[2];
      ogg_int32_t          dydi[2];
      ogg_int32_t          dxdj[2];
      ogg_int32_t          dydj[2];
      ogg_int32_t          ddxdidj;
      ogg_int32_t          ddydidj;
      ptrdiff_t            o0;
      int                  w0[4];
      int                  w[3];
      int                  dwdi[4];
      int                  dwdj[4];
      int                  ddwdidj[4];
      int                  i;
      int                  j;
      int                  k;
      od_mc_setup_mvc_split(x0+0,dxdi+0,dxdj+0,&ddxdidj,_mvx,
       MIDXS[5]/*0,1,0,0*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+0,dydi+0,dydj+0,&ddydidj,_mvy,
       MIDXS[5]/*0,1,0,0*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(x0+1,dxdi+1,dxdj+1,&ddxdidj,_mvx,
       MIDXS[7]/*0,1,1,1*/,r,_c,0,0,_log_blk_sz);
      od_mc_setup_mvc_split(y0+1,dydi+1,dydj+1,&ddydidj,_mvy,
       MIDXS[7]/*0,1,1,1*/,r,_c,0,0,_log_blk_sz);
      for(k=0;k<2;k++){
        mvp[k]=_src+(_mvx[2+k+r&3]>>18)+(_mvy[2+k+r&3]>>18)*_systride;
        mvxf[k]=_mvx[2+k+r&3]&0x3FFFF;
        mvyf[k]=_mvy[2+k+r&3]&0x3FFFF;
      }
      o0=0;
      od_mc_setup_w_split(w0,dwdi,dwdj,ddwdidj,r,
       _c,_s1&&(_c+1-r&3)!=1,_s3&&(_c+3-r&2)!=0,_log_blk_sz);
      for(k=0;k<2;k++){
        x[k]=x0[k];
        y[k]=y0[k];
      }
      for(k=0;k<4;k++)w[k]=w0[k];
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        ptrdiff_t      o;
        o=o0;
        dst=_dst;
        for(i=0;i<=blk_sz;i++){
          const unsigned char *p;
          ogg_int32_t          c[4];
          ogg_int32_t          xf;
          ogg_int32_t          yf;
          for(k=0;k<2;k++){
            p=_src+o+(x[k]>>18)+(y[k]>>18)*_systride;
            xf=x[k]&0x3FFFF;
            yf=y[k]&0x3FFFF;
            c[k]=(p[0]*(0x40000-xf)+p[1]*xf>>18)*(0x40000-yf)+
             (p[_systride]*(0x40000-xf)+p[_systride+1]*xf>>18)*yf>>18;
            x[k]+=dxdi[k];
            y[k]+=dydi[k];
            p=mvp[k]+o;
            c[k+2]=(p[0]*(0x40000-mvxf[k])+p[1]*mvxf[k]>>18)*
             (0x40000-mvyf[k])+(p[_systride]*(0x40000-mvxf[k])+
             p[_systride+1]*mvxf[k]>>18)*mvyf[k]>>18;
          }
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           (x[0]-dxdi[0])/(double)0x40000,(y[0]-dydi[0])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           (x[1]-dxdi[1])/(double)0x40000,(y[1]-dydi[1])/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[2],
           _mvx[2+r&3]/(double)0x40000,_mvy[2+r&3]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[3],
           _mvx[3+r&3]/(double)0x40000,_mvy[3+r&3]/(double)0x40000,
           i<blk_sz?"::":"\n");
          dst[0]=c[0]*w[0]+c[1]*w[1]+c[2]*w[2]+c[3]*w[3]>>log_blk_sz2+1;
          o++;
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        for(k=0;k<2;k++){
          x0[k]+=dxdj[k];
          y0[k]+=dydj[k];
          dxdi[k]+=ddxdidj;
          dydi[k]+=ddydidj;
          x[k]=x0[k];
          y[k]=y0[k];
        }
        o0+=_systride;
        _dst+=_dystride;
        for(k=0;k<4;k++){
          w0[k]+=dwdj[k];
          dwdi[k]+=ddwdidj[k];
          w[k]=w0[k];
        }
      }
    }break;
    case OD_MC_INTERP_BBBB:{
      const unsigned char *mvp[4];
      ogg_int32_t          mvxf[4];
      ogg_int32_t          mvyf[4];
      ptrdiff_t            o0;
      int                  w0[4];
      int                  w[4];
      int                  dwdi[4];
      int                  dwdj[4];
      int                  ddwdidj[4];
      int                  i;
      int                  j;
      int                  k;
      for(k=0;k<4;k++){
        mvp[k]=_src+(_mvx[k]>>18)+(_mvy[k]>>18)*_systride;
        mvxf[k]=_mvx[k]&0x3FFFF;
        mvyf[k]=_mvy[k]&0x3FFFF;
      }
      o0=0;
      od_mc_setup_w_split(w0,dwdi,dwdj,ddwdidj,r,_c,_s1,_s3,_log_blk_sz);
      for(k=0;k<4;k++)w[k]=w0[k];
      for(j=0;j<=blk_sz;j++){
        unsigned char *dst;
        ptrdiff_t      o;
        o=o0;
        dst=_dst;
        for(i=0;i<=blk_sz;i++){
          ogg_int32_t c[4];
          for(k=0;k<4;k++){
            const unsigned char *p;
            p=mvp[k]+o;
            c[k]=(p[0]*(0x40000-mvxf[k])+p[1]*mvxf[k]>>18)*(0x40000-mvyf[k])+
             (p[_systride]*(0x40000-mvxf[k])+
             p[_systride+1]*mvxf[k]>>18)*mvyf[k]>>18;
          }
          printf("%3X<%8.4lf,%8.4lf>",w[0],
           _mvx[0]/(double)0x40000,_mvy[0]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[1],
           _mvx[1]/(double)0x40000,_mvy[1]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>",w[2],
           _mvx[2]/(double)0x40000,_mvy[2]/(double)0x40000);
          printf("%3X<%8.4lf,%8.4lf>%s",w[3],
           _mvx[3]/(double)0x40000,_mvy[3]/(double)0x40000,
           i<blk_sz?"::":"\n");
          dst[0]=(unsigned char)(
           c[0]*w[0]+c[1]*w[1]+c[2]*w[2]+c[3]*w[3]>>log_blk_sz2+1);
          o++;
          dst++;
          for(k=0;k<4;k++)w[k]+=dwdi[k];
        }
        o0+=_systride;
        _dst+=_dystride;
        for(k=0;k<4;k++){
          w0[k]+=dwdj[k];
          dwdi[k]+=ddwdidj[k];
          w[k]=w0[k];
        }
      }
    }
  }
}

#if 1
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
         etype&1?'B':'V',etype&2?'B':'V',
         etype&4?'B':'V',etype&8?'B':'V',etype);
        printf("<%8.4lf,%8.4lf> <%8.4lf,%8.4lf>\n",
         mvx[0]/(double)0x40000,mvy[0]/(double)0x40000,
         mvx[1]/(double)0x40000,mvy[1]/(double)0x40000);
        printf("<%8.4lf,%8.4lf> <%8.4lf,%8.4lf>\n",
         mvx[3]/(double)0x40000,mvy[3]/(double)0x40000,
         mvx[2]/(double)0x40000,mvy[2]/(double)0x40000);
        od_mc_predict8(dst[0],sizeof(dst[0]),
         img[j+1<<log_blk_sz]+(i+1<<log_blk_sz),sizeof(img[0]),mvx,mvy,
         etype,log_blk_sz);
        for(y=0;y<=blk_sz;y++){
          for(x=0;x<=blk_sz;x++){
            printf("%2X%c",dst[y][x],x<blk_sz?' ':'\n');
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
          if(0||!(c&1)){
            etype&=~(1<<(c+1&3));
            if(c==3)etype|=(etype&8)>>3;
            else etype|=(etype&1<<c)<<1;
            s1=1;
          }
          else s1=0;
          if(!s1||!(etype>>c&1)){
            mvx[c+1&3]=mvx[c]+mvx[c+1&3]>>1;
            mvy[c+1&3]=mvy[c]+mvy[c+1&3]>>1;
          }
          if(0||(c&1)){
            etype&=~(1<<(c+2&3));
            if(c==1)etype|=(etype&1)<<3;
            else etype|=(etype&1<<(c+3&3))>>1;
            s3=1;
          }
          else s3=0;
          if(!s3||!(etype<<(-c&3)&8)){
            mvx[c+3&3]=mvx[c]+mvx[c+3&3]>>1;
            mvy[c+3&3]=mvy[c]+mvy[c+3&3]>>1;
          }
          printf("Block (%i.%i,%i.%i): size %i, "
           "interpolation type: %c%c%c%c (0x%X)\n",
           i,((c+1&3)>>1)*5,j,(c>>1)*5,1<<log_blk_sz-1,
           etype&1?'B':'V',etype&2?'B':'V',
           etype&4?'B':'V',etype&8?'B':'V',etype);
          printf("<%9.5lf,%9.5lf> <%9.5lf,%9.5lf>\n",
           mvx[0]/(double)0x40000,mvy[0]/(double)0x40000,
           mvx[1]/(double)0x40000,mvy[1]/(double)0x40000);
          printf("<%9.5lf,%9.5lf> <%9.5lf,%9.5lf>\n",
           mvx[3]/(double)0x40000,mvy[3]/(double)0x40000,
           mvx[2]/(double)0x40000,mvy[2]/(double)0x40000);
          od_mc_predict8_split(dst2[c][0],sizeof(dst2[c][0]),
           img[(j+1<<1|(c>>1))<<log_blk_sz-1]+
           ((i+1<<1|((c+1&3)>>1))<<log_blk_sz-1),
           sizeof(img[0]),mvx,mvy,etype,c,s1,s3,log_blk_sz-1);
          memset(mismatch[c][0],0,sizeof(mismatch[c]));
          switch(c){
            case 0:{
              for(x=0;x<=blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst[0][x])mismatch[c][0][x]++;
              }
              /*for(y=1;y<=blk_sz>>1;y++){
                if(dst2[c][y][0]!=dst[y][0])mismatch[c][y][0]++;
              }*/
            }break;
            case 1:{
              for(x=0;x<=blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst[0][x+(blk_sz>>1)])mismatch[c][0][x]++;
              }
              /*for(y=1;y<=blk_sz>>1;y++){
                if(dst2[c][y][blk_sz>>1]!=dst[y][blk_sz]){
                  mismatch[c][y][blk_sz>>1]++;
                }
              }*/
              for(y=0;y<=blk_sz>>1;y++){
                if(dst2[c][y][0]!=dst2[0][y][blk_sz>>1])mismatch[c][y][0]++;
              }
            }break;
            case 2:{
              for(x=0;x<=blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst2[1][blk_sz>>1][x])mismatch[c][0][x]++;
              }
              /*for(y=0;y<=blk_sz>>1;y++){
                if(dst2[c][y][blk_sz>>1]!=dst[y+(blk_sz>>1)][blk_sz]){
                  mismatch[c][y][blk_sz>>1]++;
                }
              }*/
              for(x=0;x<blk_sz>>1;x++){
                if(dst2[c][blk_sz>>1][x]!=dst[blk_sz][x+(blk_sz>>1)]){
                  mismatch[c][blk_sz>>1][x]++;
                }
              }
            }break;
            case 3:{
              for(x=0;x<=blk_sz>>1;x++){
                if(dst2[c][0][x]!=dst2[0][blk_sz>>1][x])mismatch[c][0][x]++;
              }
              for(y=1;y<=blk_sz>>1;y++){
                if(dst2[c][y][blk_sz>>1]!=dst2[2][y][0]){
                  mismatch[c][y][blk_sz>>1]++;
                }
              }
              for(x=0;x<=blk_sz>>1;x++){
                if(dst2[c][blk_sz>>1][x]!=dst[blk_sz][x]){
                  mismatch[c][blk_sz>>1][x]++;
                }
              }
              /*for(y=0;y<blk_sz>>1;y++){
                if(dst2[c][y][0]!=dst[y+(blk_sz>>1)][0])mismatch[c][y][0]++;
              }*/
            }break;
          }
          for(y=0;y<=blk_sz>>1;y++){
            for(x=0;x<=blk_sz>>1;x++){
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
