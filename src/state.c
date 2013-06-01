/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "state.h"
#if defined(OD_X86ASM)
# include "x86/x86int.h"
#endif



/*Initializes the buffers used for reference frames.
  These buffers are padded with 16 extra pixels on each side, to allow
   (relatively) unrestricted motion vectors without special casing reading
   outside the image boundary.
  If chroma is decimated in either direction, the padding is reduced by an
   appropriate factor on the appropriate sides.*/
static void od_state_ref_imgs_init(od_state *_state,int _nrefs,int _nio){
  daala_info    *info;
  od_img        *img;
  od_img_plane  *iplane;
  unsigned char *ref_img_data;
  size_t         data_sz;
  int            plane_buf_width;
  int            plane_buf_height;
  int            imgi;
  int            pli;
  int            y;
  OD_ASSERT(_nrefs>=3);
  OD_ASSERT(_nrefs<=4);
  OD_ASSERT(_nio==2);
  info=&_state->info;
  data_sz=0;
  /*TODO: Check for overflow before allocating.*/
  for(pli=0;pli<info->nplanes;pli++){
    /*Reserve space for this plane in _nrefs reference images.*/
    plane_buf_width=_state->frame_width+(OD_UMV_PADDING<<1)
     >>info->plane_info[pli].xdec;
    plane_buf_height=_state->frame_height+(OD_UMV_PADDING<<1)
     >>info->plane_info[pli].ydec;
    data_sz+=plane_buf_width*plane_buf_height*_nrefs<<2;
#if defined(OD_DUMP_IMAGES)
    /*Reserve space for this plane in 1 visualization image.*/
    data_sz+=plane_buf_width*plane_buf_height<<2;
#endif
    /*Reserve space for this plane in _nio input/output images.*/
    data_sz+=plane_buf_width*plane_buf_height*_nio;
  }
  /*Reserve space for the line buffer in the up-sampler.*/
  data_sz+=(_state->frame_width+(OD_UMV_PADDING<<1)<<1)*8;
  _state->ref_img_data=ref_img_data=(unsigned char *)_ogg_malloc(data_sz);
  /*Fill in the reference image structures.*/
  for(imgi=0;imgi<_nrefs;imgi++){
    img=_state->ref_imgs+imgi;
    img->nplanes=info->nplanes;
    img->width=_state->frame_width<<1;
    img->height=_state->frame_height<<1;
    for(pli=0;pli<img->nplanes;pli++){
      plane_buf_width=(_state->frame_width+(OD_UMV_PADDING<<1)<<1)
       >>info->plane_info[pli].xdec;
      plane_buf_height=(_state->frame_height+(OD_UMV_PADDING<<1)<<1)
       >>info->plane_info[pli].ydec;
      iplane=img->planes+pli;
      iplane->data=ref_img_data
       +((OD_UMV_PADDING<<1)>>info->plane_info[pli].xdec)
       +plane_buf_width*((OD_UMV_PADDING<<1)>>info->plane_info[pli].ydec);
      ref_img_data+=plane_buf_width*plane_buf_height;
      iplane->xdec=info->plane_info[pli].xdec;
      iplane->ydec=info->plane_info[pli].ydec;
      iplane->xstride=1;
      iplane->ystride=plane_buf_width;
    }
  }
  /*Fill in the reconstruction image structure.*/
  for(imgi=0;imgi<_nio;imgi++){
    img=_state->io_imgs+imgi;
    img->nplanes=info->nplanes;
    img->width=_state->frame_width;
    img->height=_state->frame_height;
    for(pli=0;pli<img->nplanes;pli++){
      plane_buf_width=_state->frame_width+(OD_UMV_PADDING<<1)
       >>info->plane_info[pli].xdec;
      plane_buf_height=_state->frame_height+(OD_UMV_PADDING<<1)
       >>info->plane_info[pli].ydec;
      iplane=img->planes+pli;
      iplane->data=ref_img_data
       +(OD_UMV_PADDING>>info->plane_info[pli].xdec)
       +plane_buf_width*(OD_UMV_PADDING>>info->plane_info[pli].ydec);
      ref_img_data+=plane_buf_width*plane_buf_height;
      iplane->xdec=info->plane_info[pli].xdec;
      iplane->ydec=info->plane_info[pli].ydec;
      iplane->xstride=1;
      iplane->ystride=plane_buf_width;
    }
  }
  /*Fill in the line buffers.*/
  for(y=0;y<8;y++){
    _state->ref_line_buf[y]=ref_img_data+(OD_UMV_PADDING<<1);
    ref_img_data+=_state->frame_width+(OD_UMV_PADDING<<1)<<1;
  }
  /*Mark all of the reference image buffers available.*/
  for(imgi=0;imgi<_nrefs;imgi++)_state->ref_imgi[imgi]=-1;
#if defined(OD_DUMP_IMAGES)
  /*Fill in the visualization image structure.*/
  img=&_state->vis_img;
  img->nplanes=info->nplanes;
  img->width=_state->frame_width+(OD_UMV_PADDING<<1)<<1;
  img->height=_state->frame_height+(OD_UMV_PADDING<<1)<<1;
  for(pli=0;pli<img->nplanes;pli++){
    iplane=img->planes+pli;
    plane_buf_width=img->width>>info->plane_info[pli].xdec;
    plane_buf_height=img->height>>info->plane_info[pli].ydec;
    iplane->data=ref_img_data;
    ref_img_data+=plane_buf_width*plane_buf_height;
    iplane->xdec=info->plane_info[pli].xdec;
    iplane->ydec=info->plane_info[pli].ydec;
    iplane->xstride=1;
    iplane->ystride=plane_buf_width;
  }
#endif
}

static void od_state_mvs_init(od_state *_state){
  int nhmvbs;
  int nvmvbs;
  nhmvbs=_state->nhmbs+1<<2;
  nvmvbs=_state->nvmbs+1<<2;
  _state->mv_grid=(od_mv_grid_pt **)od_calloc_2d(nvmvbs+1,nhmvbs+1,
   sizeof(_state->mv_grid[0][0]));
}

void od_state_opt_vtbl_init_c(od_state *_state){
  _state->opt_vtbl.mc_predict1imv8=od_mc_predict1imv8_c;
  _state->opt_vtbl.mc_predict1fmv8=od_mc_predict1fmv8_c;
  _state->opt_vtbl.mc_blend_full8=od_mc_blend_full8_c;
  _state->opt_vtbl.mc_blend_full_split8=od_mc_blend_full_split8_c;
}

static void od_state_opt_vtbl_init(od_state *_state){
#if defined(OD_X86ASM)
  od_state_opt_vtbl_init_x86(_state);
#else
  od_state_opt_vtbl_init_c(_state);
#endif
}

int od_state_init(od_state *_state,const daala_info *_info){
  int nplanes;
  int pli;
  /*First validate the parameters.*/
  if(_info==NULL)return OD_EFAULT;
  nplanes=_info->nplanes;
  if(nplanes<=0||nplanes>OD_NPLANES_MAX)return OD_EINVAL;
  /*The first plane (the luma plane) must not be subsampled.*/
  if(_info->plane_info[0].xdec||_info->plane_info[0].ydec)return OD_EINVAL;
  memset(_state,0,sizeof(*_state));
  memcpy(&_state->info,_info,sizeof(*_info));
  /*Frame size is a multiple of a super block.*/
  _state->frame_width=_info->pic_width+(OD_SUPERBLOCK_SIZE-1)&~(OD_SUPERBLOCK_SIZE-1);
  _state->frame_height=_info->pic_height+(OD_SUPERBLOCK_SIZE-1)&~(OD_SUPERBLOCK_SIZE-1);
  _state->nhmbs=_state->frame_width>>4;
  _state->nvmbs=_state->frame_height>>4;
  od_state_opt_vtbl_init(_state);
  od_state_ref_imgs_init(_state,4,2);
  od_state_mvs_init(_state);
  for(pli=0;pli<nplanes;pli++){
    _state->adapt_row[pli].ctx=(od_adapt_ctx *)_ogg_malloc(
     _state->nhmbs*sizeof(*_state->adapt_row[pli].ctx));
    _state->adapt_row[pli].nhmbs = _state->nhmbs;
  }
  _state->nhsb=(_state->frame_width>>5);
  _state->nvsb=(_state->frame_height>>5);
  _state->bsize=(unsigned char *)_ogg_malloc(
      (_state->nhsb+2)*4 *
      (_state->nvsb+2)*4);
  _state->bstride = (_state->nhsb+2)*4;
  _state->bsize += 4*_state->bstride+4;
  return 0;
}

void od_state_clear(od_state *_state){
  int nplanes;
  int pli;
  nplanes=_state->info.nplanes;
  for(pli=nplanes;pli-->0;)_ogg_free(_state->adapt_row[pli].ctx);
  od_free_2d(_state->mv_grid);
  _ogg_free(_state->ref_img_data);
  _state->bsize -= 4*_state->bstride+4;
  _ogg_free(_state->bsize);
}

#if 0
/*Upsamples the reconstructed image to a reference image.
  TODO: Pipeline with reconstruction.*/
void od_state_upsample8(od_state *_state,int _refi){
  int pli;
  for(pli=0;pli<_state->io_imgs[OD_FRAME_REC].nplanes;pli++){
    od_img_plane  *siplane;
    od_img_plane  *diplane;
    unsigned char *src;
    unsigned char *dst;
    int            xpad;
    int            ypad;
    int            w;
    int            h;
    int            x;
    int            y;
    siplane=_state->io_imgs[OD_FRAME_REC].planes+pli;
    diplane=_state->ref_imgs[_refi].planes+pli;
    xpad=OD_UMV_PADDING>>siplane->xdec;
    ypad=OD_UMV_PADDING>>siplane->ydec;
    w=_state->io_imgs[OD_FRAME_REC].width>>siplane->xdec;
    h=_state->io_imgs[OD_FRAME_REC].height>>siplane->ydec;
    src=siplane->data;
    dst=diplane->data-(diplane->ystride<<1)*ypad;
    for(y=-ypad;y<h+ypad+2;y++){
      /*Horizontal filtering:*/
      if(y<h+ypad){
        unsigned char *buf;
        buf=_state->ref_line_buf[y&3];
        for(x=-xpad;x<-1;x++){
          *(buf+(x<<1))=src[0];
          *(buf+(x<<1|1))=src[0];
        }
        *(buf-2)=src[0];
        *(buf-1)=OD_CLAMP255(132*src[0]-4*src[1]+64>>7);
        buf[0]=OD_CLAMP255(121*src[0]+7*src[1]+64>>7);
        /*buf[0]=src[0];*/
        buf[1]=OD_CLAMP255(68*(src[0]+src[1])-4*(src[0]+src[2])+64>>7);
        for(x=1;x<w-2;x++){
          buf[x<<1]=OD_CLAMP255(7*src[x-1]+114*src[x]+7*src[x+1]+64>>7);
          /*buf[x<<1]=src[x];*/
          buf[x<<1|1]=OD_CLAMP255(68*(src[x]+src[x+1])-
           4*(src[x-1]+src[x+2])+64>>7);
        }
        buf[x<<1]=OD_CLAMP255(7*src[x-1]+114*src[x]+7*src[x+1]+64>>7);
        /*buf[x<<1]=src[x];*/
        buf[x<<1|1]=OD_CLAMP255(68*(src[x]+src[x+1])-
         4*(src[x-1]+src[x+1])+64>>7);
        x++;
        buf[x<<1]=OD_CLAMP255(7*src[x-1]+114*src[x]+7*src[x]+64>>7);
        /*buf[x<<1]=src[x];*/
        buf[x<<1|1]=OD_CLAMP255(132*src[x]-4*src[x-1]+64>>7);
        for(x++;x<w+xpad;x++){
          buf[x<<1]=src[w-1];
          buf[x<<1|1]=src[w-1];
        }
        if(y>=0&&y+1<h)src+=siplane->ystride;
      }
      /*Vertical filtering:*/
      if(y>=-ypad+2){
        if(y<1||y>h+1){
          memcpy(dst-(xpad<<1),_state->ref_line_buf[y-2&3]-(xpad<<1),
           w+(xpad<<1)<<1);
          /*fprintf(stderr,"%3i: ",y-2<<1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
          memcpy(dst-(xpad<<1),_state->ref_line_buf[y-2&3]-(xpad<<1),
           w+(xpad<<1)<<1);
          /*fprintf(stderr,"%3i: ",y-2<<1|1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
        }
        else{
          unsigned char *buf[4];
          buf[0]=_state->ref_line_buf[y-3&3];
          buf[1]=_state->ref_line_buf[y-2&3];
          buf[2]=_state->ref_line_buf[y-1&3];
          buf[3]=_state->ref_line_buf[y-0&3];
          for(x=-xpad<<1;x<w+xpad<<1;x++){
            *(dst+x)=OD_CLAMP255(7**(buf[0]+x)+114**(buf[1]+x)+
             7**(buf[2]+x)+64>>7);
            /**(dst+x)=*(buf[1]+x);*/
          }
          /*fprintf(stderr,"%3i: ",y-2<<1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
          for(x=-xpad<<1;x<w+xpad<<1;x++){
            *(dst+x)=OD_CLAMP255(68*(*(buf[1]+x)+*(buf[2]+x))-
             4*(*(buf[0]+x)+*(buf[3]+x))+64>>7);
          }
          /*fprintf(stderr,"%3i: ",y-2<<1|1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
        }
      }
    }
  }
}
#else
/*Upsamples the reconstructed image to a reference image.
  TODO: Pipeline with reconstruction.*/
void od_state_upsample8(od_state *_state,od_img *_dst,const od_img *_src){
  int pli;
  for(pli=0;pli<_state->io_imgs[OD_FRAME_REC].nplanes;pli++){
    const od_img_plane  *siplane;
    od_img_plane        *diplane;
    const unsigned char *src;
    unsigned char       *dst;
    int                  xpad;
    int                  ypad;
    int                  w;
    int                  h;
    int                  x;
    int                  y;
    siplane=_src->planes+pli;
    diplane=_dst->planes+pli;
    xpad=OD_UMV_PADDING>>siplane->xdec;
    ypad=OD_UMV_PADDING>>siplane->ydec;
    w=_src->width>>siplane->xdec;
    h=_src->height>>siplane->ydec;
    src=siplane->data;
    dst=diplane->data-(diplane->ystride<<1)*ypad;
    for(y=-ypad;y<h+ypad+3;y++){
      /*Horizontal filtering:*/
      if(y<h+ypad){
        unsigned char *buf;
        buf=_state->ref_line_buf[y&7];
        memset(buf-(xpad<<1),src[0],xpad-2<<1);
        /*for(x=-xpad;x<-2;x++){
          *(buf+(x<<1))=src[0];
          *(buf+(x<<1|1))=src[0];
        }*/
        *(buf-4)=src[0];
        *(buf-3)=OD_CLAMP255(31*src[0]+src[1]+16>>5);
        *(buf-2)=src[0];
        *(buf-1)=OD_CLAMP255(36*src[0]-5*src[1]+src[1]+16>>5);
        buf[0]=src[0];
        buf[1]=OD_CLAMP255(20*(src[0]+src[1])-
         5*(src[0]+src[2])+src[0]+src[3]+16>>5);
        buf[2]=src[1];
        buf[3]=OD_CLAMP255(20*(src[1]+src[2])-
         5*(src[0]+src[3])+src[0]+src[4]+16>>5);
        for(x=2;x<w-3;x++){
          buf[x<<1]=src[x];
          buf[x<<1|1]=OD_CLAMP255(20*(src[x]+src[x+1])-
           5*(src[x-1]+src[x+2])+src[x-2]+src[x+3]+16>>5);
        }
        buf[x<<1]=src[x];
        buf[x<<1|1]=OD_CLAMP255(20*(src[x]+src[x+1])-
         5*(src[x-1]+src[x+2])+src[x-2]+src[x+2]+16>>5);
        x++;
        buf[x<<1]=src[x];
        buf[x<<1|1]=OD_CLAMP255(20*(src[x]+src[x+1])-
         5*(src[x-1]+src[x+1])+src[x-2]+src[x+1]+16>>5);
        x++;
        buf[x<<1]=src[x];
        buf[x<<1|1]=OD_CLAMP255(36*src[x]-5*src[x-1]+src[x-2]+16>>5);
        x++;
        buf[x<<1]=src[w-1];
        buf[x<<1|1]=OD_CLAMP255(31*src[w-1]+src[w-2]+16>>5);
        memset(buf+(++x<<1),src[w-1],xpad-1<<1);
        /*for(x++;x<w+xpad;x++){
          buf[x<<1]=src[w-1];
          buf[x<<1|1]=src[w-1];
        }*/
        if(y>=0&&y+1<h)src+=siplane->ystride;
      }
      /*Vertical filtering:*/
      if(y>=-ypad+3){
        if(y<1||y>h+3){
          memcpy(dst-(xpad<<1),_state->ref_line_buf[y-3&7]-(xpad<<1),
           w+(xpad<<1)<<1);
          /*fprintf(stderr,"%3i: ",y-3<<1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
          memcpy(dst-(xpad<<1),_state->ref_line_buf[y-3&7]-(xpad<<1),
           w+(xpad<<1)<<1);
          /*fprintf(stderr,"%3i: ",y-3<<1|1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
        }
        else{
          unsigned char *buf[6];
          buf[0]=_state->ref_line_buf[y-5&7];
          buf[1]=_state->ref_line_buf[y-4&7];
          buf[2]=_state->ref_line_buf[y-3&7];
          buf[3]=_state->ref_line_buf[y-2&7];
          buf[4]=_state->ref_line_buf[y-1&7];
          buf[5]=_state->ref_line_buf[y-0&7];
          memcpy(dst-(xpad<<1),_state->ref_line_buf[y-3&7]-(xpad<<1),
           w+(xpad<<1)<<1);
          /*fprintf(stderr,"%3i: ",y-3<<1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
          for(x=-xpad<<1;x<w+xpad<<1;x++){
            *(dst+x)=OD_CLAMP255(20*(*(buf[2]+x)+*(buf[3]+x))-
             5*(*(buf[1]+x)+*(buf[4]+x))+
             *(buf[0]+x)+*(buf[5]+x)+16>>5);
          }
          /*fprintf(stderr,"%3i: ",y-3<<1|1);
          for(x=-xpad<<1;x<w+xpad<<1;x++)fprintf(stderr,"%02X",*(dst+x));
          fprintf(stderr,"\n");*/
          dst+=diplane->ystride;
        }
      }
    }
  }
}
#endif

/*The data used to build the following two arrays.*/
const int OD_VERT_D[22]={
/*0 1     4 5     8 9     12 13      17 18*/
  0,0,1,1,0,0,1,2,0,0,2,1,0,-1,1,1,0,-1,0,1,1,-1
};

/*The vector offsets in the X direction for each vertex from the upper-left,
   indexed by [exterior corner][split state][vertex].*/
const int *const OD_VERT_SETUP_DX[4][4]={
  {
    OD_VERT_D+ 9,OD_VERT_D+ 1,OD_VERT_D+ 9,OD_VERT_D+ 1,
  },{
    OD_VERT_D+13,OD_VERT_D+13,OD_VERT_D+ 1,OD_VERT_D+ 1,
  },{
    OD_VERT_D+18,OD_VERT_D+ 1,OD_VERT_D+18,OD_VERT_D+ 1,
  },{
    OD_VERT_D+ 5,OD_VERT_D+ 5,OD_VERT_D+ 1,OD_VERT_D+ 1,
  }
};


/*The vector offsets in the Y direction for each vertex from the upper-left,
   indexed by [exterior corner][split state][vertex].*/
const int *const OD_VERT_SETUP_DY[4][4]={
  {
    OD_VERT_DY+ 4,OD_VERT_DY+ 4,OD_VERT_DY+ 0,OD_VERT_DY +0,
  },{
    OD_VERT_DY+ 8,OD_VERT_DY+ 0,OD_VERT_DY+ 8,OD_VERT_DY +0,
  },{
    OD_VERT_DY+12,OD_VERT_DY+12,OD_VERT_DY+ 0,OD_VERT_DY +0,
  },{
    OD_VERT_DY+17,OD_VERT_DY+ 0,OD_VERT_DY+17,OD_VERT_DY +0,
  }
};


void od_state_pred_block_from_setup(od_state *_state,unsigned char *_buf,
 int _ystride,int _ref,int _pli,int _vx,int _vy,int _c,int _s,int _log_mvb_sz){
  od_img_plane  *iplane;
  od_mv_grid_pt *grid[4];
  ogg_int32_t    mvx[4];
  ogg_int32_t    mvy[4];
  const int     *dxp;
  const int     *dyp;
  int            etype;
  int            k;
  iplane=_state->ref_imgs[_state->ref_imgi[_ref]].planes+_pli;
  dxp=OD_VERT_SETUP_DX[_c][_s];
  dyp=OD_VERT_SETUP_DY[_c][_s];
  for(k=0;k<4;k++){
    grid[k]=_state->mv_grid[_vy+(dyp[k]<<_log_mvb_sz)]+
     _vx+(dxp[k]<<_log_mvb_sz);
    mvx[k]=(ogg_int32_t)grid[k]->mv[0]<<14-iplane->xdec;
    mvy[k]=(ogg_int32_t)grid[k]->mv[1]<<14-iplane->ydec;
  }
  /*Note: We pull the edge types from the grid points _after_ they've been
    offset to account for unsplit edges.
   This means that we do not rely on all the flags that would be implicitly set
    to the edge type of the unsplit edge to be up to date.
   This might also cause some other flags to be set incorrectly, but those are
    precisely the ones we override below.*/
  etype=grid[0]->right|grid[1]->down<<1|grid[3]->right<<2|grid[0]->down<<3;
  if(!(_s&1)){
    etype&=~(1<<(_c+1&3));
    /*if(_c==3)etype|=(etype&8)>>3;
    else etype|=(etype&1<<_c)<<1;*/
    etype|=(etype|etype<<4)>>3&1<<(_c+1&3);
    if(etype>>_c&1){
      mvx[_c+1&3]=mvx[_c]+mvx[_c+1&3]>>1;
      mvy[_c+1&3]=mvy[_c]+mvy[_c+1&3]>>1;
    }
  }
  if(!(_s&2)){
    etype&=~(1<<(_c+2&3));
    /*if(_c==1)etype|=(etype&1)<<3;
    else etype|=(etype&1<<(_c+3&3))>>1;*/
    etype|=(etype|etype<<4)>>1&1<<(_c+2&3);
    if(etype<<(-_c&3)&8){
      mvx[_c+3&3]=mvx[_c]+mvx[_c+3&3]>>1;
      mvy[_c+3&3]=mvy[_c]+mvy[_c+3&3]>>1;
    }
  }
  /*fprintf(stderr,"interpolation type: %c%c%c%c (0x%X) ",
   etype&1?'V':'B',etype&2?'V':'B',etype&4?'V':'B',etype&8?'V':'B',etype);*/
  od_mc_predict8(_state,_buf,_ystride,iplane->data+
   (_vy-2<<3-iplane->ydec)*iplane->ystride+(_vx-2<<3-iplane->xdec),
   iplane->ystride,mvx,mvy,etype,_c,_s,_log_mvb_sz+2-iplane->xdec,
   _log_mvb_sz+2-iplane->ydec);
}

void od_state_pred_block(od_state *_state,unsigned char *_buf,int _ystride,
 int _ref,int _pli,int _vx,int _vy,int _log_mvb_sz){
  int half_mvb_sz;
  half_mvb_sz=1<<_log_mvb_sz-1;
  if(_log_mvb_sz>0&&_state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
    od_img_plane *iplane;
    int           half_xblk_sz;
    int           half_yblk_sz;
    iplane=_state->ref_imgs[_state->ref_imgi[_ref]].planes+_pli;
    half_xblk_sz=1<<_log_mvb_sz+1-iplane->xdec;
    half_yblk_sz=1<<_log_mvb_sz+1-iplane->ydec;
    od_state_pred_block(_state,_buf,
     _ystride,_ref,_pli,_vx,_vy,_log_mvb_sz-1);
    od_state_pred_block(_state,_buf+half_xblk_sz,
     _ystride,_ref,_pli,_vx+half_mvb_sz,_vy,_log_mvb_sz-1);
    od_state_pred_block(_state,_buf+half_yblk_sz*_ystride,
     _ystride,_ref,_pli,_vx,_vy+half_mvb_sz,_log_mvb_sz-1);
    od_state_pred_block(_state,_buf+half_yblk_sz*_ystride+half_xblk_sz,
     _ystride,_ref,_pli,_vx+half_mvb_sz,_vy+half_mvb_sz,_log_mvb_sz-1);
  }
  else{
    int c;
    int s;
    if(_log_mvb_sz<2){
      int mask;
      mask=(1<<_log_mvb_sz+1)-1;
      c=!!(_vx&mask);
      if(_vy&mask)c=3-c;
      s=_state->mv_grid[_vy+(OD_VERT_DY[c+1&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+1&3]<<_log_mvb_sz)].valid|
       _state->mv_grid[_vy+(OD_VERT_DY[c+3&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+3&3]<<_log_mvb_sz)].valid<<1;
    }
    else{
      c=0;
      s=3;
    }
    od_state_pred_block_from_setup(_state,_buf,_ystride,_ref,_pli,
     _vx,_vy,c,s,_log_mvb_sz);
  }
}

int od_state_dump_yuv(od_state *_state,od_img *_img,const char *_suf){
  char  fname[128];
  FILE *fp;
  int   pic_width;
  int   pic_height;
  int   y;
  int   pli;
  static const char *CHROMA_TAGS[4]={" C420jpeg",""," C422jpeg"," C444"};
  sprintf(fname,"%08i%s.y4m",
   (int)daala_granule_basetime(_state,_state->cur_time),_suf);
  fp=fopen(fname,"wb");
  pic_width=_state->info.pic_width;
  pic_height=_state->info.pic_height;
  OD_ASSERT(_img->nplanes>=3);
  fprintf(fp,"YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   pic_width,pic_height,1,1,0,0,
   CHROMA_TAGS[(_img->planes[1].xdec==0)+(_img->planes[1].ydec==0)*2]);
  fprintf(fp,"FRAME\n");
  for(pli=0;pli<3;pli++){
    int xdec;
    int ydec;
    int ystride;
    xdec=_img->planes[pli].xdec;
    ydec=_img->planes[pli].ydec;
    ystride=_img->planes[pli].ystride;
    for(y=0;y<pic_height+ydec>>ydec;y++){
      if(fwrite(_img->planes[pli].data+ystride*y,
       (pic_width+xdec>>xdec),1,fp)<1){
        fprintf(stderr,"Error writing to \"%s\".\n","fixme");
        return EXIT_FAILURE;
      }
    }
  }
  fclose(fp);
  return 0;
}

#if defined(OD_DUMP_IMAGES)
# include <png.h>
#include <zlib.h>

/*State visualization.
  None of this is particularly fast.*/

/*static const unsigned char OD_YCbCr_BORDER[3]={ 49,109,184};*/
static const unsigned char OD_YCbCr_BORDER[3]={113, 72,137};
static const unsigned char OD_YCbCr_BEDGE[3]= { 41,240,110};
static const unsigned char OD_YCbCr_VEDGE[3]= {145, 54, 34};
static const unsigned char OD_YCbCr_MV[3]=    { 81, 90,240};

void od_img_draw_point(od_img *_img,int _x,int _y,
 const unsigned char _ycbcr[3]){
  if(_x<0||_y<0||_x>=_img->width||_y>=_img->height)return;
  *(_img->planes[0].data+_img->planes[0].ystride*(_y>>_img->planes[0].ydec)+
   (_x>>_img->planes[0].xdec))=_ycbcr[0];
  if(_img->nplanes>=3){
    *(_img->planes[1].data+_img->planes[1].ystride*(_y>>_img->planes[1].ydec)+
     (_x>>_img->planes[1].xdec))=_ycbcr[1];
    *(_img->planes[2].data+_img->planes[2].ystride*(_y>>_img->planes[2].ydec)+
     (_x>>_img->planes[2].xdec))=_ycbcr[2];
  }
}

void od_img_draw_line(od_img *_img,int _x0,int _y0,int _x1,int _y1,
 const unsigned char _ycbcr[3]){
  int x0[2];
  int x1[2];
  int dx[2];
  int step[2];
  int steep;
  int err;
  int derr;
  steep=abs(_y1-_y0)>abs(_x1-_x0);
  x0[0]=_x0;
  x0[1]=_y0;
  x1[0]=_x1;
  x1[1]=_y1;
  dx[0]=abs(x1[0]-x0[0]);
  dx[1]=abs(x1[1]-x0[1]);
  err=0;
  derr=dx[1-steep];
  step[0]=((x0[0]<x1[0])<<1)-1;
  step[1]=((x0[1]<x1[1])<<1)-1;
  od_img_draw_point(_img,x0[0],x0[1],_ycbcr);
  while(x0[steep]!=x1[steep]){
    x0[steep]+=step[steep];
    err+=derr;
    if(err<<1>dx[steep]){
      x0[1-steep]+=step[1-steep];
      err-=dx[steep];
    }
    od_img_draw_point(_img,x0[0],x0[1],_ycbcr);
  }
}

static void od_state_draw_mv_grid_block(od_state *_state,int _vx,int _vy,
 int _log_mvb_sz){
  int half_mvb_sz;
  half_mvb_sz=1<<_log_mvb_sz-1;
  if(_log_mvb_sz>0&&_state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
    od_state_draw_mv_grid_block(_state,_vx,_vy,_log_mvb_sz-1);
    od_state_draw_mv_grid_block(_state,_vx+half_mvb_sz,_vy,_log_mvb_sz-1);
    od_state_draw_mv_grid_block(_state,_vx,_vy+half_mvb_sz,_log_mvb_sz-1);
    od_state_draw_mv_grid_block(_state,_vx+half_mvb_sz,_vy+half_mvb_sz,
     _log_mvb_sz-1);
  }
  else{
    od_mv_grid_pt *grid[4];
    ogg_int32_t    mvx[4];
    ogg_int32_t    mvy[4];
    const int     *dxp;
    const int     *dyp;
    int            mvb_sz;
    int            etype;
    int            x0;
    int            y0;
    int            k;
    int            c;
    int            s;
    mvb_sz=1<<_log_mvb_sz;
    if(_log_mvb_sz<2){
      int mask;
      mask=(1<<_log_mvb_sz+1)-1;
      c=!!(_vx&mask);
      if(_vy&mask)c=3-c;
      s=_state->mv_grid[_vy+(OD_VERT_DY[c+1&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+1&3]<<_log_mvb_sz)].valid|
       _state->mv_grid[_vy+(OD_VERT_DY[c+3&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+3&3]<<_log_mvb_sz)].valid<<1;
    }
    else{
      c=0;
      s=3;
    }
    dxp=OD_VERT_SETUP_DX[c][s];
    dyp=OD_VERT_SETUP_DY[c][s];
    for(k=0;k<4;k++){
      grid[k]=_state->mv_grid[_vy+(dyp[k]<<_log_mvb_sz)]+
       _vx+(dxp[k]<<_log_mvb_sz);
      mvx[k]=grid[k]->mv[0];
      mvy[k]=grid[k]->mv[1];
    }
    etype=grid[0]->right|grid[1]->down<<1|grid[3]->right<<2|grid[0]->down<<3;
    if(!(s&1)){
      etype&=~(1<<(c+1&3));
      etype|=(etype|etype<<4)>>3&1<<(c+1&3);
      if(etype>>c&1){
        mvx[c+1&3]=mvx[c]+mvx[c+1&3]>>1;
        mvy[c+1&3]=mvy[c]+mvy[c+1&3]>>1;
      }
    }
    if(!(s&2)){
      etype&=~(1<<(c+2&3));
      etype|=(etype|etype<<4)>>1&1<<(c+2&3);
      if(etype<<(-c&3)&8){
        mvx[c+3&3]=mvx[c]+mvx[c+3&3]>>1;
        mvy[c+3&3]=mvy[c]+mvy[c+3&3]>>1;
      }
    }
    x0=(_vx-2<<3)+(OD_UMV_PADDING<<1);
    y0=(_vy-2<<3)+(OD_UMV_PADDING<<1);
    od_img_draw_line(&_state->vis_img,x0,y0,
     x0+(mvb_sz<<3),y0,etype&1?OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
    od_img_draw_line(&_state->vis_img,x0+(mvb_sz<<3),y0,
     x0+(mvb_sz<<3),y0+(mvb_sz<<3),etype&2?OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
    od_img_draw_line(&_state->vis_img,x0,y0+(mvb_sz<<3),
     x0+(mvb_sz<<3),y0+(mvb_sz<<3),etype&4?OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
    od_img_draw_line(&_state->vis_img,x0,y0,
     x0,y0+(mvb_sz<<3),etype&8?OD_YCbCr_VEDGE:OD_YCbCr_BEDGE);
  }
}

void od_state_draw_mv_grid(od_state *_state){
  int     vx;
  int     vy;
  int     nhmvbs;
  int     nvmvbs;
  nhmvbs=_state->nhmbs+1<<2;
  nvmvbs=_state->nvmbs+1<<2;
  for(vy=0;vy<nvmvbs;vy+=4)for(vx=0;vx<nhmvbs;vx+=4){
    od_state_draw_mv_grid_block(_state,vx,vy,2);
  }
}

static void od_state_draw_mvs_block(od_state *_state,int _vx,int _vy,
 int _log_mvb_sz){
  int half_mvb_sz;
  half_mvb_sz=1<<_log_mvb_sz-1;
  if(_log_mvb_sz>0&&_state->mv_grid[_vy+half_mvb_sz][_vx+half_mvb_sz].valid){
    od_state_draw_mvs_block(_state,_vx,_vy,_log_mvb_sz-1);
    od_state_draw_mvs_block(_state,_vx+half_mvb_sz,_vy,_log_mvb_sz-1);
    od_state_draw_mvs_block(_state,_vx,_vy+half_mvb_sz,_log_mvb_sz-1);
    od_state_draw_mvs_block(_state,_vx+half_mvb_sz,_vy+half_mvb_sz,
     _log_mvb_sz-1);
  }
  else{
    od_mv_grid_pt *grid[4];
    ogg_int32_t    mvx[4];
    ogg_int32_t    mvy[4];
    const int     *dxp;
    const int     *dyp;
    int            etype;
    int            x0;
    int            y0;
    int            k;
    int            c;
    int            s;
    if(_log_mvb_sz<2){
      int mask;
      mask=(1<<_log_mvb_sz+1)-1;
      c=!!(_vx&mask);
      if(_vy&mask)c=3-c;
      s=_state->mv_grid[_vy+(OD_VERT_DY[c+1&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+1&3]<<_log_mvb_sz)].valid|
       _state->mv_grid[_vy+(OD_VERT_DY[c+3&3]<<_log_mvb_sz)][
       _vx+(OD_VERT_DX[c+3&3]<<_log_mvb_sz)].valid<<1;
    }
    else{
      c=0;
      s=3;
    }
    dxp=OD_VERT_SETUP_DX[c][s];
    dyp=OD_VERT_SETUP_DY[c][s];
    for(k=0;k<4;k++){
      grid[k]=_state->mv_grid[_vy+(dyp[k]<<_log_mvb_sz)]+
       _vx+(dxp[k]<<_log_mvb_sz);
      mvx[k]=grid[k]->mv[0];
      mvy[k]=grid[k]->mv[1];
    }
    etype=grid[0]->right|grid[1]->down<<1|grid[3]->right<<2|grid[0]->down<<3;
    if(!(s&1)){
      etype&=~(1<<(c+1&3));
      etype|=(etype|etype<<4)>>3&1<<(c+1&3);
      if(etype>>c&1){
        mvx[c+1&3]=mvx[c]+mvx[c+1&3]>>1;
        mvy[c+1&3]=mvy[c]+mvy[c+1&3]>>1;
      }
    }
    if(!(s&2)){
      etype&=~(1<<(c+2&3));
      etype|=(etype|etype<<4)>>1&1<<(c+2&3);
      if(etype<<(-c&3)&8){
        mvx[c+3&3]=mvx[c]+mvx[c+3&3]>>1;
        mvy[c+3&3]=mvy[c]+mvy[c+3&3]>>1;
      }
    }
    x0=(_vx-2<<3)+(OD_UMV_PADDING<<1);
    y0=(_vy-2<<3)+(OD_UMV_PADDING<<1);
    for(k=0;k<4;k++){
      x0=(_vx-2+(dxp[k]<<_log_mvb_sz)<<3)+(OD_UMV_PADDING<<1);
      y0=(_vy-2+(dyp[k]<<_log_mvb_sz)<<3)+(OD_UMV_PADDING<<1);
      /*od_img_draw_point(&_state->vis_img,x0,y0,OD_YCbCr_MV);*/
      od_img_draw_line(&_state->vis_img,x0,y0,
       x0+OD_DIV_ROUND_POW2(grid[k]->mv[0],2,2),
       y0+OD_DIV_ROUND_POW2(grid[k]->mv[1],2,2),OD_YCbCr_MV);
    }
  }
}

void od_state_draw_mvs(od_state *_state){
  int     vx;
  int     vy;
  int     nhmvbs;
  int     nvmvbs;
  nhmvbs=_state->nhmbs+1<<2;
  nvmvbs=_state->nvmbs+1<<2;
  for(vy=0;vy<nvmvbs;vy+=4)for(vx=0;vx<nhmvbs;vx+=4){
    od_state_draw_mvs_block(_state,vx,vy,2);
  }
}

void od_state_fill_vis(od_state *_state){
  od_img *img;
  od_img *ref_img;
  int     pli;
  int     x;
  int     y;
  img=&_state->vis_img;
  /*Upsample the reconstructed image for better quality.*/
  /*Adjust the data pointers so that the padding works like the reference
     images.*/
  for(pli=0;pli<img->nplanes;pli++){
    img->planes[pli].data+=(OD_UMV_PADDING<<1>>img->planes[pli].xdec)+
     img->planes[pli].ystride*(OD_UMV_PADDING<<1>>img->planes[pli].ydec);
  }
  od_state_upsample8(_state,img,_state->io_imgs+OD_FRAME_REC);
  /*Upsample the input image, as well, and subtract it to get a difference
     image.*/
  ref_img=_state->ref_imgs+_state->ref_imgi[OD_FRAME_SELF];
  od_state_upsample8(_state,ref_img,&_state->input);
  for(y=0;y<ref_img->height;y++){
    for(x=0;x<ref_img->width;x++){
      int diff;
      diff=*(img->planes[0].data+
       img->planes[0].ystride*(y>>img->planes[0].ydec)+
       (x>>img->planes[0].xdec))-
       *(ref_img->planes[0].data+
       ref_img->planes[0].ystride*(y>>ref_img->planes[0].ydec)+
       (x>>ref_img->planes[0].xdec));
      /*Scale the differences by 2 to make them visible.*/
      diff=OD_CLAMP255((diff<<1)+128);
      *(img->planes[0].data+img->planes[0].ystride*(y>>img->planes[0].ydec)+
       (x>>img->planes[0].xdec))=(unsigned char)diff;
    }
  }
  /*Undo the adjustment.*/
  for(pli=0;pli<img->nplanes;pli++){
    img->planes[pli].data-=(OD_UMV_PADDING<<1>>img->planes[pli].xdec)+
     img->planes[pli].ystride*(OD_UMV_PADDING<<1>>img->planes[pli].ydec);
  }
  /*Clear the border region.*/
  for(y=0;y<(OD_UMV_PADDING<<1>>img->planes[0].ydec);y++){
    memset(img->planes[0].data+(img->planes[0].ystride)*y,0,
     img->width>>img->planes[0].xdec);
  }
  for(;y<img->height-(OD_UMV_PADDING<<1)>>img->planes[0].ydec;y++){
    memset(img->planes[0].data+img->planes[0].ystride*y,0,
     OD_UMV_PADDING<<1>>img->planes[0].xdec);
    memset(img->planes[0].data+img->planes[0].ystride*y+
     (img->width-(OD_UMV_PADDING<<1)>>img->planes[0].xdec),0,
     OD_UMV_PADDING<<1>>img->planes[0].xdec);
  }
  for(;y<img->height>>img->planes[0].ydec;y++){
    memset(img->planes[0].data+(img->planes[0].ystride)*y,0,
     img->width>>img->planes[0].xdec);
  }
  /*Clear the chroma planes.*/
  for(pli=1;pli<img->nplanes;pli++){
    memset(img->planes[pli].data,128,
     (img->height>>img->planes[pli].ydec)*(img->width>>img->planes[pli].xdec));
  }
  od_img_draw_line(img,(OD_UMV_PADDING<<1)-1,(OD_UMV_PADDING<<1)-1,
   img->width-(OD_UMV_PADDING<<1),(OD_UMV_PADDING<<1)-1,OD_YCbCr_BORDER);
  od_img_draw_line(img,(OD_UMV_PADDING<<1)-2,(OD_UMV_PADDING<<1)-2,
   img->width-(OD_UMV_PADDING<<1)+1,(OD_UMV_PADDING<<1)-2,OD_YCbCr_BORDER);
  od_img_draw_line(img,img->width-(OD_UMV_PADDING<<1),(OD_UMV_PADDING<<1)-1,
   img->width-(OD_UMV_PADDING<<1),img->height-(OD_UMV_PADDING<<1),
   OD_YCbCr_BORDER);
  od_img_draw_line(img,img->width-(OD_UMV_PADDING<<1)+1,(OD_UMV_PADDING<<1)-2,
   img->width-(OD_UMV_PADDING<<1)+1,img->height-(OD_UMV_PADDING<<1)+1,
   OD_YCbCr_BORDER);
  od_img_draw_line(img,(OD_UMV_PADDING<<1)-1,img->height-(OD_UMV_PADDING<<1),
   img->width-(OD_UMV_PADDING<<1),img->height-(OD_UMV_PADDING<<1),
   OD_YCbCr_BORDER);
  od_img_draw_line(img,(OD_UMV_PADDING<<1)-2,img->height-(OD_UMV_PADDING<<1)+1,
   img->width-(OD_UMV_PADDING<<1)+1,img->height-(OD_UMV_PADDING<<1)+1,
   OD_YCbCr_BORDER);
  od_img_draw_line(img,(OD_UMV_PADDING<<1)-1,(OD_UMV_PADDING<<1)-1,
   (OD_UMV_PADDING<<1)-1,img->height-(OD_UMV_PADDING<<1),OD_YCbCr_BORDER);
  od_img_draw_line(img,(OD_UMV_PADDING<<1)-2,(OD_UMV_PADDING<<1)-2,
   (OD_UMV_PADDING<<1)-2,img->height-(OD_UMV_PADDING<<1)+1,OD_YCbCr_BORDER);
  od_state_draw_mv_grid(_state);
  od_state_draw_mvs(_state);
}

/*Dump a PNG of the reconstructed image, or a reference frame.*/
int od_state_dump_img(od_state *_state,od_img *_img,const char *_suf){
  png_structp    png;
  png_infop      info;
  png_bytep     *data;
  FILE          *fp;
  char           fname[128];
  unsigned char *p_rows[3];
  unsigned char *p[3];
  int            nplanes;
  int            pli;
  int            x;
  int            y;
  sprintf(fname,"%08i%s.png",
   (int)daala_granule_basetime(_state,_state->cur_time),_suf);
  fp=fopen(fname,"wb");
  png=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
  if(png==NULL){
    fclose(fp);
    return OD_EFAULT;
  }
  info=png_create_info_struct(png);
  if(info==NULL){
    png_destroy_write_struct(&png,NULL);
    fclose(fp);
    return OD_EFAULT;
  }
  if(setjmp(png_jmpbuf(png))){
    png_destroy_write_struct(&png,&info);
    fclose(fp);
    return OD_EFAULT;
  }
  data=(png_bytep *)od_malloc_2d(_img->height,6*_img->width,sizeof(data[0][0]));
  if(_img->nplanes<3)nplanes=1;
  else nplanes=3;
  for(pli=0;pli<nplanes;pli++)p_rows[pli]=_img->planes[pli].data;
  /*Chroma up-sampling is just done with a box filter.
    This is very likely what will actually be used in practice on a real
     display, and also removes one more layer to search in for the source of
     artifacts.
    As an added bonus, it's dead simple.*/
  for(y=0;y<_img->height;y++){
    int mask;
    /*LOOP VECTORIZES.*/
    for(pli=0;pli<nplanes;pli++)p[pli]=p_rows[pli];
    for(x=0;x<_img->width;x++){
      float    yval;
      float    cbval;
      float    crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      /*This is intentionally slow and very accurate.*/
      yval=(p[0][0]-16)*(1.0F/219);
      if(nplanes>=3){
        cbval=(p[1][0]-128)*(2*(1-0.114F)/224);
        crval=(p[2][0]-128)*(2*(1-0.299F)/224);
      }
      else cbval=crval=0;
      rval=OD_CLAMPI(0,(int)(65535*(yval+crval)+0.5F),65535);
      gval=OD_CLAMPI(0,(int)(65535*
       (yval-cbval*(0.114F/0.587F)-crval*(0.299F/0.587F))+0.5F),65535);
      bval=OD_CLAMPI(0,(int)(65535*(yval+cbval)+0.5F),65535);
      data[y][6*x+0]=(unsigned char)(rval>>8);
      data[y][6*x+1]=(unsigned char)(rval&0xFF);
      data[y][6*x+2]=(unsigned char)(gval>>8);
      data[y][6*x+3]=(unsigned char)(gval&0xFF);
      data[y][6*x+4]=(unsigned char)(bval>>8);
      data[y][6*x+5]=(unsigned char)(bval&0xFF);
      for(pli=0;pli<nplanes;pli++){
        mask=(1<<_img->planes[pli].xdec)-1;
        p[pli]+=(x&mask)==mask;
      }
    }
    for(pli=0;pli<nplanes;pli++){
      mask=(1<<_img->planes[pli].ydec)-1;
      p_rows[pli]+=((y&mask)==mask)*_img->planes[pli].ystride;
    }
  }
  png_init_io(png,fp);
  png_set_compression_level(png,Z_BEST_COMPRESSION);
  png_set_IHDR(png,info,_img->width,_img->height,16,PNG_COLOR_TYPE_RGB,
   PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
  /*TODO: Define real colorspace.*/
  /*png_set_gAMA(png,info,2.2);
  png_set_cHRM_fixed(png,info,31006,31616,67000,32000,21000,71000,14000,8000);*/
  png_set_pHYs(png,info,_state->info.pixel_aspect_numerator,
   _state->info.pixel_aspect_denominator,0);
  png_set_rows(png,info,data);
  png_write_png(png,info,PNG_TRANSFORM_IDENTITY,NULL);
  png_write_end(png,info);
  png_destroy_write_struct(&png,&info);
  od_free_2d(data);
  fclose(fp);
  return 0;
}
#endif

void od_state_mc_predict(od_state *state, int ref) {
  unsigned char  __attribute__((aligned(16))) buf[16][16];
  od_img *img;
  int nhmvbs;
  int nvmvbs;
  int pli;
  int vx;
  int vy;
  nhmvbs = (state->nhmbs + 1) << 2;
  nvmvbs = (state->nvmbs + 1) << 2;
  img = state->io_imgs + OD_FRAME_REC;
  for (vy = 0; vy < nvmvbs; vy += 4) {
    for (vx = 0; vx < nhmvbs; vx += 4) {
      for (pli = 0; pli < img->nplanes; pli++) {
        od_img_plane *iplane;
        unsigned char *p;
        int blk_w;
        int blk_h;
        int blk_x;
        int blk_y;
        int y;
        od_state_pred_block(state, buf[0], sizeof(buf[0]), ref, pli, vx, vy,
         2);
        /*Copy the predictor into the image, with clipping.*/
        iplane = img->planes + pli;
        blk_w = 16 >> iplane->xdec;
        blk_h = 16 >> iplane->ydec;
        blk_x = (vx - 2) << (2 - iplane->xdec);
        blk_y = (vy - 2) << (2 - iplane->ydec);
        p = buf[0];
        if (blk_x < 0) {
          blk_w += blk_x;
          p -= blk_x;
          blk_x = 0;
        }
        if (blk_y < 0) {
          blk_h += blk_y;
          p -= blk_y*sizeof(buf[0]);
          blk_y = 0;
        }
        if (blk_x + blk_w > img->width >> iplane->xdec) {
          blk_w = (img->width >> iplane->xdec) - blk_x;
        }
        if (blk_y + blk_h > img->height >> iplane->ydec) {
          blk_h = (img->height >> iplane->ydec) - blk_y;
        }
        for (y = blk_y; y < blk_y + blk_h; y++) {
          memcpy(iplane->data + y*iplane->ystride + blk_x,
           p, blk_w);
          p += sizeof(buf[0]);
        }
      }
    }
  }
}


ogg_int64_t daala_granule_basetime(void *_encdec,ogg_int64_t _granpos){
  od_state *state;
  state=(od_state *)_encdec;
  if(_granpos>=0){
    ogg_int64_t key_time;
    ogg_int64_t delta_time;
    key_time=_granpos>>state->info.keyframe_granule_shift;
    delta_time=_granpos-(key_time<<state->info.keyframe_granule_shift);
    return key_time+delta_time;
  }
  return -1;
}

double daala_granule_time(void *_encdec,ogg_int64_t _granpos){
  od_state    *state;
  ogg_int64_t  base_time;
  state=(od_state *)_encdec;
  base_time=daala_granule_basetime(_encdec,_granpos);
  if(base_time>=0){
    return base_time*(double)state->info.timebase_denominator/
     state->info.timebase_numerator;
  }
  return -1;
}

void od_extract_bsize(unsigned char *_bsize_out,int _bstride_out,
 const unsigned char *_bsize_in,int _bstride_in,int _dec){
  int j;
  int i;
  for (j=0;j<1<<(2-_dec);j++) {
    _bsize_out[_bstride_out*(j+1)]=
     OD_CLAMPI(0,_bsize_in[_bstride_in*(1<<_dec)*j-1]-_dec,2);
    _bsize_out[_bstride_out*(j+1)+(1<<(2-_dec))+1]=
     OD_CLAMPI(0,_bsize_in[_bstride_in*(1<<_dec)*j+4]-_dec,2);
  }
  for (i=0;i<1<<(2-_dec);i++) {
    _bsize_out[i+1]=
     OD_CLAMPI(0,_bsize_in[_bstride_in*(0-1)+(1<<_dec)*i]-_dec,2);
    _bsize_out[_bstride_out*((1<<(2-_dec))+1)+i+1]=
     OD_CLAMPI(0,_bsize_in[_bstride_in*(0+4)+(1<<_dec)*i]-_dec,2);
  }
  for (j=0;j<1<<(2-_dec);j++) {
    for (i=0;i<1<<(2-_dec);i++) {
      _bsize_out[_bstride_out*(j+1)+i+1]=
       OD_CLAMPI(0,_bsize_in[_bstride_in*(1<<_dec)*j+(1<<_dec)*i]-_dec,2);
    }
  }
}
