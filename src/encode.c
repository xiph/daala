/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

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

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "encint.h"
#include "generic_code.h"
#include "filter.h"
#include "dct.h"
#include "intra.h"
#include "pvq.h"
#include "pvq_code.h"
#include "tf.h"

static int od_enc_init(od_enc_ctx *_enc,const daala_info *_info){
  int ret;
  ret=od_state_init(&_enc->state,_info);
  if(ret<0)return ret;
  oggbyte_writeinit(&_enc->obb);
  od_ec_enc_init(&_enc->ec,65025);
  _enc->packet_state=OD_PACKET_INFO_HDR;
  _enc->mvest=od_mv_est_alloc(_enc);
  return 0;
}

static void od_enc_clear(od_enc_ctx *_enc){
  od_mv_est_free(_enc->mvest);
  od_ec_enc_clear(&_enc->ec);
  oggbyte_writeclear(&_enc->obb);
  od_state_clear(&_enc->state);
}

daala_enc_ctx *daala_encode_alloc(const daala_info *_info){
  od_enc_ctx *enc;
  if(_info==NULL)return NULL;
  enc=(od_enc_ctx *)_ogg_malloc(sizeof(*enc));
  if(od_enc_init(enc,_info)<0){
    _ogg_free(enc);
    return NULL;
  }
  return enc;
}

void daala_encode_free(daala_enc_ctx *_enc){
  if(_enc!=NULL){
    od_enc_clear(_enc);
    _ogg_free(_enc);
  }
}

int daala_encode_ctl(daala_enc_ctx *_enc,int _req,void *_buf,size_t _buf_sz){
  (void)_enc;
  (void)_buf;
  (void)_buf_sz;
  switch(_req){
    default:return OD_EIMPL;
  }
}

void od_state_mc_predict(od_state *_state,int _ref){
  unsigned char  __attribute__((aligned(16))) buf[16][16];
  od_img        *img;
  int            nhmvbs;
  int            nvmvbs;
  int            pli;
  int            vx;
  int            vy;
  nhmvbs=_state->nhmbs+1<<2;
  nvmvbs=_state->nvmbs+1<<2;
  img=_state->io_imgs+OD_FRAME_REC;
  for(vy=0;vy<nvmvbs;vy+=4){
    for(vx=0;vx<nhmvbs;vx+=4){
      for(pli=0;pli<img->nplanes;pli++){
        od_img_plane  *iplane;
        unsigned char *p;
        int            blk_w;
        int            blk_h;
        int            blk_x;
        int            blk_y;
        int            y;
        od_state_pred_block(_state,buf[0],sizeof(buf[0]),_ref,pli,vx,vy,2);
        /*Copy the predictor into the image, with clipping.*/
        iplane=img->planes+pli;
        blk_w=16>>iplane->xdec;
        blk_h=16>>iplane->ydec;
        blk_x=vx-2<<2-iplane->xdec;
        blk_y=vy-2<<2-iplane->ydec;
        p=buf[0];
        if(blk_x<0){
          blk_w+=blk_x;
          p-=blk_x;
          blk_x=0;
        }
        if(blk_y<0){
          blk_h+=blk_y;
          p-=blk_y*sizeof(buf[0]);
          blk_y=0;
        }
        if(blk_x+blk_w>img->width>>iplane->xdec){
          blk_w=(img->width>>iplane->xdec)-blk_x;
        }
        if(blk_y+blk_h>img->height>>iplane->ydec){
          blk_h=(img->height>>iplane->ydec)-blk_y;
        }
        for(y=blk_y;y<blk_y+blk_h;y++){
          memcpy(iplane->data+y*iplane->ystride+blk_x,p,blk_w);
          p+=sizeof(buf[0]);
        }
      }
    }
  }
}

static void od_img_plane_copy_pad8(od_img_plane *_dst,
 int _plane_width,int _plane_height,od_img_plane *_src,
 int _pic_x,int _pic_y,int _pic_width,int _pic_height){
  unsigned char *dst_data;
  ptrdiff_t      dstride;
  int            y;
  dstride=_dst->ystride;
  /*If we have _no_ data, just encode a dull green.*/
  if(_pic_width==0||_pic_height==0){
    dst_data=_dst->data;
    for(y=0;y<_plane_height;y++){
      memset(dst_data,0,_plane_width*sizeof(*dst_data));
      dst_data+=dstride;
    }
  }
  /*Otherwise, copy what we do have, and add our own padding.*/
  else{
    unsigned char *src_data;
    unsigned char *src;
    unsigned char *dst;
    ptrdiff_t      sxstride;
    ptrdiff_t      systride;
    int            x;
    /*Step 1: Copy the data we do have.*/
    sxstride=_src->xstride;
    systride=_src->ystride;
    dst_data=_dst->data;
    src_data=_src->data;
    dst=dst_data+dstride*_pic_y+_pic_x;
    src=src_data;
    for(y=0;y<_pic_height;y++){
      if(sxstride==1)memcpy(dst,src,_pic_width);
      else for(x=0;x<_pic_width;x++)dst[x]=*(src+sxstride*x);
      dst+=dstride;
      src+=systride;
    }
    /*Step 2: Perform a low-pass extension into the padding region.*/
    /*Left side.*/
    for(x=_pic_x;x-->0;){
      dst=dst_data+dstride*_pic_y+x;
      for(y=0;y<_pic_height;y++){
        dst[0]=2*dst[1]+(dst-(dstride&-(y>0)))[1]
         +(dst+(dstride&-(y+1<_pic_height)))[1]+2>>2;
        dst+=dstride;
      }
    }
    /*Right side.*/
    for(x=_pic_x+_pic_width;x<_plane_width;x++){
      dst=dst_data+dstride*_pic_y+x-1;
      for(y=0;y<_pic_height;y++){
        dst[1]=2*dst[0]+(dst-(dstride&-(y>0)))[0]
         +(dst+(dstride&-(y+1<_pic_height)))[0]+2>>2;
        dst+=dstride;
      }
    }
    /*Top.*/
    dst=dst_data+dstride*_pic_y;
    for(y=_pic_y;y-->0;){
      for(x=0;x<_plane_width;x++){
        (dst-dstride)[x]=2*dst[x]+dst[x-(x>0)]+dst[x+(x+1<_plane_width)]+2>>2;
      }
      dst-=dstride;
    }
    /*Bottom.*/
    dst=dst_data+dstride*(_pic_y+_pic_height);
    for(y=_pic_y+_pic_height;y<_plane_height;y++){
      for(x=0;x<_plane_width;x++){
        dst[x]=2*(dst-dstride)[x]+(dst-dstride)[x-(x>0)]
         +(dst-dstride)[x+(x+1<_plane_width)]+2>>2;
      }
      dst+=dstride;
    }
  }
}

static void od_resample_luma_coeffs(od_coeff *_l,int _lstride,
 const od_coeff *_c,int _cstride,int _xdec,int _ydec,int _n){
  if(_xdec){
    if(_ydec)od_tf_up_hv_lp(_l,_lstride,_c,_cstride,_n,_n,_n);
    else od_tf_up_h_lp(_l,_lstride,_c,_cstride,_n,_n);
  }
  else{
    OD_ASSERT(_ydec);
    if(_ydec)od_tf_up_v_lp(_l,_lstride,_c,_cstride,_n,_n);
  }
}

static double mode_bits=0;
static double mode_count=0;

int daala_encode_img_in(daala_enc_ctx *_enc,od_img *_img,int _duration){
  int refi;
  int nplanes;
  int pli;
  int scale;
  int frame_width;
  int frame_height;
  int pic_x;
  int pic_y;
  int pic_width;
  int pic_height;
  if(_enc==NULL||_img==NULL)return OD_EFAULT;
  if(_enc->packet_state==OD_PACKET_DONE)return OD_EINVAL;
  /*Check the input image dimensions to make sure they're compatible with the
     declared video size.*/
  nplanes=_enc->state.info.nplanes;
  if(_img->nplanes!=nplanes)return OD_EINVAL;
  for(pli=0;pli<nplanes;pli++){
    if(_img->planes[pli].xdec!=_enc->state.info.plane_info[pli].xdec||
     _img->planes[pli].ydec!=_enc->state.info.plane_info[pli].ydec){
      return OD_EINVAL;
    }
  }
  frame_width=_enc->state.info.frame_width;
  frame_height=_enc->state.info.frame_height;
  pic_x=_enc->state.info.pic_x;
  pic_y=_enc->state.info.pic_y;
  pic_width=_enc->state.info.pic_width;
  pic_height=_enc->state.info.pic_height;
  if(_img->width!=frame_width||_img->height!=frame_height){
    /*The buffer does not match the frame size.
      Check to see if it matches the picture size.*/
    if(_img->width!=pic_width||_img->height!=pic_height){
      /*It doesn't; we don't know how to handle it yet.*/
      return OD_EINVAL;
    }
  }
  /*Copy and pad the image.*/
  for(pli=0;pli<nplanes;pli++){
    od_img_plane plane;
    int          plane_x;
    int          plane_y;
    int          plane_width;
    int          plane_height;
    *&plane=*(_img->planes+pli);
    plane_x=pic_x>>plane.xdec;
    plane_y=pic_y>>plane.ydec;
    plane_width=(pic_x+pic_width+(1<<plane.xdec)-1>>plane.xdec)-plane_x;
    plane_height=(pic_y+pic_height+(1<<plane.ydec)-1>>plane.ydec)-plane_y;
    if(_img->width==frame_width&&_img->height==frame_height){
      plane.data-=plane.ystride*(ptrdiff_t)plane_y
       +plane.xstride*(ptrdiff_t)plane_x;
    }
    od_img_plane_copy_pad8(_enc->state.io_imgs[OD_FRAME_INPUT].planes+pli,
     frame_width>>plane.xdec,frame_height>>plane.ydec,
     &plane,plane_x,plane_y,plane_width,plane_height);
  }
  /*Initialize the entropy coder.*/
  od_ec_enc_reset(&_enc->ec);
  /*Update the buffer state.*/
  if(_enc->state.ref_imgi[OD_FRAME_SELF]>=0){
    _enc->state.ref_imgi[OD_FRAME_PREV]=
     _enc->state.ref_imgi[OD_FRAME_SELF];
    /*TODO: Update golden frame.*/
    if(_enc->state.ref_imgi[OD_FRAME_GOLD]<0){
      _enc->state.ref_imgi[OD_FRAME_GOLD]=_enc->state.ref_imgi[OD_FRAME_SELF];
      /*TODO: Mark keyframe timebase.*/
    }
  }
  /*Select a free buffer to use for this reference frame.*/
  for(refi=0;refi==_enc->state.ref_imgi[OD_FRAME_GOLD]||
   refi==_enc->state.ref_imgi[OD_FRAME_PREV]||
   refi==_enc->state.ref_imgi[OD_FRAME_NEXT];refi++);
  _enc->state.ref_imgi[OD_FRAME_SELF]=refi;
  memcpy(&_enc->state.input,_img,sizeof(_enc->state.input));
  /*TODO: Incrment frame count.*/
  if(_enc->state.ref_imgi[OD_FRAME_PREV]>=0/*&&
   daala_granule_basetime(_enc,_enc->state.cur_time)>=19*/){
#if defined(OD_DUMP_IMAGES)&&defined(OD_ANIMATE)
    _enc->state.ani_iter=0;
#endif
    fprintf(stderr,"Predicting frame %i:\n",
     (int)daala_granule_basetime(_enc,_enc->state.cur_time));
#if 0
    od_mv_est(_enc->mvest,OD_FRAME_PREV,452/*118*/);
#endif
    od_state_mc_predict(&_enc->state,OD_FRAME_PREV);
#if defined(OD_DUMP_IMAGES)
    /*Dump reconstructed frame.*/
/*    od_state_dump_img(&_enc->state,_enc->state.io_imgs+OD_FRAME_REC,"rec");*/
    od_state_fill_vis(&_enc->state);
    od_state_dump_img(&_enc->state,&_enc->state.vis_img,"vis");
#endif
  }
  scale=10;/*atoi(getenv("QUANT"));*/
  /*TODO: Encode image.*/
  {
    GenericEncoder model_dc[OD_NPLANES_MAX];
    GenericEncoder model_g[OD_NPLANES_MAX];
    GenericEncoder model_ym[OD_NPLANES_MAX];
    ogg_uint16_t   mode_p0[OD_INTRA_NMODES];
    int            ex_dc[OD_NPLANES_MAX];
    int            ex_g[OD_NPLANES_MAX];
    od_coeff      *ctmp[OD_NPLANES_MAX];
    od_coeff      *dtmp[OD_NPLANES_MAX];
    od_coeff      *ltmp[OD_NPLANES_MAX];
    od_coeff      *lbuf[OD_NPLANES_MAX];
    signed char   *modes;
    int            xdec;
    int            ydec;
    int            nvmbs;
    int            nhmbs;
    int            mby;
    int            mbx;
    int            iyfill;
    int            oyfill;
    int            mi;
    int            h;
    int            w;
    int            y;
    int            x;
    nhmbs=_enc->state.nhmbs;
    nvmbs=_enc->state.nvmbs;
    /*Initialize the data needed for each plane.*/
    modes=_ogg_calloc((frame_width>>2)*(frame_height>>2),sizeof(*modes));
    for(mi=0;mi<OD_INTRA_NMODES;mi++)mode_p0[mi]=32768/OD_INTRA_NMODES;
    for(pli=0;pli<nplanes;pli++){
      od_pvq_adapt_ctx *pvq_adapt_row;
      generic_model_init(model_dc+pli);
      generic_model_init(model_g+pli);
      generic_model_init(model_ym+pli);
      ex_dc[pli]=pli>0?8:32768;
      ex_g[pli]=8;
      xdec=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w=frame_width>>xdec;
      h=frame_height>>ydec;
      ctmp[pli]=_ogg_calloc(w*h,sizeof(*ctmp[pli]));
      dtmp[pli]=_ogg_calloc(w*h,sizeof(*dtmp[pli]));
      /*We predict chroma planes from the luma plane.
        Since chroma can be subsampled, we cache subsampled versions of the
         luma plane in the frequency domain.
        We can share buffers with the same subsampling.*/
      if(pli>0){
        int plj;
        if(xdec||ydec){
          for(plj=1;plj<pli;plj++){
            if(xdec==_enc->state.io_imgs[OD_FRAME_INPUT].planes[plj].xdec
             &&ydec==_enc->state.io_imgs[OD_FRAME_INPUT].planes[plj].ydec){
              ltmp[pli]=NULL;
              lbuf[pli]=ltmp[plj];
            }
          }
          if(plj>=pli)lbuf[pli]=ltmp[pli]=_ogg_calloc(w*h,sizeof(*ltmp));
        }
        else{
          ltmp[pli]=NULL;
          lbuf[pli]=ctmp[pli];
        }
      }
      else lbuf[pli]=ltmp[pli]=NULL;
      pvq_adapt_row=_enc->state.pvq_adapt_row[pli];
      for(mbx=0;mbx<nhmbs;mbx++){
        pvq_adapt_row[mbx].mean_k_q8=OD_K_ROW_INIT_Q8;
        pvq_adapt_row[mbx].mean_sum_ex_q8=OD_SUM_EX_ROW_INIT_Q8;
        pvq_adapt_row[mbx].mean_count_q8=OD_COUNT_ROW_INIT_Q8;
        pvq_adapt_row[mbx].mean_count_ex_q8=OD_COUNT_EX_ROW_INIT_Q8;
#if !defined(OD_DISABLE_PVQ_CODE1)
        pvq_adapt_row[mbx].mean_pos_q4=OD_POS_ROW_INIT_Q4;
#endif
      }
      od_adapt_row_init(&_enc->state.adapt_row[pli]);
    }
    iyfill=0;
    oyfill=0;
    for(mby=0;mby<nvmbs;mby++){
      int next_iyfill;
      int next_oyfill;
      int ixfill;
      int oxfill;
      od_adapt_ctx adapt_hmean[OD_NPLANES_MAX];
      int hmean_k_q7[OD_NPLANES_MAX];
      int hmean_sum_ex_q7[OD_NPLANES_MAX];
      int hmean_count_q7[OD_NPLANES_MAX];
      int hmean_count_ex_q7[OD_NPLANES_MAX];
#if !defined(OD_DISABLE_PVQ_CODE1)
      int hmean_pos_q3[OD_NPLANES_MAX];
#endif
      for(pli=0;pli<nplanes;pli++){
        hmean_k_q7[pli]=OD_K_INIT_Q7;
        hmean_sum_ex_q7[pli]=OD_SUM_EX_INIT_Q7;
        hmean_count_q7[pli]=OD_COUNT_INIT_Q7;
        hmean_count_ex_q7[pli]=OD_COUNT_EX_INIT_Q7;
#if !defined(OD_DISABLE_PVQ_CODE1)
        hmean_pos_q3[pli]=OD_POS_INIT_Q3;
#endif
        od_adapt_hmean_init(&adapt_hmean[pli]);
      }
      next_iyfill=mby+1<nvmbs?(mby+1<<4)+8:frame_height;
      next_oyfill=mby+1<nvmbs?(mby+1<<4)-8:frame_height;
      ixfill=0;
      oxfill=0;
      for(mbx=0;mbx<nhmbs;mbx++){
        int      next_ixfill;
        int      next_oxfill;
        next_ixfill=mbx+1<nhmbs?(mbx+1<<4)+8:frame_width;
        next_oxfill=mbx+1<nhmbs?(mbx+1<<4)-8:frame_width;
        for(pli=0;pli<nplanes;pli++){
          od_pvq_adapt_ctx *pvq_adapt_row;
          od_pvq_adapt_ctx  pvq_adapt;
          od_adapt_row_ctx *adapt_row;
          od_adapt_ctx      adapt;
          od_coeff         *c;
          od_coeff         *d;
          od_coeff         *l;
          unsigned char    *data;
          int               ystride;
          int               by;
          int               bx;
          int               nk;
          int               k_total;
          int               sum_ex_total_q8;
          int               ncount;
          int               count_total_q8;
          int               count_ex_total_q8;
#if !defined(OD_DISABLE_PVQ_CODE1)
          int               npos;
          int               pos_total;
#endif
          c=ctmp[pli];
          d=dtmp[pli];
          l=lbuf[pli];
          xdec=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
          ydec=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
          w=frame_width>>xdec;
          h=frame_height>>ydec;
          /*Construct the luma predictors for chroma planes.*/
          if(ltmp[pli]!=NULL){
            OD_ASSERT(pli>0);
            OD_ASSERT(l==ltmp[pli]);
            for(by=mby<<2-ydec;by<mby+1<<2-ydec;by++){
              for(bx=mbx<<2-xdec;bx<mbx+1<<2-xdec;bx++){
                od_resample_luma_coeffs(l+(by<<2)*w+(bx<<2),w,
                 dtmp[0]+(by<<2+ydec)*frame_width+(bx<<2+xdec),frame_width,
                 xdec,ydec,4);
              }
            }
          }
          /*Collect the image data needed for this macro block.*/
          data=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].data;
          ystride=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
          for(y=iyfill>>ydec;y<next_iyfill>>ydec;y++){
            for(x=ixfill>>xdec;x<next_ixfill>>xdec;x++){
              c[y*w+x]=data[ystride*y+x]-128;
            }
          }
          /*Apply the prefilter across the bottom block edges.*/
          for(by=mby<<2-ydec;by<(mby+1<<2-ydec)-(mby+1>=nvmbs);by++){
            for(x=ixfill>>xdec;x<next_ixfill>>xdec;x++){
              od_coeff p[4];
              for(y=0;y<4;y++)p[y]=c[((by<<2)+y+2)*w+x];
              od_pre_filter4(p,p);
              for(y=0;y<4;y++)c[((by<<2)+y+2)*w+x]=p[y];
            }
          }
          /*Apply the prefilter across the right block edges.*/
          for(y=iyfill>>ydec;y<next_iyfill>>ydec;y++){
            for(bx=mbx<<2-xdec;bx<(mbx+1<<2-xdec)-(mbx+1>=nhmbs);bx++){
              od_pre_filter4(c+y*w+(bx<<2)+2,c+y*w+(bx<<2)+2);
            }
          }
          pvq_adapt_row=_enc->state.pvq_adapt_row[pli]+mbx;
          pvq_adapt.mean_k_q8=(hmean_k_q7[pli]>>OD_K_ADAPT_SPEED)
           +pvq_adapt_row->mean_k_q8;
          pvq_adapt.mean_sum_ex_q8=
           (hmean_sum_ex_q7[pli]>>OD_SUM_EX_ADAPT_SPEED)
           +pvq_adapt_row->mean_sum_ex_q8;
          nk=k_total=sum_ex_total_q8=0;
          pvq_adapt.mean_count_q8=(hmean_count_q7[pli]>>OD_DELTA_ADAPT_SPEED)
           +pvq_adapt_row->mean_count_q8;
          pvq_adapt.mean_count_ex_q8=
           (hmean_count_ex_q7[pli]>>OD_DELTA_ADAPT_SPEED)
           +pvq_adapt_row->mean_count_ex_q8;
          ncount=count_total_q8=count_ex_total_q8=0;
#if !defined(OD_DISABLE_PVQ_CODE1)
          pvq_adapt.mean_pos_q4=OD_MAXI((hmean_pos_q3>>OD_POS_ADAPT_SPEED)
           +pvq_adapt_row->mean_pos_q4,1);
          npos=pos_total=0;
#endif
          adapt_row=&_enc->state.adapt_row[pli];
          od_adapt_update_stats(adapt_row, mbx, &adapt_hmean[pli], &adapt);
          for(by=mby<<2-ydec;by<mby+1<<2-ydec;by++){
            for(bx=mbx<<2-xdec;bx<mbx+1<<2-xdec;bx++){
              od_coeff pred[4*4];
              od_coeff predt[4*4];
              ogg_int16_t pvq_scale[4*4];
              int sgn;
              int qg;
              int cblock[16];
              int zzi;
              int vk;
              vk=0;
              /*fDCT a 4x4 block.*/
              od_bin_fdct4x4(d+(by<<2)*w+(bx<<2),w,c+(by<<2)*w+(bx<<2),w);
              for(zzi=0;zzi<16;zzi++)pvq_scale[zzi]=0;
              if(bx>0&&by>0){
                if(pli==0){
                  ogg_uint16_t mode_cdf[OD_INTRA_NMODES];
                  ogg_uint32_t mode_dist[OD_INTRA_NMODES];
                  int          m_l;
                  int          m_ul;
                  int          m_u;
                  int          mode;
                  m_l=modes[by*(w>>2)+bx-1];
                  m_ul=modes[(by-1)*(w>>2)+bx-1];
                  m_u=modes[(by-1)*(w>>2)+bx];
                  od_intra_pred_cdf(mode_cdf,OD_INTRA_PRED_PROB_4x4[pli],
                   mode_p0,OD_INTRA_NMODES,m_l,m_ul,m_u);
                  od_intra_pred4x4_dist(mode_dist,d+(by<<2)*w+(bx<<2),w,pli);
                  /*Lambda = 1*/
                  mode=od_intra_pred_search(mode_p0,mode_cdf,mode_dist,
                   OD_INTRA_NMODES,128,m_l,m_ul,m_u);
                  od_intra_pred4x4_get(pred,d+(by<<2)*w+(bx<<2),w,mode);
                  od_ec_encode_cdf_unscaled(&_enc->ec,mode,mode_cdf,
                   OD_INTRA_NMODES);
                  mode_bits-=M_LOG2E*log(
                   (mode_cdf[mode]-(mode==0?0:mode_cdf[mode-1]))/
                   (float)mode_cdf[OD_INTRA_NMODES-1]);
                  mode_count++;
                  modes[by*(w>>2)+bx]=mode;
                }
                else{
                  int chroma_weights_q8[3];
                  int mode;
                  mode=modes[(by<<ydec)*(frame_width>>2)+(bx<<xdec)];
                  chroma_weights_q8[0]=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                  chroma_weights_q8[1]=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                  chroma_weights_q8[2]=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                  mode=modes[(by<<ydec)*(frame_width>>2)+((bx<<xdec)+xdec)];
                  chroma_weights_q8[0]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                  chroma_weights_q8[1]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                  chroma_weights_q8[2]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                  mode=modes[((by<<ydec)+ydec)*(frame_width>>2)+(bx<<xdec)];
                  chroma_weights_q8[0]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                  chroma_weights_q8[1]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                  chroma_weights_q8[2]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                  mode=modes[((by<<ydec)+ydec)*(frame_width>>2)
                   +((bx<<xdec)+xdec)];
                  chroma_weights_q8[0]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                  chroma_weights_q8[1]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                  chroma_weights_q8[2]+=OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                  od_chroma_pred4x4(pred,d+(by<<2)*w+(bx<<2),
                   l+(by<<2)*w+(bx<<2),w,chroma_weights_q8);
                }
              }
              else{
                for(zzi=0;zzi<16;zzi++)pred[zzi]=0;
                if(bx>0)pred[0]=d[(by<<2)*w+(bx-1<<2)];
                else if(by>0)pred[0]=d[(by-1<<2)*w+(bx<<2)];
                if(pli==0)modes[by*(w>>2)+bx]=0;
              }
              /*Quantize*/
              for(y=0;y<4;y++){
                for(x=0;x<4;x++){
                  d[((by<<2)+y)*w+(bx<<2)+x]=cblock[OD_ZIG4[y*4+x]]=
                   d[((by<<2)+y)*w+(bx<<2)+x];
                   /*/scale;*//*OD_DIV_ROUND((p[k]-pred[1]),scale);*/
                  predt[OD_ZIG4[y*4+x]]=pred[y*4+x];
                }
              }
              sgn=(cblock[0]-predt[0])<0;
              cblock[0]=(int)floor(pow(fabs(cblock[0]-predt[0])/scale,0.75));
              generic_encode(&_enc->ec,model_dc+pli,cblock[0],ex_dc+pli,0);
              if(cblock[0])od_ec_enc_bits(&_enc->ec,sgn,1);
              cblock[0]=(int)(pow(cblock[0],4.0/3)*scale);
              cblock[0]*=sgn?-1:1;
              cblock[0]+=predt[0];
#if 1
              quant_pvq(cblock+1,predt+1,pvq_scale,pred+1,15,scale,&qg);
              generic_encode(&_enc->ec,model_g+pli,abs(qg),ex_g+pli,0);
              if(qg)od_ec_enc_bits(&_enc->ec,qg<0,1);
              vk=0;
              for(zzi=0;zzi<15;zzi++)vk+=abs(pred[zzi+1]);
              /*No need to code vk because we can get it from qg.*/
              /*Expectation is that half the pulses will go in y[m].*/
              if(vk!=0){
                int ex_ym;
                ex_ym=(65536/2)*vk;
                generic_encode(&_enc->ec,model_ym+pli,vk-pred[1],&ex_ym,0);
              }
              pvq_encoder(&_enc->ec,pred+2,14,vk-abs(pred[1]),&pvq_adapt);
#else
              vk=quant_scalar(cblock+1,predt+1,pvq_scale,pred+1,15,11,
               &pvq_adapt);
              {
                /*This stupid crap is jmspeex's fault.
                  Don't ever enable this code.*/
                static int ex_k=32768;
                generic_encode(&_enc->ec,&model_ym,&ex_k,4);
              }
              /*Treat first component (y[m]) like all others.*/
              pvq_encoder(&_enc->ec,pred+1,15,vk,&pvq_adapt);
#endif
              if(pvq_adapt.k>=0){
                nk++;
                k_total+=pvq_adapt.k;
                sum_ex_total_q8+=pvq_adapt.sum_ex_q8;
              }
              if(pvq_adapt.count_q8>=0){
                ncount++;
                count_total_q8+=pvq_adapt.count_q8;
                count_ex_total_q8+=pvq_adapt.count_ex_q8;
              }
#if !defined(OD_DISABLE_PVQ_CODE1)
              if(pvq_adapt.pos>=0){
                npos++;
                pos_total+=pvq_adapt.pos;
              }
#endif
              /*Dequantize*/
              for(y=0;y<4;y++){
                for(x=0;x<4;x++){
                  d[((by<<2)+y)*w+(bx<<2)+x]=cblock[OD_ZIG4[y*4+x]];
                }
              }
              /*iDCT the 4x4 block.*/
              od_bin_idct4x4(c+(by<<2)*w+(bx<<2),w,d+(by<<2)*w+(bx<<2),w);
            }
          }
          if(nk>0){
            pvq_adapt_row->k=OD_DIVU_SMALL(k_total<<7,nk);
            hmean_k_q7[pli]+=
             pvq_adapt_row->k-hmean_k_q7[pli]>>OD_K_ADAPT_SPEED;
            pvq_adapt_row->mean_k_q8+=
             hmean_k_q7[pli]-pvq_adapt_row->mean_k_q8>>OD_SUM_EX_ADAPT_SPEED;
            pvq_adapt_row->sum_ex_q8=OD_DIVU_SMALL(sum_ex_total_q8,nk)>>1;
            hmean_sum_ex_q7[pli]+=pvq_adapt_row->sum_ex_q8
             -hmean_sum_ex_q7[pli]>>OD_SUM_EX_ADAPT_SPEED;
            pvq_adapt_row->mean_sum_ex_q8+=hmean_sum_ex_q7[pli]
             -pvq_adapt_row->mean_sum_ex_q8>>OD_K_ADAPT_SPEED;
          }
          else pvq_adapt_row->sum_ex_q8=pvq_adapt_row->k=-1;
          if(ncount>0){
            pvq_adapt_row->count_q8=
             OD_DIVU_SMALL(count_total_q8,ncount)>>1;
            hmean_count_q7[pli]+=pvq_adapt_row->count_q8
             -hmean_count_q7[pli]>>OD_DELTA_ADAPT_SPEED;
            pvq_adapt_row->mean_count_q8+=hmean_count_q7[pli]
             -pvq_adapt_row->mean_count_q8>>OD_DELTA_ADAPT_SPEED;
            pvq_adapt_row->count_ex_q8=
             OD_DIVU_SMALL(count_ex_total_q8,ncount)>>1;
            hmean_count_ex_q7[pli]+=pvq_adapt_row->count_ex_q8
             -hmean_count_ex_q7[pli]>>OD_DELTA_ADAPT_SPEED;
            pvq_adapt_row->mean_count_ex_q8+=hmean_count_ex_q7[pli]
             -pvq_adapt_row->mean_count_ex_q8>>OD_DELTA_ADAPT_SPEED;
          }
          else pvq_adapt_row->count_ex_q8=pvq_adapt_row->count_q8=-1;
#if !defined(OD_DISABLE_PVQ_CODE1)
          if(npos>0){
            pvq_adapt_row->pos=OD_DIVU_SMALL(pos_total<<3,npos);
            hmean_pos_q3[pli]+=
             pvq_adapt_row->pos-hmean_pos_q3[pli]>>OD_POS_ADAPT_SPEED;
            pvq_adapt_row->mean_pos_q4+=
             hmean_pos_q3[pli]-pvq_adapt_row->mean_pos_q4>>OD_POS_ADAPT_SPEED;
          }
          else pvq_adapt_row->pos=-1;
#endif
          if(nk>0){
            adapt.curr[OD_ADAPT_K_Q8] = OD_DIVU_SMALL(k_total<<8,nk);
            adapt.curr[OD_ADAPT_SUM_EX_Q8] = OD_DIVU_SMALL(sum_ex_total_q8,nk);
          } else {
            adapt.curr[OD_ADAPT_K_Q8] = OD_ADAPT_NO_VALUE;
            adapt.curr[OD_ADAPT_SUM_EX_Q8] = OD_ADAPT_NO_VALUE;
          }
          if (ncount>0)
          {
            adapt.curr[OD_ADAPT_COUNT_Q8] = OD_DIVU_SMALL(count_total_q8,ncount);
            adapt.curr[OD_ADAPT_COUNT_EX_Q8] = OD_DIVU_SMALL(count_ex_total_q8,ncount);
          } else {
            adapt.curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
            adapt.curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
          }
          od_adapt_mb(adapt_row, mbx, &adapt_hmean[pli], &adapt);

          /*Apply the postfilter across the left block edges.*/
          for(y=oyfill>>ydec;y<next_oyfill>>ydec;y++){
            for(bx=(mbx<<2-xdec)+(mbx<=0);bx<mbx+1<<2-xdec;bx++){
              od_post_filter4(c+y*w+(bx<<2)-2,c+y*w+(bx<<2)-2);
            }
          }
          /*Apply the postfilter across the top block edges.*/
          for(by=(mby<<2-ydec)+(mby<=0);by<mby+1<<2-ydec;by++){
            for(x=oxfill>>xdec;x<next_oxfill>>xdec;x++){
              od_coeff p[4];
              for(y=0;y<4;y++)p[y]=c[((by<<2)+y-2)*w+x];
              od_post_filter4(p,p);
              for(y=0;y<4;y++)c[((by<<2)+y-2)*w+x]=p[y];
            }
          }
          data=_enc->state.io_imgs[OD_FRAME_REC].planes[pli].data;
          for(y=oyfill>>ydec;y<next_oyfill>>ydec;y++){
            for(x=oxfill>>xdec;x<next_oxfill>>xdec;x++){
              data[ystride*y+x]=OD_CLAMP255(c[y*w+x]+128);
            }
          }
        }
        ixfill=next_ixfill;
        oxfill=next_oxfill;
      }
      for(pli=0;pli<nplanes;pli++){
        od_pvq_adapt_ctx *pvq_adapt_row;
        pvq_adapt_row=_enc->state.pvq_adapt_row[pli];
        hmean_k_q7[pli]=OD_K_INIT_Q7;
        hmean_sum_ex_q7[pli]=OD_SUM_EX_INIT_Q7;
        hmean_count_q7[pli]=OD_COUNT_INIT_Q7;
        hmean_count_ex_q7[pli]=OD_COUNT_EX_INIT_Q7;
#if !defined(OD_DISABLE_PVQ_CODE1)
        hmean_pos_q3[pli]=OD_POS_INIT_Q3;
#endif
        for(mbx=nhmbs;mbx-->0;){
          if(pvq_adapt_row[mbx].k>=0){
            pvq_adapt_row[mbx].mean_k_q8+=(hmean_k_q7[pli]>>OD_K_ADAPT_SPEED)
             -(hmean_k_q7[pli]>>2*OD_K_ADAPT_SPEED);
            hmean_k_q7[pli]+=
             pvq_adapt_row[mbx].k-hmean_k_q7[pli]>>OD_K_ADAPT_SPEED;
            pvq_adapt_row[mbx].mean_sum_ex_q8+=
             (hmean_sum_ex_q7[pli]>>OD_SUM_EX_ADAPT_SPEED)
             -(hmean_sum_ex_q7[pli]>>2*OD_SUM_EX_ADAPT_SPEED);
            hmean_sum_ex_q7[pli]+=pvq_adapt_row[mbx].sum_ex_q8
             -hmean_sum_ex_q7[pli]>>OD_SUM_EX_ADAPT_SPEED;
          }
          if(pvq_adapt_row[mbx].count_q8>=0){
            pvq_adapt_row[mbx].mean_count_q8+=(hmean_count_q7[pli]>>OD_DELTA_ADAPT_SPEED)
             -(hmean_count_q7[pli]>>2*OD_DELTA_ADAPT_SPEED);
            hmean_count_q7[pli]+=
             pvq_adapt_row[mbx].count_q8-hmean_count_q7[pli]>>OD_DELTA_ADAPT_SPEED;
            pvq_adapt_row[mbx].mean_count_ex_q8+=
             (hmean_count_ex_q7[pli]>>OD_DELTA_ADAPT_SPEED)
             -(hmean_count_ex_q7[pli]>>2*OD_DELTA_ADAPT_SPEED);
            hmean_count_ex_q7[pli]+=pvq_adapt_row[mbx].count_ex_q8
             -hmean_count_ex_q7[pli]>>OD_DELTA_ADAPT_SPEED;
          }
#if !defined(OD_DISABLE_PVQ_CODE1)
          if(pvq_adapt_row[mbx].pos>=0){
            pvq_adapt_row[mbx].mean_pos_q4+=
             (hmean_pos_q3[pli]>>OD_POS_ADAPT_SPEED)
             -(hmean_pos_q3[pli]>>2*OD_POS_ADAPT_SPEED);
            hmean_pos_q3[pli]+=
             pvq_adapt_row[mbx].pos-hmean_pos_q3[pli]>>OD_POS_ADAPT_SPEED;
          }
#endif
        }
        od_adapt_row(&_enc->state.adapt_row[pli], &adapt_hmean[pli]);
      }
      iyfill=next_iyfill;
      oyfill=next_oyfill;
    }
    for(pli=nplanes;pli-->0;){
      _ogg_free(ltmp[pli]);
      _ogg_free(dtmp[pli]);
      _ogg_free(ctmp[pli]);
    }
    _ogg_free(modes);
  }
  for(pli=0;pli<nplanes;pli++){
    unsigned char *data;
    ogg_int64_t    mc_sqerr;
    ogg_int64_t    enc_sqerr;
    ogg_uint32_t   npixels;
    int            ystride;
    int            xdec;
    int            ydec;
    int            w;
    int            h;
    int            x;
    int            y;
#ifdef OD_DPCM
    int            err_accum;
    err_accum=0;
#endif
    mc_sqerr=0;
    enc_sqerr=0;
    data=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].data;
    ystride=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
    xdec=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
    ydec=_enc->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
    w=frame_width>>xdec;
    h=frame_height>>ydec;
    npixels=w*h;
    for(y=0;y<h;y++){
      unsigned char *prev_rec_row;
      unsigned char *rec_row;
      unsigned char *inp_row;
      rec_row=_enc->state.io_imgs[OD_FRAME_REC].planes[pli].data+
       _enc->state.io_imgs[OD_FRAME_REC].planes[pli].ystride*y;
      prev_rec_row=rec_row
       -_enc->state.io_imgs[OD_FRAME_REC].planes[pli].ystride;
      inp_row=data+ystride*y;
      memcpy(_enc->state.ref_line_buf[1],rec_row,w);
      for(x=0;x<w;x++){
        int rec_val;
        int inp_val;
        int diff;
        rec_val=rec_row[x];
        inp_val=inp_row[x];
        diff=inp_val-rec_val;
        mc_sqerr+=diff*diff;
#ifdef OD_DPCM
        {
          int pred_diff;
          int qdiff;
          /*DPCM code the residual with uniform quantization.
            This provides simulated residual coding errors, without introducing
             blocking artifacts.*/
          if(x>0)pred_diff=rec_row[x-1]-_enc->state.ref_line_buf[1][x-1];
          else pred_diff=0;
          if(y>0){
            if(x>0){
              pred_diff+=prev_rec_row[x-1]-_enc->state.ref_line_buf[0][x-1];
            }
            pred_diff+=prev_rec_row[x]-_enc->state.ref_line_buf[0][x];
            if(x+1<w){
              pred_diff+=prev_rec_row[x+1]-_enc->state.ref_line_buf[0][x+1];
            }
          }
          pred_diff=OD_DIV_ROUND_POW2(pred_diff,2,2);
          qdiff=(((diff-pred_diff)+(diff-pred_diff>>31)+(5+err_accum))/10)*10+
           pred_diff;
          /*qdiff=(OD_DIV_ROUND_POW2(diff-pred_diff,3,4+err_accum)<<3)+
           pred_diff;*/
          /*fprintf(stderr,"d-p_d: %3i  e_a: %3i  qd-p_d: %3i  e_a: %i\n",
           diff-pred_diff,err_accum,qdiff-pred_diff,diff-qdiff);*/
          err_accum+=diff-qdiff;
          rec_row[x]=OD_CLAMP255(rec_val+qdiff);
        }
#else
/*        rec_row[x]=inp_val;*/
#endif
        diff=inp_val-rec_row[x];
        enc_sqerr+=diff*diff;
      }
      prev_rec_row=_enc->state.ref_line_buf[0];
      _enc->state.ref_line_buf[0]=_enc->state.ref_line_buf[1];
      _enc->state.ref_line_buf[1]=prev_rec_row;
    }
    /*printf("Bytes: %d  ex_dc: %d ex_g: %d ex_k: %d\n",(od_ec_enc_tell(&_enc->ec)+7)>>3,ex_dc,ex_g,ex_k);*/
    if(_enc->state.ref_imgi[OD_FRAME_PREV]>=0){
      fprintf(stderr,
       "Plane %i, Squared Error: %12lli  Pixels: %6u  PSNR:  %5.2f\n",
       pli,(long long)mc_sqerr,npixels,10*log10(255*255.0*npixels/mc_sqerr));
    }
    fprintf(stderr,
     "Encoded Plane %i, Squared Error: %12lli  Pixels: %6u  PSNR:  %5.2f\n",
     pli,(long long)enc_sqerr,npixels,10*log10(255*255.0*npixels/enc_sqerr));
  }
  fprintf(stderr, "mode bits: %f/%f=%f\n", mode_bits,mode_count,mode_bits/mode_count);
  /*Dump YUV*/
  od_state_dump_yuv(&_enc->state,_enc->state.io_imgs+OD_FRAME_REC,"out");
  _enc->packet_state=OD_PACKET_READY;
  od_state_upsample8(&_enc->state,
   _enc->state.ref_imgs+_enc->state.ref_imgi[OD_FRAME_SELF],
   _enc->state.io_imgs+OD_FRAME_REC);
#if defined(OD_DUMP_IMAGES)
  /*Dump reference frame.*/
  /*od_state_dump_img(&_enc->state,
   _enc->state.ref_img+_enc->state.ref_imigi[OD_FRAME_SELF],"ref");*/
#endif
  if(_enc->state.info.frame_duration==0)_enc->state.cur_time+=_duration;
  else _enc->state.cur_time+=_enc->state.info.frame_duration;
  return 0;
}

int daala_encode_packet_out(daala_enc_ctx *_enc,int _last,ogg_packet *_op){
  ogg_uint32_t nbytes;
  if(_enc==NULL||_op==NULL)return OD_EFAULT;
  else if(_enc->packet_state<=0||_enc->packet_state==OD_PACKET_DONE)return 0;
  _op->packet=od_ec_enc_done(&_enc->ec,&nbytes);
  _op->bytes=nbytes;
  fprintf(stderr,"Output Bytes: %ld\n",_op->bytes);
  _op->b_o_s=0;
  _op->e_o_s=_last;
  _op->packetno=0;
  _op->granulepos=_enc->state.cur_time;
  if(_last)_enc->packet_state=OD_PACKET_DONE;
  else _enc->packet_state=OD_PACKET_EMPTY;
  return 1;
}
