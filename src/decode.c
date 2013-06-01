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
#include <math.h>
#include <string.h>
#include "decint.h"
#include "generic_code.h"
#include "filter.h"
#include "dct.h"
#include "intra.h"
#include "pvq.h"
#include "pvq_code.h"
#include "block_size_dec.h"

static int od_dec_init(od_dec_ctx *dec, const daala_info *info,
 const daala_setup_info *setup) {
  int ret;
  (void)setup;
  ret = od_state_init(&dec->state, info);
  if (ret < 0) return ret;
  dec->packet_state = OD_PACKET_DATA;
  return 0;
}

static void od_dec_clear(od_dec_ctx *dec) {
  od_state_clear(&dec->state);
}

daala_dec_ctx *daala_decode_alloc(const daala_info *info,
 const daala_setup_info *setup) {
  od_dec_ctx *dec;
  if (info == NULL) return NULL;
  dec = (od_dec_ctx *)_ogg_malloc(sizeof(*dec));
  if (od_dec_init(dec, info, setup) < 0) {
    _ogg_free(dec);
    return NULL;
  }
  return dec;
}

void daala_decode_free(daala_dec_ctx *dec) {
  if (dec != NULL) {
    od_dec_clear(dec);
    _ogg_free(dec);
  }
}

int daala_decode_ctl(daala_dec_ctx *dec, int req, void *buf, size_t buf_sz) {
  (void)dec;
  (void)buf;
  (void)buf_sz;
  switch(req) {
    default: return OD_EIMPL;
  }
}

static void od_decode_mv(daala_dec_ctx *dec, od_mv_grid_pt *mvg,
 int mv_res, int width, int height) {
  int ox;
  int oy;
  ox = laplace_decode(&dec->ec, 2269 >> mv_res, width << 1);
  oy = laplace_decode(&dec->ec, 569 >> mv_res, height << 1);
  /*Deinterleave positive and negative values.*/
  ox = (ox >> 1) ^ -(ox & 1);
  oy = (oy >> 1) ^ -(oy & 1);
  mvg->mv[0] = ox << mv_res;
  mvg->mv[1] = oy << mv_res;
}

int daala_decode_packet_in(daala_dec_ctx *dec, od_img *img,
 const ogg_packet *op) {
  int nplanes;
  int pli;
  int scale[OD_NPLANES_MAX];
  int frame_width;
  int frame_height;
  int pic_width;
  int pic_height;
  int nvsb;
  int nhsb;
  int refi;
  int i;
  int j;
  int is_keyframe;
  if (dec == NULL || img == NULL || op == NULL) return OD_EFAULT;
  if (dec->packet_state != OD_PACKET_DATA) return OD_EINVAL;
  if (op->e_o_s) dec->packet_state = OD_PACKET_DONE;
  od_ec_dec_init(&dec->ec, op->packet, op->bytes);
  /*Read the packet type bit.*/
  if (od_ec_decode_bool_q15(&dec->ec, 16384)) return OD_EBADPACKET;
  is_keyframe = od_ec_decode_bool_q15(&dec->ec, 16384);
  /*Update the buffer state.*/
  if (dec->state.ref_imgi[OD_FRAME_SELF] >= 0) {
    dec->state.ref_imgi[OD_FRAME_PREV] =
     dec->state.ref_imgi[OD_FRAME_SELF];
    /*TODO: Update golden frame.*/
    if (dec->state.ref_imgi[OD_FRAME_GOLD] < 0) {
      dec->state.ref_imgi[OD_FRAME_GOLD] =
       dec->state.ref_imgi[OD_FRAME_SELF];
      /*TODO: Mark keyframe timebase.*/
    }
  }
  /*Select a free buffer to use for this reference frame.*/
  for (refi = 0; refi == dec->state.ref_imgi[OD_FRAME_GOLD]
   || refi == dec->state.ref_imgi[OD_FRAME_PREV]
   || refi == dec->state.ref_imgi[OD_FRAME_NEXT]; refi++);
  dec->state.ref_imgi[OD_FRAME_SELF] = refi;
  nhsb = dec->state.nhsb;
  nvsb = dec->state.nvsb;
  for(i = -4; i < (nhsb+1)*4; i++) {
    for(j = -4; j < 0; j++) {
      dec->state.bsize[(j*dec->state.bstride) + i] = 3;
    }
    for(j = nvsb*4; j < (nvsb+1)*4; j++) {
      dec->state.bsize[(j*dec->state.bstride) + i] = 3;
    }
  }
  for(j = -4; j < (nvsb+1)*4; j++) {
    for(i = -4; i < 0; i++) {
      dec->state.bsize[(j*dec->state.bstride) + i] = 3;
    }
    for(i = nhsb*4; i < (nhsb+1)*4; i++) {
      dec->state.bsize[(j*dec->state.bstride) + i] = 3;
    }
  }
  for(i = 0; i < nvsb; i++) {
    for(j = 0; j < nhsb; j++) {
      od_block_size_decode(&dec->ec,
       &dec->state.bsize[4*dec->state.bstride*i + 4*j], dec->state.bstride);
    }
  }
  if(dec->state.ref_imgi[OD_FRAME_PREV] >= 0 && !is_keyframe){
    /* Input the motion vectors. */
    int nhmvbs;
    int nvmvbs;
    int vx;
    int vy;
    od_img *img;
    int width;
    int height;
    int mv_res;
    od_mv_grid_pt *mvp;
    od_mv_grid_pt **grid;
    od_state_mvs_clear(&dec->state);
    nhmvbs = (dec->state.nhmbs + 1) << 2;
    nvmvbs = (dec->state.nvmbs + 1) << 2;
    img = dec->state.io_imgs + OD_FRAME_REC;
    mv_res = dec->state.mv_res = od_ec_dec_uint(&dec->ec, 3);
    width = (img->width + 32) << (3 - mv_res);
    height = (img->height + 32) << (3 - mv_res);
    grid = dec->state.mv_grid;
    /*Level 0.*/
    for (vy = 0; vy <= nvmvbs; vy += 4) {
      for (vx = 0; vx <= nhmvbs; vx += 4) {
        mvp = &grid[vy][vx];
        mvp->valid = 1;
        od_decode_mv(dec, mvp, mv_res, width, height);
      }
    }
    /*Level 1.*/
    for (vy = 2; vy <= nvmvbs; vy += 4) {
      for (vx = 2; vx <= nhmvbs; vx += 4) {
        mvp = &grid[vy][vx];
        mvp->valid = od_ec_decode_bool_q15(&dec->ec, 16384);
        if (mvp->valid) od_decode_mv(dec, mvp, mv_res, width, height);
      }
    }
    /*Level 2.*/
    for (vy = 0; vy <= nvmvbs; vy += 2) {
      for (vx = 2*((vy & 3) == 0); vx <= nhmvbs; vx += 4) {
        mvp = &grid[vy][vx];
        if (vy-2 >= 0 && grid[vy-2][vx].valid
         && vx-2 >= 0 && grid[vy][vx-2].valid
         && vy+2 <= nvmvbs && grid[vy+2][vx].valid
         && vx+2 <= nhmvbs && grid[vy][vx+2].valid) {
          mvp->valid = od_ec_decode_bool_q15(&dec->ec, 16384);
          if (mvp->valid) od_decode_mv(dec, mvp, mv_res, width, height);
        }
      }
    }
    od_state_mc_predict(&dec->state, OD_FRAME_PREV);
  }
  frame_width = dec->state.frame_width;
  frame_height = dec->state.frame_height;
  pic_width = dec->state.info.pic_width;
  pic_height = dec->state.info.pic_height;
  {
    GenericEncoder model_dc[OD_NPLANES_MAX];
    GenericEncoder model_g[OD_NPLANES_MAX];
    GenericEncoder model_ym[OD_NPLANES_MAX];
    ogg_uint16_t mode_p0[OD_INTRA_NMODES];
    int ex_dc[OD_NPLANES_MAX];
    int ex_g[OD_NPLANES_MAX];
    od_coeff *ctmp[OD_NPLANES_MAX];
    od_coeff *dtmp[OD_NPLANES_MAX];
    od_coeff *mctmp[OD_NPLANES_MAX];
    od_coeff *mdtmp[OD_NPLANES_MAX];
    od_coeff *ltmp[OD_NPLANES_MAX];
    od_coeff *lbuf[OD_NPLANES_MAX];
    signed char *modes;
    int xdec;
    int ydec;
    int nvmbs;
    int nhmbs;
    int mby;
    int mbx;
    int mi;
    int h;
    int w;
    int y;
    int x;
    nhmbs = dec->state.nhmbs;
    nvmbs = dec->state.nvmbs;
    /*Initialize the data needed for each plane.*/
    modes = _ogg_calloc((frame_width >> 2)*(frame_height >> 2),
     sizeof(*modes));
    for (mi = 0; mi < OD_INTRA_NMODES; mi++) {
      mode_p0[mi] = 32768/OD_INTRA_NMODES;
    }
    nplanes = dec->state.info.nplanes;
    /*Apply the prefilter to the motion-compensated reference.*/
    if (!is_keyframe) {
      for (pli = 0; pli < nplanes; pli++) {
        xdec = dec->state.io_imgs[OD_FRAME_REC].planes[pli].xdec;
        ydec = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ydec;
        w = frame_width >> xdec;
        h = frame_height >> ydec;
        mctmp[pli] = _ogg_calloc(w*h, sizeof(*mctmp[pli]));
        mdtmp[pli] = _ogg_calloc(w*h, sizeof(*mdtmp[pli]));
        /*Collect the image data needed for this plane.*/
        {
          unsigned char *mdata;
          int ystride;
          mdata = dec->state.io_imgs[OD_FRAME_REC].planes[pli].data;
          ystride = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ystride;
          for (y=0;y<h;y++) {
            for (x=0;x<w;x++) mctmp[pli][y*w+x]=mdata[ystride*y+x]-128;
          }
        }
        /*Apply the prefilter across the entire image.*/
        {
          int sby;
          int sbx;
          /* This code assumes 4:4:4 or 4:2:0 input. */
          OD_ASSERT(xdec==ydec);
          /*Apply the prefilter down the bottom block edge columns.*/
          for (sby = 0; sby < nvsb; sby++) {
            for (sbx = 0; sbx < nhsb; sbx++) {
              unsigned char btmp[6*6];
              od_extract_bsize(btmp,6,&dec->state.bsize[dec->state.bstride*(sby<<2)+(sbx<<2)],dec->state.bstride,xdec);
              od_apply_filter(&mctmp[pli][(sby<<(5-ydec))*w+(sbx<<(5-xdec))],w,0,0,3-xdec,
               &btmp[6*1+1],6,OD_BOTTOM_EDGE,sby<nvsb-1?OD_BOTTOM_EDGE:0,0);
            }
          }
          /*Apply the prefilter across the right block edge rows.*/
          for (sby = 0; sby < nvsb; sby++) {
            for (sbx = 0; sbx < nhsb; sbx++) {
              unsigned char btmp[6*6];
              od_extract_bsize(btmp,6,&dec->state.bsize[dec->state.bstride*(sby<<2)+(sbx<<2)],dec->state.bstride,xdec);
              od_apply_filter(&mctmp[pli][(sby<<(5-ydec))*w+(sbx<<(5-xdec))],w,0,0,3-xdec,
               &btmp[6*1+1],6,OD_RIGHT_EDGE,sbx<nhsb-1?OD_RIGHT_EDGE:0,0);
            }
          }
        }
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      generic_model_init(model_dc + pli);
      generic_model_init(model_g + pli);
      generic_model_init(model_ym + pli);
      ex_dc[pli] = pli > 0 ? 8 : 32768;
      ex_g[pli] = 8;
      xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      scale[pli] = od_ec_dec_uint(&dec->ec, 512);
      ctmp[pli] = _ogg_calloc(w*h, sizeof(*ctmp[pli]));
      dtmp[pli] = _ogg_calloc(w*h, sizeof(*dtmp[pli]));
      /*We predict chroma planes from the luma plane.
        Since chroma can be subsampled, we cache subsampled versions of the
         luma plane in the frequency domain.
        We can share buffers with the same subsampling.*/
      if (pli > 0) {
        int plj;
        if (xdec || ydec) {
          for (plj = 1; plj < pli; plj++) {
            if (xdec == dec->state.io_imgs[OD_FRAME_INPUT].planes[plj].xdec
             && ydec == dec->state.io_imgs[OD_FRAME_INPUT].planes[plj].ydec) {
              ltmp[pli] = NULL;
              lbuf[pli] = ltmp[plj];
            }
          }
          if (plj >= pli) {
            lbuf[pli] = ltmp[pli] = _ogg_calloc(w*h, sizeof(*ltmp));
          }
        }
        else{
          ltmp[pli] = NULL;
          lbuf[pli] = ctmp[pli];
        }
      }
      else lbuf[pli] = ltmp[pli] = NULL;
      od_adapt_row_init(&dec->state.adapt_row[pli]);
    }
    for (mby = 0; mby < nvmbs; mby++) {
      od_adapt_ctx adapt_hmean[OD_NPLANES_MAX];
      for (pli = 0; pli < nplanes; pli++) {
        od_adapt_hmean_init(&adapt_hmean[pli]);
      }
      for (mbx = 0; mbx < nhmbs; mbx++) {
        for (pli = 0; pli < nplanes; pli++) {
          od_adapt_row_ctx *adapt_row;
          od_adapt_ctx adapt;
          od_coeff *c;
          od_coeff *d;
          od_coeff *mc;
          od_coeff *md;
          od_coeff *l;
          unsigned char *data;
          int ystride;
          int by;
          int bx;
          int nk;
          int k_total;
          int sum_ex_total_q8;
          int ncount;
          int count_total_q8;
          int count_ex_total_q8;
          unsigned char* mdata;
          c = ctmp[pli];
          d = dtmp[pli];
          if (!is_keyframe) {
            mc = mctmp[pli];
            md = mdtmp[pli];
          }
          l = lbuf[pli];
          xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
          ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
          w = frame_width >> xdec;
          h = frame_height >> ydec;
          /*Construct the luma predictors for chroma planes.*/
          if (ltmp[pli] != NULL) {
            OD_ASSERT(pli > 0);
            OD_ASSERT(l == ltmp[pli]);
            for (by = mby << (2 - ydec); by < (mby + 1) << (2 - ydec); by++) {
              for (bx = mbx << (2 - xdec); bx < (mbx + 1) << (2 - xdec);
               bx++) {
                od_resample_luma_coeffs(l + (by << 2)*w + (bx<<2), w,
                 dtmp[0] + (by << (2 + ydec))*frame_width + (bx<<(2 + xdec)),
                 frame_width, xdec, ydec, 4);
              }
            }
          }
          ystride = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ystride;
          data = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].data;
          mdata = dec->state.io_imgs[OD_FRAME_REC].planes[pli].data;
          nk = k_total = sum_ex_total_q8 = 0;
          ncount = count_total_q8 = count_ex_total_q8 = 0;
          adapt_row = &dec->state.adapt_row[pli];
          od_adapt_update_stats(adapt_row, mbx, &adapt_hmean[pli], &adapt);
          for (by = mby << (2 - ydec); by < (mby + 1) << (2 - ydec); by++) {
            for (bx = mbx << (2 - xdec); bx < (mbx + 1) << (2 - xdec); bx++) {
              od_coeff pred[4*4];
              od_coeff predt[4*4];
              ogg_int16_t pvq_scale[4*4];
              int sgn;
              int qg;
              int zzi;
              int vk;
#ifdef OD_LOLOSSLESS
              od_coeff backup[4*4];
#endif
              vk = 0;
              if (!is_keyframe) {
                od_bin_fdct4x4(md + (by << 2)*w + (bx << 2), w,
                    mc + (by << 2)*w + (bx << 2), w);
              }
              for (zzi = 0; zzi < 16; zzi++) pvq_scale[zzi] = 0;
              if (is_keyframe) {
                if (bx > 0 && by > 0) {
                  if (pli == 0) {
                    ogg_uint16_t mode_cdf[OD_INTRA_NMODES];
                    int m_l;
                    int m_ul;
                    int m_u;
                    int mode;
                    od_coeff *ur;
                    od_coeff *coeffs[4];
                    int strides[4];
                    ur = (by > 0 && (((bx + 1) < (mbx + 1) << (2 - xdec))
                     || (by == mby << (2 - ydec)))) ?
                     d + ((by - 1) << 2)*w + ((bx + 1) << 2) :
                     d + ((by - 1) << 2)*w + (bx << 2);
                    m_l = modes[by*(w >> 2) + bx - 1];
                    m_ul = modes[(by - 1)*(w >> 2) + bx - 1];
                    m_u = modes[(by - 1)*(w >> 2) + bx];
                    coeffs[0] = d + ((by - 1) << 2)*w + ((bx - 1) << 2);
                    coeffs[1] = d + ((by - 1) << 2)*w + (bx << 2);
                    coeffs[2] = ur;
                    coeffs[3] = d + (by << 2)*w + ((bx - 1) << 2);
                    strides[0] = w;
                    strides[1] = w;
                    strides[2] = w;
                    strides[3] = w;
                    od_intra_pred_cdf(mode_cdf, OD_INTRA_PRED_PROB_4x4[pli],
                     mode_p0, OD_INTRA_NMODES, m_l, m_ul, m_u);
                    mode = od_ec_decode_cdf_unscaled(&dec->ec, mode_cdf,
                     OD_INTRA_NMODES);
                    od_intra_pred4x4_get(pred, coeffs, strides, mode);
                    modes[by*(w >> 2) + bx] = mode;
                    od_intra_pred_update(mode_p0, OD_INTRA_NMODES, mode,
                     m_l, m_ul, m_u);
                  }
                  else{
                    int chroma_weights_q8[3];
                    int mode;
                    mode = modes[(by << ydec)*(frame_width >> 2) + (bx << xdec)];
                    chroma_weights_q8[0] = OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                    chroma_weights_q8[1] = OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                    chroma_weights_q8[2] = OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                    mode = modes[(by << ydec)*(frame_width >> 2)
                     + ((bx << xdec) + xdec)];
                    chroma_weights_q8[0] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                    chroma_weights_q8[1] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                    chroma_weights_q8[2] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                    mode = modes[((by << ydec) + ydec)*(frame_width >> 2)
                     + (bx << xdec)];
                    chroma_weights_q8[0] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                    chroma_weights_q8[1] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                    chroma_weights_q8[2] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                    mode = modes[((by << ydec) + ydec)*(frame_width >> 2)
                     + ((bx << xdec) + xdec)];
                    chroma_weights_q8[0] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][0];
                    chroma_weights_q8[1] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][1];
                    chroma_weights_q8[2] +=
                     OD_INTRA_CHROMA_WEIGHTS_Q6[mode][2];
                    od_chroma_pred4x4(pred, d + (by << 2)*w + (bx << 2),
                     l + (by << 2)*w + (bx << 2), w, chroma_weights_q8);
                  }
                }
                else{
                  for (zzi = 0; zzi < 16; zzi++) pred[zzi] = 0;
                  if (bx > 0) pred[0] = d[(by << 2)*w + ((bx - 1) << 2)];
                  else if (by > 0) pred[0] = d[((by - 1) << 2)*w + (bx << 2)];
                  if (pli == 0) modes[by*(w >> 2) + bx] = 0;
                }
              }
              else {
                int x;
                int y;
                int i;
                i = 0;
                for( y=0; y<4; y++ ) {
                  for( x=0; x<4; x++ ) {
                    pred[i++] = md[(y + (by << 2))*w + (x + (bx << 2))];
                  }
                }
              }
              /*Zig-zag*/
              for (y = 0; y < 4; y++) {
                for (x = 0; x < 4; x++) {
                  predt[OD_ZIG4[y*4 + x]] = pred[y*4 + x];
                }
              }
#ifdef OD_LOLOSSLESS
              for (zzi = 0; zzi < 16; zzi++) {
                backup[zzi] = od_ec_dec_uint(&dec->ec, 65536);
              }
#endif
              sgn = 0;
              pred[0] = generic_decode(&dec->ec, model_dc + pli, ex_dc + pli,
               0);
              if (pred[0]) sgn = od_ec_dec_bits(&dec->ec,1);
              pred[0] = (int)(pow(pred[0],4.0/3)*scale[pli]);
              pred[0] *= sgn ? -1 : 1;
              pred[0] += predt[0];
              qg = generic_decode(&dec->ec, model_g + pli, ex_g + pli, 0);
              if (qg) qg *= od_ec_dec_bits(&dec->ec, 1) ? -1 : 1;
              vk = pvq_unquant_k(&predt[1], 15, qg, scale[pli]);
              pred[1] = 0;
              if (vk != 0) {
                int ex_ym;
                ex_ym = (65536/2)*vk;
                pred[1] = vk - generic_decode(&dec->ec, model_ym + pli,
                 &ex_ym, 0);
              }
              pvq_decoder(&dec->ec, pred + 2, 14, vk - abs(pred[1]), &adapt);
              dequant_pvq(pred + 1, predt + 1, pvq_scale, 15, scale[pli], qg);
              if (adapt.curr[OD_ADAPT_K_Q8] >= 0) {
                nk++;
                k_total += adapt.curr[OD_ADAPT_K_Q8];
                sum_ex_total_q8 += adapt.curr[OD_ADAPT_SUM_EX_Q8];
              }
              if (adapt.curr[OD_ADAPT_COUNT_Q8] >= 0) {
                ncount++;
                count_total_q8 += adapt.curr[OD_ADAPT_COUNT_Q8];
                count_ex_total_q8 += adapt.curr[OD_ADAPT_COUNT_EX_Q8];
              }
#ifdef OD_LOLOSSLESS
              for (zzi = 0; zzi < 16; zzi++) {
                pred[zzi] = backup[zzi] - 32768;
              }
#endif
              /*Dequantize*/
              for (y = 0; y < 4; y++) {
                for (x = 0; x < 4; x++) {
                  d[((by << 2) + y)*w + (bx << 2) + x] =
                   pred[OD_ZIG4[y*4 + x]];
                }
              }
              /*iDCT the 4x4 block.*/
              od_bin_idct4x4(c + (by << 2)*w + (bx << 2), w, d + (by << 2)*w
               + (bx << 2), w);
            }
          }
          if (nk > 0) {
            adapt.curr[OD_ADAPT_K_Q8] = OD_DIVU_SMALL(k_total << 8, nk);
            adapt.curr[OD_ADAPT_SUM_EX_Q8] =
             OD_DIVU_SMALL(sum_ex_total_q8, nk);
          } else {
            adapt.curr[OD_ADAPT_K_Q8] = OD_ADAPT_NO_VALUE;
            adapt.curr[OD_ADAPT_SUM_EX_Q8] = OD_ADAPT_NO_VALUE;
          }
          if (ncount > 0)
          {
            adapt.curr[OD_ADAPT_COUNT_Q8] =
             OD_DIVU_SMALL(count_total_q8, ncount);
            adapt.curr[OD_ADAPT_COUNT_EX_Q8] =
             OD_DIVU_SMALL(count_ex_total_q8, ncount);
          } else {
            adapt.curr[OD_ADAPT_COUNT_Q8] = OD_ADAPT_NO_VALUE;
            adapt.curr[OD_ADAPT_COUNT_EX_Q8] = OD_ADAPT_NO_VALUE;
          }
          od_adapt_mb(adapt_row, mbx, &adapt_hmean[pli], &adapt);
        }
      }
      for (pli = 0; pli < nplanes; pli++) {
        od_adapt_row(&dec->state.adapt_row[pli], &adapt_hmean[pli]);
      }
    }
    for (pli = 0; pli < nplanes; pli++) {
      xdec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].xdec;
      ydec = dec->state.io_imgs[OD_FRAME_INPUT].planes[pli].ydec;
      w = frame_width >> xdec;
      h = frame_height >> ydec;
      /*Apply the postfilter across the entire image.*/
      {
        int sby;
        int sbx;
        /* This code assumes 4:4:4 or 4:2:0 input. */
        OD_ASSERT(xdec==ydec);
        /*Apply the postfilter across the right block edge rows.*/
        for (sby = 0; sby < nvsb; sby++) {
          for (sbx = 0; sbx < nhsb; sbx++) {
            unsigned char btmp[6*6];
            od_extract_bsize(btmp,6,&dec->state.bsize[dec->state.bstride*(sby<<2)+(sbx<<2)],dec->state.bstride,xdec);
            od_apply_filter(&ctmp[pli][(sby<<(5-ydec))*w+(sbx<<(5-xdec))],w,0,0,3-xdec,
             &btmp[6*1+1],6,OD_RIGHT_EDGE,sbx<nhsb-1?OD_RIGHT_EDGE:0,1);
          }
        }
        /*Apply the postfilter down the bottom block edge columns.*/
        for (sby = 0; sby < nvsb; sby++) {
          for (sbx = 0; sbx < nhsb; sbx++) {
            unsigned char btmp[6*6];
            od_extract_bsize(btmp,6,&dec->state.bsize[dec->state.bstride*(sby<<2)+(sbx<<2)],dec->state.bstride,xdec);
            od_apply_filter(&ctmp[pli][(sby<<(5-ydec))*w+(sbx<<(5-xdec))],w,0,0,3-xdec,
             &btmp[6*1+1],6,OD_BOTTOM_EDGE,sby<nvsb-1?OD_BOTTOM_EDGE:0,1);
          }
        }
      }
      {
        unsigned char *data;
        int ystride;
        data = dec->state.io_imgs[OD_FRAME_REC].planes[pli].data;
        ystride = dec->state.io_imgs[OD_FRAME_REC].planes[pli].ystride;
        for (y=0;y<h;y++) {
          for (x=0;x<w;x++) {
            data[ystride*y+x]=OD_CLAMP255(ctmp[pli][y*w+x]+128);
          }
        }
      }
    }
    for (pli = nplanes; pli-- > 0;) {
      _ogg_free(ltmp[pli]);
      _ogg_free(dtmp[pli]);
      _ogg_free(ctmp[pli]);
      if (!is_keyframe) {
        _ogg_free(mdtmp[pli]);
        _ogg_free(mctmp[pli]);
      }
    }
    _ogg_free(modes);
  }
  /*Dump YUV*/
  od_state_dump_yuv(&dec->state, dec->state.io_imgs + OD_FRAME_REC, "decout");
  od_state_upsample8(&dec->state,
   dec->state.ref_imgs + dec->state.ref_imgi[OD_FRAME_SELF],
   dec->state.io_imgs + OD_FRAME_REC);
  /*Return decoded frame.*/
  *img = dec->state.io_imgs[OD_FRAME_REC];
  img->width = dec->state.info.pic_width;
  img->height = dec->state.info.pic_height;
  return 0;
}
