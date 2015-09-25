/*Daala video codec
Copyright (c) 2002-2015 Daala project contributors.  All rights reserved.

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
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vidinput.h"
#include "../src/encint.h"
#if defined(_WIN32)
# include <io.h>
# include <fcntl.h>
#endif
#include "getopt.h"
#include <math.h>

static void usage(char **_argv) {
  fprintf(stderr, "Usage: %s <input> <output>\n"
   "    <input> must be a YUV4MPEG file.\n"
   "    <output> will be <input> upsampled using od_img_upsample8.\n"
   , _argv[0]);
}

static const char *CHROMA_TAGS[4] = { " C420jpeg", "", " C422jpeg", " C444" };

/*Upsamples the reconstructed image to a reference image.
  TODO: Pipeline with reconstruction.*/
void od_img_upsample8(od_state *state, od_img *dimg, const od_img *simg) {
  int pli;
  for (pli = 0; pli < simg->nplanes; pli++) {
    const od_img_plane *siplane;
    od_img_plane *diplane;
    const unsigned char *src;
    unsigned char *dst;
    int xpad;
    int ypad;
    int w;
    int h;
    int x;
    int y;
    siplane = simg->planes + pli;
    diplane = dimg->planes + pli;
    xpad = OD_BUFFER_PADDING >> siplane->xdec;
    ypad = OD_BUFFER_PADDING >> siplane->ydec;
    w = simg->width >> siplane->xdec;
    h = simg->height >> siplane->ydec;
    src = siplane->data;
    dst = diplane->data - (diplane->ystride << 1)*ypad;
    for (y = -ypad; y < h + ypad + 3; y++) {
      /*Horizontal filtering:*/
      if (y < h + ypad) {
        unsigned char *buf;
        buf = state->ref_line_buf[y & 7];
        memset(buf - (xpad << 1), src[0], (xpad - 2) << 1);
        /*for (x = -xpad; x < -2; x++) {
          *(buf + (x << 1)) = src[0];
          *(buf + (x << 1 | 1)) = src[0];
        }*/
        *(buf - 4) = src[0];
        *(buf - 3) = OD_CLAMP255((31*src[0] + src[1] + 16) >> 5);
        *(buf - 2) = src[0];
        *(buf - 1) = OD_CLAMP255((36*src[0] - 5*src[1] + src[1] + 16) >> 5);
        buf[0] = src[0];
        buf[1] = OD_CLAMP255((20*(src[0] + src[1])
         - 5*(src[0] + src[2]) + src[0] + src[3] + 16) >> 5);
        buf[2] = src[1];
        buf[3] = OD_CLAMP255((20*(src[1] + src[2])
         - 5*(src[0] + src[3]) + src[0] + src[4] + 16) >> 5);
        for (x = 2; x < w - 3; x++) {
          buf[x << 1] = src[x];
          buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
           - 5*(src[x - 1] + src[x + 2]) + src[x - 2] + src[x + 3] + 16) >> 5);
        }
        buf[x << 1] = src[x];
        buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
         - 5*(src[x - 1] + src[x + 2]) + src[x - 2] + src[x + 2] + 16) >> 5);
        x++;
        buf[x << 1] = src[x];
        buf[x << 1 | 1] = OD_CLAMP255((20*(src[x] + src[x + 1])
         - 5*(src[x - 1] + src[x + 1]) + src[x - 2] + src[x + 1] + 16) >> 5);
        x++;
        buf[x << 1] = src[x];
        buf[x << 1 | 1] =
         OD_CLAMP255((36*src[x] - 5*src[x - 1] + src[x - 2] + 16) >> 5);
        x++;
        buf[x << 1] = src[w - 1];
        buf[x << 1 | 1] = OD_CLAMP255((31*src[w - 1] + src[w - 2] + 16) >> 5);
        memset(buf + (++x << 1), src[w - 1], (xpad - 1) << 1);
        /*for (x++; x < w + xpad; x++) {
          buf[x << 1] = src[w - 1];
          buf[x << 1 | 1]=src[w - 1];
        }*/
        if (y >= 0 && y + 1 < h) src += siplane->ystride;
      }
      /*Vertical filtering:*/
      if (y >= -ypad + 3) {
        if (y < 1 || y > h + 3) {
          OD_COPY(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          OD_COPY(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
        else {
          unsigned char *buf[6];
          buf[0] = state->ref_line_buf[(y - 5) & 7];
          buf[1] = state->ref_line_buf[(y - 4) & 7];
          buf[2] = state->ref_line_buf[(y - 3) & 7];
          buf[3] = state->ref_line_buf[(y - 2) & 7];
          buf[4] = state->ref_line_buf[(y - 1) & 7];
          buf[5] = state->ref_line_buf[(y - 0) & 7];
          OD_COPY(dst - (xpad << 1),
           state->ref_line_buf[(y - 3) & 7] - (xpad << 1),
           (w + (xpad << 1)) << 1);
          /*fprintf(stderr, "%3i: ", (y - 3) << 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            *(dst + x) = OD_CLAMP255((20*(*(buf[2] + x) + *(buf[3] + x))
             - 5*(*(buf[1] + x) + *(buf[4] + x))
             + *(buf[0] + x) + *(buf[5] + x) + 16) >> 5);
          }
          /*fprintf(stderr, "%3i: ", (y - 3) << 1 | 1);
          for (x = -xpad << 1; x < (w + xpad) << 1; x++) {
            fprintf(stderr, "%02X", *(dst + x));
          }
          fprintf(stderr, "\n");*/
          dst += diplane->ystride;
        }
      }
    }
  }
}


int main(int _argc, char **_argv) {
  const char *optstring = "";
  const struct option long_options[] = {
    { NULL, 0, NULL, 0 }
  };
  FILE *fin;
  FILE *fout;
  video_input vid;
  video_input_info info;
  int frameno;
  int pli;
  int xdec[3];
  int ydec[3];
  int w[3];
  int h[3];
  int long_option_index;
  int c;
  od_state state;
  daala_info dinfo;

  while ((c = getopt_long(_argc, _argv, optstring, long_options,
   &long_option_index)) != EOF) {
    switch (c) {
      default: {
        usage(_argv);
        exit(EXIT_FAILURE);
      }
    }
  }
  if (optind+2 != _argc) {
    usage(_argv);
    exit(EXIT_FAILURE);
  }
  fin = strcmp(_argv[optind], "-") == 0 ? stdin : fopen(_argv[optind], "rb");
  if (fin == NULL) {
    fprintf(stderr, "Unable to open '%s' for extraction.\n", _argv[optind]);
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "Opening %s as input...\n", _argv[optind]);
  if (video_input_open(&vid, fin) < 0) exit(EXIT_FAILURE);
  video_input_get_info(&vid, &info);

  daala_info_init(&dinfo);
  for (pli = 0; pli < 3; pli++) {
    xdec[pli] = pli && !(info.pixel_fmt&1);
    ydec[pli] = pli && !(info.pixel_fmt&2);
    h[pli] = info.pic_h >> ydec[pli];
    w[pli] = info.pic_w >> xdec[pli];

    dinfo.plane_info[pli].xdec = xdec[pli];
    dinfo.plane_info[pli].ydec = ydec[pli];
  }
  dinfo.nplanes = 3;
  dinfo.pic_height = h[0];
  dinfo.pic_width = w[0];

  od_state_init(&state, &dinfo);

  fout = strcmp(_argv[optind+1], "-") == 0 ? stdout : fopen(_argv[optind+1],
   "wb");
  if (fout == NULL) {
    fprintf(stderr, "Error opening output file \"%s\".\n", _argv[optind+1]);
    return 1;
  }
  fprintf(fout, "YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   info.pic_w*2, info.pic_h*2, (unsigned)info.fps_n,
   (unsigned)info.fps_d, info.par_n, info.par_d,
   CHROMA_TAGS[ydec[1] ? xdec[1] ? 0 : 2 : 3]);
  for (frameno = 0;; frameno++) {
    video_input_ycbcr in;
    int ret = 0;
    char tag[5];
    od_img *simg = &state.ref_imgs[1];
    od_img *dimg = &state.ref_imgs[0];
    int x, y;
    ret = video_input_fetch_frame(&vid, in, tag);
    if (ret == 0) break;
    for (pli = 0; pli < 3; pli++) {
      od_img_plane *siplane = simg->planes + pli;
      unsigned char *src = siplane->data;
      int src_stride = siplane->ystride;
      int plane_width = simg->width >> xdec[pli];
      int plane_height = simg->height >> ydec[pli];
      for (y = 0; y < h[pli]; y++) {
        for (x = 0; x < w[pli]; x++) {
          int cy = y + (int)(info.pic_y >> ydec[pli]);
          int cx = x + (int)(info.pic_x >> xdec[pli]);
          src[y*src_stride + x] = in[pli].data[cy*in[pli].stride + cx];
        }
      }
      /*From od_img_plane_copy_pad8*/
      /*Right side.*/
      for (x = w[pli]; x < plane_width; x++) {
        src = siplane->data + x - 1;
        for (y = 0; y < h[pli]; y++) {
          src[1] = (2*src[0] + (src - (src_stride & -(y > 0)))[0]
           + (src + (src_stride & -(y + 1 < h[pli])))[0] + 2) >> 2;
          src += src_stride;
        }
      }
      /*Bottom.*/
      src = siplane->data + src_stride*h[pli];
      for (y = h[pli]; y < plane_height; y++) {
        for (x = 0; x < plane_width; x++) {
          src[x] = (2*(src - src_stride)[x] + (src - src_stride)[x - (x > 0)]
           + (src - src_stride)[x + (x + 1 < plane_width)] + 2) >> 2;
        }
        src += src_stride;
      }
    }
    od_img_upsample8(&state, dimg, simg);
    fprintf(fout, "FRAME\n");
    for (pli = 0; pli < 3; pli++) {
      od_img_plane *diplane = dimg->planes + pli;
      unsigned char *dst = diplane->data;
      for (y = 0; y < 2*h[pli]; y++) {
        if (fwrite(dst + diplane->ystride*y, 2*w[pli], 1, fout) < 1) {
          fprintf(stderr, "Error writing to output.\n");
          return EXIT_FAILURE;
        }
      }
    }
    fprintf(stderr, "Completed frame %d.\n", frameno);
  }
  video_input_close(&vid);
  if (fout != stdout) fclose(fout);
  return EXIT_SUCCESS;
}
