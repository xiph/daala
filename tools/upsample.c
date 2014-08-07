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
#include "../src/state.h"
#if defined(_WIN32)
# include <io.h>
# include <fcntl.h>
#endif
#include "getopt.h"
#include <math.h>

static void usage(char **_argv) {
  fprintf(stderr, "Usage: %s <input> <output>\n"
   "    <input> must be a YUV4MPEG file.\n"
   "    <output> will be <input> upsampled using od_state_upsample8.\n"
   , _argv[0]);
}

static const char *CHROMA_TAGS[4] = { " C420jpeg", "", " C422jpeg", " C444" };

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
    od_img *simg = &state.io_imgs[OD_FRAME_REC];
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
    od_state_upsample8(&state, dimg, simg);
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
