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
#if defined(_WIN32)
# include <io.h>
# include <fcntl.h>
#endif
#include "getopt.h"

static void usage(char **_argv) {
  fprintf(stderr, "Usage: %s [options] <input> <output>\n"
   "    <input> must be a YUV4MPEG file.\n\n"
   "    Options:\n\n"
   "      --ds-algo <n>  Downsampling algorithm:\n"
   "                     0 (default) => 2x2 box average,\n"
   "                     1 => discard every other row and every other column.\n",
   _argv[0]);
}

static const char *CHROMA_TAGS[4] = { " C420jpeg", "", " C422jpeg", " C444" };

enum downsamplers {
  BOX_DOWNSAMPLER = 0,
  DISCARD_DOWNSAMPLER,
  DOWNSAMPLERS_MAX
};

int main(int _argc, char **_argv) {
  const char *optstring = "";
  const struct option long_options[] = {
    { "ds-algo", required_argument, NULL, 0 },
    { NULL, 0, NULL, 0 }
  };
  FILE *fin;
  FILE *fout;
  video_input vid;
  video_input_info info;
  int frameno;
  int pli;
  unsigned char *buf[3];
  int buf_stride[3];
  int xdec[3];
  int ydec[3];
  int w[3];
  int h[3];
  int long_option_index;
  int c;
  int ds_algo;

  ds_algo = BOX_DOWNSAMPLER;
  while ((c = getopt_long(_argc, _argv, optstring, long_options,
   &long_option_index)) != EOF) {
    switch (c) {
      case 0: {
        if (strcmp(long_options[long_option_index].name,"ds-algo") == 0) {
          ds_algo = atoi(optarg);
          if (ds_algo < 0 || ds_algo >= DOWNSAMPLERS_MAX) {
            fprintf(stderr, "Invalid value for --ds-algo");
            exit(1);
          }
        }
        break;
      }
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

  for (pli = 0; pli < 3; pli++) {
    xdec[pli] = pli && !(info.pixel_fmt&1);
    ydec[pli] = pli && !(info.pixel_fmt&2);
    h[pli] = info.pic_h >> ydec[pli];
    w[pli] = info.pic_w >> xdec[pli];

    if (h[pli] & 1) {
      fprintf(stderr, "Plane %d has odd height: %d\n", pli, h[pli]);
      return 1;
    }
    if (w[pli] & 1) {
      fprintf(stderr, "Plane %d has odd width: %d\n", pli, w[pli]);
      return 1;
    }

    buf[pli]=malloc(w[pli]/2*h[pli]/2*sizeof(**buf));
    buf_stride[pli]=w[pli]/2;
  }
  fout = strcmp(_argv[optind+1], "-") == 0 ? stdout : fopen(_argv[optind+1],
   "wb");
  if (fout == NULL) {
    fprintf(stderr, "Error opening output file \"%s\".\n", _argv[optind+1]);
    return 1;
  }
  fprintf(fout, "YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   info.pic_w/2, info.pic_h/2, (unsigned)info.fps_n,
   (unsigned)info.fps_d, info.par_n, info.par_d,
   CHROMA_TAGS[ydec[1] ? xdec[1] ? 0 : 2 : 3]);
  for (frameno = 0;; frameno++) {
    video_input_ycbcr in;
    int ret = 0;
    char tag[5];
    int x, y;
    ret = video_input_fetch_frame(&vid, in, tag);
    if (ret == 0) break;
    for (pli = 0; pli < 3; pli++) {
      unsigned char *src = in[pli].data;
      int src_stride = in[pli].stride;
      for (y = 0; y < h[pli]; y += 2) {
        for (x = 0; x < w[pli]; x += 2) {
          int cy = y + (info.pic_y >> ydec[pli]);
          int cx = x + (info.pic_x >> xdec[pli]);
          if (ds_algo == BOX_DOWNSAMPLER) {
            buf[pli][y/2*buf_stride[pli] + x/2] = (src[cy*src_stride + cx] +
              src[cy*src_stride + cx + 1] + src[(cy + 1)*src_stride + cx] +
              src[(cy + 1)*src_stride + cx + 1])/4;
          } else {
            buf[pli][y/2*buf_stride[pli] + x/2] = src[cy*src_stride + cx];
          }
        }
      }
    }
    fprintf(fout, "FRAME\n");
    for (pli = 0; pli < 3; pli++) {
      if (fwrite(buf[pli], w[pli]/2*h[pli]/2, 1, fout) < 1) {
        fprintf(stderr, "Error writing to output.\n");
        return EXIT_FAILURE;
      }
    }
    fprintf(stderr, "Completed frame %d.\n", frameno);
  }
  video_input_close(&vid);
  if (fout != stdout) fclose(fout);
  return EXIT_SUCCESS;
}
