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
   "      --plane <n>        Plane to preserve, all other planes will\n"
   "                         be blanked, default = 0\n"
   "      --output-mono-y4m  Output as a 4:0:0 YUV4MPEG instead of keeping\n"
   "                         the original subsampling and blanking the\n"
   "                         other planes (only implemented for --plane 0).\n",
   _argv[0]);
}

static const char *CHROMA_TAGS[4] = { " C420jpeg", "", " C422jpeg", " C444" };

int main(int _argc, char **_argv) {
  const char *optstring = "";
  const struct option long_options[] = {
    { "plane", required_argument, NULL, 0 },
    { "output-mono-y4m", no_argument, NULL, 0 },
    { NULL, 0, NULL, 0 }
  };
  FILE *fin;
  FILE *fout;
  video_input vid;
  video_input_info info;
  int frameno;
  int pli;
  int w[3];
  int h[3];
  int xdec[3];
  int ydec[3];
  unsigned char *blank;
  int long_option_index;
  int c;
  int output_mono_y4m;
  int preserved_plane;

  output_mono_y4m = 0;
  preserved_plane = 0;
  while ((c = getopt_long(_argc, _argv, optstring, long_options,
   &long_option_index)) != EOF) {
    switch (c) {
      case 0: {
        if (strcmp(long_options[long_option_index].name,
         "output-mono-y4m") == 0) {
          output_mono_y4m = 1;
        } else if (strcmp(long_options[long_option_index].name,
         "plane") == 0) {
          preserved_plane = atoi(optarg);
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
  if (output_mono_y4m && preserved_plane != 0) {
    fprintf(stderr, "--output-mono-y4m is only supported for plane 0.\n");
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
  fout = strcmp(_argv[optind+1], "-") == 0 ? stdout : fopen(_argv[optind+1],
   "wb");
  if (fout == NULL) {
    fprintf(stderr, "Error opening output file \"%s\".\n", _argv[optind+1]);
    return 1;
  }
  for (pli = 0; pli < 3; pli++) {
    xdec[pli] = pli && !(info.pixel_fmt&1);
    ydec[pli] = pli && !(info.pixel_fmt&2);
    w[pli] = info.pic_w >> xdec[pli];
    h[pli] = info.pic_h >> ydec[pli];
  }
  blank = malloc(w[0]*h[0]);
  memset(blank, 128, w[0]*h[0]);
  fprintf(fout, "YUV4MPEG2 W%i H%i F%i:%i Ip A%i:%i%s\n",
   info.pic_w, info.pic_h, (unsigned)info.fps_n,
   (unsigned)info.fps_d, info.par_n, info.par_d,
   output_mono_y4m ? " Cmono" : CHROMA_TAGS[ydec[1] ? xdec[1] ? 0 : 2 : 3]);
  for (frameno = 0;; frameno++) {
    video_input_ycbcr in;
    int ret = 0;
    char tag[5];
    unsigned char *src;
    int src_stride;
    int i;
    ret = video_input_fetch_frame(&vid, in, tag);
    if (ret == 0) break;
    fprintf(fout, "FRAME\n");

    for (pli = 0; pli < 3; pli++) {
      src = in[pli].data;
      src_stride = in[pli].stride;
      if (pli == preserved_plane) {
        for (i = 0; i < h[pli]; i++) {
          if (fwrite(src + (i + (info.pic_y >> ydec[pli]))*src_stride +
           (info.pic_x >> xdec[pli]), w[pli], 1, fout) < 1) {
            fprintf(stderr, "Error writing to output.\n");
            return EXIT_FAILURE;
          }
        }
      } else if (!output_mono_y4m) {
        if (fwrite(blank, w[pli]*h[pli], 1, fout) < 1) {
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
