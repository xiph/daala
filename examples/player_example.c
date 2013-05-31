/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

#include <stdio.h>
#include <assert.h>

#include <daala/codec.h>
#include <daala/daaladec.h>

#include "SDL.h"

#define PACKET_START 0
#define PACKET_HEADER 1
#define PACKET_DATA 2
#define PACKET_DONE 3

#define OD_MINI(_a,_b)      ((_a)<(_b)?(_a):(_b))
#define OD_MAXI(_a,_b)      ((_a)>(_b)?(_a):(_b))
#define OD_CLAMPI(_a,_b,_c) (OD_MAXI(_a,OD_MINI(_b,_c)))
#define OD_SIGNMASK(_a)     (-((_a)<0))
#define OD_FLIPSIGNI(_a,_b) ((_a)+OD_SIGNMASK(_b)^OD_SIGNMASK(_b))
#define OD_DIV_ROUND(_x,_y) (((_x)+OD_FLIPSIGNI((_y)>>1,_x))/(_y))
#define OD_CLAMP255(_x)     ((unsigned char)((((_x)<0)-1)&((_x)|-((_x)>255))))

void img_to_rgb(SDL_Surface *surf, const od_img *img);

int main(int argc, char *argv[]) {
  SDL_Surface *screen;
  SDL_Event event;
  int done;
  daala_info di;
  daala_comment dc;
  daala_setup_info *dsi;
  ogg_sync_state oy;
  ogg_page page;
  int ret;
  FILE *input;
  size_t bytes;
  int eof;
  int packet_state;
  ogg_stream_state os;
  ogg_packet packet;
  daala_dec_ctx *dctx;
  SDL_Surface *surf;
  od_img img;
  int frame;
  int paused;
  int step;

  if (argc != 2) {
    fprintf(stderr, "usage: %s input.ogg\n", argv[0]);
    exit(1);
  }

  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    fprintf(stderr, "error: unable to init SDL");
    exit(1);
  }
  atexit(SDL_Quit);

  screen = NULL;
  surf = NULL;
  paused = 1;

restart:
  ret = ogg_sync_init(&oy);
  assert(ret == 0);
  daala_info_init(&di);
  daala_comment_init(&dc);
  dsi = NULL;
  packet_state = PACKET_START;
  frame = 0;
  done = 0;
  step = 0;

  input = fopen(argv[1], "rb");
  assert(input != NULL);
  
  eof = 0;
  while (!eof) {
    while (!eof && ogg_sync_pageout(&oy, &page) != 1) {
      char *buffer = ogg_sync_buffer(&oy, 4096);
      assert(buffer != NULL);
      bytes = fread(buffer, 1, 4096, input);
      if (bytes > 0) {
        ret = ogg_sync_wrote(&oy, bytes);
        assert(ret == 0);
      } else {
        eof = 1;
      }
    }

    if (packet_state == PACKET_START) {
      ret = ogg_stream_init(&os, ogg_page_serialno(&page));
      assert(ret == 0);
      packet_state = PACKET_HEADER;
    }

    ret = ogg_stream_pagein(&os, &page);
    assert(ret == 0);

    while (!done && ogg_stream_packetout(&os, &packet) == 1) {
      switch (packet_state) {
      case PACKET_HEADER:
        ret = daala_decode_header_in(&di, &dc, &dsi, &packet);
        assert(ret >= 0);
        if (ret != 0) break;
        dctx = daala_decode_alloc(&di, dsi);
        assert(dctx != NULL);
        /*daala_setup_free(dsi);*/
        packet_state = PACKET_DATA;
      case PACKET_DATA:
        if (screen == NULL) {
          screen = SDL_SetVideoMode(di.pic_width, di.pic_height, 24, 
           SDL_SWSURFACE | SDL_DOUBLEBUF);
          assert(screen != NULL);
          surf = SDL_GetVideoSurface();
          assert(surf != NULL);
        }

        ret = daala_decode_packet_in(dctx, &img, &packet);
        assert(ret == 0);

        SDL_LockSurface(surf);
        img_to_rgb(surf, &img);
        SDL_UnlockSurface(surf);
        SDL_Flip(screen);

        if (step) {
          step = 0;
          paused = 1;
        }

        do {
          while (SDL_PollEvent(&event)) {
            switch (event.type) {
            case SDL_QUIT:
              done = 1;
              break;
            case SDL_KEYDOWN:
              switch (event.key.keysym.sym) {
              case SDLK_ESCAPE:
                done = 1;
                break;
              case SDLK_SPACE:
                paused = !paused;
                break;
              case SDLK_PERIOD:
                step = 1;
                paused = 0;
                break;
              case SDLK_r:
                done = 0;
                step = 0;
                goto restart;
              default:
                break;
              }
            }
          }
        } while (paused && !done);
      }
    }
  }

  while (!done) {
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
      case SDL_QUIT:
        done = 1;
        break;
      case SDL_KEYDOWN:
        if (event.key.keysym.sym == SDLK_ESCAPE) done = 1;
        if (event.key.keysym.sym == SDLK_r) {
          paused = 0;
          done = 0;
          step = 0;
          goto restart;
        }
        break;
      }
    }
  }

  return 0;
}

void img_to_rgb(SDL_Surface *surf, const od_img *img) {
  unsigned char *y_row;
  unsigned char *cb_row;
  unsigned char *cr_row;
  unsigned char *y;
  unsigned char *cb;
  unsigned char *cr;
  int            y_stride;
  int            cb_stride;
  int            cr_stride;
  int            width;
  int            height;
  int            xdec;
  int            ydec;
  int            i;
  int            j;
  unsigned char *pixels;
  int pitch;
  pixels = (unsigned char *)surf->pixels;
  pitch = surf->pitch;
  width = img->width;
  height = img->height;
  /*Assume both C planes are decimated.*/
  xdec = img->planes[1].xdec;
  ydec = img->planes[1].ydec;
  y_stride = img->planes[0].ystride;
  cb_stride = img->planes[1].ystride;
  cr_stride = img->planes[2].ystride;
  y_row = img->planes[0].data;
  cb_row = img->planes[1].data;
  cr_row = img->planes[2].data;
  /*Chroma up-sampling is just done with a box filter.
    This is very likely what will actually be used in practice on a real
     display, and also removes one more layer to search in for the source of
     artifacts.
    As an added bonus, it's dead simple.*/
  for (j = 0; j < height; j++) {
    int dc;
    y = y_row;
    cb = cb_row;
    cr = cr_row;
    for (i = 0; i < 3 * width;) {
      int64_t  yval;
      int64_t  cbval;
      int64_t  crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      yval=*y-16;
      cbval=*cb-128;
      crval=*cr-128;
      /*This is intentionally slow and very accurate.*/
      rval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval+4490222169144LL*crval,9745792000LL),65535);
      gval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval-534117096223LL*cbval-1334761232047LL*crval,
       9745792000LL),65535);
      bval=OD_CLAMPI(0,(int32_t)OD_DIV_ROUND(
       2916394880000LL*yval+5290866304968LL*cbval,9745792000LL),65535);
      *(pixels + pitch*j + i++)=(unsigned char)(bval>>8);
      *(pixels + pitch*j + i++)=(unsigned char)(gval>>8);
      *(pixels + pitch*j + i++)=(unsigned char)(rval>>8);
      dc = y - y_row & 1 | 1 - xdec;
      y++;
      cb += dc;
      cr += dc;
    }
    y_row += y_stride;
    dc = -(j & 1 | 1 - ydec);
    cb_row += dc & cb_stride;
    cr_row += dc & cr_stride;
  }
}
