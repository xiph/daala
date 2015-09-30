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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include <daala/codec.h>
#include <daala/daaladec.h>

#include <ogg/ogg.h>

#include "SDL.h"

#define ODS_NONE 0
#define ODS_HEADER 1
#define ODS_DATA 2

typedef struct {
  SDL_Window *screen;
  SDL_Renderer *renderer;
  daala_info di;
  daala_comment dc;
  ogg_sync_state oy;
  FILE *input;
  const char *input_path;
  ogg_stream_state os;
  daala_dec_ctx *dctx;
  SDL_Texture *texture;
  od_img img;
  int width;
  int height;
  int done;
  int od_state;
  int paused;
  int restart;
  int slow;
  int loop;
  int step;
  int fullscreen;
  int valid;
  int plane_mask;
} player_example;

static void ogg_to_daala_packet(daala_packet *dp, ogg_packet *op) {
  dp->packet     = op->packet;
  dp->bytes      = op->bytes;

  dp->b_o_s      = op->b_o_s;
  dp->e_o_s      = op->e_o_s;

  dp->granulepos = op->granulepos;
  dp->packetno   = op->packetno;
}

enum {
  OD_LUMA_MASK = 1 << 0,
  OD_CB_MASK = 1 << 1,
  OD_CR_MASK = 1 << 2,
  OD_ALL_MASK = OD_LUMA_MASK | OD_CB_MASK | OD_CR_MASK
};

static void img_to_rgb(SDL_Texture *tex, const od_img *img, int plane_mask);
static int next_plane(int plane_mask);
static void wait_to_refresh(uint32_t *previous_ticks, uint32_t ms_per_frame);

int player_example_init(player_example *player);
void build_yuv_to_rgb_table();
player_example *player_example_create();
int player_example_clear(player_example *player);
int player_example_free(player_example *player);
int player_example_reset(player_example *player);
int player_example_play(player_example *player);
int player_example_restart(player_example *player);
int player_example_open(player_example *player, char *path);

int player_example_daala_stream_init(player_example *player, int serial) {
  int ret;
  if (player == NULL) return -1;
  ret = ogg_stream_init(&player->os, serial);
  if (ret != 0) return -1;
  daala_info_init(&player->di);
  daala_comment_init(&player->dc);
  player->od_state = ODS_HEADER;
  return 0;
}

int player_example_daala_stream_clear(player_example *player) {
  if (player == NULL) return -1;
  daala_info_clear(&player->di);
  daala_comment_clear(&player->dc);
  if (player->dctx != NULL) {
    daala_decode_free(player->dctx);
    player->dctx = NULL;
  }
  ogg_stream_clear(&player->os);
  player->od_state = ODS_NONE;
  return 0;
}

int player_example_init(player_example *player) {
  if (player == NULL) return -1;
  player->screen = NULL;
  player->renderer = NULL;
  player->texture = NULL;
  player->width = 0;
  player->height = 0;
  player->paused = 0;
  player->restart = 0;
  player->slow = 0;
  player->loop = 0;
  player->done = 0;
  player->step = 0;
  player->fullscreen = 0;
  player->valid = 0;
  player->od_state = ODS_NONE;
  player->plane_mask = OD_ALL_MASK;
  return 0;
}

int player_example_clear(player_example *player) {
  if (player == NULL) return -1;
  memset(player, 0, sizeof(player_example));
  return 0;
}

player_example *player_example_create() {
  int ret;
  player_example *player;
  player = (player_example *)malloc(sizeof(player_example));
  if (player == NULL) return NULL;
  ret = player_example_init(player);
  if (ret != 0) {
    free(player);
    return NULL;
  }
  return player;
}

int player_example_free(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  ret = player_example_clear(player);
  if (ret == 0) {
    free(player);
    return 0;
  }
  return -1;
}

int player_example_reset(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  ret = player_example_clear(player);
  if (ret != 0) return -1;
  ret = player_example_init(player);
  if (ret != 0) return -1;
  return 0;
}

int player_example_open_input(player_example *player, const char *path) {
  if ((player == NULL) || ((path == NULL) || (path[0] == '\0'))) return -1;
  if ((path[0] == '-') && (path[1] == '\0')) {
    player->input = stdin;
  }
  else {
    player->input = fopen(path, "rb");
  }
  if (player->input == NULL) {
    player->input_path = "";
    return -1;
  }
  player->input_path = path;
  ogg_sync_init(&player->oy);
  return 0;
}

int player_example_close_input(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  if ((player->input == stdin) || (player->input == NULL)) return -1;
  ret = fclose(player->input);
  player->input = NULL;
  ogg_sync_clear(&player->oy);
  if (ret != 0) return -1;
  return 0;
}

int player_example_input_restart(player_example *player) {
  int ret;
  if (player == NULL) return -1;
  if (player->input == stdin) return -1;
  ret = player_example_close_input(player);
  if (ret != 0) return -1;
  ret = player_example_open_input(player, player->input_path);
  return ret;
}

void player_example_display_frame(player_example *player) {
  SDL_RenderClear(player->renderer);
  SDL_RenderCopy(player->renderer, player->texture, NULL, NULL);
  SDL_RenderPresent(player->renderer);
}

void player_example_handle_event(player_example *player, SDL_Event *event) {
  switch (event->type) {
    case SDL_QUIT: {
      player->done = 1;
      break;
    }
    case SDL_KEYDOWN: {
      switch (event->key.keysym.sym) {
        case SDLK_q: {
          player->done = 1;
          break;
        }
        case SDLK_s: {
          player->slow = !player->slow;
          break;
        }
        case SDLK_p: {
          player->plane_mask = next_plane(player->plane_mask);
          break;
        }
        case SDLK_l: {
          player->loop = !player->loop;
          break;
        }
        case SDLK_ESCAPE: {
          player->done = 1;
          break;
        }
        case SDLK_SPACE: {
          player->paused = !player->paused;
          player->step = 0;
          break;
        }
        case SDLK_RIGHT:
        case SDLK_PERIOD: {
          player->step = 1;
          player->paused = 1;
          break;
        }
        case SDLK_HOME:
        case SDLK_r: {
          player->restart = 1;
          if (player->paused) {
            player->step = 1;
          }
          break;
        }
        case SDLK_f: {
          player->fullscreen = !player->fullscreen;
          SDL_SetWindowFullscreen(player->screen,
              player->fullscreen ? SDL_WINDOW_FULLSCREEN_DESKTOP : 0);
          img_to_rgb(player->texture, &player->img, player->plane_mask);
          player_example_display_frame(player);
          break;
        }
        default: break;
      }
      break;
    }
    case SDL_WINDOWEVENT: {
      switch (event->window.event) {
        case SDL_WINDOWEVENT_RESIZED: {
          player_example_display_frame(player);
          break;
        }
        default: break;
      }
    }
  }
}

void player_example_wait_user_input(player_example *player) {
  SDL_Event event;
  if (SDL_WaitEvent(&event)) {
    player_example_handle_event(player, &event);
  }
}

void player_example_check_user_input(player_example *player) {
  SDL_Event event;
  while (SDL_PollEvent(&event)) {
    player_example_handle_event(player, &event);
  }
}

int player_example_play(player_example *player) {
  size_t bytes;
  char *buffer;
  int ret;
  ogg_page page;
  ogg_packet op;
  daala_packet dp;
  daala_setup_info *dsi;
  uint32_t ms_per_frame;
  uint32_t ticks = 0;
  ms_per_frame = 0;
  dsi = NULL;
  while (!player->done) {
    while (ogg_sync_pageout(&player->oy, &page) != 1) {
      buffer = ogg_sync_buffer(&player->oy, 4096);
      if (buffer == NULL) return -1;
      bytes = fread(buffer, 1, 4096, player->input);
      if (bytes > 0) {
        ret = ogg_sync_wrote(&player->oy, bytes);
        if (ret != 0) return -1;
      }
      else {
        if (!player->valid) {
          fprintf(stderr, "Invalid Ogg\n");
          exit(1);
        }
        if (player->od_state != ODS_NONE) {
          ret = player_example_daala_stream_clear(player);
          if (ret != 0) return -1;
        }
        if (player->input == stdin) {
          return 0;
        }
        if (player->loop == 1) {
          player_example_input_restart(player);
          continue;
        }
        for (;;) {
          player_example_wait_user_input(player);
          if (player->restart) {
            ret = player_example_input_restart(player);
            player->restart = 0;
            if (ret != 0) return -1;
            break;
          }
          if (player->done) {
            return 0;
          }
        }
      }
    }
    if (ogg_page_bos(&page)) {
      ret = player_example_daala_stream_init(player,
       ogg_page_serialno(&page));
      if (ret != 0) return -1;
    }
    ret = ogg_stream_pagein(&player->os, &page);
    if (ret != 0) return -1;
    while (ogg_stream_packetout(&player->os, &op) == 1) {
      ogg_to_daala_packet(&dp, &op);
      switch (player->od_state) {
        case ODS_HEADER: {
          ret =
           daala_decode_header_in(&player->di, &player->dc, &dsi, &dp);
          if (ret < 0) {
            if (memcmp(dp.packet, "fishead", dp.bytes)) {
              fprintf(stderr, "Ogg Skeleton streams not supported\n");
            }
            return -1;
          }
          if (ret != 0) break;
          player->dctx = daala_decode_alloc(&player->di, dsi);
          if (player->dctx == NULL) return -1;
          daala_setup_free(dsi);
          dsi = NULL;
          player->od_state = ODS_DATA;
          if (player->di.timebase_numerator
           && player->di.timebase_denominator) {
            ms_per_frame = 1000 /
             (player->di.timebase_numerator / player->di.timebase_denominator);
            ticks = SDL_GetTicks();
          }
          break;
        }
        case ODS_DATA: {
          if ((player->screen == NULL)
           || (player->width != player->di.pic_width)
           || (player->height != player->di.pic_height)) {
            player->width = player->di.pic_width;
            player->height = player->di.pic_height;
            player->screen = SDL_CreateWindow("Daala example player",
                SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                player->width, player->height,
                SDL_WINDOW_ALLOW_HIGHDPI);
            if (player->screen == NULL) return -1;
            player->renderer = SDL_CreateRenderer(player->screen, -1, 0);
            if (player->renderer == NULL) return -1;
            player->texture = SDL_CreateTexture(player->renderer,
                SDL_PIXELFORMAT_ARGB8888,
                SDL_TEXTUREACCESS_STREAMING,
                player->width, player->height);
            if (player->texture == NULL) return -1;
          }
          ret = daala_decode_packet_in(player->dctx, &player->img, &dp);
          if (ret != 0) return -1;
          player->valid = 1;
          if ((player->slow) && (!player->step)) {
            SDL_Delay(420);
          }
          player_example_check_user_input(player);
          while ((player->paused) && (!player->done)) {
            if (player->restart) {
              break;
            }
            if (player->step) {
              player->step = 0;
              break;
            }
            player_example_wait_user_input(player);
          }
          if ((!player->restart) && (!player->done)) {
            wait_to_refresh(&ticks, ms_per_frame);
            img_to_rgb(player->texture, &player->img, player->plane_mask);
            player_example_display_frame(player);
          }
          break;
        }
      }
    }
    if ((player->restart) || (ogg_page_eos(&page))) {
      ret = player_example_daala_stream_clear(player);
      if (ret != 0) return -1;
    }
    if (player->restart) {
      ret = player_example_input_restart(player);
      player->restart = 0;
      if (ret != 0) return -1;
    }
  }
  if (player->od_state != ODS_NONE) {
    ret = player_example_daala_stream_clear(player);
    if (ret != 0) return -1;
  }
  return 0;
}

int main(int argc, char *argv[]) {
  int ret;
  char *input;
  int start_paused;
  player_example *player;
  build_yuv_to_rgb_table();
  daala_log_init();
  if ((argc == 3) && (memcmp(argv[1], "-p", 2) == 0)) {
    start_paused = 1;
    input = argv[2];
  }
  else {
    if ((argc != 2) || ((argc == 2)
     && ((memcmp(argv[1], "-h", 2) == 0)
     || (memcmp(argv[1] + 1, "-h", 2) == 0)))) {
      fprintf(stderr, "usage: %s input.ogg\n%s\n", argv[0],
       "\nProgram Options:\n-p to start paused\n- to read from stdin\n\n"
       "Playback Control: \n"
       "r to restart\nl to loop\ns for slow\n. to step\nspace to pause\n"
       "p to switch planes\nq to quit");
      exit(1);
    } else if ((argc == 2)
            && memcmp(argv[1], "--version", 9) == 0
            && strlen(argv[1]) == strlen("--version")) {
      fprintf(stderr, "%s\n", daala_version_string());
      exit(0);
    }
    start_paused = 0;
    input = argv[1];
  }
  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    fprintf(stderr, "error: unable to init SDL\n");
    fprintf(stderr, "SDL_Error: %s\n", SDL_GetError());
    exit(1);
  }
  atexit(SDL_Quit);

  player = player_example_create();
  if (player == NULL) {
    fprintf(stderr, "player example error: create player\n");
    return -1;
  }
  ret = player_example_open_input(player, input);
  if (ret != 0) {
    fprintf(stderr, "player example error: could not open: %s\n", input);
    player_example_free(player);
    return -1;
  }
  if (start_paused == 1) {
    player->step = 1;
    player->paused = 1;
  }
  ret = player_example_play(player);
  if (ret != 0) {
    fprintf(stderr, "player example error: playback error\n");
    exit(1);
  }
  ret = player_example_free(player);

  return ret;
}

#define OD_MINI(a, b) ((a) < (b) ? (a) : (b))
#define OD_MAXI(a, b) ((a) > (b) ? (a) : (b))
#define OD_CLAMPI(a, b, c) (OD_MAXI(a, OD_MINI(b, c)))
#define OD_SIGNMASK(a) (-((a) < 0))
#define OD_FLIPSIGNI(a, b) (((a) + OD_SIGNMASK(b)) ^ OD_SIGNMASK(b))
#define OD_DIV_ROUND(x, y) (((x) + OD_FLIPSIGNI((y)  >>  1, x))/(y))
#define OD_CLAMP255(x) \
  ((unsigned char)((((x) < 0) - 1) & ((x) | -((x) > 255))))

/*Lookup Tables for YUV-to-RGB Conversion*/
static uint8_t rvalTable [256*256];
static uint8_t gvalTable [256*256*256];
static uint8_t bvalTable [256*256];

void build_yuv_to_rgb_table() {
  int y;
  int cr;
  int cb;
  int64_t yval;
  int64_t cbval;
  int64_t crval;
  unsigned rval;
  unsigned gval;
  unsigned bval;
  int plane_mask = OD_ALL_MASK;
  for (y = 0; y < 256; y ++) {
    yval = (plane_mask & OD_LUMA_MASK) * (y - 16)
           + (((plane_mask & OD_LUMA_MASK) ^ OD_LUMA_MASK) << 7);
    for (cr = 0; cr < 256; cr ++) {
      crval = ((plane_mask & OD_CR_MASK) >> 2) * (cr - 128);
      rval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 4490222169144LL*crval, 9745792000LL), 65535);
      rvalTable[y*256+cr] = rval >> 8;
    }
    for (cb = 0; cb < 256; cb ++) {
      cbval = ((plane_mask & OD_CB_MASK) >> 1) * (cb - 128);
      bval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 5290866304968LL*cbval, 9745792000LL), 65535);
      bvalTable[y*256+cb] = bval >> 8;
    }
    for (cb = 0; cb < 256; cb ++) {
      for (cr = 0; cr < 256; cr ++) {
        crval = ((plane_mask & OD_CR_MASK) >> 2) * (cr - 128);
        cbval = ((plane_mask & OD_CB_MASK) >> 1) * (cb - 128);
        gval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
         2916394880000LL*yval - 534117096223LL*cbval - 1334761232047LL*crval,
         9745792000LL), 65535);
        gvalTable[y*256*256+cb*256+cr] = gval >> 8;
      }
    }
  }
}

void img_to_rgb(SDL_Texture *texture, const od_img *img, int plane_mask) {
  unsigned char *y_row;
  unsigned char *cb_row;
  unsigned char *cr_row;
  unsigned char *y;
  unsigned char *cb;
  unsigned char *cr;
  int y_stride;
  int cb_stride;
  int cr_stride;
  int width;
  int height;
  int xdec;
  int ydec;
  int i;
  int j;
  unsigned char *pixels;
  int pitch;
  /*Assume both C planes are decimated.*/
  xdec = img->planes[1].xdec;
  ydec = img->planes[1].ydec;
  y_stride = img->planes[0].ystride;
  cb_stride = img->planes[1].ystride;
  cr_stride = img->planes[2].ystride;
  y_row = img->planes[0].data;
  cb_row = img->planes[1].data;
  cr_row = img->planes[2].data;
  /*Lock the texture in video memory for update.*/
  if (SDL_LockTexture(texture, NULL, (void**)&pixels, &pitch)) {
    fprintf(stderr, "Couldn't lock video texture!");
    exit(1);
  }
  /*The texture memory is only allocated for the cropped frame.  The
    od_img is rounded up to superblock. */
  if(SDL_QueryTexture(texture, NULL, NULL, &width, &height)){
    fprintf(stderr, "Couldn't query video texture!");
    exit(1);
  }
  width = OD_MINI(img->width, width);
  height = OD_MINI(img->height, height);
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
    for (i = 0; i < 4 * width;) {
      int64_t yval;
      int64_t cbval;
      int64_t crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      if (plane_mask == OD_ALL_MASK) {
        /*Uuse precomputed lookup table, only available for OD_ALL_MASK*/
        rval = rvalTable[(*y)*256+(*cr)];
        gval = gvalTable[(*y)*256*256+(*cb)*256+(*cr)];
        bval = bvalTable[(*y)*256+(*cb)];
      } else {
        yval = (plane_mask & OD_LUMA_MASK) * (*y - 16)
         + (((plane_mask & OD_LUMA_MASK) ^ OD_LUMA_MASK) << 7);
        cbval = ((plane_mask & OD_CB_MASK) >> 1) * (*cb - 128);
        crval = ((plane_mask & OD_CR_MASK) >> 2) * (*cr - 128);
        /*This is intentionally slow and very accurate.*/
        rval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
         2916394880000LL*yval + 4490222169144LL*crval, 9745792000LL),
         65535) >> 8;
        gval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
         2916394880000LL*yval - 534117096223LL*cbval - 1334761232047LL*crval,
         9745792000LL), 65535) >> 8;
        bval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
         2916394880000LL*yval + 5290866304968LL*cbval, 9745792000LL),
         65535) >> 8;
      }
      *(pixels + pitch*j + i++) = (unsigned char)(bval);
      *(pixels + pitch*j + i++) = (unsigned char)(gval);
      *(pixels + pitch*j + i++) = (unsigned char)(rval);
      *(pixels + pitch*j + i++) = 0;
      dc = ((y - y_row) & 1) | (1 - xdec);
      y++;
      cb += dc;
      cr += dc;
    }
    y_row += y_stride;
    dc = -((j & 1) | (1 - ydec));
    cb_row += dc & cb_stride;
    cr_row += dc & cr_stride;
  }
  SDL_UnlockTexture(texture);
}

int next_plane(int plane_mask) {
  return OD_MINI(plane_mask << 1, OD_ALL_MASK) >>
   ((plane_mask == OD_ALL_MASK) << 1);
}

void wait_to_refresh(uint32_t *previous_ticks, uint32_t ms_per_frame)
{
  uint32_t tmp;
  /* play dumb until we've parsed the frame rate */
  if (!*previous_ticks)
      return;

  tmp = SDL_GetTicks();
  while (tmp - *previous_ticks < ms_per_frame) {
      SDL_Delay(ms_per_frame - (tmp - *previous_ticks));
      tmp = SDL_GetTicks();
  }
  *previous_ticks = tmp;
}
