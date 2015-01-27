/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

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
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../block_size.h"
#include "../filter.h"

#define IMAGES (10)

#define ROWS (32)
#define COLS (32)

#define PRINT_BSIZE  (0)
#define PRINT_LUMA   (0)
#define PRINT_CHROMA (0)
#define ZERO_FILTERS (0)

#define BX (0)
#define BY (0)

void randomize_bsize(unsigned char *bsize, int stride, int w, int h) {
  int i;
  int j;
  for (i=0; i<h; i++) {
    for (j=0; j<w; j++) {
      bsize[i*stride+j]=rand()&1;
    }
  }
  for (i = 0; i<h; i+=2) {
    for (j = 0; j<w; j+=2) {
      if ((rand() & 3) == 0) {
        int x;
        int y;
        for (x = 0;x < 2; x++) {
          for (y = 0; y < 2; y++) {
            bsize[(i + x)*stride + j + y]=OD_BLOCK_16X16;
          }
        }
      }
    }
  }
  for (i = 0; i < h; i += 4) {
    for (j = 0; j < w; j += 4) {
      if ((rand() & 7) == 0) {
        int x;
        int y;
        for (x = 0;x < 4; x++) {
          for (y = 0; y< 4; y++) {
            bsize[(i + x)*stride + j + y]=OD_BLOCK_32X32;
          }
        }
      }
    }
  }
}

#if PRINT_LUMA|PRINT_CHROMA
static void print_block(const char *label, od_coeff *_c,int _stride,int _sz,
 int _x,int _y) {
  int j;
  int i;
  printf("%s\n",label);
  for (j=0;j<_sz;j++) {
    for (i=0;i<_sz;i++) {
      printf("%s%4i",i>0?",":"",_c[_stride*(_sz*_y+j)+_sz*_x+i]);
    }
    printf("\n");
  }
}
#endif

#define OD_BLOCK_SIZE4x4_DEC(bsize, bstride, bx, by, dec) \
 OD_MAXI(0, \
  OD_BLOCK_SIZE4x4(bsize, bstride, (bx) << (dec), (by) << (dec)) - (dec))

static void filter_cols(od_coeff *c, unsigned char *bsize, int bstride,
 int inv, int dec) {
  int i;
  int j;
  int k;
  int l;
  int x;
  int f;
  unsigned char d;
  unsigned char n;
  od_coeff s[4 << OD_NBSIZES];
  od_coeff t[4 << OD_NBSIZES];
  /*Apply the column filter.*/
  for (i = 0; i < COLS << 3 - dec; i++) {
    for (j = 1; j < ROWS << 3 - dec; j++) {
      d = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, i, j, dec);
      if (!(j & ((1 << d) - 1))) {
        f = d;
        n = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, i, j - 1, dec);
        /*If current block is larger, then search *my* neighbors.*/
        if (d > n) {
          x = i & ~((1 << d) - 1);
          for (k = 0; k < 1 << d; k++) {
            f = OD_MINI(f,
             OD_BLOCK_SIZE4x4_DEC(bsize, bstride, x + k, j - 1, dec));
          }
        }
        /*If neighbor is larger, then search *its* neighbors.*/
        if (d < n) {
          x = i & ~((1 << n) - 1);
          for (k = 0; k < 1 << n; k++) {
            f = OD_MINI(f,
             OD_BLOCK_SIZE4x4_DEC(bsize, bstride, x + k, j, dec));
          }
        }
        /*Apply the filter 4 columns at a time.*/
        for (k = 0; k < 4; k++) {
          for (l = 0; l < 4 << f; l++) {
            s[l] = c[(COLS << 5 - dec)*(4*j - (2 << f) + l) + 4*i + k];
          }
          (*(inv ? OD_POST_FILTER : OD_PRE_FILTER)[f])(t, s);
          for (l = 0; l < 4 << f; l++) {
            c[(COLS << 5 - dec)*(4*j - (2 << f) + l) + 4*i + k] = t[l];
#if ZERO_FILTERS
            c[(COLS << 5 - dec)*(4*j - (2 << f) + l) + 4*i + k] = 0;
#endif
          }
        }
      }
    }
  }
}

static void filter_rows(od_coeff *c, unsigned char *bsize, int bstride,
 int inv, int dec) {
  int i;
  int j;
  int k;
  int l;
  int y;
  int f;
  unsigned char d;
  unsigned char n;
  od_coeff t[4 << OD_NBSIZES];
  /*Apply the row filter.*/
  for (j = 0; j < ROWS << 3 - dec; j++) {
    for (i = 1; i < COLS << 3 - dec; i++) {
      d = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, i, j, dec);
      if (!(i & ((1 << d) - 1))) {
        f = d;
        n = OD_BLOCK_SIZE4x4_DEC(bsize, bstride, i - 1, j, dec);
        /*If current block is larger, then search *my* neighbors.*/
        if (d > n) {
          y = j & ~((1 << d) - 1);
          for (k = 0; k < 1 << d; k++) {
            f = OD_MINI(f,
             OD_BLOCK_SIZE4x4_DEC(bsize, bstride, i - 1, y + k, dec));
          }
        }
        /*If neighbor is larger, then search *its* neighbors.*/
        if (d < n) {
          y = j & ~((1 << n) - 1);
          for (k = 0; k < 1 << n; k++) {
            f = OD_MINI(f,
             OD_BLOCK_SIZE4x4_DEC(bsize, bstride, i, y + k, dec));
          }
        }
        /*Apply the filter 4 rows at a time.*/
        for (k = 0; k < 4; k++) {
          od_coeff *s;
          s = &c[(COLS << 5 - dec)*(4*j + k) + 4*i - (2 << f)];
          (*(inv ? OD_POST_FILTER : OD_PRE_FILTER)[f])(t, s);
          for (l = 0; l < 4 << f; l++) {
            s[l] = t[l];
#if ZERO_FILTERS
            s[l] = 0;
#endif
          }
        }
      }
    }
  }
}

static int test_plane(const char *p, od_coeff *ref, od_coeff *luma, int dec) {
  int j;
  int i;
  for (j = 0; j < ROWS << 5 - dec; j++) {
    for (i = 0; i < COLS << 5 - dec; i++) {
      if (luma[(COLS << 5 - dec)*j + i] != ref[(COLS << 5 - dec)*j + i]) {
        printf("Failed!\n  Error on %s plane, %i not equal %i at %i %i\n",
         p, luma[(COLS << 5 - dec)*j + i], ref[(COLS << 5 - dec)*j + i], i, j);
        return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
  unsigned int   seed;
  int ret;
  od_coeff *luma;
  od_coeff *chroma;
  unsigned char *bsize;
  int bstride;
  od_coeff *luma_ref;
  od_coeff *chroma_ref;
  od_coeff *luma_test;
  od_coeff *chroma_test;
  int i;
  int j;
  int img;

  if (argc > 2) {
    fprintf(stderr, "Usage: %s [<seed>]\n", argv[0]);
    return EXIT_FAILURE;
  }
  if (argc > 1) {
    seed = atoi(argv[1]);
  }
  else {
    const char *env_seed;
    env_seed = getenv("SEED");
    if (env_seed) {
      seed = atoi(env_seed);
    }
    else {
      seed = time(NULL);
    }
  }
  srand(seed);
  fprintf(stderr, "Testing pre/post-filters... Random seed: %u (%.4X).\n",
   seed, rand() & 65535);

  ret = EXIT_SUCCESS;

  /*Allocate memory for reference data.*/
  luma = (od_coeff *)_ogg_malloc(COLS*32*ROWS*32*sizeof(*luma));
  chroma = (od_coeff *)_ogg_malloc(COLS*16*ROWS*16*sizeof(*chroma));
  luma_ref = (od_coeff *)_ogg_malloc(COLS*32*ROWS*32*sizeof(*luma_ref));
  chroma_ref = (od_coeff *)_ogg_malloc(COLS*16*ROWS*16*sizeof(*chroma_ref));
  luma_test = (od_coeff *)_ogg_malloc(COLS*32*ROWS*32*sizeof(*luma_test));
  chroma_test = (od_coeff *)_ogg_malloc(COLS*16*ROWS*16*sizeof(*chroma_test));
  bsize =
   (unsigned char *)_ogg_malloc((COLS*4 + 2)*(ROWS*4 + 2)*sizeof(*bsize));
  bstride = COLS*4 + 2;

  /*Set the block size decision of the image padding to 32x32.*/
  for (j = 0; j < ROWS*4 + 2; j++){
    bsize[bstride*j + 0] = OD_BLOCK_32X32;
    bsize[bstride*j + (COLS*4 + 2) - 1] = OD_BLOCK_32X32;
  }
  for (i = 0; i < COLS*4 + 2; i++) {
    bsize[bstride*0 + i] = OD_BLOCK_32X32;
    bsize[bstride*((ROWS*4 + 2) - 1) + i] = OD_BLOCK_32X32;
  }
  bsize += bstride*1 + 1;

  for (img = 1; img <= IMAGES; img++) {
    printf("Test #%2i: Generating random %ix%i luma, %ix%i chroma... ",
     img, COLS*32, ROWS*32, COLS*16, ROWS*16);
    fflush(stdout);

    /*Generate a random block size decision.*/
    randomize_bsize(bsize, bstride, COLS*4, ROWS*4);

#if PRINT_BSIZE
    /* Print the decisions. */
    for (j = -1; j < ROWS*4 + 1; j++) {
      for (i = -1; i < COLS*4 + 1; i++) {
        printf("%s%i", i > 0 ? "," : "", bsize[bstride*j + i]);
      }
      printf("\n");
    }
#endif

    /* Generate a random luma plane. */
    for (j = 0; j < ROWS*32; j++) {
      for (i = 0; i < COLS*32; i++) {
        luma[COLS*32*j + i] = (od_coeff)rand()%1024 - 512;
        luma_ref[COLS*32*j + i] = luma[COLS*32*j + i];
        luma_test[COLS*32*j + i] = luma[COLS*32*j + i];
      }
    }
#if PRINT_LUMA
    print_block("Original luma", luma, COLS*32, 32, BX, BY);
#endif
    /*Apply reference prefilter.*/
    filter_cols(luma_ref, bsize, bstride, 0, 0);
    filter_rows(luma_ref, bsize, bstride, 0, 0);
#if PRINT_LUMA
    print_block("Pre-filtered luma (ref)",luma_ref, COLS*32, 32, BX, BY);
#endif
    /*Apply per-block prefilter.*/
    for (j = 0; j < ROWS; j++) {
      for (i = 0; i < COLS; i++) {
        od_apply_prefilter(luma_test, COLS*32, i, j, 3, bsize, bstride, 0, 0,
         (i > 0 ? OD_LEFT_EDGE : 0) | (j < ROWS - 1 ? OD_BOTTOM_EDGE : 0));
      }
    }
#if PRINT_LUMA
    print_block("Pre-filtered luma (test)", luma_test, COLS*32, 32, BX, BY);
#endif
    /*Test the prefiltered luma.*/
    if (test_plane("pre-filtered luma", luma_ref, luma_test, 0)) {
      ret = EXIT_FAILURE;
      break;
    }
    /*Apply reference postfilter.*/
    filter_rows(luma_ref, bsize, bstride, 1, 0);
    filter_cols(luma_ref, bsize, bstride, 1, 0);
#if PRINT_LUMA
    print_block("Post-filtered luma (ref)", luma_ref, COLS*32, 32, BX, BY);
#endif
    /*Apply per-block postfilter.*/
    for (j = 0; j < ROWS; j++) {
      for (i = 0; i < COLS; i++) {
        od_apply_postfilter(luma_test, COLS*32, i, j, 3, bsize, bstride, 0, 0,
         (i < COLS - 1 ? OD_RIGHT_EDGE : 0) | (j > 0 ? OD_TOP_EDGE : 0));
      }
    }
#if PRINT_LUMA
    print_block("Post-filtered luma (test)", luma_test, COLS*32, 32, BX, BY);
#endif
    /*Test the post filtered luma.*/
    if (test_plane("post-filtered luma", luma, luma_test, 0)) {
      ret = EXIT_FAILURE;
      break;
    }
    /*Test that the reference filter algorithm is correct.*/
    if (test_plane("reference luma", luma, luma_ref, 0)) {
      ret = EXIT_FAILURE;
      break;
    }

    /* Generate a random chroma plane. */
    for (j = 0; j < ROWS*16; j++) {
      for (i = 0; i < COLS*16; i++) {
        chroma[COLS*16*j + i] = (od_coeff)rand()%1024 - 512;
        chroma_ref[COLS*16*j + i] = chroma[COLS*16*j + i];
        chroma_test[COLS*16*j + i] = chroma[COLS*16*j + i];
      }
    }
#if PRINT_CHROMA
    print_block("Original chroma", chroma, COLS*16, 16, BX, BY);
#endif
    /*Apply reference prefilter.*/
    filter_cols(chroma_ref, bsize, bstride, 0, 1);
    filter_rows(chroma_ref, bsize, bstride, 0, 1);
#if PRINT_CHROMA
    print_block("Pre-filtered chroma (ref)", chroma_ref, COLS*16, 16, BX, BY);
#endif
    /*Apply per-block prefilter.*/
    for (j = 0; j < ROWS; j++) {
      for (i = 0; i < COLS; i++) {
        od_apply_prefilter(chroma_test, COLS*16, i, j, 3, bsize, bstride, 1, 1,
         (i > 0 ? OD_LEFT_EDGE : 0) | (j < ROWS - 1 ? OD_BOTTOM_EDGE : 0));
      }
    }
#if PRINT_CHROMA
    print_block("Pre-filtered chroma (test)", chroma_test, COLS*16, 16, BX, BY);
#endif
    /*Test the prefiltered chroma.*/
    if (test_plane("pre-filtered chroma", chroma_ref, chroma_test, 1)) {
      ret = EXIT_FAILURE;
      break;
    }
    /*Apply reference postfilter.*/
    filter_rows(chroma_ref, bsize, bstride, 1, 1);
    filter_cols(chroma_ref, bsize, bstride, 1, 1);
#if PRINT_CHROMA
    print_block("Post-filtered chroma (ref)", chroma_ref, COLS*16, 16, BX, BY);
#endif
    /*Apply per-block postfilter.*/
    for (j = 0; j < ROWS; j++) {
      for (i = 0; i < COLS; i++) {
        od_apply_postfilter(chroma_test, COLS*16, i, j, 3, bsize, bstride, 1, 1,
         (i < COLS - 1 ? OD_RIGHT_EDGE : 0) | (j > 0 ? OD_TOP_EDGE : 0));
      }
    }
#if PRINT_CHROMA
    print_block("Post-filtered chroma (test)", chroma_test, COLS*16, 16, BX, BY);
#endif
    /*Test the post filtered luma.*/
    if (test_plane("post-filtered chroma", chroma, chroma_test, 1)) {
      ret = EXIT_FAILURE;
      break;
    }
    /*Test that the reference filter algorithm is correct.*/
    if (test_plane("reference chroma", chroma, chroma_ref, 1)) {
      ret = EXIT_FAILURE;
      break;
    }
    printf("Passed!\n");
  }

  /*Free allocated memory.*/
  bsize -= bstride*1 + 1;
  _ogg_free(luma);
  _ogg_free(chroma);
  _ogg_free(luma_ref);
  _ogg_free(chroma_ref);
  _ogg_free(luma_test);
  _ogg_free(chroma_test);
  _ogg_free(bsize);

  return ret;
}
