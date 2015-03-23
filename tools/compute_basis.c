
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../src/dct.h"
#include "../src/filter.h"

void usage(char *str) {
  fprintf(stderr,
   "usage: %s <log size - 2> <coeff | coeff420 | mag | mag420>\n", str);
}

#define OD_BASIS_SIZE (3*OD_BSIZE_MAX)
#define OD_BASIS_PULSE 1024

/* Computes the synthesis basis functions and their magnitudes. The lapping
   filter is selected exactly as it would in the codec, including the 4x4
   blocks that have 8x8 lapping on one side and 4x4 on the other. The exact
   lapping filters are controlled by the "left" and "right" variables. */
int main(int argc, char **argv) {
  int i;
  int j;
  int n;
  int ln;
  int left;
  int right;
  int magnitude;
  int dec;
  od_coeff x0[OD_BASIS_SIZE];
  od_coeff y0[OD_BASIS_SIZE];
  od_coeff *x;
  od_coeff *y;
  if (argc != 3) {
    usage(argv[0]);
    return 1;
  }
  ln = atoi(argv[1]);
  dec = 0;
  if (strcmp("mag", argv[2]) == 0) {
    magnitude = 1;
  }
  else if (strcmp("coeff", argv[2]) == 0) {
    magnitude = 0;
  }
  else if (strcmp("mag420", argv[2]) == 0) {
    magnitude = 1;
    dec = 1;
  }
  else if (strcmp("coeff420", argv[2]) == 0) {
    magnitude = 0;
    dec = 1;
  }
  else {
    usage(argv[0]);
    return 1;
  }
  n = 4 << ln;
  /* The first lapping filter is applied based on a larger (unsplit)
     block size. */
  left = OD_FILT_SIZE(OD_MINI(OD_NBSIZES - 1, ln + 1), dec);
  right = OD_FILT_SIZE(ln, dec);
  for (i = 0; i < n; i++) {
    OD_CLEAR(x0, OD_BASIS_SIZE);
    OD_CLEAR(y0, OD_BASIS_SIZE);
    x = &x0[n];
    y = &y0[n];
    x[i] = OD_BASIS_PULSE;
    OD_IDCT_1D[ln](y, 1, x);
    /* We need to apply left before right for 4x4 because the wider lapping
       is always applied first. */
    OD_POST_FILTER[left](y - (2 << left), y - (2 << left));
    OD_POST_FILTER[right](y + n - (2 << right), y + n - (2 << right));
    if (magnitude) {
      double sum;
      sum = 0;
      for (j = 0; j < n + (2 << left) + (2 << right); j++) {
        sum += y[j - (2 << left)]*y[j - (2 << left)];
      }
      printf("%f ", sqrt(sum)/OD_BASIS_PULSE);
    }
    else {
      for (j = 0; j < n + (2 << left) + (2 << right); j++) {
        printf("%d ", y[j - (2 << left)]);
      }
      printf("\n");
    }
  }
  if (magnitude) printf("\n");
  return 0;
}
