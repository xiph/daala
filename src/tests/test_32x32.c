#include <stdlib.h>
#include <stdio.h>
#include "../dct.c"
#include "../tf.c"
#include "../internal.c"

int main() {
  od_coeff v[32*32];
  od_coeff t[32*32];
  od_coeff o[32*32];
  int i;
  int j;
  printf("Input block:\n");
  for (i = 0; i < 32; i++) {
    printf("    ");
    for (j = 0; j < 32; j++) {
      v[i*32 + j] = random();
      if (random() & 1) v[i*32 + j] = -v[i*32 + j];
      printf("%d ", v[i*32 + j]);
    }
    printf("\n");
  }
  printf("\n");
  od_bin_fxform32x32(t, 32, v, 32);
  printf("Transformed block:\n");
  for (i = 0; i < 32; i++) {
    printf("    ");
    for (j = 0; j < 32; j++) {
      printf("%d ", t[i*32 + j]);
    }
    printf("\n");
  }
  printf("\n");
  od_bin_ixform32x32(o, 32, t, 32);
  printf("Reconstructed block:\n");
  for (i = 0; i < 32; i++) {
    printf("    ");
    for (j = 0; j < 32; j++) {
      printf("%d ", o[i*32 + j]);
    }
    printf("\n");
  }
  printf("\n");
  printf("Residual:\n");
  for (i = 0; i < 32; i++) {
    printf("    ");
    for (j = 0; j < 32; j++) {
      printf("%d ", o[i*32 + j] - v[i*32 + j]);
    }
    printf("\n");
  }
  return 0;
}
