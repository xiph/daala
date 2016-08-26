#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "../src/odintrin.h"
#define OD_SQRT_TBL_SHIFT (10)

int main()
{
  int i;
  static int n[13] = {0, 0, 0, 0, 8, 15, 32, 0, 128, 0, 512, 0, 2048};
  printf("static const od_val16 table[2][13] = {\n");
  printf("{");
  for (i = 0; i < 13; i++) {
    if (n[i]) {
      printf("%d, ", OD_MINI(32767, (int32_t)floor(.5
       + (1 << OD_SQRT_TBL_SHIFT)*sqrt((n[i] + 2)/2.))));
    }
    else {
      printf("0, ");
    }
  }
  printf("},\n");
  printf("{");
  for (i = 0; i < 13; i++) {
    if (n[i]) {
      printf("%d, ", OD_MINI(32767, (int32_t)floor(.5
       + (1 << OD_SQRT_TBL_SHIFT)*sqrt((n[i] + 3)/2.))));
    }
    else {
      printf("0, ");
    }
  }
  printf("}};\n");
  return 0;
}
