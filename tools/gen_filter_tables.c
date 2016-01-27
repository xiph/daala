/*Daala video codec
Copyright (c) 2015 Daala project contributors.  All rights reserved.

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void od_filter_dering_direction_offsets(int *offset, int dir, int bstride) {
  int k;
  int f;
  if (dir <= 4) {
    f = dir - 2;
    for (k = 1; k <= 3; k++) offset[k - 1] = f*k/2*bstride + k;
  }
  else {
    f = 6 - dir;
    for (k = 1; k <= 3; k++) offset[k - 1] = k*bstride + f*k/2;
  }
}

int main(int argc, char **argv) {
  printf("int OD_DIRECTION_OFFSETS_TABLE[8][3] = {\n");
  int offset[3];
  for (int dir = 0; dir < 8; dir++) {
    int k;
    int f;
    printf("  {");
    if (dir <= 4) {
      f = dir - 2;
      for (int k = 1; k <= 3; k++)  {
        printf("%2i*OD_FILT_BSTRIDE + %2i", f*k/2, k);
        if (k < 3) {
          printf(",");
        }
      }
    } else {
      f = 6 - dir;
      for (int k = 1; k <= 3; k++)  {
        printf("%2i*OD_FILT_BSTRIDE + %2i", k, f*k/2);
        if (k < 3) {
          printf(",");
        }
      }
    }
    printf(" }");
    if (dir < 15) {
      printf(",\n");
    }
  }
  printf("\n};\n");
}
