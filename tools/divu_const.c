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
#include "../src/odintrin.h"

#define NBITS (32)

/*This program generates the constants for use with OD_DIVU_SMALL().*/
int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("usage: %s <max>\n", argv[0]);
    return EXIT_FAILURE;
  }
  else {
    int max;
    int d;
    max = atoi(argv[1]);
    printf("ogg_uint32_t OD_DIVU_SMALL_CONSTS[OD_DIVU_DMAX][2]={\n ");
    for (d = 1; d <= max; d++) {
      if ((d & (d - 1)) == 0) {
        printf(" {0xFFFFFFFF,0xFFFFFFFF},");
      }
      else {
        unsigned long long t;
        unsigned long long r;
        int m;
        m = OD_LOG2(d);
        t = (1UL << m + NBITS)/d;
        r = (t*d + d) & ((1UL << NBITS) - 1);
        if (r <= 1UL << m) {
          printf(" {0x%llX,         0},", t + 1);
        }
        else {
          printf(" {0x%llX,0x%llX},", t, t);
        }
        if (d % 3 == 0) {
          printf("\n ");
        }
      }
    }
    printf("\n};\n");
    return EXIT_SUCCESS;
  }
}
