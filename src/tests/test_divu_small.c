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
#include "../internal.h"

int main(int argc, char *argv[]) {
  unsigned int seed;
  int          d;
  unsigned int x;
  int          i;
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
  fprintf(stderr, "Random seed: %u (%.4X).\n", seed, rand() & 65535);
  fprintf(stderr, "Testing all divisions up to 1,000,000... ");
  for (d = 1; d <= OD_DIVU_DMAX; d++) {
    for (x = 1; x <= 1000000; x++) {
      if (x/d != OD_DIVU_SMALL(x, d)) {
        fprintf(stderr, "Failure! x=%i d=%i x/d=%i != %i\n",
         x, d, x/d, OD_DIVU_SMALL(x, d));
        return EXIT_FAILURE;
      }
    }
  }
  fprintf(stderr, "Passed!\n");
  fprintf(stderr, "Testing 1,000,000 random divisions... ");
  for (d = 1; d < OD_DIVU_DMAX; d++) {
    for (i = 0; i < 1000000; i++) {
      x = rand();
      if (x/d != OD_DIVU_SMALL(x, d)) {
        fprintf(stderr,"Failure! x=%i d=%i x/d=%i != %i\n",
         x, d, x/d, OD_DIVU_SMALL(x, d));
        return EXIT_FAILURE;
      }
    }
  }
  fprintf(stderr, "Passed!\n");
  return EXIT_SUCCESS;
}
