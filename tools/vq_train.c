/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "../src/dct.h"

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MAX_ENTRIES (4096)
#define MAX_DIMS    (128)

double inner_prod(const double *x, const double *y, int n) {
  double sum;
  int i;
  sum = 0;
  for (i = 0; i < n; i++) sum += x[i]*y[i];
  return sum;
}

void normalize(double *x, int n) {
  int i;
  double sum;
  sum = 1e-30;
  for (i = 0; i < n; i++) sum += x[i]*x[i];
  sum = 1./sqrt(sum);
  for (i = 0; i < n; i++) x[i] *= sum;
}

/* Returns the distance to the closest K=2 codeword. We can take a shortcut
 because there are only two possibilities: both pulses at the position with
 largest magnitude, or one pulse at each of the two largest magnitudes. */
double pvq_dist_k2(const double *data, int n) {
  double xbest1, xbest2;
  int i;
  xbest1 = xbest2 = -1;
  for (i = 0; i < n; i++) {
    if (fabs(data[i]) > xbest2) {
      if (fabs(data[i]) > xbest1) {
        xbest2 = xbest1;
        xbest1 = fabs(data[i]);
      }
      else {
        xbest2 = fabs(data[i]);
      }
    }
  }
  return 2 - 2*MAX(xbest1, M_SQRT1_2*(xbest1 + xbest2));
}


int find_nearest(const double *data, const double *codebook, int nb_entries,
 int n, double *sign, double *err) {
  double best_dist;
  double best_sign;
  int best_id;
  int i;
  best_dist = -1;
  best_id = 0;
  best_sign = 1;
  for (i = 0; i < nb_entries; i++) {
    double dist;
    dist = inner_prod(data, &codebook[i*n], n);
    if (fabs(dist) > best_dist) {
      best_dist = fabs(dist);
      best_sign = dist > 0 ? 1 : -1;
      best_id = i;
    }
  }
  if (sign) *sign = best_sign;
  if (err) *err = 2 - 2*best_dist;
  return best_id;
}

double vq_train(const double *data, int nb_vectors, double *codebook,
 int nb_entries, int n, int nb_iter, int exclude_pvq) {
  int i, j;
  int iter;
  double rms;
  double accum[MAX_ENTRIES*MAX_DIMS];
  /* Start with a codebook made of randomly selected vectors. */
  for (i = 0; i < nb_entries; i++) {
    int id;
    id = rand()%nb_vectors;
    for (j = 0; j < n; j++) {
      /* Add some noise just in case we pick the same vector twice. */
      codebook[i*n + j] = data[id*n + j] + .01*(rand()%3 - 1);
    }
    normalize(&codebook[i*n], n);
  }

  rms = 0;
  for (iter = 0; iter < nb_iter; iter++) {
    double pvq_err;
    rms = 0;
    for (i = 0; i < nb_entries*n; i++) accum[i] = 0;
    for (i = 0; i < nb_vectors; i++) {
      int id;
      double sign;
      double err;
      id = find_nearest(&data[i*n], codebook, nb_entries, n, &sign, &err);
      pvq_err = pvq_dist_k2(&data[i*n], n);
      /*printf("%f ", err);*/
      if (!exclude_pvq || err < pvq_err) {
        rms += err;
        for (j = 0; j < n; j++) accum[id*n + j] += sign*data[i*n + j];
      }
      else rms += pvq_err;
    }
    for (i = 0; i < nb_entries; i++) normalize(&accum[i*n], n);
    for (i = 0; i < nb_entries*n; i++) codebook[i] = accum[i];
    rms = sqrt(rms/nb_vectors);
    fprintf(stderr, "RMS: %f\n", rms);
  }
  return rms;
}

int main(int argc, char **argv)
{
  int i,j;
  int nb_vectors, nb_entries, ndim;
  double *data, *codebook;
  double rms;
  unsigned seed;
  seed = time(NULL);
  srand(seed);
  
  if (argc != 4) {
    fprintf(stderr, "usage: vq_train <dimensions> <max vectors> <bits>\n");
    return 1;
  }
  ndim = atoi(argv[1]);
  nb_vectors = atoi(argv[2]);
  nb_entries = 1<<atoi(argv[3]);
  
  data = malloc(nb_vectors*ndim*sizeof(*data));
  codebook = malloc(nb_entries*ndim*sizeof(*codebook));
  if (data == NULL || codebook == NULL) {
    fprintf(stderr, "malloc() failed, giving up.\n");
    return 1;
  }
  for (i = 0;i < nb_vectors; i++) {
    if (feof(stdin))
      break;
    for (j = 0; j < ndim; j++) {
      scanf("%lf ", &data[i*ndim + j]);
    }
    normalize(&data[i*ndim], ndim);
  }
  nb_vectors = i;
  fprintf(stderr, "read %d vectors\n", nb_vectors);

  rms = vq_train(data, nb_vectors, codebook, nb_entries, ndim, 100, 1);
#if 0
  for (i = 0; i < nb_vectors; i++)
  {
    double sign;
    int nearest;
    nearest = find_nearest(&data[i*ndim], codebook, nb_entries, ndim, &sign,
     NULL);
    printf("%d %f\n", nearest, sign);
  }
#endif

  printf("/* Automatically generated by vq_train. */\n");
  printf("/* Seed was %u. */\n", seed);
  printf("/* RMS training error is %f. */\n", rms);
  printf("const double codebook[%d*%d] = {\n", nb_entries, ndim);
  for (i = 0; i < nb_entries; i++) {
    for(j = 0; j < ndim; j++) printf("%f, ", codebook[i*ndim + j]);
    printf("\n");
  }
  printf("};\n");
  free(data);
  free(codebook);
  return 0;
}
