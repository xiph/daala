#include <stdio.h>
#include <stdlib.h>
#include "../src/dct.h"
#include "../src/filter.h"
#include "../src/tf.h"

#define OD_BASIS_SIZE (4*OD_BSIZE_MAX)
#define OD_BASIS_PULSE 1024

int main(int argc, char *argv[]) {
  int ln;
  int n;
  od_coeff x0[OD_BASIS_SIZE*OD_BASIS_SIZE];
  od_coeff *x1;
  od_coeff *x2;
  int i;
  int j;
  int k;
  int l;
  double sum[OD_NBSIZES - 1][4];
  (void)argc;
  (void)argv;
  for (ln = 0; ln < OD_NBSIZES - 1; ln++) {
    n = 4 << ln;
    for (j = 0; j < 2; j++) {
      for (i = 0; i < 2; i++) {
        /*printf("%i %i\n",i,j);*/
        OD_CLEAR(x0, OD_BASIS_SIZE*OD_BASIS_SIZE);
        x0[n*(j + 1)*OD_BASIS_SIZE + n*(i + 1)] = OD_BASIS_PULSE;
        OD_HAAR_KERNEL(x0[n*1*OD_BASIS_SIZE + n*1], x0[n*1*OD_BASIS_SIZE + n*2],
                       x0[n*2*OD_BASIS_SIZE + n*1], x0[n*2*OD_BASIS_SIZE + n*2]);
        for (l = 0; l < 2; l++) {
          for (k = 0; k < 2; k++) {
            x1 = &x0[n*(l + 1)*OD_BASIS_SIZE + n*(k + 1)];
            OD_IDCT_2D_C[ln](x1, OD_BASIS_SIZE, x1, OD_BASIS_SIZE);
          }
        }
        od_postfilter_split(&x0[n*OD_BASIS_SIZE + n], OD_BASIS_SIZE, ln + 1, 0, 0, 0, 0, 1, 1);
        x1 = &x0[n - 2];
        x2 = &x0[3*n - 2];
        for (l = 0; l < 4*n; l++) {
          (*OD_POST_FILTER[0])(x1 + l*OD_BASIS_SIZE, x1 + l*OD_BASIS_SIZE);
          (*OD_POST_FILTER[0])(x2 + l*OD_BASIS_SIZE, x2 + l*OD_BASIS_SIZE);
        }
        x1 = &x0[(n - 2)*OD_BASIS_SIZE];
        x2 = &x0[(3*n - 2)*OD_BASIS_SIZE];
        for (k = 0; k < 4*n; k++) {
          od_coeff t1[4];
          od_coeff t2[4];
          for (l = 0; l < 4; l++) {
            t1[l] = x1[OD_BASIS_SIZE*l + k];
            t2[l] = x2[OD_BASIS_SIZE*l + k];
          }
          (*OD_POST_FILTER[0])(t1, t1);
          (*OD_POST_FILTER[0])(t2, t2);
          for (l = 0; l < 4; l++) {
            x1[OD_BASIS_SIZE*l + k] = t1[l];
            x2[OD_BASIS_SIZE*l + k] = t2[l];
          }
        }
        /*x1 = &x0[(n - 2)*OD_BASIS_SIZE + n - 2];
        for (l = 0; l < 2*n + 4; l++) {
          for (k = 0; k < 2*n + 4; k++) {
            printf("%5i ", x1[l*OD_BASIS_SIZE + k]);
          }
          printf("\n");
        }*/
        x1 = &x0[(n - 2)*OD_BASIS_SIZE + n - 2];
        sum[ln][j*2 + i] = 0;
        for (l = 0; l < 2*n + 4; l++) {
          for (k = 0; k < 2*n + 4; k++) {
            sum[ln][j*2 + i] += x1[OD_BASIS_SIZE*l + k]*x1[OD_BASIS_SIZE*l + k];
          }
        }
      }
    }
  }
  for (ln = 0; ln < OD_NBSIZES - 1; ln++) {
    for (i = 0; i < 4; i++) {
      printf("%s%f", i > 0 ? ", " : "", sqrt(sum[ln][i])/OD_BASIS_PULSE);
    }
    printf("\n");
  }
  printf("\n");
  printf("const od_coeff OD_DC_QM[OD_NBSIZES - 1][2] = {\n");
  for (ln = 0; ln < OD_NBSIZES - 1; ln++) {
    int hv_quant;
    int diag_quant;
    hv_quant = (int)(OD_BASIS_PULSE/sqrt(sum[ln][1])*64 + 0.5);
    diag_quant = (int)(OD_BASIS_PULSE/sqrt(sum[ln][3])*64 + 0.5);
    printf("%s{%i, %i}", ln > 0 ? ", " : "  ", hv_quant, diag_quant);
  }
  printf("\n};\n");
  return EXIT_SUCCESS;
}
