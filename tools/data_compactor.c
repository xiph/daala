/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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

/*This generates intradata.c or od_intra_data.c based on whether TRAIN is set
   or not. It uses the original non-sparse versions of these to generate
   smaller ones.*/

#include <stdio.h>

#if defined(TRAIN)
#define PREFIX NE
#define CONST ""
#include "tools/od_intra_data.c"
#else
#define PREFIX OD
#define CONST "const "
#include "src/intradata.c"
#endif

#define M2S(x) #x
#define MM2S(x) M2S(x)
#define EXPAND(prefix,stem,n) EXPAND_(prefix,stem,n)
#define EXPAND_(prefix,stem,n) \
  prefix ## _ ## stem ## _ ## n ## x ## n

#define PRINT_FORMAT_double "%.17g"
#define PRINT_FORMAT_int "%d"

#define PRINT_ARRAY_MODE_NXN(typ,n,data)                                \
  {                                                                     \
    int m;                                                              \
    int i;                                                              \
    int j;                                                              \
    printf(CONST #typ " " MM2S(PREFIX) "_" #data "_" #n "x" #n "[OD_INTRA_NMODES][" #n "][" #n "]={\n"); \
    for (m = 0; m < OD_INTRA_NMODES; m++) {                             \
      printf("/* Mode %d */\n", m);                                     \
      printf("  {\n");                                                  \
      for (i = 0; i < n; i++) {                                         \
        printf("    { ");                                               \
        for (j = 0; j < n; j++) {                                       \
          printf(PRINT_FORMAT_##typ, EXPAND(PREFIX,data,n)[m][i][j]);   \
          if (j < n-1) printf(", ");                                    \
        }                                                               \
        if (i < n-1) {                                                  \
          printf(" },\n");                                              \
        } else {                                                        \
          printf(" }\n");                                               \
        }                                                               \
      }                                                                 \
      if (m < OD_INTRA_NMODES-1) {                                      \
        printf("  },\n");                                               \
      } else {                                                          \
        printf("  }\n");                                                \
      }                                                                 \
    }                                                                   \
    printf("};\n");                                                     \
  }

#define PRINT_COMPACT_MODE_NXN(typ,n,mults,data)                        \
  {                                                                     \
    int sum;                                                            \
    int m;                                                              \
    int i;                                                              \
    int j;                                                              \
    int k;                                                              \
    int s;                                                              \
    sum = 0;                                                            \
    for (m = 0; m < OD_INTRA_NMODES; m++) {                             \
      for (i = 0; i < n; i++) {                                         \
        for (j = 0; j < n; j++) {                                       \
          sum += EXPAND(PREFIX,mults,n)[m][i][j];                       \
        }                                                               \
      }                                                                 \
    }                                                                   \
    printf(CONST #typ " " MM2S(PREFIX) "_" #data "_" #n "x" #n "_DATA[%d]={\n", sum); \
    s = 0;                                                              \
    for (m = 0; m < OD_INTRA_NMODES; m++) {                             \
      for (i = 0; i < n; i++) {                                         \
        for (j = 0; j < n; j++) {                                       \
          int nvals;                                                    \
          nvals = EXPAND(PREFIX,mults,n)[m][i][j];                      \
          for (k = 0; k < nvals; k++) {                                 \
            printf(PRINT_FORMAT_##typ, EXPAND(PREFIX,data,n)[m][i][j][k]); \
            if (s < sum-1) printf(",");                                 \
            printf("\n");                                               \
            s++;                                                        \
          }                                                             \
        }                                                               \
      }                                                                 \
    }                                                                   \
    printf("};\n");                                                     \
    printf(CONST #typ " *" MM2S(PREFIX) "_" #data "_" #n "x" #n "[OD_INTRA_NMODES][" #n "][" #n "]={\n"); \
    sum = 0;                                                            \
    for (m = 0; m < OD_INTRA_NMODES; m++) {                             \
      printf("/* Mode %d */\n", m);                                     \
      printf("  {\n");                                                  \
      for (i = 0; i < n; i++) {                                         \
        printf("    { ");                                               \
        for (j = 0; j < n; j++) {                                       \
          int nvals;                                                    \
          nvals = EXPAND(PREFIX,mults,n)[m][i][j]; \
            if (nvals > 0) {                                            \
              printf(MM2S(PREFIX) "_" #data "_" #n "x" #n "_DATA+%d", sum); \
            } else {                                                    \
              printf("NULL");                                           \
            }                                                           \
            if (j < n-1) printf(", ");                                  \
            sum += nvals;                                               \
        }                                                               \
        printf("    }");                                                \
        if (i < n-1) printf(",");                                       \
        printf("\n");                                                   \
      }                                                                 \
      printf("  }");                                                    \
      if (m < OD_INTRA_NMODES-1) printf(",");                           \
      printf("\n");                                                     \
    }                                                                   \
    printf("};\n");                                                     \
  }

#if defined(TRAIN)
#define PRINT_DATA_NXN(n) \
  PRINT_ARRAY_MODE_NXN(double, n, PRED_OFFSETS);    \
  PRINT_ARRAY_MODE_NXN(int, n, PRED_MULTS);                \
  PRINT_COMPACT_MODE_NXN(double, n, PRED_MULTS, PRED_WEIGHTS); \
  PRINT_COMPACT_MODE_NXN(int, n, PRED_MULTS, PRED_PARAMX); \
  PRINT_COMPACT_MODE_NXN(int, n, PRED_MULTS, PRED_PARAMY);
#else
#define PRINT_DATA_NXN(n) \
  PRINT_ARRAY_MODE_NXN(int, n, PRED_MULTS);                \
  PRINT_COMPACT_MODE_NXN(double, n, PRED_MULTS, PRED_WEIGHTS); \
  PRINT_COMPACT_MODE_NXN(int, n, PRED_MULTS, PRED_PARAMX); \
  PRINT_COMPACT_MODE_NXN(int, n, PRED_MULTS, PRED_PARAMY);
#endif

int main() {
#if defined(TRAIN)
  printf("#include \"od_defs.h\"\n");
  printf("#include \"od_intra.h\"\n\n");
#else
  printf("#include \"intra.h\"\n\n");
#endif

  PRINT_DATA_NXN(4);
  PRINT_DATA_NXN(8);
  PRINT_DATA_NXN(16);

#if !defined(TRAIN)
  {
    int i;
    int m;
    int c;
    printf("const unsigned char "
     "OD_INTRA_PRED_PROB_4x4[3][OD_INTRA_NMODES][OD_INTRA_NCONTEXTS]={\n");
    for (i = 0; i < 3; i++) {
      printf("{");
      for (m = 0; m < OD_INTRA_NMODES; m++) {
        printf("{");
        for (c = 0; c < OD_INTRA_NCONTEXTS; c++) {
          printf("%u", OD_INTRA_PRED_PROB_4x4[i][m][c]);
          if (c < OD_INTRA_NCONTEXTS-1) printf(", ");
        }
        printf("}");
        if (m < OD_INTRA_NMODES-1) printf(",");
        printf("\n");
      }
      printf("}");
      if (i < 2) printf(",");
      printf("\n");
    }
    printf("};\n");
  }
#endif

  return 0;
}
