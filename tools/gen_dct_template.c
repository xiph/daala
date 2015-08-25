#include <stdio.h>
#include <stdlib.h>

const char SUF[]=
 "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ";

int br(int v, int ln) {
  int i;
  int r;
  r = 0;
  for (i = 0; i < ln; i++) {
    r <<= 1;
    r |= v&1;
    v >>= 1;
  }
  return r;
}

void printDCT(FILE *fp, int ln) {
  int n;
  int i;
  char t;
  n = 1 << ln;
  t = 'o' + ln;
  fprintf(fp, "#define OD_FDCT_%i(", n);
  for (i = 0; i < n; i++) {
    fprintf(fp, "%s%c%c", i > 0 ? ", " : "", t, SUF[br(i, ln)]);
  }
  fprintf(fp, ") \\\n");
  fprintf(fp, "  /* Embedded %i-point orthonormal Type-II fDCT. */ \\\n", n);
  fprintf(fp, "  do { \\\n");
  for (i = n/2; i < n; i++) {
    fprintf(fp, "    int %c%ch; \\\n", t, SUF[i]);
  }
  for (i = 0; i < n/2; i += 2) {
    char a;
    char b;
    char c;
    char d;
    a = SUF[i];
    b = SUF[i + 1];
    c = SUF[n - 1 - i - 1];
    d = SUF[n - 1 - i];
    fprintf(fp, "    %c%c = %c%c - %c%c; \\\n", t, d, t, a, t, d);
    fprintf(fp, "    %c%ch = OD_DCT_RSHIFT(%c%c, 1); \\\n", t, d, t, d);
    fprintf(fp, "    %c%c -= %c%ch; \\\n", t, a, t, d);
    fprintf(fp, "    %c%c += %c%c; \\\n", t, c, t, b);
    fprintf(fp, "    %c%ch = OD_DCT_RSHIFT(%c%c, 1); \\\n", t, c, t, c);
    fprintf(fp, "    %c%c = %c%ch - %c%c; \\\n", t, b, t, c, t, b);
  }
  fprintf(fp, "    OD_FDCT_%i_ASYM(", 1 << ln - 1);
  for (i = 0; i < n/2; i += 2) {
    char a;
    char b;
    a = SUF[2*br(i, ln - 1)];
    b = SUF[2*br(i + 1, ln - 1)];
    fprintf(fp, "%s%c%c, %c%c, %c%ch", i > 0 ? ", " : "", t, a, t, b, t, b);
  }
  fprintf(fp, "); \\\n");
  fprintf(fp, "    OD_FDST_%i_ASYM(", 1 << ln - 1);
  for (i = n/2; i-- > 0; ) {
    char a;
    a = SUF[2*br(i, ln - 1) + 1];
    fprintf(fp, "%s%c%c", i < n/2 - 1 ? ", " : "", t, a);
  }
  fprintf(fp, "); \\\n");
  fprintf(fp, "  } while (0)\n\n");
  fprintf(fp, "#define OD_IDCT_%i(", n);
  for (i = 0; i < n; i++) {
    fprintf(fp, "%s%c%c", i > 0 ? ", " : "", t, SUF[br(i, ln)]);
  }
  fprintf(fp, ") \\\n");
  fprintf(fp, "  /* Embedded %i-point orthonormal Type-II fDCT. */ \\\n", n);
  fprintf(fp, "  do { \\\n");
  for (i = 0; i < n/2; i++) {
    fprintf(fp, "    int %c%ch; \\\n", t, SUF[2*i + 1]);
  }
  fprintf(fp, "    OD_IDST_%i_ASYM(", n/2);
  for (i = n/2; i-- > 0; ) {
    char a;
    a = SUF[br(i, ln - 1) + n/2];
    fprintf(fp, "%s%c%c", i < n/2 - 1 ? ", " : "", t, a);
  }
  fprintf(fp, "); \\\n");
  fprintf(fp, "    OD_IDCT_%i_ASYM(", n/2);
  for (i = 0; i < n/4; i++) {
    fprintf(fp, "%s%c%c", i > 0 ? ", " : "", t, SUF[br(i, ln - 1)]);
  }
  for (i = n/4; i < n/2; i++) {
    char a;
    a = SUF[br(i, ln - 1)];
    fprintf(fp, ", %c%c, %c%ch", t, a, t, a);
  }
  fprintf(fp, "); \\\n");
  for (i = 0; i < n/2; i += 2) {
    char a;
    char b;
    char c;
    char d;
    a = SUF[i];
    b = SUF[i + 1];
    c = SUF[n - 1 - i - 1];
    d = SUF[n - 1 - i];
    fprintf(fp, "    %c%ch = OD_DCT_RSHIFT(%c%c, 1); \\\n", t, d, t, d);
    fprintf(fp, "    %c%c += %c%ch; \\\n", t, a, t, d);
    fprintf(fp, "    %c%c = %c%c - %c%c; \\\n", t, d, t, a, t, d);
    fprintf(fp, "    %c%c = %c%ch - %c%c; \\\n", t, c, t, b, t, c);
    fprintf(fp, "    %c%c -= %c%c; \\\n", t, b, t, c);
  }
  fprintf(fp, "  } while (0)\n\n");
}

void printAsymDCT(FILE *fp, int ln) {
  int n;
  int i;
  char t;
  n = 1 << ln;
  t = 'o' + ln;
  fprintf(fp, "#define OD_FDCT_%i_ASYM(", n);
  for (i = 0; i < n; i += 2) {
    char a;
    char b;
    a = SUF[br(i, ln)];
    b = SUF[br(i + 1, ln)];
    fprintf(fp, "%s%c%c, %c%c, %c%ch", i > 0 ? ", " : "", t, a, t, b, t, b);
  }
  fprintf(fp, ") \\\n");
  fprintf(fp, "  /* Embedded %i-point asymmetric Type-II fDCT. */ \\\n", n);
  fprintf(fp, "  do { \\\n");
  for (i = 0; i < n/2; i += 2) {
    char a;
    char b;
    char c;
    char d;
    a = SUF[i];
    b = SUF[i + 1];
    c = SUF[n - 1 - i - 1];
    d = SUF[n - 1 - i];
    fprintf(fp, "    %c%c += %c%ch; \\\n", t, a, t, d);
    fprintf(fp, "    %c%c = %c%c - %c%c; \\\n", t, d, t, a, t, d);
    if (i + 1 < n/2) {
      fprintf(fp, "    %c%c = %c%ch - %c%c; \\\n", t, b, t, c, t, b);
      fprintf(fp, "    %c%c -= %c%c; \\\n", t, c, t, b);
    }
  }
  fprintf(fp, "    OD_FDCT_%i(", 1 << ln - 1);
  for (i = 0; i < n/2; i++) {
    fprintf(fp, "%s%c%c", i > 0 ? ", " : "", t, SUF[2*br(i, ln - 1)]);
  }
  fprintf(fp, "); \\\n");
  fprintf(fp, "    OD_FDST_%i(", 1 << ln - 1);
  for (i = n/2; i-- > 0; ) {
    fprintf(fp, "%s%c%c", i < n/2 - 1 ? ", " : "", t, SUF[2*br(i, ln - 1) + 1]);
  }
  fprintf(fp, "); \\\n");
  fprintf(fp, "  } while (0)\n\n");
  fprintf(fp, "#define OD_IDCT_%i_ASYM(", n);
  for (i = 0; i < n/2; i++) {
    fprintf(fp, "%s%c%c", i > 0 ? ", " : "", t, SUF[2*br(i, ln - 1)]);
  }
  for (i = 0; i < n/2; i++) {
    char a;
    a = SUF[2*br(i, ln - 1) + 1];
    fprintf(fp, ", %c%c, %c%ch", t, a, t, a);
  }
  fprintf(fp, ") \\\n");
  fprintf(fp, "  /* Embedded %i-point asymmetric Type-II iDCT. */ \\\n", n);
  fprintf(fp, "  do { \\\n");
  fprintf(fp, "    OD_IDST_%i(", n/2);
  for (i = n/2; i-- > 0; ) {
    char a;
    a = SUF[br(i, ln - 1) + n/2];
    fprintf(fp, "%s%c%c", i < n/2 - 1 ? ", " : "", t, a);
  }
  fprintf(fp, "); \\\n");
  fprintf(fp, "    OD_IDCT_%i(", n/2);
  for (i = 0; i < n/2; i++) {
    fprintf(fp, "%s%c%c", i > 0 ? ", " : "", t, SUF[br(i, ln - 1)]);
  }
  fprintf(fp, "); \\\n");
  for (i = 0; i < n/2; i += 2) {
    char a;
    char b;
    char c;
    char d;
    a = SUF[i];
    b = SUF[i + 1];
    c = SUF[n - 1 - i - 1];
    d = SUF[n - 1 - i];
    fprintf(fp, "    %c%c = %c%c - %c%c; \\\n", t, d, t, a, t, d);
    fprintf(fp, "    %c%ch = OD_DCT_RSHIFT(%c%c, 1); \\\n", t, d, t, d);
    fprintf(fp, "    %c%c -= %c%ch; \\\n", t, a, t, d);
    if (i + 1 < n/2) {
      fprintf(fp, "    %c%c += %c%c; \\\n", t, b, t, c);
      fprintf(fp, "    %c%ch = OD_DCT_RSHIFT(%c%c, 1); \\\n", t, b, t, b);
      fprintf(fp, "    %c%c = %c%ch - %c%c; \\\n", t, c, t, b, t, c);
    }
  }
  fprintf(fp, "  } while (0)\n");
}

int main(int _argc,char *_argv[]) {
  printAsymDCT(stdout, 4);
}
