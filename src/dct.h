/*Daala video codec
Copyright (c) 2002-2010 Daala project contributors.  All rights reserved.

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

#if !defined(_dct_H)
# define _dct_H (1)
# include <math.h>
# include "filter.h"

/*These are approximations of a true type-II DCT that can be computed entirely
   with integer operations.
  The integerization is done with lifting steps to allow them to be perfectly
   reversible using only small precision outputs.
  If no quantization is performed, the inputs can be recovered exactly.
  This is in contrast to most linear, perfect reconstruction transforms, such
   as a true DCT, which might easily have two implementations which return
   different answers, due to the limited precision of their calculations.
  The drawback, of course, is that these transforms are not truly linear, since
   they use rounding and truncation.

  The 4-point and 16-point transforms are scaled to match the orthonormal
   version of the DCT (where the L^2 norm of each basis function is 1).
  The 8-point transform requires a scale factor of sqrt(0.5) to be applied
   to match the orthonormal version.

  The inverse version of this transform does not come close to passing the
   IEEE1180-1990 conformance test for an iDCT, with peak errors as large as 11
   in the 8x8 case.
  However, that test uses evenly distributed random inputs.
  For natural images, samples are much more correlated, and a common model is
   an AR(0) process with a correlation coefficient of 0.95.
  For this model, the mean squared error is very small, less then 2E-5 in the
   8x8 case.

  The 4-point transform adds one bit of dynamic range to the outputs (e.g.,
   inputs in the range [-256,255) will produce outputs in the range [-512,511).
  The 8-point and 16-point transforms both add 2 bits of dynamic range (e.g.,
   inputs in the range [-256,255) will produce outputs in the range
   [-1024,1023).*/

typedef void (*od_dct_func_2d)(od_coeff *out, int out_stride,
 const od_coeff *in, int in_stride);

typedef void (*od_fdct_func_1d)(od_coeff *out, const od_coeff *in,
 int in_stride);

typedef void (*od_idct_func_1d)(od_coeff *out, int out_stride,
 const od_coeff *in);

extern const od_dct_func_2d OD_FDCT_2D_C[OD_NBSIZES + 1];
extern const od_dct_func_2d OD_IDCT_2D_C[OD_NBSIZES + 1];

extern const od_fdct_func_1d OD_FDCT_1D[OD_NBSIZES + 1];

extern const od_idct_func_1d OD_IDCT_1D[OD_NBSIZES + 1];

/*A reversible integer approximation of a 4-point type-II DCT.
  y:       The destination vector (of size 4).
           This may be the same as the source.
  x:       The source vector (of size 4).
  xstride: The stride of the source.*/
void od_bin_fdct4(od_coeff y[4], const od_coeff *x, int xstride);

/*An inverse of the reversible integer approximation of a 4-point type-II DCT
   above.
  x:       The destination vector (of size 4).
           This may be the same as the source.
  xstride: The stride of the destination.
  y:       The source vector (of size 4).*/
void od_bin_idct4(od_coeff *x, int xstride, const od_coeff y[4]);

void od_bin_fdct4x4(od_coeff *y, int ystride, const od_coeff *x, int xstride);
void od_bin_idct4x4(od_coeff *x, int xstride, const od_coeff *y, int ystride);

/*A reversible integer approximation of an 8-point type-II DCT based on
   Loeffler's factorization.
  The final output should be scaled by M_SQRT1_2 in order to approximate a true
   orthonormal DCT.
  y:       The destination vector (of size 8).
           This may be the same as the source.
  x:       The source vector (of size 8).
  xstride: The stride of the source.*/
void od_bin_fdct8(od_coeff y[8], const od_coeff *x, int xstride);

/*An inverse of the reversible integer approximation of an 8-point type-II DCT
   above.
  The final output should be scaled by M_SQRT2 in order to approximate a true
   orthonormal DCT.
  x:       The destination vector (of size 8).
           This may be the same as the source.
  xstride: The stride of the destination.
  y:       The source vector (of size 8).*/
void od_bin_idct8(od_coeff x[8], int xstride, const od_coeff *y);

void od_bin_fdct8x8(od_coeff *y, int ystride, const od_coeff *x, int xstride);
void od_bin_idct8x8(od_coeff *x, int xstride, const od_coeff *y, int ystride);

/*A reversible integer approximation of a 16-point type-II DCT based on
   Loeffler's factorization.
  y:       The destination vector (of size 16).
           This may be the same as the source.
  x:       The source vector (of size 16).
  xstride: The stride of the source.*/
void od_bin_fdct16(od_coeff y[16], const od_coeff *x, int xstride);

/*An inverse of the reversible integer approximation of a 16-point type-II DCT
   above.
  x:       The destination vector (of size 16).
           This may be the same as the source.
  xstride: The stride of the destination.
  y:       The source vector (of size 16).*/
void od_bin_idct16(od_coeff *x, int xstride, const od_coeff y[16]);

void od_bin_fdct16x16(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct16x16(od_coeff *x, int xstride,
 const od_coeff *y, int ystride);

/*A reversible integer approximation of a 32-point type-II DCT based on
   Loeffler's factorization.
  y:       The destination vector (of size 32).
           This may be the same as the source.
  x:       The source vector (of size 32).
  xstride: The stride of the source.*/
void od_bin_fdct32(od_coeff y[32], const od_coeff *x, int xstride);

/*An inverse of the reversible integer approximation of a 32-point type-II DCT
   above.
  x:       The destination vector (of size 32).
           This may be the same as the source.
  xstride: The stride of the destination.
  y:       The source vector (of size 32).*/
void od_bin_idct32(od_coeff *x, int xstride, const od_coeff y[32]);

void od_bin_fdct32x32(od_coeff *y, int ystride,
 const od_coeff *x, int xstride);
void od_bin_idct32x32(od_coeff *x, int xstride,
 const od_coeff *y, int ystride);

void od_haar(od_coeff *y, int ystride, const od_coeff *x, int xstride, int ln);
void od_haar_inv(od_coeff *x, int xstride, const od_coeff *y, int ystride,
 int ln);

# if defined(OD_CHECKASM)
void od_dct_check(int ln, const od_coeff *ref, const od_coeff *x,
 int xstride);
# endif

#endif
