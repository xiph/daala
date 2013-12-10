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

# if !defined(M_PI)
#  define M_PI      (3.1415926535897932384626433832795)
# endif

# if !defined(M_SQRT2)
#  define M_SQRT2 (1.41421356237309504880168872420970)
# endif

# if !defined(M_SQRT1_2)
#  define M_SQRT1_2 (0.70710678118654752440084436210485)
# endif

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

typedef void (*od_dct_func_2d)(od_coeff *_out, int _out_stride,
 const od_coeff *_in, int _in_stride);

typedef void (*od_fdct_func_1d)(od_coeff *_out, const od_coeff *_in,
 int _in_stride);

typedef void (*od_idct_func_1d)(od_coeff *_out, int _out_stride,
 const od_coeff *_in);

extern const od_dct_func_2d OD_FDCT_2D[OD_NBSIZES + 1];
extern const od_dct_func_2d OD_IDCT_2D[OD_NBSIZES + 1];

extern const od_fdct_func_1d OD_FDCT_1D[OD_NBSIZES + 1];

extern const od_idct_func_1d OD_IDCT_1D[OD_NBSIZES + 1];

extern const int OD_TRANS_QUANT_ADJ[3];

extern const unsigned char OD_ZIG4[16];
extern const unsigned char OD_ZIG8[64];
extern const unsigned char OD_ZIG16[256];
extern const unsigned char *OD_DCT_ZIGS[OD_NBSIZES + 1];

/*A reversible integer approximation of a 4-point type-II DCT.
  _y:       The destination vector (of size 4).
            This may be the same as the source.
  _x:       The source vector (of size 4).
  _xstride: The stride of the source.*/
void od_bin_fdct4(od_coeff _y[4], const od_coeff *_x, int _xstride);

/*An inverse of the reversible integer approximation of a 4-point type-II DCT
   above.
  _x:       The destination vector (of size 4).
            This may be the same as the source.
  _xstride: The stride of the destination.
  _y:       The source vector (of size 4).*/
void od_bin_idct4(od_coeff *_x, int _xstride, const od_coeff _y[4]);

void od_bin_fdct4x4(od_coeff *_y, int _ystride, const od_coeff *_x,
 int _xstride);
void od_bin_idct4x4(od_coeff *_x, int _xstride, const od_coeff *_y,
 int _ystride);

/*A reversible integer approximation of an 8-point type-II DCT based on
   Loeffler's factorization.
  The final output should be scaled by M_SQRT1_2 in order to approximate a true
   orthonormal DCT.
  _y:       The destination vector (of size 8).
            This may be the same as the source.
  _x:       The source vector (of size 8).
  _xstride: The stride of the source.*/
void od_bin_fdct8(od_coeff _y[8], const od_coeff *_x, int _xstride);

/*An inverse of the reversible integer approximation of an 8-point type-II DCT
   above.
  The final output should be scaled by M_SQRT2 in order to approximate a true
   orthonormal DCT.
  _x:       The destination vector (of size 8).
            This may be the same as the source.
  _xstride: The stride of the destination.
  _y:       The source vector (of size 8).*/
void od_bin_idct8(od_coeff _x[8], int _xstride, const od_coeff *_y);

void od_bin_fdct8x8(od_coeff *_y, int _ystride, const od_coeff *_x,
 int _xstride);
void od_bin_idct8x8(od_coeff *_x, int _xstride, const od_coeff *_y,
 int _ystride);

/*A reversible integer approximation of a 16-point type-II DCT based on
   Loeffler's factorization.
  _y:       The destination vector (of size 16).
            This may be the same as the source.
  _x:       The source vector (of size 16).
  _xstride: The stride of the source.*/
void od_bin_fdct16(od_coeff _y[16], const od_coeff *_x, int _xstride);

/*An inverse of the reversible integer approximation of a 16-point type-II DCT
   above.
  _x:       The destination vector (of size 16).
            This may be the same as the source.
  _xstride: The stride of the destination.
  _y:       The source vector (of size 16).*/
void od_bin_idct16(od_coeff *_x, int _xstride, const od_coeff _y[16]);

void od_bin_fdct16x16(od_coeff *_y, int _ystride,
 const od_coeff *_x, int _xstride);
void od_bin_idct16x16(od_coeff *_x, int _xstride,
 const od_coeff *_y, int _ystride);

#endif
