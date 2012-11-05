/*Daala video codec
Copyright (c) 2002-2010 Daala project contributors.  All rights reserved.

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

#if !defined(_dct_H)
# define _dct_H (1)
# include "filter.h"

#if !defined(M_PI)
# define M_PI      (3.14159265358979323846264383832795)
#endif

#if !defined(M_SQRT2)
# define M_SQRT2 (1.41421356237309504880168872420970)
#endif

#if !defined(M_SQRT1_2)
# define M_SQRT1_2 (0.70710678118654752440084436210485)
#endif

#define SIN_3PI_8 (0.92387953251128675612818318939679)

/*These are approximations of a true type-II DCT that can be computed entirely
   with integer operations.
  The integerization is done with lifting steps to allow perfect
   reconstruction using only small precision outputs.
  If no quantization is performed the inputs can be recovered exactly.
  This is in contrast to most linear, perfect reconstruction transforms, such
   as a true DCT, which might easily have two implementations which return
   different answers, due to the limited precision of their calculations.
  The drawback, of course, is that these transforms are not truly linear, since
   they use rounding and truncation.

  Of course, to match the actual DCT, a scaling factor must be applied to each
   output.
  This can be combined with the quantization step.
  This allows hardware implementations of a true DCT to be used as an
   approximation of this transform for speed.
  This naturally prevents lossless decoding, for which the integer transform
   must be used exactly.

  The inverse version of this transform does not come close to passing the
   IEEE1180-1990 conformance test for an iDCT, with peak errors as large as 10
   in the 8x8 case.
  However, that test uses evenly distributed random inputs.
  For natural images, samples are much more correlated, and a common model is
   an AR(0) process with a correlation coefficient of 0.95.
  For this model, the mean squared error is very small, less then 2E-5 in the
   8x8 case.
  Thus it is conceivable that such hardware acceleration could produce
   acceptable output quality.
  It should not be used for encoding, however.

  For each DCT, for input of one bit-width, the range of each output is within
   a (larger) constant bit-width, except for the DC component, which requires
   one more bit.
  The DC coefficient is always the sum of the input coefficients.*/

/*The (1-D) scaling factors that make a true DCT approximation out of the
   integer transform.*/
static const double DCT4_FSCALE[4]={
  0.5,
  M_SQRT2/SIN_3PI_8,
  1,
  M_SQRT2*SIN_3PI_8
};


/*A perfect-reconstruction integerization of a 4-point type-II DCT.
  The final output should be scaled by
   {
    1/2,
    2/sqrt(1+sqrt(1/2)),
    1,
    sqrt(1+sqrt(1/2))
   }
   in order to approximate a true DCT.
  y: The destination vector (of size 4).
     This may be the same as the source.
  x: The source vector (of size 4).*/
void od_bin_fdct4(od_coeff _y[],const od_coeff _x[]);
/*The (1-D) scaling factors that make a true DCT approximation out of the
   integer transform.*/
static const double DCT8_FSCALE[8]={
  M_SQRT2,
  M_SQRT1_2,
  1/SIN_3PI_8,
  1,
  M_SQRT1_2,
  1,
  SIN_3PI_8,
  M_SQRT1_2
};
/*A perfect-reconstruction integerization of an 8-point type-II DCT based on
   Loeffler's factorization.
  The final output should be scaled by
   {
    sqrt(1/8),
    sqrt(1/2),
    2/sqrt(2+sqrt(2)),
    1,
    sqrt(1/2),
    1,
    sqrt(2+sqrt(2))/2,
    sqrt(1/2)
   }
   in order to approximate a true DCT.
  y: The destination vector (of size 8).
     This may be the same as the source.
  x: The source vector (of size 8).*/
void od_bin_fdct8(od_coeff _y[],const od_coeff _x[]);
/*The (1-D) scaling factors that make a true DCT approximation out of the
   integer transform.*/
static const double DCT16_FSCALE[16]={
  1.0,
  1.0,
  1.0,
  M_SQRT1_2/SIN_3PI_8,
  M_SQRT1_2/SIN_3PI_8,
  M_SQRT1_2*SIN_3PI_8,
  M_SQRT1_2,
  M_SQRT1_2,
  1.0,
  M_SQRT1_2,
  M_SQRT1_2,
  M_SQRT1_2/SIN_3PI_8,
  M_SQRT1_2*SIN_3PI_8,
  M_SQRT1_2*SIN_3PI_8,
  1.0,
  0.5
};
/*A perfect-reconstruction integerization of a 16-point type-II DCT based on
   Loeffler's factorization.
  The final output should be scaled by
   {
    1/4,
    1/2
    1/2
    sqrt(2/(2+sqrt(2))),
    sqrt(2/(2+sqrt(2))),
    sqrt((2+sqrt(2)/8),
    sqrt(1/2),
    sqrt(1/2),
    1/2,
    sqrt(1/2),
    sqrt(1/2),
    sqrt(2/(2+sqrt(2))),
    sqrt((2+sqrt(2)/8),
    sqrt((2+sqrt(2)/8),
    1/2,
    1/2
   }
   in order to approximate a true DCT.
  y: The destination vector (of size 16).
     This may be the same as the source.
  x: The source vector (of size 16).*/
void od_bin_fdct16(od_coeff _y[],const od_coeff _x[]);

/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT4_ISCALE[4]={
  2,
  M_SQRT1_2*SIN_3PI_8,
  1,
  M_SQRT1_2/SIN_3PI_8
};
/*An inverse of the perfect-reconstruction integerization of a 4-point type-II
   DCT above.
  The input should be scaled by
   {
    2,
    sqrt(1+sqrt(1/2))/2,
    1,
    1/sqrt(1+sqrt(1/2))
   }
   in order to approximate a true inverse DCT.
  x: The destination vector (of size 4).
     This may be the same as the source.
  y: The source vector (of size 4).*/
void od_bin_idct4(od_coeff _x[],const od_coeff _y[]);
/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT8_ISCALE[8]={
  M_SQRT1_2,
  M_SQRT2,
  SIN_3PI_8,
  1,
  M_SQRT2,
  1,
  1/SIN_3PI_8,
  M_SQRT2
};
/*An inverse of the perfect-reconstruction integerization of an 8-point type-II
   DCT based on Loeffler's factorization above.
  The input should be scaled by
   {
    sqrt(8),
    sqrt(2),
    sqrt(2+sqrt(2))/2,
    1,
    sqrt(2),
    1,
    2/sqrt(2+sqrt(2)),
    sqrt(2)
   }
   in order to approximate a true DCT.
  x: The destination vector (of size 8).
     This may be the same as the source.
  y: The source vector (of size 8).*/
void od_bin_idct8(od_coeff _x[],const od_coeff _y[]);
/*The (1-D) scaling factors that make a true iDCT approximation out of the
   integer transform.*/
static const double DCT16_ISCALE[16]={
  1,
  1,
  1,
  M_SQRT2*SIN_3PI_8,
  M_SQRT2*SIN_3PI_8,
  M_SQRT2/SIN_3PI_8,
  M_SQRT2,
  M_SQRT2,
  1,
  M_SQRT2,
  M_SQRT2,
  M_SQRT2*SIN_3PI_8,
  M_SQRT2/SIN_3PI_8,
  M_SQRT2/SIN_3PI_8,
  1,
  2
};
/*An inverse of the perfect-reconstruction integerization of a 16-point type-II
   DCT based on Loeffler's factorization above.
  The input should be scaled by
   {
    4,
    2
    2
    sqrt((2+sqrt(2))/2),
    sqrt((2+sqrt(2))/2),
    2sqrt(2/(2+sqrt(2)),
    sqrt(2),
    sqrt(2),
    2,
    sqrt(2),
    sqrt(2),
    sqrt((2+sqrt(2))/2),
    2sqrt(2/(2+sqrt(2)),
    2sqrt(2/(2+sqrt(2)),
    2,
    2
   }
   in order to approximate a true DCT.
  x: The destination vector (of size 16).
     This may be the same as the source.
  y: The source vector (of size 16).*/
void od_bin_idct16(od_coeff _x[],const od_coeff _y[]);

static const unsigned char od_zig4[16]={0,1,5,6,2,4,7,12,3,8,11,13,9,10,14,15};

#endif
