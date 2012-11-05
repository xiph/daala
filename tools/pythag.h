/*Daala video codec
Copyright (c) 2005-2008 Daala project contributors.  All rights reserved.

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

#if !defined(_cpack_linalg_pythag_H)
# define _cpack_linalg_pythag_H (1)
# include <math.h>

/*Computes the sqrt(a*a+b*b) without destructive underflow.*/
# if defined(__GNUC__)&& \
 (defined(__USE_MISC)||defined(__USE_XOPEN)||defined(__USE_ISOC99))
#  define pythag hypot
# else
static double pythag(double _a,double _b){
  double absa;
  double absb;
  double s;
  int    sb;
  absa=fabs(_a);
  absb=fabs(_b);
  /*frexp() should be good enough (instead of using the recriprocal), and is
     finally properly optimized in gcc 4.3, so it might actually be faster.
    It also doesn't introduce rounding error due to the division and reciprocal
     multiplication.*/
  s=frexp(absa>absb?absa:absb,&sb);
  /*Technically <1 is sufficient, but <= is perhaps safer against
     non-conforming implementations of frexp().*/
  if(s>0&&s<=1){
    double a2;
    double b2;
    a2=ldexp(absa,-sb);
    b2=ldexp(absb,-sb);
    return ldexp(sqrt(a2*a2+b2*b2),sb);
  }
  /*frexp() returns 0 when its argument is 0, and +/- Inf when its argument is
     +/- Inf.*/
  else return s;
}
# endif

#endif
