/*Some common macros for potential platform-specific optimization.*/
#include <math.h>
#include <limits.h>
#if !defined(_odintrin_H)
# define _odintrin_H (1)

/*Some specific platforms may have optimized intrinsic or inline assembly
   versions of these functions which can substantially improve performance.
  We define macros for them to allow easy incorporation of these non-ANSI
   features.*/

/*Note that we do not provide a macro for abs(), because it is provided as a
   library function, which we assume is translated into an intrinsic to avoid
   the function call overhead and then implemented in the smartest way for the
   target platform.
  With modern gcc (4.x), this is true: it uses cmov instructions if the
   architecture supports it and branchless bit-twiddling if it does not (the
   speed difference between the two approaches is not measurable).
  Interestingly, the bit-twiddling method was patented in 2000 (US 6,073,150)
   by Sun Microsystems, despite prior art dating back to at least 1996:
   http://web.archive.org/web/19961201174141/www.x86.org/ftp/articles/pentopt/PENTOPT.TXT
  On gcc 3.x, however, our assumption is not true, as abs() is translated to a
   conditional jump, which is horrible on deeply piplined architectures (e.g.,
   all consumer architectures for the past decade or more).
  Also be warned that -C*abs(x) where C is a constant is mis-optimized as
   abs(C*x) on every gcc release before 4.2.3.
  See bug http://gcc.gnu.org/bugzilla/show_bug.cgi?id=34130 */

/*Modern gcc (4.x) can compile the naive versions of min and max with cmov if
   given an appropriate architecture, but the branchless bit-twiddling versions
   are just as fast, and do not require any special target architecture.
  Earlier gcc versions (3.x) compiled both code to the same assembly
   instructions, because of the way they represented ((_b)>(_a)) internally.*/
/*#define OD_MAXI(_a,_b)      ((_a)<(_b)?(_b):(_a))*/
#define OD_MAXI(_a,_b)      ((_a)-((_a)-(_b)&-((_b)>(_a))))
/*#define OD_MINI(_a,_b)      ((_a)>(_b)?(_b):(_a))*/
#define OD_MINI(_a,_b)      ((_a)+((_b)-(_a)&-((_b)<(_a))))
/*This has a chance of compiling branchless, and is just as fast as the
   bit-twiddling method, which is slightly less portable, since it relies on a
   sign-extended rightshift, which is not guaranteed by ANSI (but present on
   every relevant platform).*/
#define OD_SIGNI(_a)        (((_a)>0)-((_a)<0))
/*Slightly more portable than relying on a sign-extended right-shift (which is
   not guaranteed by ANSI), and just as fast, since gcc (3.x and 4.x both)
   compile it into the right-shift anyway.*/
#define OD_SIGNMASK(_a)     (-((_a)<0))
/*Clamps an integer into the given range.
  If _a>_c, then the lower bound _a is respected over the upper bound _c (this
   behavior is required to meet our documented API behavior).
  _a: The lower bound.
  _b: The value to clamp.
  _c: The upper boud.*/
#define OD_CLAMPI(_a,_b,_c) (OD_MAXI(_a,OD_MINI(_b,_c)))
/*Clamps a signed integer between 0 and 255, returning an unsigned char.
  This assumes a char is 8 bits.*/
#define OD_CLAMP255(_x)     ((unsigned char)((((_x)<0)-1)&((_x)|-((_x)>255))))
/*Divides an integer by a power of two, truncating towards 0.
  _dividend: The integer to divide.
  _shift:    The non-negative power of two to divide by.
  _rmask:    (1<<_shift)-1*/
#define OD_DIV_POW2(_dividend,_shift,_rmask)\
  ((_dividend)+(OD_SIGNMASK(_dividend)&(_rmask))>>(_shift))
/*Divides _x by 65536, truncating towards 0.*/
#define OD_DIV2_16(_x) OD_DIV_POW2(_x,16,0xFFFF)
/*Divides _x by 2, truncating towards 0.*/
#define OD_DIV2(_x) OD_DIV_POW2(_x,1,0x1)
/*Divides _x by 8, truncating towards 0.*/
#define OD_DIV8(_x) OD_DIV_POW2(_x,3,0x7)
/*Divides _x by 16, truncating towards 0.*/
#define OD_DIV16(_x) OD_DIV_POW2(_x,4,0xF)
/*Right shifts _dividend by _shift, adding _rval, and subtracting one for
   negative dividends first.
  When _rval is (1<<_shift-1), this is equivalent to division with rounding
   ties away from zero.*/
#define OD_DIV_ROUND_POW2(_dividend,_shift,_rval)\
  ((_dividend)+OD_SIGNMASK(_dividend)+(_rval)>>(_shift))
/*Divides a _x by 2, rounding towards even numbers.*/
#define OD_DIV2_RE(_x) ((_x)+((_x)>>1&1)>>1)
/*Divides a _x by (1<<(_shift)), rounding towards even numbers.*/
#define OD_DIV_POW2_RE(_x,_shift) \
 ((_x)+((_x)>>(_shift)&1)+((1<<(_shift))-1>>1)>>(_shift))
/*Count leading zeros.
  This macro should only be used for implementing od_ilog(), if it is defined.
  All other code should use OD_ILOG() instead.*/
#if defined(__GNUC__)
# if __GNUC_PREREQ(3,4)
#  if INT_MAX>=2147483647
#   define OD_CLZ0 sizeof(unsigned)*CHAR_BIT
#   define OD_CLZ(_x) (__builtin_clz(_x))
#  elif LONG_MAX>=2147483647L
#   define OD_CLZ0 sizeof(unsigned long)*CHAR_BIT
#   define OD_CLZ(_x) (__builtin_clzl(_x))
#  endif
# endif
#endif
#if defined(OD_CLZ)
# define OD_ILOG(_x)  (OD_CLZ0-OD_CLZ(_x))
/*Note that __builtin_clz is not defined when _x==0, according to the gcc
   documentation (and that of the x86 BSR instruction that implements it), so
   we have to special-case it.
  We define a special version of the macro to use when _x can be zero.*/
# define OD_ILOG0(_x) (OD_ILOG(_x)&-!!(_x))
#else
# define OD_ILOG(_x)  (od_ilog(_x))
# define OD_ILOG0(_x) (od_ilog(_x))
#endif

/*Swaps two integers _a and _b if _a>_b.*/
/*#define OD_SORT2I(_a,_b)\
  if((_a)>(_b)){\
    int t__;\
    t__=(_a);\
    (_a)=(_b);\
    (_b)=t__;\
  }*/
/*This branchless version is significantly faster than the above
   straightforward implementation on modern processors.*/
#define OD_SORT2I(_a,_b)\
  do{\
    int t__;\
    t__=((_a)^(_b))&-((_b)<(_a));\
    (_a)^=t__;\
    (_b)^=t__;\
  }\
  while(0)



/*All of these macros should expect floats as arguments.*/
/*These two should compile as a single SSE instruction.*/
#define OD_MINF(_a,_b)      ((_a)<(_b)?(_a):(_b))
#define OD_MAXF(_a,_b)      ((_a)>(_b)?(_a):(_b))
#define OD_CLAMPF(_a,_b,_c) (OD_MAXF(_a,OD_MINF(_b,_c)))
#if defined(__GNUC__)
# define OD_FABSF(_f)        (fabsf(_f))
# define OD_SQRTF(_f)        (sqrtf(_f))
# define OD_POWF(_b,_e)      (powf(_b,_e))
# define OD_LOGF(_f)         (logf(_f))
# define OD_IFLOORF(_f)      (floorf(_f))
# define OD_ICEILF(_f)       (ceilf(_f))
#else
# define OD_FABSF(_f)        ((float)fabs(_f))
# define OD_SQRTF(_f)        ((float)sqrt(_f))
# define OD_POWF(_b,_e)      ((float)pow(_b,_e))
# define OD_LOGF(_f)         ((float)log(_f))
# define OD_IFLOORF(_f)      ((int)floor(_f))
# define OD_ICEILF(_f)       ((int)ceil(_f))
#endif

#endif
