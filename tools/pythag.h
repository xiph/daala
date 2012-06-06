#if !defined(_cpack_linalg_pythag_H)
# define _cpack_linalg_pythag_H (1)
# include <math.h>

/*Computes the sqrt(a*a+b*b) without destructive underflow.*/
# if defined(__GNUC__)
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
