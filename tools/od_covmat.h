#if !defined(_od_covmat_H)
# define _od_covmat_H (0)

#include <stdio.h>

typedef struct od_covmat od_covmat;

struct od_covmat{
  int     sz;
  double  w;
  double *mean;
  double *cov;
  double *work;
};

void od_covmat_init(od_covmat *_this,int _sz);
void od_covmat_clear(od_covmat *_this);
void od_covmat_reset(od_covmat *_this);
void od_covmat_add(od_covmat *_this,const double *_data,double _w);
void od_covmat_combine(od_covmat *_a,const od_covmat *_b);
void od_covmat_correct(od_covmat *_this);
void od_covmat_normalize(od_covmat *_this);
void od_covmat_collapse(od_covmat *_this,int _n,double *_r);
void od_covmat_expand(od_covmat *_this,int _n,const double *_r);
void od_covmat_print(od_covmat *_this,FILE *_fp);

#endif
