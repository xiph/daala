#if !defined(_tools_tran_tools_H)
# define _tools_tran_tools_H (0)

#include <stdio.h>
#include "intra_fit_tools.h"
#include "../src/internal.h"

#define NE_BITS (6)
#define NE_DISABLE_FILTER (0)

typedef struct image_ctx image_ctx;

struct image_ctx{
  const char *name;
  int         nxblocks;
  int         nyblocks;
};

void image_ctx_init(image_ctx *_this,const char *_name,int _nxblocks,
 int _nyblocks);

typedef struct trans_data trans_data;

struct trans_data{
  int     sz;
  int     n;
  double *mean;
  double *cov;
  double *work;
};

void trans_data_init(trans_data *_this,int _sz);
void trans_data_clear(trans_data *_this);
void trans_data_add(trans_data *_this,const unsigned char *_data);
void trans_data_combine(trans_data *_a,const trans_data *_b);
void trans_data_correct(trans_data *_this);
void trans_data_normalize(trans_data *_this);
void trans_data_collapse(trans_data *_this,int _n,double *_r);
void trans_data_expand(trans_data *_this,int _n,const double *_r);
void trans_data_print(trans_data *_this,FILE *_fp);

void covariance_collapse(const double *_cov,int _sz,int _n,double *_r,
 double *_work);
void covariance_expand(double *_cov,int _sz,int _n,const double *_r);

typedef struct trans_ctx trans_ctx;

struct trans_ctx{
  image_ctx  img;
  trans_data td;
};

typedef struct param_data param_data;

struct param_data{
  int    p[B_SZ/2-1];
  int    q[B_SZ/2-1];
  int    s[B_SZ/2];
  double grgt[B_SZ*B_SZ];
  double hi[2*B_SZ*B_SZ];
  double cg;
  double sba;
  double w;
  double bo;
};

typedef void (*ne_filter_func_double)(double _out[],const double _in[]);

extern const ne_filter_func_double NE_PRE_FILTER_DOUBLE[OD_NBSIZES];
extern const ne_filter_func_double NE_POST_FILTER_DOUBLE[OD_NBSIZES];

void auto_regressive_collapsed(double *_out,int _sz,int _n,double _r);

void analysis(double *_out,int _out_stride,const double *_in,int _in_stride,
 int _n);
void synthesis(double *_out,int _out_stride,const double *_in,int _in_stride,
 int _n);

double coding_gain_1d(const double _r[2*B_SZ*2*B_SZ]);
double coding_gain_1d_collapsed(const double _r[2*B_SZ]);
double coding_gain_2d(const double _r[2*B_SZ*2*B_SZ*2*B_SZ*2*B_SZ]);
double coding_gain_2d_collapsed(const double _r[2*B_SZ*2*B_SZ]);

extern const double *SUBSET1_1D[OD_NBSIZES];
extern const double *SUBSET3_1D[OD_NBSIZES];

extern const double *SUBSET1_2D[OD_NBSIZES];
extern const double *SUBSET3_2D[OD_NBSIZES];

#endif
