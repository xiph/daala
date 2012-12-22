#if !defined(_tools_stats_tools_H)
# define _tools_stats_tools_H (1)

#include "intra_fit_tools.h"
#include "../src/intra.h"

#define INPUT_SCALE (16)

typedef struct mode_data mode_data;

struct mode_data{
  int    n;
  double satd_avg[B_SZ*B_SZ];
  double mean;
  double var;
  double ref_mean[B_SZ*B_SZ];
  double ref_cov[B_SZ*B_SZ][B_SZ*B_SZ];
  double pred_mean[B_SZ*B_SZ];
  double pred_cov[B_SZ*B_SZ][B_SZ*B_SZ];
};

void mode_data_init(mode_data *_md);
void mode_data_add_input(mode_data *_md,const unsigned char *_data,int _stride);
void mode_data_add_block(mode_data *_md,const od_coeff *_block,int _stride,
 int _ref);
void mode_data_correct(mode_data *_md);
void mode_data_print(mode_data *_md,const char *_label,double *_scale);
void mode_data_combine(mode_data *_a,const mode_data *_b);

typedef struct intra_stats intra_stats;

struct intra_stats{
  mode_data fr;
  mode_data md[OD_INTRA_NMODES];
};

void intra_stats_init(intra_stats *_this);
void intra_stats_correct(intra_stats *_this);
void intra_stats_print(intra_stats *_this,const char *_label,double *_scale);
void intra_stats_combine(intra_stats *_this,const intra_stats *_that);

/* Is there a good way to have a static initializer in C? */
extern double VP8_SCALE[B_SZ];
extern double OD_SCALE[B_SZ];

void vp8_scale_init(double _vp8_scale[B_SZ]);
void od_scale_init(double _od_scale[B_SZ]);

typedef struct image_data image_data;

struct image_data{
  const char      *name;
  int              nxblocks;
  int              nyblocks;
  int             *mode;
  od_coeff        *pre;
  int              pre_stride;
  od_coeff        *fdct;
  int              fdct_stride;
  od_coeff        *pred;
  int              pred_stride;
  od_coeff        *idct;
  int              idct_stride;
  od_coeff        *post;
  int              post_stride;
};

void image_data_init(image_data *_this,const char *_name,int _nxblocks,
 int _nyblocks);
void image_data_clear(image_data *_this);

void image_data_pre_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj);
void image_data_fdct_block(image_data *_this,int _bi,int _bj);
void image_data_mode_block(image_data *_this,int _bi,int _bj);
void image_data_pred_block(image_data *_this,int _bi,int _bj);
void image_data_idct_block(image_data *_this,int _bi,int _bj);
void image_data_post_block(image_data *_this,int _bi,int _bj);

#endif
