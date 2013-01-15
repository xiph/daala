#if !defined(_tools_stats_tools_H)
# define _tools_stats_tools_H (1)

#include "intra_fit_tools.h"
#include "image.h"
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

int vp8_select_mode(const unsigned char *_data,int _stride,double *_weight);
int od_select_mode(const od_coeff *_block,int _stride,double *_weight); 

extern od_rgba16_pixel COLORS[OD_INTRA_NMODES];

void image_draw_block(od_rgba16_image *_image,int _x,int _y,
 const unsigned char *_block,int _stride);
int image_write_png(od_rgba16_image *_image,const char *_name);

typedef struct image_files image_files;

struct image_files{
  od_rgba16_image raw;
  od_rgba16_image map;
  od_rgba16_image pred;
  od_rgba16_image res;
};

void image_files_init(image_files *_this,int _nxblocks,int _nyblocks);
void image_files_clear(image_files *_this);

void image_files_write(image_files *_this,const char *_name,const char *_suf);

typedef struct image_data image_data;

struct image_data{
  const char      *name;
  int              nxblocks;
  int              nyblocks;
  unsigned char   *mode;
  double          *weight;
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
void image_data_pred_block(image_data *_this,int _bi,int _bj);
void image_data_stats_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,intra_stats *_stats);
void image_data_idct_block(image_data *_this,int _bi,int _bj);
void image_data_post_block(image_data *_this,int _bi,int _bj);
void image_data_files_block(image_data *_this,const unsigned char *_data,
 int _stride,int _bi,int _bj,image_files *_files);

int image_data_save_map(image_data *_this);
int image_data_load_map(image_data *_this);

extern int NE_FILTER_PARAMS4[4];

extern const od_filter_func NE_PRE_FILTER[OD_NBSIZES];
extern const od_filter_func NE_POST_FILTER[OD_NBSIZES];

extern double NE_PRED_OFFSETS_4x4[OD_INTRA_NMODES][4][4];
extern double NE_PRED_WEIGHTS_4x4[OD_INTRA_NMODES][4][4][2*4][3*4];

extern double NE_PRED_OFFSETS_8x8[OD_INTRA_NMODES][8][8];
extern double NE_PRED_WEIGHTS_8x8[OD_INTRA_NMODES][8][8][2*8][3*8];

extern double NE_PRED_OFFSETS_16x16[OD_INTRA_NMODES][16][16];
extern double NE_PRED_WEIGHTS_16x16[OD_INTRA_NMODES][16][16][2*16][3*16];

void ne_intra_pred4x4_mult(const od_coeff *_c,int _stride,int _mode,double *_p);
void print_betas(FILE *_fp);

#endif
