#if !defined(_tools_init_intra_tools_H)
# define _tools_init_intra_tools_H (1)
#include "vidinput.h"
#include "../src/odintrin.h"

#define B_SZ_LOG (2)
#define B_SZ     (1<<B_SZ_LOG)

#define OD_INTRA_DC (0)
#define OD_INTRA_TM (1)
#define OD_INTRA_HU (2)
#define OD_INTRA_HE (3)
#define OD_INTRA_HD (4)
#define OD_INTRA_RD (5)
#define OD_INTRA_VR (6)
#define OD_INTRA_VE (7)
#define OD_INTRA_VL (8)
#define OD_INTRA_LD (9)
#define OD_INTRA_NMODES (10)



typedef int (*plane_start_func)(void *_ctx,const char *_name,
 const th_info *_ti,int _pli,int _nxblocks,int _nyblocks);
typedef void (*block_func)(void *_ctx,const unsigned char *_data,int _stride,
 int _bi,int _bj);
typedef int (*plane_finish_func)(void *_ctx);



/*Computes the starting offset and number of blocks which can be intra
   predicted with full context (i.e., all of their neighbors) available.*/
void get_intra_dims(const th_info *_ti,int _pli,
 int *_x0,int *_y0,int *_nxblocks,int *_nyblocks);

char *get_map_filename(const char *_name,int _pli,int _nxblocks,int _nyblocks);

int apply_to_blocks(void *_ctx,plane_start_func _start,block_func _block,
 plane_finish_func _finish,int _argc,const char **_argv);

#endif
