#if !defined(_encint_H)
# define _encint_H (1)
# include "../include/daala/daalaenc.h"
# include "state.h"



typedef struct daala_enc_ctx od_enc_ctx;
typedef struct od_mv_est_ctx od_mv_est_ctx;



/*Constants for the packet state machine specific to the encoder.*/

/*No packet currently ready to output.*/
#define OD_PACKET_EMPTY       (0)
/*A packet ready to output.*/
#define OD_PACKET_READY       (1)

/*The number of fractional bits of precision in our \lambda values.*/
#define OD_LAMBDA_SCALE       (5)



struct daala_enc_ctx{
  od_state        state;
  oggbyte_buffer  obb;
  int             packet_state;
  od_mv_est_ctx  *mvest;
};

od_mv_est_ctx *od_mv_est_alloc(od_enc_ctx *_enc);
void           od_mv_est_free(od_mv_est_ctx *_est);
void           od_mv_est(od_mv_est_ctx *_est,int _ref,int _lambda);


#endif
