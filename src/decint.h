#if !defined(_encint_H)
# define _encint_H (1)
# include "../include/daala/daaladec.h"
# include "state.h"



typedef struct daala_dec_ctx od_dec_ctx;



/*Constants for the packet state machine specific to the decoder.*/

/*Next packet to read: Data packet.*/
#define OD_PACKET_DATA        (0)



struct daala_dec_ctx{
  od_state        state;
  oggbyte_buffer  obb;
  int             packet_state;
};

#endif
