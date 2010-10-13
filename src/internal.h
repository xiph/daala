#if !defined(_internal_H)
# define _internal_H (1)
# include <limits.h>
# include "../include/daala/codec.h"
# include "odintrin.h"



# define OD_VERSION_MAJOR (0)
# define OD_VERSION_MINOR (0)
# define OD_VERSION_SUB   (0)

# define OD_VENDOR_STRING "derf's experimental encoder library " __DATE__

/*Constants for the packet state machine common between encoder and decoder.*/

/*Next packet to emit/read: Codec info header.*/
#define OD_PACKET_INFO_HDR    (-3)
/*Next packet to emit/read: Comment header.*/
#define OD_PACKET_COMMENT_HDR (-2)
/*Next packet to emit/read: Codec setup header.*/
#define OD_PACKET_SETUP_HDR   (-1)
/*Next more packets to emit/read.*/
#define OD_PACKET_DONE        (INT_MAX)



#if 1
/*Currently this structure is only in Tremor, and is read-only.*/
typedef struct oggbyte_buffer oggbyte_buffer;

/*Simple libogg1-style buffer.*/
struct oggbyte_buffer{
  unsigned char *buf;
  unsigned char *ptr;
  long           storage;
};

/*Encoding functions.*/
void oggbyte_writeinit(oggbyte_buffer *_b);
void oggbyte_writetrunc(oggbyte_buffer *_b,long _bytes);
void oggbyte_write1(oggbyte_buffer *_b,unsigned _value);
void oggbyte_write4(oggbyte_buffer *_b,ogg_uint32_t _value);
void oggbyte_writecopy(oggbyte_buffer *_b,void *_source,long _bytes);
void oggbyte_writeclear(oggbyte_buffer *_b);
/*Decoding functions.*/
void oggbyte_readinit(oggbyte_buffer *_b,unsigned char *_buf,long _bytes);
int oggbyte_look1(oggbyte_buffer *_b);
int oggbyte_look4(oggbyte_buffer *_b,ogg_uint32_t *_val);
void oggbyte_adv1(oggbyte_buffer *_b);
void oggbyte_adv4(oggbyte_buffer *_b);
int oggbyte_read1(oggbyte_buffer *_b);
int oggbyte_read4(oggbyte_buffer *_b,ogg_uint32_t *_val);
/*Shared functions.*/
void oggbyte_reset(oggbyte_buffer *_b);
long oggbyte_bytes(oggbyte_buffer *_b);
unsigned char *oggbyte_get_buffer(oggbyte_buffer *_b);

#endif

int od_ilog(ogg_uint32_t _v);
void **od_malloc_2d(size_t _height,size_t _width,size_t _sz);
void **od_calloc_2d(size_t _height,size_t _width,size_t _sz);
void od_free_2d(void *_ptr);

#endif
