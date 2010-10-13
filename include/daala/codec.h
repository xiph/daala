/**\mainpage
 *
 * \section intro Introduction
 *
 * This is the documentation of the <tt>libdaala</tt> C API.
 * <tt>libdaala</tt> is the reference implementation for Daala, a free video
 *  codec.
 *
 * \subsection Organization
 * <tt>libdaala</tt> is actually subdivided into three separate libraries.
 * - <tt>libdaalaenc<tt>, containing routines exclusive to the encoder.
 *   You must link to this if you use any of the functions listed in
 *    \ref encfuncs.
 * - <tt>libdaaladec</tt>, containing routines exclusive to the decoder.
 *   You must link to this if you use any of the functions listed in
 *    \ref decfuncs.
 * - <tt>libdaalabase</tt>, containing routines shared by the encoder and the
 *    decoder.
 *   You must link to this if you use any of the functions in this API, not
 *    just those listed in \ref basefuncs.*/

/**\file
 * The shared <tt>libdaala<tt/> C API.*/
#if !defined(_O_DAALA_CODEC_H_)
# define _O_DAALA_CODEC_H_ (1)
/*Pick up typedefs.*/
#include <ogg/ogg.h>

#if defined(__cplusplus)
extern "C" {
#endif
#if __GNUC_PREREQ(4,0)
# pragma GCC visibility push(default)
#endif

/*TODO: remove this ugliness*/
# if defined(_MSC_VER)
#  pragma warning(disable:4100 4115 4125 4127 4152 4505 4554 4711)
# endif

/**\name Error codes*/
/*@{*/
/**An invalid pointer was provided.*/
#define OD_EFAULT     (-1)
/**An invalid argument was provided.*/
#define OD_EINVAL     (-10)
/**The contents of the header were incomplete, invalid, or unexpected.*/
#define OD_EBADHEADER (-20)
/**The header does not belong to a Daala stream.*/
#define OD_ENOTFORMAT (-21)
/**The bitstream version is too high.*/
#define OD_EVERSION   (-22)
/**The specified function is not implemented.*/
#define OD_EIMPL      (-23)
/**There were errors in the video data packet.*/
#define OD_EBADPACKET (-24)
/*@}*/

/**\name Colorspaces
 * The currently defined color space tags.*/
/*@{*/
/**The color space was not specified at the encoder.
 * It may be conveyed by an external means.*/
#define OD_CS_UNSPECIFIED   (0)
/**A Y'CbCr color space designed for NTSC content*/
#define OD_CS_ITU_REC_470M  (1)
/**A Y'CbCr color space designed for PAL/SECAM content.*/
#define OD_CS_ITU_REC_470BG (2)
/**A Y'CbCr color space designed for HD content.*/
#define OD_CS_ITU_REC_790   (3)
/**A Y'CgCo color space designed for sRGB content.*/
#define OD_CS_YCgCo         (4)
/**The total number of currently defined color spaces.*/
#define OD_CS_NSPACES       (5)
/*@}*/

typedef struct od_img_plane od_img_plane;
typedef struct od_img       od_img;
typedef struct daala_comp   daala_comp;
typedef struct daala_info   daala_info;

struct od_img_plane{
  unsigned char *data;
  char           xdec;
  char           ydec;
  char           xstride;
  int            ystride;
};

struct od_img{
  od_img_plane *planes;
  int           nplanes;
  int           width;
  int           height;
};

struct daala_comp{
  int xdec;
  int ydec;
};

struct daala_info{
  unsigned char  version_major;
  unsigned char  version_minor;
  unsigned char  version_sub;
  ogg_uint32_t   frame_width;
  ogg_uint32_t   frame_height;
  ogg_uint32_t   pic_x;
  ogg_uint32_t   pic_y;
  ogg_uint32_t   pic_width;
  ogg_uint32_t   pic_height;
  ogg_uint32_t   pixel_aspect_numerator;
  ogg_uint32_t   pixel_aspect_denominator;
  ogg_uint32_t   timebase_numerator;
  ogg_uint32_t   timebase_denominator;
  ogg_uint32_t   frame_duration;
  unsigned       keyframe_granule_shift;
  int            ncomps;
  daala_comp    *comps;
};

void daala_info_init(daala_info *_info);
void daala_info_clear(daala_info *_info);

/**The comment information.
 *
 * This structure holds in the in-stream metadata corresponding to the
 *  'comment' header packet.
 * The comment header is meant to be used much like someone jotting a quick
 *  note on the label of a disc.
 * It should be a short, to the point text note that can be more than a couple
 *  words, but not more than a short paragraph.
 *
 * The metadata is stored as a series of (tag, value) pairs, in length-encoded
 *  string vectors.
 * The first occurrence of the '=' character delimits the tag and value.
 * A particular tag may occur more than once, and order is significant.
 * The character set encoding for the strings is always UTF-8, but the tag
 *  names are limited to ASCII, and treated as case-insensitive.
 * See the Daala specification for details.
 *
 * In filling in this structure, daala_decode_header() will null-terminate the
 *  user_comment strings for safety.
 * However, the bitstream format itself treats them as 8-bit clean vectors,
 *  possibly containing null characters, and so the length array should be
 *  treated as their authoritative length.*/
typedef struct daala_comment{
  /**The array of comment string vectors.*/
  char **user_comments;
  /**An array of the corresponding lengths of each vector, in bytes.*/
  int   *comment_lengths;
  /**The total number of comment strings.*/
  int    comments;
  /**The null-terminated vendor string.
     This identifies the software used to encode the stream.*/
  char  *vendor;
}daala_comment;

void daala_comment_init(daala_comment *_dc);
void daala_comment_clear(daala_comment *_dc);

ogg_int64_t daala_granule_basetime(void *_encdec,ogg_int64_t _granpos);
double      daala_granule_time(void *_encdec,ogg_int64_t _granpos);

#if __GNUC_PREREQ(4,0)
# pragma GCC visibility pop
#endif
#if defined(__cplusplus)
}
#endif

#endif
