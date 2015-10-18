/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/*Daala video codec
Copyright (c) 2006-2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

/**\mainpage
 *
 * \section intro Introduction
 *
 * This is the documentation of the C API for the reference implementation
 *  of Daala, a free video codec.
 *
 * \subsection Organization
 * The implementation is subdivided into three separate libraries.
 * - <tt>libdaalaenc</tt> contains routines exclusive to the encoder.
 *   You must link to this if you use any of the functions listed in
 *    \ref encfuncs.
 * - <tt>libdaaladec</tt> contains routines exclusive to the decoder.
 *   You must link to this if you use any of the functions listed in
 *    \ref decfuncs.
 * - <tt>libdaalabase</tt>, containing routines shared by the encoder and the
 *    decoder.
 *   You must link to this if you use any of the functions in this API,
 *    regardless of whether you use the encoder, decoder or both.*/

/**\file
 * The shared <tt>libdaala</tt> C API.*/
#if !defined(_daala_codec_H)
# define _daala_codec_H (1)

# if defined(__cplusplus)
extern "C" {
# endif

#include "daala_integer.h"

/*Enable special features for gcc and compatible compilers.*/
# if defined(__GNUC__) && defined(__GNUC_MINOR__) && defined(__GNUC_PATCHLEVEL__)
#  define OD_GNUC_PREREQ(maj, min, pat)                             \
  ((__GNUC__ << 16) + (__GNUC_MINOR__ << 8) + __GNUC_PATCHLEVEL__ >= ((maj) << 16) + ((min) << 8) + pat)
# else
#  define OD_GNUC_PREREQ(maj, min, pat) (0)
# endif

#if OD_GNUC_PREREQ(4, 0, 0)
# pragma GCC visibility push(default)
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
# define OD_WARN_UNUSED_RESULT __attribute__((__warn_unused_result__))
#else
# define OD_WARN_UNUSED_RESULT
#endif

#if OD_GNUC_PREREQ(3, 4, 0)
# define OD_ARG_NONNULL(x) __attribute__((__nonnull__(x)))
#else
# define OD_ARG_NONNULL(x)
#endif

/*TODO: remove this ugliness*/
# if defined(_MSC_VER)
#  pragma warning(disable:4100 4115 4125 4127 4152 4505 4554 4711)
# endif

/**\name Error codes*/
/*@{*/
/**No error occurred.*/
#define OD_SUCCESS (0)
/**An invalid pointer was provided.*/
# define OD_EFAULT (-1)
/**An invalid argument was provided.*/
# define OD_EINVAL (-10)
/**The contents of the header were incomplete, invalid, or unexpected.*/
# define OD_EBADHEADER (-20)
/**The header does not belong to a Daala stream.*/
# define OD_ENOTFORMAT (-21)
/**The bitstream version is too high.*/
# define OD_EVERSION (-22)
/**The specified function is not implemented.*/
# define OD_EIMPL (-23)
/**There were errors in the video data packet.*/
# define OD_EBADPACKET (-24)
/*@}*/

/**\name Colorspaces
 * The currently defined color space tags.*/
/*@{*/
/**The color space was not specified at the encoder.
 * It may be conveyed by an external means.*/
# define OD_CS_UNSPECIFIED (0)
/**A Y'CbCr color space designed for NTSC content*/
# define OD_CS_ITU_REC_470M (1)
/**A Y'CbCr color space designed for PAL/SECAM content.*/
# define OD_CS_ITU_REC_470BG (2)
/**A Y'CbCr color space designed for HD content.*/
# define OD_CS_ITU_REC_790 (3)
/**A Y'CgCo color space designed for sRGB content.*/
# define OD_CS_YCgCo (4)
/**The total number of currently defined color spaces.*/
# define OD_CS_NSPACES (5)
/*@}*/

/**The maximum number of color planes allowed in a single frame.*/
# define OD_NPLANES_MAX (4)

typedef struct od_img_plane od_img_plane;
typedef struct od_img od_img;
typedef struct daala_plane_info daala_plane_info;
typedef struct daala_info daala_info;
typedef struct daala_comment daala_comment;

const char *daala_version_string(void);

/**Initialize the logging module.
   This should be called once before invoking any other Daala functions if
    logging support is enabled and you want to use it.
   Otherwise no logging messages will be printed.
   \retval 0 Success.
             This function always succeeds.
             It only returns a value for convenience (e.g., for use in a static
              initializer).*/
int daala_log_init(void);

/** Representation of a single component within an image or frame. */
struct od_img_plane {
  /** Image data is stored as an unsigned octet type whether it's
      actually 8 bit or a multi-byte depth. */
  unsigned char *data;
  /** The decimation factor in the x and y direction. Pixels are reduced
      by a factor of 2^xdec so 0 is none, 1 is decimated by a factor of 2.
      ( YUV420 will  have xdec of 1 and ydec also of 1. YUV444 will have
      xdec and ydec set to zero ). */
  unsigned char xdec;
  unsigned char ydec;
  /** Distance in memory between two pixels horizontally next to each other.
      The value is in bytes regardless of the 'actual' underlying depth
      (either unsigned bytes for 8 bit video or unsigned 16 bit shorts for
      high-depth video). The xstride may be larger than the actual data
      width calculated from the bitdepth; this implies packed rather than
      planar data. */
  int xstride;
  /** Distance in memory between two pixels vertically next to each other.
      As with xstride, this value is always in bytes. */
  int ystride;
  /** 8 for 'normal' video precision; data is unsigned bytes centered on 128.
      Greater-than-8 indicates high-depth video; data is unnormalized
      host-endian order unsigned signed 16-bit shorts (two octets).
      For example, 10 bit video would declare a bit depth of 10, use the
      lower 10 bits of each 16 bit short, and center on 512. */
  int bitdepth;
};

/** Representation of an image or video frame. */
struct od_img {
  /** Typical 3 planes for Y, Cb, and  Cr. Can have a 4th plane for alpha */
  od_img_plane planes[OD_NPLANES_MAX];
  /** Number of planes (1 for greyscale, 3 for YCbCr, 4 for YCbCr+Alpha ) */
  int nplanes;
  /** Width and height in pixels */
  int32_t width;
  int32_t height;
};

/** Subsampling factors for a plane as a power of 2.
 *  4:2:0 would have {0, 0} for Y and {1, 1} for Cb and Cr. */
struct daala_plane_info {
  unsigned char xdec;
  unsigned char ydec;
};

/**\name Bit Depths
 * The three video bit depths currently supported by Daala.*/
/*@{*/
/**8-bit mode.*/
#define OD_BITDEPTH_MODE_8 (1)
/**10-bit mode.*/
#define OD_BITDEPTH_MODE_10 (2)
/**12-bit mode.*/
#define OD_BITDEPTH_MODE_12 (3)
/*@}*/

/** Configuration parameters for a codec instance. */
struct daala_info {
  unsigned char version_major;
  unsigned char version_minor;
  unsigned char version_sub;
  /** pic_width,_height form a region of interest to encode */
  int32_t pic_width;
  int32_t pic_height;
  uint32_t pixel_aspect_numerator;
  uint32_t pixel_aspect_denominator;
  uint32_t timebase_numerator;
  uint32_t timebase_denominator;
  uint32_t frame_duration;
  int keyframe_granule_shift;
  /** bitdepth_mode is one of the three OD_BITDEPTH_MODE_X choices allowed
   * above. */
  int bitdepth_mode;
  int nplanes;
  daala_plane_info plane_info[OD_NPLANES_MAX];
   /** key frame rate defined how often a key frame is emitted by encoder in
    * number of frames. So 10 means every 10th frame is a keyframe.  */
  int keyframe_rate;
};

typedef struct {
  unsigned char *packet;
  long bytes;
  long b_o_s;
  long e_o_s;

  int64_t granulepos;
  int64_t packetno;
} daala_packet;

void daala_info_init(daala_info *info);
void daala_info_clear(daala_info *info);

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
 * In filling in this structure, daala_decode_header_in() will null-terminate
 *  the user_comment strings for safety.
 * However, the bitstream format itself treats them as 8-bit clean vectors,
 *  possibly containing null characters, and so the length array should be
 *  treated as their authoritative length.*/
struct daala_comment {
  /**The array of comment string vectors.*/
  char **user_comments;
  /**An array of the corresponding lengths of each vector, in bytes.*/
  int *comment_lengths;
  /**The total number of comment strings.*/
  int comments;
  /**The null-terminated vendor string.
     This identifies the software used to encode the stream.*/
  char *vendor;
};

/**Initializes a daala_comment section. Users should free the returned
   data with daala_comment_clear().
   \param dc A #daala_comment structure.*/
void daala_comment_init(daala_comment *dc);
/**Free resources allocated for metadata.
   \param dc A #daala_comment structure.*/
void daala_comment_clear(daala_comment *dc);

int64_t daala_granule_basetime(void *encdec, int64_t granpos);
double daala_granule_time(void *encdec, int64_t granpos);
/**Determines whether a Daala packet is a header or not.
   This function does no verification beyond checking the packet type bit, so
    it should not be used for bitstream identification.
   Use daala_decode_headerin() for that.
   \param dpkt A daala_packet structure.
   \retval 1 The packet is a header packet.
   \retval 0 The packet is a video data packet.*/
int daala_packet_isheader(daala_packet *dpkt);
/**Determines whether a Daala packet is a key frame or not.
   This function does no verfication beyond checking the packet type and key
    frame bits, so it should not be used for bitstream identification.
   Feed the packet to an actual decoder for that.
   \param dpkt A daala_packet structure.
   \retval 1  The packet contains a key frame.
   \retval 0  The packet contains a delta frame.*/
int daala_packet_iskeyframe(daala_packet *dpkt);

# if OD_GNUC_PREREQ(4, 0, 0)
#  pragma GCC visibility pop
# endif
# if defined(__cplusplus)
}
# endif

#endif
