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

/**\file
 * The <tt>libdaala</tt> C decoding API.*/
#if !defined(_daala_daaladec_H)
# define _daala_daaladec_H (1)
# include "codec.h"

# if defined(__cplusplus)
extern "C" {
# endif
# if OD_GNUC_PREREQ(4, 0, 0)
#  pragma GCC visibility push(default)
# endif

#define OD_DECCTL_SET_BSIZE_BUFFER (7001)
#define OD_DECCTL_SET_FLAGS_BUFFER (7003)
#define OD_DECCTL_SET_MV_BUFFER    (7005)
/** Copy the motion compensated reference into a user supplied od_img.
 * \param[in]  <tt>od_img*</tt>: Pointer to the user supplied od_img.
 *              Image must be allocated by the caller, and must be the
 *              same format as the decoder output images. */
#define OD_DECCTL_SET_MC_IMG       (7007)
#define OD_DECCTL_GET_ACCOUNTING   (7009)
#define OD_DECCTL_SET_ACCOUNTING_ENABLED (7011)
#define OD_DECCTL_SET_DERING_BUFFER (7013)


#define OD_ACCT_FRAME (10)
#define OD_ACCT_MV (11)

typedef struct {
  /** x position in units of 4x4 luma blocks for layers 0-3, or vx for
     OD_ACCT_MV. Has no meaning for OD_ACCT_FRAME.*/
  int16_t x;
  /** y position in units of 4x4 luma blocks for layers 0-3, or vy for
     OD_ACCT_MV. Has no meaning for OD_ACCT_FRAME.*/
  int16_t y;
  /** layers (0..NPLANES) for color plane coefficients, or one of
      OD_ACCT_FRAME and OD_ACCT_MV. */
  unsigned char layer;
  /** For layers 0-3, 0 means 4x4, 1, means 8x8, and so on. For OD_ACCT_MV,
     it is the motion vector level. Has no meaning for OD_ACCT_FRAME. */
  unsigned char level;
  /** Integer id in the dictionary. */
  unsigned char id;
  /** Number of bits in units of 1/8 bit. */
  unsigned char bits_q3;
} od_acct_symbol;

/* Max number of entries for symbol types in the dictionary (increase as
   necessary). */
#define MAX_SYMBOL_TYPES (256)

/** Dictionary for translating strings into id. */
typedef struct {
  char *(str[MAX_SYMBOL_TYPES]);
  int nb_str;
} od_accounting_dict;

typedef struct {
  /** All recorded symbols decoded. */
  od_acct_symbol *syms;
  /** Number of symbols actually recorded. */
  int nb_syms;
  /** Dictionary for translating strings into id. */
  od_accounting_dict dict;
} od_accounting;


/**\name Decoder state
   The following data structures are opaque, and their contents are not
    publicly defined by this API.
   Referring to their internals directly is unsupported, and may break without
    warning.*/
/*@{*/
/**The decoder context.*/
typedef struct daala_dec_ctx daala_dec_ctx;
/**Setup information.
   This contains auxiliary information decoded from the setup header by
    daala_decode_header_in() to be passed to daala_decode_alloc().
   It can be re-used to initialize any number of decoders, and can be freed
    via daala_setup_free() at any time.*/
typedef struct daala_setup_info daala_setup_info;
/*@}*/

/**\defgroup decfuncs Functions for Decoding*/
/*@{*/
/**\name Functions for decoding
 * You must link to <tt>libdaalaadec</tt> if you use any of the functions in
 *  this section.
 **/
/*@{*/

/**Parses the header packets from an Ogg Daala stream.
 * To use this function:
 * - Initialize a daala_info structure using daala_info_init().
 * - Initialize a daala_comment structure using daala_comment_init().
 * - Initialize a daala_setup_info pointer to NULL.
 * - Call this function three times, passing in pointers to the same
 *    daala_info, daala_comment, and daala_setup_info pointer each time, and
 *    the three successive header packets.
 * \param info The #daala_info structure to fill in.
 *             This must have been previously initialized with
                daala_info_init().
               The application may begin using the contents of this structure
                after the first header is decoded, though it must continue to
                be passed in unmodified on all subsequent calls.
 * \param dc The #daala_comment structure to fill in.
             This must have been previously initialized with
              daala_comment_init().
             The application may immediately begin using the contents of this
              structure after the second header is decoded, though it must
              continue to be passed in on all subsequent calls.
 * \param ds A pointer to a daala_setup_info pointer to fill in.
             The contents of this pointer must be initialized to <tt>NULL</tt>
              on the first call, and the returned value must continue to be
              passed in on all subsequent calls.
 * \param dp The current header packet to process.
 * \return A positive value indicates that a Daala header was successfully
            processed and indicates the remaining number of headers to be read.
 * \retval 0 The last header was processed and the next packet will
              contain video.
 * \retval OD_EFAULT One of \a info, \a dc, or \a ds was <tt>NULL</tt>, or
                      there was a memory allocation failure.
 * \retval OD_EBADHEADER \a op was <tt>NULL</tt>, the packet was not the next
                          header packet in the expected sequence, or the
                          format fo the header data was invalid.
 * \retval OD_EVERSION The packet data was a Daala header, but for a bitstream
                        version not decodable with this version of
                        <tt>libdaaladec</tt>.
 * \retval OD_ENOTFORMAT The packet was not a Daala header.*/
int daala_decode_header_in(daala_info *info,
 daala_comment *dc, daala_setup_info **ds, const daala_packet *dp);

/**Allocates a decoder instance.
 * \param info A #daala_info struct filled via daala_decode_header_in().
 * \param setup A #daala_setup_info handle returned via
 *               daala_decode_header_in().
 * \return The initialized #daala_dec_ctx handle.
 * \retval NULL If the decoding parameters were invalid.*/
daala_dec_ctx *daala_decode_alloc(const daala_info *info,
 const daala_setup_info *setup);
/**Releases all storage used for the decoder setup information.
 * This should be called after you no longer want to create any decoders for
 *  a stream whose headers you have parsed with daala_decode_header_in().
 * \param setup The setup information to free.
 *              This can safely be <tt>NULL</tt>.*/
void daala_setup_free(daala_setup_info *setup);
/**Decoder control function.
 * This is used to provide advanced control of the decoding process.
 * \param dec A #daala_dec_ctx handle.
 * \param req The control code to process.
 *            See \ref decctlcodes "the list of available control codes"
 *             for details.
 * \param buf The parameters for this control code.
 * \param buf_sz The size of the parameter buffer.*/
int daala_decode_ctl(daala_dec_ctx *dec,
 int req, void *buf, size_t buf_sz);
/**Frees an allocated decoder instance.
 * \param dec A #daala_dec_ctx handle.*/
void daala_decode_free(daala_dec_ctx *dec);
/**Retrieves decoded video data frames.
 * \param dec A #daala_dec_ctx handle.
 * \param img A buffer to receive the decoded image data.
 * \param dp An incoming Daala packet.*/
int daala_decode_packet_in(daala_dec_ctx *dec, od_img *img,
 const daala_packet *dp);
/*@}*/

/** \defgroup decctlcodes Configuration keys for the decoder ctl interface.
 * Decoder CTL settings.
 *
 * These defines and macros are for altering the behaviour of the decoder
 * through the \ref daala_decode_ctl interface.
 *
 * No keys are currently defined.
 */

/*@}*/

# if OD_GNUC_PREREQ(4, 0, 0)
#  pragma GCC visibility pop
# endif
# if defined(__cplusplus)
}
# endif

#endif
