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
# if OD_GNUC_PREREQ(4, 0)
#  pragma GCC visibility push(default)
# endif

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
 * \param op The current header packet to process.
 * \return A positive value indicates that a Daala header was successfully
            processed.
 * \retval 0 The first video data packet was encountered after all required
              header packets were parsed.
             The packet just passed to this call should be saved and fed to
              daaala_decode_packet_in() to begin decoding video data.
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
 daala_comment *dc, daala_setup_info **ds, const ogg_packet *op);

/**Allocates a decoder instance.
 * \param info A #daala_info struct filled via daala_decode_header_in().
 * \param setup A #daala_setup_info handle returned via
 *               daala_decode_header_in().
 * \return The initialized #daala_dec_ctx handle.
 * \retval NULL If the decoding parameters were invalid.*/
extern daala_dec_ctx *daala_decode_alloc(const daala_info *info,
 const daala_setup_info *setup);
/**Releases all storage used for the decoder setup information.
 * This should be called after you no longer want to create any decoders for
 *  a stream whose headers you have parsed with daala_decode_header_in().
 * \param setup The setup information to free.
 *              This can safely be <tt>NULL</tt>.*/
extern void daala_setup_free(daala_setup_info *setup);
/**Decoder control function.
 * This is used to provide advanced control of the decoding process.
 * \param dec A #daala_dec_ctx handle.
 * \param req The control code to process.
 *            See \ref decctlcodes "the list of available control codes"
 *             for details.
 * \param buf The parameters for this control code.
 * \param buf_sz The size of the parameter buffer.*/
extern int daala_decode_ctl(daala_dec_ctx *dec,
 int req, void *buf, size_t buf_sz);
/**Frees an allocated decoder instance.
 * \param dec A #daala_dec_ctx handle.*/
extern void daala_decode_free(daala_dec_ctx *dec);
/**Retrieves decoded video data frames.
 * \param dec A #daala_dec_ctx handle.
 * \param img A buffer to receive the decoded image data.
 * \param op An incoming Ogg packet.*/
extern int daala_decode_packet_in(daala_dec_ctx *dec, od_img *img,
 const ogg_packet *op);
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

# if OD_GNUC_PREREQ(4, 0)
#  pragma GCC visibility pop
# endif
# if defined(__cplusplus)
}
# endif

#endif
