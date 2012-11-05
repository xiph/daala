/*Daala video codec
Copyright (c) 2006-2010 Daala project contributors.  All rights reserved.

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
 * The <tt>libdaala</tt> C encoding API.*/
#if !defined(_O_DAALA_DAALAENC_H)
# define _O_DAALA_DAALAENC_H (1)
# include "codec.h"

#if defined(__cplusplus)
extern "C" {
#endif
#if __GNUC_PREREQ(4,0)
# pragma GCC visibility push(default)
#endif

/**\name Encoder state
 * The following data structure is opaque, and its contents are not publicly
 *  defined by this API.
 * Referring to its internals directly is unsupported, and may break without
 *  warning.*/
/*@{*/
/**The encoder context.*/
typedef struct daala_enc_ctx daala_enc_ctx;
/*@}*/

/**\defgroup encfuncs Functions for Encoding*/
/*@{*/
/**\name Functions for encoding
 * You must link to <tt>libdaalabase</tt> and <tt>libdaalaenc</tt> if you use
 *  any of the functions in this section.
 *
 * The functions are listed in the order they are used in a typical encode.
 * The basic steps are:
 * - Fill in a #daala_info structure with details on the format of the video
 *    you wish to encode.
 * - Allocate a #daala_enc_ctx handle with daala_encode_alloc().
 * - Perform any additional encoder configuration required with
 *    daala_encode_ctl().
 * - Repeatedly call daala_encode_flusheader() to retrieve all the header
 *    packets.
 * - For each uncompressed frame:
 *   - Submit the compressed frame via daala_encode_img_in().
 *   - Repeatedly call daala_encode_packet_out() to retrieve any video data
 *      packets that are ready.
 * - Call daala_encode_free() to release all encoder memory.*/
/*@{*/
/**Allocates an encoder instance.
 * \param _info A #daala_info struct filled with the desired encoding
 *               parameters.
 * \return The initialized #daala_enc_ctx handle.
 * \retval NULL if the encoding parameters were invalid.*/
extern daala_enc_ctx *daala_encode_alloc(const daala_info *_info);
/**Encoder control function.
 * This is used to provide advanced control of the encoding process.
 * \param _enc    A #daala_enc_ctx handle.
 * \param _req    The control code to process.
 *                See \ref encctlcodes "the list of available control codes"
 *                 for details.
 * \param _buf    The parameters for this control code.
 * \param _buf_sz The size of the parameter buffer.*/
extern int daala_encode_ctl(daala_enc_ctx *_enc,int _req,void *_buf,
 size_t _buf_sz);
/**Outputs the next header packet.
 * This should be called repeatedly after encoder initialization until it
 *  return <tt>0</tt> to get all of the header packets, in order, before
 *  encoding actual video data.
 * \param _enc      A #daala_enc_ctx handle.
 * \param _comments The metadata to place in the comment header, when it is
 *                   encoded.
 * \param _op       An <tt>ogg_packet</tt> structure to fill.
 *                  All of the elements of this structure will be set,
 *                   including a pointer to the header data.
 *                  The memory for the header data is owned by
 *                   <tt>libdaala</tt>.
 * \return A positive value indicates that a header packet was successfully
 *          produced.
 * \retval 0         No packet was produced, and no more header packets remain.
 * \retval OD_EFAULT \a _enc, \a _comments or \a _op was <tt>NULL</tt>.*/
extern int daala_encode_flush_header(daala_enc_ctx *_enc,
 daala_comment *_comments,ogg_packet *_op);
/**Submits an uncompressed frame to the encoder.
 * \param _enc      A #daala_enc_ctx handle.
 * \param _img      A buffer of image data to encode.
 * \param _duration The duration to display the frame for, in timebase units.
 *                  If a non-zero frame duration was specified in the header,
 *                   then this parameter is ignored.
 * \retval 0          Success.
 * \retval OD_EFAULT  \a _enc or \a _img was <tt>NULL</tt>.
 * \retval OD_EINVAL  The image size does not match the frame size the encoder
 *                     was initialized with, or encoding has already
 *                     completed.*/
extern int daala_encode_img_in(daala_enc_ctx *_enc,od_img *_img,int _duration);
/**Retrieves encoded video data packets.
 * This should be called repeatedly after each frame is submitted to flush any
 *  encoded packets, until it returns 0.
 * The encoder will not buffer these packets as subsequent frames are
 *  compressed, so a failure to do so will result in lost video data.
 * \note Current the encoder operates in a one-frame-in, one-packet-out
 *        manner.
 *       However, this may be changed in the future.
 * \param _enc  A #daala_enc_ctx handle.
 * \param _last Set this flag to a non-zero value if no more uncompressed
 *               frames will be submitted.
 *              This ensures that a proper EOS flag is set on the last packet.
 * \param _op   An <tt>ogg_packet</tt> structure to fill.
 *              All of the elements of this structure will be set, including a
 *               pointer to the video data.
 *              The memory for the video data is owned by <tt>libdaala</tt>.
 * \return A positive value indicates that a video data packet was successfully
 *          produced.
 * \retval 0         No packet was produced, and no more encoded video data
 *                    remains.
 * \retval OD_EFAULT \a _enc or \a _op was <tt>NULL</tt>.*/
extern int daala_encode_packet_out(daala_enc_ctx *_enc,int _last,
 ogg_packet *_op);
/**Frees an allocated encoder instance.
 * \param _enc A #daala_enc_ctx handle.*/
extern void daala_encode_free(daala_enc_ctx *_enc);
/*@}*/
/*@}*/

#if __GNUC_PREREQ(4,0)
# pragma GCC visibility pop
#endif
#if defined(__cplusplus)
}
#endif

#endif
