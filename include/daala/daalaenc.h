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
 * The <tt>libdaala</tt> C encoding API.*/
#if !defined(_daala_daalaenc_H)
# define _daala_daalaenc_H (1)
# include "codec.h"

# if defined(__cplusplus)
extern "C" {
# endif
# if OD_GNUC_PREREQ(4, 0)
#  pragma GCC visibility push(default)
# endif

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
/**Allocates and initializes an encoder instance.
 * \param info A #daala_info struct filled with the desired encoding
 *              parameters.
 * \return The initialized #daala_enc_ctx handle.
 * \retval NULL if the encoding parameters were invalid.*/
extern daala_enc_ctx *daala_encode_create(const daala_info *info);
/**Encoder control function.
 * This is used to provide advanced control of the encoding process.
 * \param enc A #daala_enc_ctx handle.
 * \param req The control code to process.
 *            See \ref encctlcodes "the list of available control codes"
 *             for details.
 * \param buf The parameters for this control code.
 * \param buf_sz The size of the parameter buffer.*/
extern int daala_encode_ctl(daala_enc_ctx *enc,
 int req, void *buf, size_t buf_sz);
/**Outputs the next header packet.
 * This should be called repeatedly after encoder initialization until it
 *  return <tt>0</tt> to get all of the header packets, in order, before
 *  encoding actual video data.
 * \param enc A #daala_enc_ctx handle.
 * \param comments The metadata to place in the comment header, when it is
 *                  encoded.
 * \param op An <tt>ogg_packet</tt> structure to fill.
 *           All of the elements of this structure will be set,
 *            including a pointer to the header data.
 *           The memory for the header data is owned by
 *            <tt>libdaala</tt>.
 * \return A positive value indicates that a header packet was successfully
 *          produced.
 * \retval 0 No packet was produced, and no more header packets remain.
 * \retval OD_EFAULT \a enc, \a comments or \a op was <tt>NULL</tt>.*/
extern int daala_encode_flush_header(daala_enc_ctx *enc,
 daala_comment *comments, ogg_packet *op);
/**Submits an uncompressed frame to the encoder.
 * \param enc A #daala_enc_ctx handle.
 * \param img A buffer of image data to encode.
 * \param duration The duration to display the frame for, in timebase units.
 *                 If a non-zero frame duration was specified in the header,
 *                  then this parameter is ignored.
 * \retval 0 Success.
 * \retval OD_EFAULT \a enc or \a img was <tt>NULL</tt>.
 * \retval OD_EINVAL The image size does not match the frame size the encoder
 *                   was initialized with, or encoding has already
 *                    completed.*/
extern int daala_encode_img_in(daala_enc_ctx *enc, od_img *img, int duration);
/**Retrieves encoded video data packets.
 * This should be called repeatedly after each frame is submitted to flush any
 *  encoded packets, until it returns 0.
 * The encoder will not buffer these packets as subsequent frames are
 *  compressed, so a failure to do so will result in lost video data.
 * \note Current the encoder operates in a one-frame-in, one-packet-out
 *        manner.
 *       However, this may be changed in the future.
 * \param enc A #daala_enc_ctx handle.
 * \param last Set this flag to a non-zero value if no more uncompressed
 *              frames will be submitted.
 *             This ensures that a proper EOS flag is set on the last packet.
 * \param op An <tt>ogg_packet</tt> structure to fill.
 *           All of the elements of this structure will be set, including a
 *            pointer to the video data.
 *           The memory for the video data is owned by <tt>libdaala</tt>.
 * \return A positive value indicates that a video data packet was successfully
 *          produced.
 * \retval 0 No packet was produced, and no more encoded video data
 *            remains.
 * \retval OD_EFAULT \a enc or \a op was <tt>NULL</tt>.*/
extern int daala_encode_packet_out(daala_enc_ctx *enc,
 int last, ogg_packet *op);
/**Frees an allocated encoder instance.
 * \param enc A #daala_enc_ctx handle.*/
extern void daala_encode_free(daala_enc_ctx *enc);
/*@}*/

/** \defgroup encctlcodes Configuration keys for the encoder ctl interface.
 * Encoder CTL settings.
 *
 * These defines and macros are for altering the behaviour of the encoder
 *  through the \ref daala_encode_ctl interface.
 * These should have even values.
 */
/** Set the quantizer scale.
 * The passed buffer is interpreted as containing a single <tt>int</tt>.
 * The valid range is 0-511. */
#define OD_SET_QUANT 4000
/** Whether the motion compensation search should use the chroma planes in
    addition to the luma plane.
 * \param[in]  _buf <tt>int</tt>: 0 to disable the use of the chroma planes,
 *                   a non-zero value otherwise (the default). */
#define OD_SET_MC_USE_CHROMA 4100
/** Minimum motion vectors resolution for the motion compensation search.
 * \param[in]  _buf <tt>int</tt>: 0 => 1/8 pel (default), 1 => 1/4 pel,
 *                   2 => 1/2 pel */
#define OD_SET_MV_RES_MIN 4102
/** Minimum motion vectors level for the motion compensation search. If this
 *  level is greater than the maximum level, the maximum level will be used
 *  instead.
 * \param[in]  _buf <tt>int</tt>: level between 0 and 4
 *                   Default: 0 */
#define OD_SET_MV_LEVEL_MIN 4104
/** Maximum motion vectors level for the motion compensation search.
 * \param[in]  _buf <tt>int</tt>: level between 0 and 4
 *                   Default: 4 */
#define OD_SET_MV_LEVEL_MAX 4106

/*@}*/

# if OD_GNUC_PREREQ(4, 0)
#  pragma GCC visibility pop
# endif
# if defined(__cplusplus)
}
# endif

#endif
