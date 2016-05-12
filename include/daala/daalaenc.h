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
# if OD_GNUC_PREREQ(4, 0, 0)
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
 * - Call daala_comment_clear() to release all comments allocated.
 * - Call daala_encode_free() to release all encoder memory.*/
/*@{*/
/**Allocates and initializes an encoder instance.
 * \param info A #daala_info struct filled with the desired encoding
 *              parameters.
 * \return The initialized #daala_enc_ctx handle.
 * \retval NULL if the encoding parameters were invalid.*/
daala_enc_ctx *daala_encode_create(const daala_info *info);
/**Encoder control function.
 * This is used to provide advanced control of the encoding process.
 * \param enc A #daala_enc_ctx handle.
 * \param req The control code to process.
 *            See \ref encctlcodes "the list of available control codes"
 *             for details.
 * \param buf The parameters for this control code.
 * \param buf_sz The size of the parameter buffer.*/
int daala_encode_ctl(daala_enc_ctx *enc,
 int req, void *buf, size_t buf_sz);
/**Outputs the next header packet.
 * This should be called repeatedly after encoder initialization until it
 *  return <tt>0</tt> to get all of the header packets, in order, before
 *  encoding actual video data.
 * \param enc A #daala_enc_ctx handle.
 * \param comments The metadata to place in the comment header, when it is
 *                  encoded. Users should free the returned data with
 *                  daala_comment_clear().
 * \param dp A <tt>daala_packet</tt> structure to fill.
 *           All of the elements of this structure will be set,
 *            including a pointer to the header data.
 *           The memory for the header data is owned by
 *            <tt>libdaala</tt>.
 * \return A positive value indicates that a header packet was successfully
 *          produced.
 * \retval 0 No packet was produced, and no more header packets remain.
 * \retval OD_EFAULT \a enc, \a comments or \a op was <tt>NULL</tt>.*/
int daala_encode_flush_header(daala_enc_ctx *enc,
 daala_comment *comments, daala_packet *dp);
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
int daala_encode_img_in(daala_enc_ctx *enc, daala_image *img, int duration);
/**Retrieves encoded video data packets.
 * This should be called repeatedly after each frame is submitted to flush any
 *  encoded packets, until it returns 0.
 * The encoder will not buffer these packets as subsequent frames are
 *  compressed, so a failure to do so will result in lost video data.
 * \note Current the encoder operates in a one-frame-in, one-packet-out
 *        manner.
 *       However, this may be changed in the future.
 * \param enc A #daala_enc_ctx handle.
 *             This ensures that a proper EOS flag is set on the last packet.
 * \param dp A <tt>daala_packet</tt> structure to fill.
 *           All of the elements of this structure will be set, including a
 *            pointer to the video data.
 *           The memory for the video data is owned by <tt>libdaala</tt>.
 * \return A positive value indicates that a video data packet was successfully
 *          produced.
 * \retval 0 No packet was produced, and no more encoded video data
 *            remains.
 * \retval OD_EFAULT \a enc or \a op was <tt>NULL</tt>.*/
int daala_encode_packet_out(daala_enc_ctx *enc, int last, daala_packet *dp);
/**Frees an allocated encoder instance.
 * \param enc A #daala_enc_ctx handle.*/
void daala_encode_free(daala_enc_ctx *enc);
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
/** Configure the encoder's computational complexity level.
 * \see OD_GET_COMPLEXITY
 * \param[in]  _buf <tt>int</tt>: The new encoder complexity level.
 *                  Values must lie in the range 0...10, inclusive, with
 *                   higher values requiring more CPU but generally producing
 *                   better quality at a given bitrate. */
#define OD_SET_COMPLEXITY 4002
/** Get the encoder's computational complexity level.
 * \see OD_SET_COMPLEXITY
 * \param[in]  _buf <tt>int</tt>: Returns a value in the range 0...10,
 *                   inclusive.
 *                  Higher values indicate higher CPU requirements, but
 *                   generally producing better quality at a given bitrate. */
#define OD_GET_COMPLEXITY 4004
/** Whether activity masking should be used or not.
 * \param[in]  _buf <tt>int</tt>: 0 to disable the use of activity masking,
 *                   a non-zero value otherwise (the default). */
#define OD_SET_ACTIVITY_MASKING 4006
/** Which quantization matrix to use.
 * \param[in]  _buf <tt>int</tt>: 0 => flat quantization matrix,
 *                   1 => HVS (the default). */
#define OD_SET_QM 4008
/** Whether the bilinear postprocessing filter should be used or not.
 * \param[in]  _buf <tt>int</tt>: 0 to disable the bilinear postprocessing filter,
 *                   a non-zero value otherwise (the default). */
#define OD_SET_DERING 4010

/** Whether the motion compensation search should use the chroma planes in
    addition to the luma plane.
 * \param[in]  _buf <tt>int</tt>: 0 to disable the use of the chroma planes,
 *                   a non-zero value otherwise (the default). */
#define OD_SET_MC_CHROMA 4100
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
/** Whether the SATD metric should be used in motion compensation or not.
 * \param[in]  _buf <tt>int</tt>: 0 to disable the use of SATD (the default),
 *                   a non-zero value otherwise. */
#define OD_SET_MC_SATD 4108
/** Number of B frames used by encoder.
 * \param[in]  _buf <tt>int</tt>: The number of B frames.
 *                  Values must lie in the range 0...OD_MAX_B_FRAMES,
 *                  inclusive. */
#define OD_SET_B_FRAMES 4110
/**Sets the current encoding bitrate.
 * Once a bitrate is set, the encoder must use a rate-controlled mode for all
 *  future frames (this restriction may be relaxed in a future version).
 * Due to the buffer delay, the exact bitrate of each section of the encode is
 *  not guaranteed.
 * The encoder may have already used more bits than allowed for the frames it
 *  has encoded, expecting to make them up in future frames, or it may have
 *  used fewer, holding the excess in reserve.
 * The exact transition between two bitrates is not well-defined by this
 *  API, but may be affected by flags set with #OD_SET_RATE_FLAGS.
 * After a number of frames equal to the buffer delay, one may expect further
 *  output to average at the target bitrate.
 *
 * \param[in] buf <tt>long</tt>: The new target bitrate, in bits per second.
 * \retval OD_SUCCESS    Success.
 * \retval OD_EFAULT     \a enc or \a buf is <tt>NULL</tt>.
 * \retval OD_EINVAL     The target bitrate was not positive.
 *                       A future version of this library may allow passing 0
 *                        to disabled rate-controlled mode and return to a
 *                        quality-based mode, in which case this function will
 *                        not return an error for that value.
 * \retval OD_EIMPL      Not supported by this implementation.*/
#define OD_SET_BITRATE 4112
/**Modifies the default bitrate management behavior.
 * Use to allow or disallow frame dropping, and to enable or disable capping
 *  bit reservoir overflows and underflows.
 * See \ref encctlcodes "the list of available flags".
 * The flags are set by default to
 *  <tt>#OD_RATECTL_DROP_FRAMES|#OD_RATECTL_CAP_OVERFLOW</tt>.
 *
 * \param[in] buf <tt>int</tt>: Any combination of
 *                  \ref ratectlflags "the available flags":
 *                 - #OD_RATECTL_DROP_FRAMES: Enable frame dropping.
 *                 - #OD_RATECTL_CAP_OVERFLOW: Don't bank excess bits for later
 *                    use.
 *                 - #OD_RATECTL_CAP_UNDERFLOW: Don't try to make up shortfalls
 *                    later.
 * \retval OD_SUCCESS    Success.
 * \retval OD_EFAULT \a enc or \a buf is <tt>NULL</tt>.
 * \retval OD_EINVAL \a buf_sz is not <tt>sizeof(int)</tt> or rate control
 *                    is not enabled.
 * \retval OD_EIMPL   Not supported by this implementation in the current
 *                    encoding mode.*/
#define OD_SET_RATE_FLAGS 4114
/**Sets the size of the bitrate management bit reservoir as a function
 *  of number of frames.
 * The reservoir size affects how quickly bitrate management reacts to
 *  instantaneous changes in the video complexity.
 * Larger reservoirs react more slowly, and provide better overall quality, but
 *  require more buffering by a client, adding more latency to live streams.
 * By default, libdaala sets the reservoir to the maximum distance between
 *  keyframes, subject to a minimum and maximum limit.
 * This call may be used to increase or decrease the reservoir, increasing or
 *  decreasing the allowed temporary variance in bitrate.
 * An implementation may impose some limits on the size of a reservoir it can
 *  handle, in which case the actual reservoir size may not be exactly what was
 *  requested.
 * The actual value set will be returned.
 *
 * \param[in]  buf <tt>int</tt>: Requested size of the reservoir measured in
 *                   frames.
 * \param[out] buf <tt>int</tt>: The actual size of the reservoir set.
 * \retval OD_EFAULT \a enc or \a buf is <tt>NULL</tt>.
 * \retval OD_EINVAL \a buf_sz is not <tt>sizeof(int)</tt>, or rate control
 *                    is not enabled.  The buffer has an implementation
 *                    defined minimum and maximum size and the value in buf
 *                    will be adjusted to match the actual value set.
 * \retval OD_EIMPL   Not supported by this implementation in the current
 *                    encoding mode.*/
#define OD_SET_RATE_BUFFER 4116
/**Enable pass 1 of two-pass encoding mode and retrieve the first pass metrics.
 * Pass 1 mode must be enabled before the first frame is encoded, and a target
 *  bitrate must have already been specified to the encoder.
 * Although this does not have to be the exact rate that will be used in the
 *  second pass, closer values may produce better results.
 * The first call returns the size of the two-pass header data, along with some
 *  placeholder content, and sets the encoder into pass 1 mode implicitly.
 * This call sets the encoder to pass 1 mode implicitly.
 * Then, a subsequent call must be made after each call to
 *  daala_encode_img_in() to retrieve the metrics for that frame.
 * An additional, final call must be made to retrieve the summary data,
 *  containing such information as the total number of frames, etc.
 * This must be stored in place of the placeholder data that was returned
 *  in the first call, before the frame metrics data.
 * All of this data must be presented back to the encoder during pass 2 using
 *  #OD_2PASS_IN.
 *
 * \param[out] <tt>char *</tt>buf: Returns a pointer to internal storage
 *              containing the two pass metrics data.
 *             This storage is only valid until the next call, or until the
 *              encoder context is freed, and must be copied by the
 *              application.
 * \retval >=0       The number of bytes of metric data available in the
 *                    returned buffer.
 * \retval OD_EFAULT \a enc or \a buf is <tt>NULL</tt>.
 * \retval OD_EINVAL \a buf_sz is not <tt>sizeof(char *)</tt>, no target
 *                    bitrate has been set, or the first call was made after
 *                    the first frame was submitted for encoding.
 * \retval OD_EIMPL   Not supported by this implementation.*/
#define OD_2PASS_OUT 4118
/**Submits two-pass encoding metric data collected the first encoding pass to
 *  the second pass.
 * The first call must be made before the first frame is encoded, and a target
 *  bitrate must have already been specified to the encoder.
 * It sets the encoder to pass 2 mode implicitly; this cannot be disabled.
 * The encoder may require reading data from some or all of the frames in
 *  advance, depending on, e.g., the reservoir size used in the second pass.
 * You must call this function repeatedly before each frame to provide data
 *  until either a) it fails to consume all of the data presented or b) all of
 *  the pass 1 data has been consumed.
 * In the first case, you must save the remaining data to be presented after
 *  the next frame.
 * You can call this function with a NULL argument to get an upper bound on
 *  the number of bytes that will be required before the next frame.
 *
 * When pass 2 is first enabled, the default bit reservoir is set to the entire
 *  file; this gives maximum flexibility but can lead to very high peak rates.
 * You can subsequently set it to another value with #OD_SET_RATE_BUFFER
 *  (e.g., to set it to the keyframe interval for non-live streaming), however,
 *  you may then need to provide more data before the next frame.
 *
 * \param[in] buf <tt>char[]</tt>: A buffer containing the data returned by
 *                  #OD_2PASS_OUT in pass 1.
 *                 You may pass <tt>NULL</tt> for \a buf to return an upper
 *                  bound on the number of additional bytes needed before the
 *                  next frame.
 *                 The summary data returned at the end of pass 1 must be at
 *                  the head of the buffer on the first call with a
 *                  non-<tt>NULL</tt> \a buf, and the placeholder data
 *                  returned at the start of pass 1 should be omitted.
 *                 After each call you should advance this buffer by the number
 *                  of bytes consumed.
 * \retval >0            The number of bytes of metric data required/consumed.
 * \retval 0             No more data is required before the next frame.
 * \retval OD_EFAULT     \a enc is <tt>NULL</tt>.
 * \retval OD_EINVAL     No target bitrate has been set, or the first call was
 *                        made after the first frame was submitted for
 *                        encoding.
 * \retval OD_ENOTFORMAT The data did not appear to be pass 1 from a compatible
 *                        implementation of this library.
 * \retval OD_EBADHEADER The data was invalid; this may be returned when
 *                        attempting to read an aborted pass 1 file that still
 *                        has the placeholder data in place of the summary
 *                        data.
 * \retval OD_EIMPL       Not supported by this implementation.*/
#define OD_2PASS_IN 4120
/*@}*/

/**\name OD_SET_RATE_FLAGS flags
 * \anchor ratectlflags
 * These are the flags available for use with #OD_SET_RATE_FLAGS.*/
/*@{*/
/**Drop frames to keep within bitrate buffer constraints.
 * This can have a severe impact on quality, but is the only way to ensure that
 *  bitrate targets are met at low rates during sudden bursts of activity.
 * It is enabled by default.*/
#define OD_RATECTL_DROP_FRAMES   (0x1)
/**Ignore bitrate buffer overflows.
 * If the encoder uses so few bits that the reservoir of available bits
 *  overflows, ignore the excess.
 * The encoder will not try to use these extra bits in future frames.
 * At high rates this may cause the result to be undersized, but allows a
 *  client to play the stream using a finite buffer; it should normally be
 *  enabled, which is the default.*/
#define OD_RATECTL_CAP_OVERFLOW  (0x2)
/**Ignore bitrate buffer underflows.
 * If the encoder uses so many bits that the reservoir of available bits
 *  underflows, ignore the deficit.
 * The encoder will not try to make up these extra bits in future frames.
 * At low rates this may cause the result to be oversized; it should normally
 *  be disabled, which is the default.*/
#define OD_RATECTL_CAP_UNDERFLOW (0x4)
/*@}*/

# if OD_GNUC_PREREQ(4, 0, 0)
#  pragma GCC visibility pop
# endif
# if defined(__cplusplus)
}
# endif

#endif
