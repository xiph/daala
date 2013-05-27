/*Daala video codec
Copyright (c) 2013 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/


#if !defined(_logging_file_H)
# define _path_file_H (1)

#include <stdarg.h>

/* Exhaustive list of all the log facilities. */
typedef enum {
  OD_LOG_GENERIC = 0,
  OD_LOG_ENTROPY_CODER,

  /* Add new facilities here. */

  OD_LOG_FACILITY_MAX
} od_log_facility;

/* The log levels.

   We list integers here to make it easier to know how to set
   logging levels externally.
*/
typedef enum {
  OD_LOG_INVALID = 0,
  OD_LOG_ERR = 1,
  OD_LOG_WARN = 2,
  OD_LOG_NOTICE = 3,
  OD_LOG_INFO = 4,
  OD_LOG_DEBUG = 5,
  OD_LOG_LEVEL_MAX
} od_log_level;

/* Initializer for logging. This MUST be called prior to engaging
    any kind of multithreaded constructs, lest it create race
    conditions.

   od_log_init() determines what levels to log at by examining
    the environment variable OD_LOG_MODULES, which is of the
    form "<facility-name>:<level>(,<facility-name>:<level>)*".

   For instance:
    OD_LOG_MODULES=generic:5

   Messages are logged if the level for a facility is >= the
    level passed to OD_LOG
 */
typedef int (*od_logger_function)(od_log_facility facility,
                                  od_log_level level,
                                  const char *fmt, va_list ap);
int od_log_init(od_logger_function logger);


/* To log a message, use OD_LOG, as follows:

   OD_LOG((OD_LOG_GENERIC, OD_LOG_WARN, "Error code = %d", 22));

   Note that the double parentheses are not an accident. They
    are necessary to compensate for C89's lack of variadic
    macros. This way, in non-debugging mode, the entire block
    is not compiled.
*/

#ifdef OD_LOGGING_ENABLED
# define OD_LOG(a)
#else
# define OD_LOG(a) od_log a

int od_log(od_log_facility fac, od_log_level level, const char *fmt, ...)
#ifdef __GNUC__
  __attribute__ ((format (printf, 3, 4)))
#endif
;

#endif

#endif
