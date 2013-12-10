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

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include "logging.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include "daala/codec.h"
#include "internal.h"

static unsigned long od_log_levels[OD_LOG_FACILITY_MAX] = { 0 };

static const char *od_log_module_names[OD_LOG_FACILITY_MAX] = {
  "generic",
  "encoder",
  "motion-estimation",
  "motion-compensation",
  "entropy-coder",
  "pvq"
};

static const char *od_log_level_names[OD_LOG_LEVEL_MAX] = {
  "INVALID",
  "ERR",
  "WARN",
  "NOTICE",
  "INFO",
  "DEBUG",
};


static int od_log_fprintf_stderr(od_log_facility facility,
 od_log_level level,
 unsigned int flags,
 const char *fmt, va_list ap);
static const char *od_log_facility_name(od_log_facility facility);
static const char *od_log_level_name(od_log_level level);

static od_logger_function od_logger = od_log_fprintf_stderr;


static const char *od_log_facility_name(od_log_facility fac) {
  /* Check for invalid input */
  OD_ASSERT(fac < OD_LOG_FACILITY_MAX);
  if (fac >= OD_LOG_FACILITY_MAX)
    return "INVALID";
  return od_log_module_names[fac];
}

static const char *od_log_level_name(od_log_level level) {
  /* Check for invalid input */
  OD_ASSERT(level < OD_LOG_LEVEL_MAX);
  if (level >= OD_LOG_LEVEL_MAX)
    return "INVALID";
  return od_log_level_names[level];
}

int daala_log_init(void) {
  return od_log_init(NULL);
}

int od_log_init(od_logger_function logger) {
  char *ptr;
  char *nextptr;
  char *comma = NULL;
  char *colon;
  char *endptr;
  unsigned long level;
  int i;
  for (i = 0; i < OD_LOG_FACILITY_MAX; ++i) {
    od_log_levels[i] = 0;
  }
  if (logger)
    od_logger = logger;
  ptr = getenv("OD_LOG_MODULES");
  if (!ptr)
    return 0;
  /* Break up the string. This is manual because strtok is evil. */
  nextptr = ptr;
  for (;;) {
    ptr = nextptr;
    if (!*ptr)
      break;
    comma = strchr(ptr, ',');
    if (comma) {
      *comma = '\0';
      nextptr = comma + 1;
    }
    else {
      nextptr = ptr + strlen(ptr);
    }
    /* At this point, ptr points to a single <facility>:<level> clause */
    colon = strchr(ptr, ':');
    if (!colon) {
      /* Generate an error and skip ahead */
      fprintf(stderr, "Bogus clause '%s'\n", ptr);
      continue;
    }
    /* At this point, we should theoretically have a valid clause */
    *colon = '\0';
    ++colon; /* Now points to the start of the number. */
    for (i = 0; i < OD_LOG_FACILITY_MAX; ++i) {
      if (od_log_module_names[i] && !strcmp(ptr, od_log_module_names[i])) {
        break;  /* Success */
      }
    }
    if (i == OD_LOG_FACILITY_MAX) {
      fprintf(stderr, "Unknown facility '%s'\n", ptr);
      continue;
    }
    if (!*colon) {
      fprintf(stderr, "Empty log level\n");
      continue;
    }
    level = strtoul(colon, &endptr, 10);
    if (*endptr) {
      fprintf(stderr, "Could not convert log level '%s'\n", colon);
      continue;
    }
    od_log_levels[i] = level;
  }
  return 0;
}

int od_logging_active_impl(od_log_facility fac, od_log_level level) {
  /* Check for invalid input */
  OD_ASSERT(fac < OD_LOG_FACILITY_MAX);
  if (fac >= OD_LOG_FACILITY_MAX)
    return 0;
  if (od_log_levels[fac] < level)
    return 0;  /* Skipped */

  return 1;
}

static int od_log_impl(od_log_facility fac, od_log_level level,
 unsigned int flags, const char *fmt, va_list ap) {
  if (!od_logging_active(fac, level))
    return 0;
  (void)od_logger(fac, level, flags, fmt, ap);

  return 0;
}

int od_log(od_log_facility fac, od_log_level level, const char *fmt, ...) {
  va_list ap;
  int rv;

  va_start(ap, fmt);
  rv = od_log_impl(fac, level, 0, fmt, ap);
  va_end(ap);

  return rv;
}

int od_log_partial(od_log_facility fac, od_log_level level,
 const char *fmt, ...) {
  va_list ap;
  int rv;

  va_start(ap, fmt);
  rv = od_log_impl(fac, level, OD_LOG_FLAG_PARTIAL, fmt, ap);
  va_end(ap);

  return rv;
}


static int od_log_fprintf_stderr(od_log_facility facility,
 od_log_level level, unsigned int flags, const char *fmt, va_list ap) {
  char fmt_buffer[1024];
  int rv;
  if (flags & OD_LOG_FLAG_PARTIAL) {
    (void)vfprintf(stderr, fmt, ap);
    return 0;
  }
  rv = snprintf(fmt_buffer, sizeof(fmt_buffer),
   "[%s/%s] %s\n", od_log_facility_name(facility),
   od_log_level_name(level), fmt);
  if ((rv < 0) || (((size_t)rv) >= sizeof(fmt_buffer))) {
    fprintf(stderr, "Error logging. Format string too long\n");
    return OD_EINVAL;
  }
  (void)vfprintf(stderr, fmt_buffer, ap);

  return 0;
}



/* Log various matrix types. Parameters are:

   T == type
   N == the name of the function
   F == fmt parameter

 */
#define DEFINE_OD_LOG_MATRIX(T, N, F) \
int od_log_matrix_##N(od_log_facility facility, \
                      od_log_level level, \
                      const char *prefix, \
                      T *values, \
                      int width, \
                      int height) \
 { \
  size_t buffer_size; \
  size_t current_size; \
  int rv; \
  char *buffer; \
  int h; \
  int w; \
 \
  if (!od_logging_active_impl(facility, level)) \
    return 0; \
 \
  /* Take an initial stab at buffer size */ \
  buffer_size = sizeof(T) * 2 * width * height; \
  current_size = 0; \
  buffer = (char *)_ogg_malloc(buffer_size); \
  if (!buffer) \
    return 0;  /* TODO: Real error value */ \
 \
  for (h=0; h<height; ++h) { \
    for (w=0; w<width; ++w) { \
re_format: \
      rv = snprintf(buffer + current_size, \
                    buffer_size - current_size, \
                    F, \
                    !w ? prefix : "", \
                    values[(width * h) + w], \
                    w != (width-1) ? ' ' : '\n'); \
      if (((size_t)rv) >= (buffer_size - current_size)) { \
        buffer_size *= 2; \
        if (!(buffer = (char *)_ogg_realloc(buffer, buffer_size))) { \
          return OD_EFAULT;  /* Out of memory */ \
        } \
        goto re_format; \
      } \
      current_size += rv; \
    } \
  } \
 \
  od_log(facility, level, "%s", buffer); \
 \
  return 0; \
}

DEFINE_OD_LOG_MATRIX(char, char, "%s%d%c")
DEFINE_OD_LOG_MATRIX(unsigned char, uchar, "%s%u%c")
DEFINE_OD_LOG_MATRIX(ogg_int16_t, int16, "%s%d%c")
DEFINE_OD_LOG_MATRIX(ogg_uint16_t, uint16, "%s%u%c")
DEFINE_OD_LOG_MATRIX(ogg_int32_t, int32, "%s%d%c")
DEFINE_OD_LOG_MATRIX(ogg_uint32_t, uint32, "%s%u%c")
DEFINE_OD_LOG_MATRIX(float, float, "%s%g%c")
