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

#include "metrics.h"

#if defined(OD_METRICS)

# include <stdio.h>
# include <stdlib.h>
# include <unistd.h>

# include "state.h"
# include "logging.h"

#endif

void write_metrics(ogg_int64_t cur_time, ogg_int64_t *metrics) {
#if defined(OD_METRICS)
  FILE *metrics_file;
  int res;
  char *runid;
  char buf[1024];
  char *filename;
  runid = getenv("OD_RUNID");
  filename = "metrics.json";
  if (runid != NULL) {
    res = snprintf(buf, sizeof(buf), "metrics.%s.json", runid);
    if (res >= 0) {
      filename = buf;
    }
  }
  if (cur_time == 1) {
    metrics_file = fopen(filename, "w");
    fprintf(metrics_file, "[\n");
  }
  else {
    metrics_file = fopen(filename, "a");
    res = ftell(metrics_file);
    fclose(metrics_file);
    truncate(filename, res - 3);
    metrics_file = fopen(filename, "a");
    fprintf(metrics_file, ",\n");
  }
  fprintf(metrics_file,
   "{\"frame\":%" OD_I64FMT ", \"total_bits\":%" OD_I64FMT ", "
   "\"mv_bits\":%" OD_I64FMT ", \"bs_bits\":%" OD_I64FMT ", "
   "\"pvq_bits\":%" OD_I64FMT ", \"dc_bits\":%" OD_I64FMT ", "
   "\"intra_bits\":%" OD_I64FMT "}\n",
   (long long)cur_time, (long long)metrics[OD_METRIC_TOTAL],
   (long long)metrics[OD_METRIC_MV],
   (long long)metrics[OD_METRIC_BLOCK_SWITCHING],
   (long long)metrics[OD_METRIC_PVQ], (long long)metrics[OD_METRIC_DC],
   (long long)metrics[OD_METRIC_INTRA]);
  fprintf(metrics_file, "]\n");
  fclose(metrics_file);
#else
  (void)cur_time;
  (void)metrics;
#endif
}
