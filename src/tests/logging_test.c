#include "../logging.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int emitted = 0;
char tmp_buf[5000];
char bogus_fmt_string[5000];

int failed = 0;

static int od_logging_test_emit(od_log_facility facility,
                                od_log_level level,
                                const char *fmt, va_list ap) {
  int rv;

  (void)facility;
  (void)level;

  emitted = 1;
  rv = vsnprintf(tmp_buf, sizeof(tmp_buf), fmt, ap);
  fprintf(stderr, "Emitted: %s\n", tmp_buf);

  return rv;
}

static void reset_result() {
  emitted = 0;
  tmp_buf[0] = '\0';
}

static void expected_something() {
  if (emitted)
    return;

  fprintf(stderr, "ERROR: Log not called\n");
  failed = 1;
}

static void expected_nothing() {
  if (!emitted)
    return;

  fprintf(stderr, "ERROR: Log called with '%s'\n", tmp_buf);
  failed = 1;
}

static void expected_result(const char *expected) {
  if (!emitted) {
    fprintf(stderr, "ERROR: Log not called\n");
    failed = 1;
    return;
  }

  if (strcmp(expected, tmp_buf)) {
    fprintf(stderr, "ERROR: Mismatch: '%s' != '%s'\n", tmp_buf, expected);
    failed = 1;
  }
}

int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

  /* Test the basic functionality. */
  setenv("OD_LOG_MODULES", "generic:3", 1);
  od_log_init(od_logging_test_emit);

  /* This should log. */
  reset_result();
  OD_LOG((OD_LOG_GENERIC, OD_LOG_ERR, "Blah blah %s:%d", "XXX", 9));
  expected_result("Blah blah XXX:9");

  /* This should not log (level too low) */
  reset_result();
  OD_LOG((OD_LOG_GENERIC, OD_LOG_DEBUG, "Blah blah %s:%d", "XXX", 9));
  expected_nothing();

  /* This should not log (facility not set) */
  reset_result();
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_ERR, "Blah blah %s:%d", "XXX", 9));
  expected_nothing();

  /* Test multiple modules */
  setenv("OD_LOG_MODULES", "generic:3,entropy-coder:5", 1);
  od_log_init(od_logging_test_emit);
  reset_result();
  OD_LOG((OD_LOG_GENERIC, OD_LOG_DEBUG, "Blah blah %s:%d", "XXX", 9));
  expected_nothing();

  reset_result();
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG, "Blah blah %s:%d", "XXX", 9));
  expected_something();

  /* Test multiple modules in the other order*/
  setenv("OD_LOG_MODULES", "entropy-coder:5,generic:3", 1);
  od_log_init(od_logging_test_emit);
  reset_result();
  OD_LOG((OD_LOG_GENERIC, OD_LOG_DEBUG, "Blah blah %s:%d", "XXX", 9));
  expected_nothing();

  reset_result();
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG, "Blah blah %s:%d", "XXX", 9));
  expected_something();

  /* Test bogus module string */
  setenv("OD_LOG_MODULES", "generic:XXX,blahblah:9,entropy-coder:5", 1);
  od_log_init(od_logging_test_emit);
  reset_result();
  OD_LOG((OD_LOG_GENERIC, OD_LOG_ERR, "Blah blah %s:%d", "XXX", 9));
  expected_nothing();

  reset_result();
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG, "Blah blah %s:%d", "XXX", 9));
  expected_something();

  /* Test a ridiculous fmt string */
  memset(bogus_fmt_string, 'X', sizeof(bogus_fmt_string));
  bogus_fmt_string[sizeof(bogus_fmt_string) - 2] = 'Y';
  bogus_fmt_string[sizeof(bogus_fmt_string) - 1] = '\0';
  OD_LOG((OD_LOG_ENTROPY_CODER, OD_LOG_DEBUG, bogus_fmt_string, "XXX", 9));

  if (failed)
    exit(1);

  exit(0);
}
