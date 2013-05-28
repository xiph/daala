noinst_PROGRAMS = ectest test_coef_coder logging_test

AM_CFLAGS = $(OGG_CFLAGS)

TESTS = ectest test_coef_coder logging_test

ectest_SOURCES = ectest.c
ectest_LDADD = -lm
test_coef_coder_SOURCES = test_coef_coder.c
test_coef_coder_LDADD = $(DAALA_ENC_LIBS) $(OGG_LIBS)
logging_test_SOURCES = logging_test.c
logging_test_LDADD = $(DAALA_ENC_LIBS) $(OGG_LIBS)
