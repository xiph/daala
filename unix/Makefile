#Daala video codec
#Copyright (c) 2002-2013 Daala project contributors.  All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#- Redistributions of source code must retain the above copyright notice, this
#  list of conditions and the following disclaimer.
#
#- Redistributions in binary form must reproduce the above copyright notice,
#  this list of conditions and the following disclaimer in the documentation
#  and/or other materials provided with the distribution.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# NOTE: This Makefile requires GNU make
# Location to put the targets.
TARGETBINDIR = .
TESTBINDIR = tests
TARGETLIBDIR = .
# Name of the targets
LIBDAALABASE_TARGET = libdaalabase.a
LIBDAALAENC_TARGET = libdaalaenc.a
LIBDAALADEC_TARGET = libdaaladec.a
DUMP_VIDEO_TARGET = dump_video
ENCODER_EXAMPLE_TARGET = encoder_example
PLAYER_EXAMPLE_TARGET = player_example
# Test targets
DCTTEST_TARGET = dcttest
ECTEST_TARGET = ectest
TEST_CHECK_INITIAL_TARGET = check_initial
TEST_COEF_CODER_TARGET = test_coef_coder
TEST_HEADER_TARGET = check_tests
TEST_LOGGING_TARGET = logging_test
TEST_DIVU_SMALL_TARGET = test_divu_small
TEST_FILTERS_TARGET = test_filters

# The command to use to generate dependency information
MAKEDEPEND = $(CC) -MM
#MAKEDEPEND = makedepend -f- -Y --
# The location of include files.
# Modify these to point to your Ogg and Vorbis include directories if they are
#  not installed in a standard location.
CINCLUDE = `sdl-config --cflags` `pkg-config ogg --cflags`
# These are gcc-only, but not actually critical.
# Extra compilation flags.
# You may get speed increases by including flags such as -O2 or -O3 or
#  -ffast-math, or additional flags, depending on your system and compiler.
# The -g flag will generally include debugging information.
CFLAGS := -g $(CFLAGS)
CFLAGS := -DOD_LOGGING_ENABLED $(CFLAGS)
CFLAGS := -DOD_ENABLE_ASSERTIONS $(CFLAGS)
#CFLAGS := -DOD_CHECKASM $(CFLAGS)
#CFLAGS := -DOD_DUMP_IMAGES $(CFLAGS)
#CFLAGS := -DOD_ENCODER_CHECK $(CFLAGS)
#CFLAGS := -DOD_ANIMATE $(CFLAGS)
#CFLAGS := -DOD_LOGGING_ENABLED $(CFLAGS)
CFLAGS := -DOD_ACCOUNTING $(CFLAGS)
CFLAGS := -fPIC $(CFLAGS)
CFLAGS := -std=c89 -pedantic $(CFLAGS)
CFLAGS := -fvisibility=hidden $(CFLAGS)
CFLAGS := -Wextra -Wno-parentheses -Wno-long-long -Wno-overlength-strings $(CFLAGS)
CFLAGS := -Wall $(CFLAGS)
CFLAGS := -I/opt/local/include $(CFLAGS)
# Set machine-specific flags
HOST := $(shell uname -m)
ifneq (,$(findstring 86,${HOST}))
CFLAGS := -DOD_X86ASM $(CFLAGS)
else ifneq (,$(findstring arm,${HOST}))
CFLAGS := -DOD_ARMASM $(CFLAGS)
#CFLAGS := -march=armv7 $(CFLAGS)
endif

# Libraries to link with, and the location of library files.
# Add -lpng -lz if you want to use -DOD_DUMP_IMAGES.
LIBS = `pkg-config ogg --libs` `sdl-config --libs` -lm
ifeq ($(findstring -DOD_DUMP_IMAGES,${CFLAGS}),-DOD_DUMP_IMAGES)
    LIBS += -lpng -lz
endif

LIBDAALAENC_LIBS = ${LIBDAALAENC_TARGET}
ifeq ($(findstring -DOD_ENCODER_CHECK,${CFLAGS}),-DOD_ENCODER_CHECK)
    LIBDAALAENC_LIBS += ${LIBDAALADEC_TARGET}
endif


CHECK_LIBS := -pthread $(shell pkg-config check --libs)

DCTTEST_LIBS = -lm
ECTEST_LIBS = -lm
TEST_COEF_CODER_LIBS =
TEST_HEADER_LIBS = ${CHECK_LIBS}
TEST_LOGGING_LIBS =
TEST_CHECK_INITIAL_LIBS = ${CHECK_LIBS}
TEST_DIVU_SMALL_LIBS =
TEST_FILTER_LIBS =

# ANYTHING BELOW THIS LINE PROBABLY DOES NOT NEED EDITING
CINCLUDE := -I../include ${CINCLUDE}
LIBSRCDIR = ../src
BINSRCDIR = ../examples
TESTSRCDIR = ${LIBSRCDIR}
WORKDIR = objs

# C source file lists
LIBDAALABASE_CSOURCES = \
adapt.c \
block_size.c \
entcode.c \
entdec.c \
filter.c \
generic_code.c \
generic_decoder.c \
info.c \
internal.c \
intra.c \
intradata.c \
laplace_decoder.c \
laplace_tables.c \
logging.c \
mc.c \
dct.c \
partition.c \
pvq.c \
pvq_code.c \
pvq_decoder.c \
state.c \
switch_table.c \
tf.c \
zigzag4.c \
zigzag8.c \
zigzag16.c \
$(if $(findstring -DOD_X86ASM,${CFLAGS}), \
x86/cpu.c \
x86/sse2mc.c \
x86/x86state.c \
) \
$(if $(findstring -DOD_LOGGING_ENABLED,${CFLAGS}), \
logging.c \
)

LIBDAALABASE_CHEADERS = \
adapt.h \
block_size.h \
entcode.h \
entdec.h \
filter.h \
generic_code.h \
laplace_code.h \
partition.h \
pvq_code.h \
internal.h \
intra.h \
logging.h \
mc.h \
mcenc.h \
odintrin.h \
pvq.h \
state.h \
tf.h \
../include/daala/codec.h \

LIBDAALADEC_CSOURCES = \
block_size_dec.c \
decode.c \
infodec.c \

LIBDAALADEC_CHEADERS =   \
${LIBDAALABASE_CHEADERS} \
decint.h \
../include/daala/daaladec.h \

LIBDAALAENC_CSOURCES = \
block_size_dec.c \
block_size_enc.c \
encode.c \
entenc.c \
generic_encoder.c \
infoenc.c \
laplace_encoder.c \
mcenc.c \
$(if $(findstring -DOD_ACCOUNTING,${CFLAGS}), \
accounting.c \
) \
pvq_encoder.c \
$(if $(findstring -DOD_X86ASM,${CFLAGS}), \
x86/x86enc.c \
x86/x86mcenc.c \
) \

LIBDAALAENC_CHEADERS = \
${LIBDAALABASE_CHEADERS} \
encint.h \
entenc.h \
laplace_encoder.h \
../include/daala/daalaenc.h \

DUMP_VIDEO_CSOURCES = dump_video.c
ENCODER_EXAMPLE_CSOURCES = encoder_example.c
PLAYER_EXAMPLE_CSOURCES = player_example.c

DCTTEST_CSOURCES = \
dct.c \
filter.c \
$(if $(findstring -DOD_ENABLE_ASSERTIONS,${CFLAGS}), \
internal.c \
)
ECTEST_CSOURCES = tests/ectest.c
TEST_COEF_CODER_CSOURCES = \
tests/test_coef_coder.c \
$(if $(findstring -DOD_LOGGING_ENABLED,${CFLAGS}), \
logging.c \
)

TEST_CHECK_INITIAL_CSOURCES = tests/check_initial.c
TEST_HEADER_CSOURCES=tests/check_main.c tests/headerencode_test.c
TEST_LOGGING_CSOURCES=tests/logging_test.c
TEST_DIVU_SMALL_CSOURCES=tests/test_divu_small.c
TEST_FILTERS_CSOURCES=tests/test_filters.c

# Create object file list.
LIBDAALABASE_OBJS:= ${LIBDAALABASE_CSOURCES:%.c=${WORKDIR}/%.o}
LIBDAALABASE_ASMS:= ${LIBDAALABASE_OBJS:%.o=%.s}
LIBDAALABASE_DEPS:= ${LIBDAALABASE_OBJS:%.o=%.d}
LIBDAALADEC_OBJS:= ${LIBDAALADEC_CSOURCES:%.c=${WORKDIR}/%.o}
LIBDAALADEC_ASMS:= ${LIBDAALADEC_OBJS:%.o=%.s}
LIBDAALADEC_DEPS:= ${LIBDAALADEC_OBJS:%.o=%.d}
LIBDAALAENC_OBJS:= ${LIBDAALAENC_CSOURCES:%.c=${WORKDIR}/%.o}
LIBDAALAENC_ASMS:= ${LIBDAALAENC_OBJS:%.o=%.s}
LIBDAALAENC_DEPS:= ${LIBDAALAENC_OBJS:%.o=%.d}
DUMP_VIDEO_OBJS:= ${DUMP_VIDEO_CSOURCES:%.c=${WORKDIR}/%.o}
ENCODER_EXAMPLE_OBJS:= ${ENCODER_EXAMPLE_CSOURCES:%.c=${WORKDIR}/%.o}
PLAYER_EXAMPLE_OBJS:= ${PLAYER_EXAMPLE_CSOURCES:%.c=${WORKDIR}/%.o}
ECTEST_OBJS:= ${ECTEST_CSOURCES:%.c=${WORKDIR}/%.o}
TEST_CHECK_INITIAL_OBJS:= ${TEST_CHECK_INITIAL_CSOURCES:%.c=${WORKDIR}/%.o}
TEST_COEF_CODER_OBJS:= ${TEST_COEF_CODER_CSOURCES:%.c=${WORKDIR}/%.o}
TEST_HEADER_OBJS:= ${TEST_HEADER_CSOURCES:%.c=${WORKDIR}/%.o}
TEST_LOGGING_OBJS:= ${TEST_LOGGING_CSOURCES:%.c=${WORKDIR}/%.o}
TEST_DIVU_SMALL_OBJS:= ${TEST_DIVU_SMALL_CSOURCES:%.c=${WORKDIR}/%.o}
TEST_FILTERS_OBJS:= ${TEST_FILTERS_CSOURCES:%.c=${WORKDIR}/%.o}
ALL_OBJS:= ${LIBDAALABASE_OBJS} ${LIBDAALADEC_OBJS} ${LIBDAALAENC_OBJS} \
 ${DUMP_VIDEO_OBJS} ${ENCODER_EXAMPLE_OBJS} ${PLAYER_EXAMPLE_OBJS} \
 ${ECTEST_OBJS} ${TEST_CHECK_INITIAL_OBJS} ${TEST_COEF_CODER_OBJS} \
 ${TEST_HEADER_OBJS} ${TEST_LOGGING_OBJS} ${TEST_DIVU_SMALL_OBJS} \
 ${TEST_FILTERS_OBJS}
# Create the dependency file list
ALL_DEPS:= ${ALL_OBJS:%.o=%.d}
# Prepend source path to file names.
LIBDAALABASE_CSOURCES:= ${LIBDAALABASE_CSOURCES:%=${LIBSRCDIR}/%}
LIBDAALABASE_CHEADERS:= ${LIBDAALABASE_CHEADERS:%=${LIBSRCDIR}/%}
LIBDAALADEC_CSOURCES:= ${LIBDAALADEC_CSOURCES:%=${LIBSRCDIR}/%}
LIBDAALADEC_CHEADERS:= ${LIBDAALADEC_CHEADERS:%=${LIBSRCDIR}/%}
LIBDAALAENC_CSOURCES:= ${LIBDAALAENC_CSOURCES:%=${LIBSRCDIR}/%}
LIBDAALAENC_CHEADERS:= ${LIBDAALAENC_CHEADERS:%=${LIBSRCDIR}/%}
DUMP_VIDEO_CSOURCES:= ${DUMP_VIDEO_CSOURCES:%=${BINSRCDIR}/%}
ENCODER_EXAMPLE_CSOURCES:= ${ENCODER_EXAMPLE_CSOURCES:%=${BINSRCDIR}/%}
PLAYER_EXAMPLE_CSOURCES:= ${PLAYER_EXAMPLE_CSOURCES:%=${BINSRCDIR}/%}
DCTTEST_CSOURCES:= ${DCTTEST_CSOURCES:%=${LIBSRCDIR}/%}
ECTEST_CSOURCES:= ${ECTEST_CSOURCES:%=${TESTSRCDIR}/%}
TEST_CHECK_INITIAL_CSOURCES:= ${TEST_CHECK_INITIAL_CSOURCES:%=${TESTSRCDIR}/%}
TEST_COEF_CODER_CSOURCES:= ${TEST_COEF_CODER_CSOURCES:%=${TESTSRCDIR}/%}
ALL_CSOURCES:= ${LIBDAALABASE_CSOURCES} ${LIBDAALADEC_CSOURCES} \
 ${LIBDAALAENC_CSOURCES} ${DUMP_VIDEO_CSOURCES} ${ENCODER_EXAMPLE_CSOURCES} \
 ${PLAYER_EXAMPLE_CSOURCES} ${DCTTEST_CSOURCES} ${ECTEST_CSOURCES} \
 ${TEST_COEF_CODER_CSOURCES} ${TEST_CHECK_INITIAL_CSOURCES}
# Prepand target path to file names.
LIBDAALABASE_TARGET:= ${TARGETLIBDIR}/${LIBDAALABASE_TARGET}
LIBDAALADEC_TARGET:= ${TARGETLIBDIR}/${LIBDAALADEC_TARGET}
LIBDAALAENC_TARGET:= ${TARGETLIBDIR}/${LIBDAALAENC_TARGET}
DUMP_VIDEO_TARGET:= ${TARGETBINDIR}/${DUMP_VIDEO_TARGET}
ENCODER_EXAMPLE_TARGET:= ${TARGETBINDIR}/${ENCODER_EXAMPLE_TARGET}
PLAYER_EXAMPLE_TARGET:= ${TARGETBINDIR}/${PLAYER_EXAMPLE_TARGET}
# Prepend test path to file names.
DCTTEST_TARGET:= ${TESTBINDIR}/${DCTTEST_TARGET}
ECTEST_TARGET:= ${TESTBINDIR}/${ECTEST_TARGET}
TEST_CHECK_INITIAL_TARGET:= ${TESTBINDIR}/${TEST_CHECK_INITIAL_TARGET}
TEST_COEF_CODER_TARGET:= ${TESTBINDIR}/${TEST_COEF_CODER_TARGET}
TEST_HEADER_TARGET:= ${TESTBINDIR}/${TEST_HEADER_TARGET}
TEST_LOGGING_TARGET:= ${TESTBINDIR}/${TEST_LOGGING_TARGET}
TEST_DIVU_SMALL_TARGET:=${TESTBINDIR}/${TEST_DIVU_SMALL_TARGET}
TEST_FILTERS_TARGET:=${TESTBINDIR}/${TEST_FILTERS_TARGET}

# Complete set of targets
ALL_TARGETS:= ${LIBDAALABASE_TARGET} ${LIBDAALADEC_TARGET} \
 ${LIBDAALAENC_TARGET} ${DUMP_VIDEO_TARGET} ${ENCODER_EXAMPLE_TARGET} \
 ${PLAYER_EXAMPLE_TARGET} ${DCTTEST_TARGET} ${ECTEST_TARGET} \
 ${TEST_COEF_CODER_TARGET} ${TEST_HEADER_TARGET} ${TEST_LOGGING_TARGET} \
 ${TEST_CHECK_INITIAL_TARGET} ${TEST_DIVU_SMALL_TARGET} ${TEST_FILTERS_TARGET}

# Targets:
# Everything (default)
all: ${ALL_TARGETS}

# libdaalabase
${LIBDAALABASE_TARGET}: ${LIBDAALABASE_OBJS}
	mkdir -p ${TARGETLIBDIR}
	ar crs $@ ${LIBDAALABASE_OBJS}

# libdaaladec
${LIBDAALADEC_TARGET}: ${LIBDAALADEC_OBJS}
	mkdir -p ${TARGETLIBDIR}
	ar crs $@ ${LIBDAALADEC_OBJS}

# libdaalaenc
${LIBDAALAENC_TARGET}: ${LIBDAALAENC_OBJS}
	mkdir -p ${TARGETLIBDIR}
	ar crs $@ ${LIBDAALAENC_OBJS}

# dump_video
${DUMP_VIDEO_TARGET}: ${DUMP_VIDEO_OBJS} ${LIBDAALABASE_TARGET} \
                       ${LIBDAALADEC_TARGET}
	mkdir -p ${TARGETBINDIR}
	${CC} ${CFLAGS} -o $@ ${DUMP_VIDEO_OBJS} \
         ${LIBDAALADEC_TARGET} ${LIBDAALABASE_TARGET} ${LIBS}

# encoder_example
${ENCODER_EXAMPLE_TARGET}: ${ENCODER_EXAMPLE_OBJS} ${LIBDAALABASE_TARGET} \
                            ${LIBDAALAENC_TARGET}
	mkdir -p ${TARGETBINDIR}
	${CC} ${CFLAGS} -o $@ ${ENCODER_EXAMPLE_OBJS} \
        ${LIBDAALAENC_LIBS} ${LIBDAALABASE_TARGET} ${LIBS}

# player_example
${PLAYER_EXAMPLE_TARGET}: ${PLAYER_EXAMPLE_OBJS} ${LIBDAALABASE_TARGET} \
                           ${LIBDAALADEC_TARGET}
	mkdir -p ${TARGETBINDIR}
	${CC} ${CFLAGS} -o $@ ${PLAYER_EXAMPLE_OBJS} \
         ${LIBDAALADEC_TARGET} ${LIBDAALABASE_TARGET} ${LIBS}

# check_initial
${TEST_CHECK_INITIAL_TARGET}: ${TEST_CHECK_INITIAL_OBJS}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${TEST_CHECK_INITIAL_OBJS} ${TEST_CHECK_INITIAL_LIBS} -o $@

# dcttest
${DCTTEST_TARGET}: ${DCTTEST_CSOURCES}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} -DOD_DCT_CHECK_OVERFLOW -DOD_DCT_TEST ${DCTTEST_CSOURCES} \
         ${DCTTEST_LIBS} -o $@

# ectest
${ECTEST_TARGET}: ${ECTEST_OBJS}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${ECTEST_OBJS} ${ECTEST_LIBS} -o $@

# test_coef_coder
${TEST_COEF_CODER_TARGET}: ${TEST_COEF_CODER_OBJS} ${LIBDAALABASE_TARGET}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${TEST_COEF_CODER_OBJS} ${TEST_COEF_CODER_LIBS} -o $@ \
	  ${LIBDAALABASE_TARGET} -lm

# test_header
${TEST_HEADER_TARGET}: ${TEST_HEADER_OBJS} ${LIBDAALADEC_TARGET} \
 ${LIBDAALAENC_TARGET} ${LIBDAALABASE_TARGET}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${TEST_HEADER_OBJS} ${TEST_HEADER_LIBS} -o $@ \
	  ${LIBDAALAENC_TARGET} ${LIBDAALADEC_TARGET} ${LIBDAALABASE_TARGET} \
	  -lm

# logging_test
${TEST_LOGGING_TARGET}: ${TEST_LOGGING_OBJS} ${LIBDAALABASE_TARGET}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${TEST_LOGGING_OBJS} ${TEST_LOGGING_LIBS} -o $@ \
	  ${LIBDAALABASE_TARGET} -lm

# test_divu_small
${TEST_DIVU_SMALL_TARGET}: ${TEST_DIVU_SMALL_OBJS} ${LIBDAALABASE_TARGET}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${TEST_DIVU_SMALL_OBJS} ${TEST_DIVU_SMALL_LIBS} -o $@ \
	  ${LIBDAALABASE_TARGET} -lm

# test_filters
${TEST_FILTERS_TARGET}: ${TEST_FILTERS_OBJS} ${LIBDAALABASE_TARGET}
	mkdir -p ${TESTBINDIR}
	${CC} ${CFLAGS} ${TEST_FILTERS_OBJS} ${TEST_FILTERS_LIBS} -o $@ \
	  ${LIBDAALABASE_TARGET} -lm

# Assembly listing
ALL_ASM := ${ALL_OBJS:%.o=%.s}
asm: ${ALL_ASM}

# check that build is complete
check: all
	${DCTTEST_TARGET}
	${ECTEST_TARGET}
	${TEST_CHECK_INITIAL_TARGET}
	${TEST_COEF_CODER_TARGET}
	${TEST_HEADER_TARGET}
	${TEST_LOGGING_TARGET}
	${TEST_DIVU_SMALL_TARGET}
	${TEST_FILTERS_TARGET}

# Remove all targets.
clean:
	${RM} ${ALL_ASM} ${ALL_OBJS} ${ALL_DEPS}
	${RM} ${ALL_TARGETS}
	-rmdir ${TESTBINDIR} ${WORKDIR}/tests ${WORKDIR}/x86 ${WORKDIR}

# Make everything depend on changes in the Makefile
${ALL_ASM} ${ALL_OBJS} ${ALL_DEPS} ${ALL_TARGETS} : Makefile

# Specify which targets are phony for GNU make
.PHONY : all clean check

# Rules
# These are set up to compile each library as a single unit, to allow gcc to do
#  Inter Module Analysis (IMA) when optimizing.
${WORKDIR}/%.d: ${LIBSRCDIR}/%.c
	mkdir -p ${dir $@}
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.o} > $@
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.s} >> $@
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT $@ >> $@
${WORKDIR}/%.s: ${LIBSRCDIR}/%.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -S -o $@ $<
${WORKDIR}/%.o: ${LIBSRCDIR}/%.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -c -o $@ $<

${WORKDIR}/%.d : ${BINSRCDIR}/%.c
	mkdir -p ${dir $@}
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.o} > $@
${WORKDIR}/%.s : ${BINSRCDIR}/%.c ${WORKDIR}/%.o
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -S -o $@ $<
${WORKDIR}/%.o : ${BINSRCDIR}/%.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -c -o $@ $<

# Include header file dependencies
include ${ALL_DEPS}
