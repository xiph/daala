# NOTE: This Makefile requires GNU make
# Location to put the targets.
TARGETBINDIR = .
TARGETLIBDIR = .
# Name of the targets
LIBDAALABASE_TARGET = libdaalabase.a
LIBDAALAENC_TARGET = libdaalaenc.a
LIBDAALADEC_TARGET = libdaaladec.a
DUMP_VIDEO_TARGET = dump_video
PLAYER_EXAMPLE_TARGET = player_example
ENCODER_EXAMPLE_TARGET = encoder_example
# The compiler to use
CC = gcc
# The command to use to generate dependency information
MAKEDEPEND = gcc -MM
#MAKEDEPEND = makedepend -f- -Y --
# The location of include files.
# Modify these to point to your Ogg and Vorbis include directories if they are
#  not installed in a standard location.
CINCLUDE = `sdl-config --cflags`
# These are gcc-only, but not actually critical.
CFLAGS = -Wall -Wno-parentheses -Wno-overlength-strings
CFLAGS += -fvisibility=hidden
# Extra compilation flags.
# You may get speed increases by including flags such as -O2 or -O3 or
#  -ffast-math, or additional flags, depending on your system and compiler.
# The -g flag will generally include debugging information.
CFLAGS += -g
#CFLAGS += -DOD_CHECKASM
#CFLAGS += -DOD_DUMP_IMAGES
#CFLAGS += -DOD_ANIMATE
#CFLAGS += -flto
#CFLAGS += -floop-interchange -floop-strip-mine -floop-block
#CFLAGS += -ftree-vectorize -ftree-vectorizer-verbose=0
#CFLAGS += -fomit-frame-pointer
CFLAGS += -fPIC
CFLAGS += -std=c89 -pedantic
# Set machine-specific flags
HOST := $(shell uname -m)
ifneq (,$(findstring 86,${HOST}))
CFLAGS += -DOD_X86ASM
#CFLAGS += -march=native
#CFLAGS += -O3 -funroll-loops -fforce-addr -msse3 -Dinline=
CFLAGS += -O0 -Dinline=
else ifneq (,$(findstring arm,${HOST}))
CFLAGS += -DOD_ARMASM
#CFLAGS += -march=armv7
endif

# Libraries to link with, and the location of library files.
# Add -lpng -lz if you want to use -DOD_DUMP_IMAGES.
LIBS = -logg -lvorbis -lvorbisenc `sdl-config --libs` -lm
ifeq ($(findstring -DOD_DUMP_IMAGES,${CFLAGS}),-DOD_DUMP_IMAGES)
    LIBS += -lpng -lz
endif

# ANYTHING BELOW THIS LINE PROBABLY DOES NOT NEED EDITING
CINCLUDE := -I../include ${CINCLUDE}
LIBSRCDIR = ../src
BINSRCDIR = ../examples
WORKDIR = objs

# C source file lists
LIBDAALABASE_CSOURCES = \
newdct.c \
entcode.c \
entdec.c \
entenc.c \
icdf_table.c \
pvq_encoder.c \
filter.c \
info.c \
internal.c \
intra.c \
intradata.c \
mc.c \
pvq.c \
state.c \
$(if $(findstring -DOD_X86ASM,${CFLAGS}), \
x86/cpu.c \
x86/sse2mc.c \
x86/x86state.c \
)

LIBDAALABASE_CHEADERS = \
internal.h \
mc.h \
entcode.h \
entdec.h \
entenc.h \
mfrngcod.h \
ecintrin.h \
odintrin.h \
pvq.h \
state.h \
../include/daala/codec.h \

LIBDAALADEC_CSOURCES = \
decode.c \

LIBDAALADEC_CHEADERS =   \
${LIBDAALABASE_CHEADERS} \
decint.h \
../include/daala/daaladec.h \

LIBDAALAENC_CSOURCES = \
encode.c \
infoenc.c \
mcenc.c \

LIBDAALAENC_CHEADERS = \
${LIBDAALABASE_CHEADERS} \
encint.h \
../include/daala/daalaenc.h \

DUMP_VIDEO_CSOURCES = dump_video.c
ENCODER_EXAMPLE_CSOURCES = encoder_example.c
PLAYER_EXAMPLE_CSOURCES = player_example.c

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
ALL_OBJS:= ${LIBDAALABASE_OBJS} ${LIBDAALADEC_OBJS} ${LIBDAALAENC_OBJS} \
 ${DUMP_VIDEO_OBJS} ${ENCODER_EXAMPLE_OBJS} ${PLAYER_EXAMPLE_OBJS}
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
ALL_CSOURCES:= ${LIBDAALABASE_CSOURCES} ${LIBDAALADEC_CSOURCES} \
 ${LIBDAALAENC_CSOURCES} ${DUMP_VIDEO_CSOURCES} ${ENCODER_EXAMPLE_CSOURCES} \
 ${PLAYER_EXAMPLE_CSOURCES}
# Prepand target path to file names.
LIBDAALABASE_TARGET:= ${TARGETLIBDIR}/${LIBDAALABASE_TARGET}
LIBDAALADEC_TARGET:= ${TARGETLIBDIR}/${LIBDAALADEC_TARGET}
LIBDAALAENC_TARGET:= ${TARGETLIBDIR}/${LIBDAALAENC_TARGET}
DUMP_VIDEO_TARGET:= ${TARGETBINDIR}/${DUMP_VIDEO_TARGET}
ENCODER_EXAMPLE_TARGET:= ${TARGETBINDIR}/${ENCODER_EXAMPLE_TARGET}
PLAYER_EXAMPLE_TARGET:= ${TARGETBINDIR}/${PLAYER_EXAMPLE_TARGET}
ALL_TARGETS:= ${LIBDAALABASE_TARGET} ${LIBDAALADEC_TARGET} \
 ${LIBDAALAENC_TARGET} ${DUMP_VIDEO_TARGET} ${ENCODER_EXAMPLE_TARGET} \
 ${PLAYER_EXAMPLE_TARGET}

# Targets:
# Everything (default)
all: ${ALL_TARGETS}

check: all

# libdaalabase
${LIBDAALABASE_TARGET}: ${LIBDAALABASE_OBJS}
	mkdir -p ${TARGETLIBDIR}
	ar cqs $@ ${LIBDAALABASE_OBJS}

# libdaaladec
${LIBDAALADEC_TARGET}: ${LIBDAALADEC_OBJS}
	mkdir -p ${TARGETLIBDIR}
	ar cqs $@ ${LIBDAALADEC_OBJS}

# libdaalaenc
${LIBDAALAENC_TARGET}: ${LIBDAALAENC_OBJS}
	mkdir -p ${TARGETLIBDIR}
	ar cqs $@ ${LIBDAALAENC_OBJS}

# dump_video
${DUMP_VIDEO_TARGET}: ${DUMP_VIDEO_OBJS} ${LIBDAALABASE_TARGET} \
                       ${LIBDAALADEC_TARGET}
	mkdir -p ${TARGETBINDIR}
	${CC} ${CFLAGS} -o $@ ${DUMP_VIDEO_OBJS} ${LIBS} \
         ${LIBDAALADEC_TARGET} ${LIBDAALABASE_TARGET}

# encoder_example
${ENCODER_EXAMPLE_TARGET}: ${ENCODER_EXAMPLE_OBJS} ${LIBDAALABASE_TARGET} \
                            ${LIBDAALAENC_TARGET}
	mkdir -p ${TARGETBINDIR}
	${CC} ${CFLAGS} -o $@ ${ENCODER_EXAMPLE_OBJS} ${LIBS} \
         ${LIBDAALAENC_TARGET} ${LIBDAALABASE_TARGET}

# player_example
${PLAYER_EXAMPLE_TARGET}: ${PLAYER_EXAMPLE_OBJS} ${LIBDAALABASE_TARGET} \
                           ${LIBDAALADEC_TARGET}
	mkdir -p ${TARGETBINDIR}
	${CC} ${CFLAGS} -o $@ ${PLAYER_EXAMPLE_OBJS} ${LIBS} \
         ${LIBDAALADEC_TARGET} ${LIBDAALABASE_TARGET}

# Assembly listing
ALL_ASM := ${ALL_OBJS:%.o=%.s}
asm: ${ALL_ASM}

# check that build is complete
check: all

# Remove all targets.
clean:
	${RM} ${ALL_ASM} ${ALL_OBJS} ${ALL_DEPS}
	${RM} ${ALL_TARGETS}
	-rmdir ${WORKDIR}/x86 ${WORKDIR}

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
	${CC} ${CINCLUDE} ${CFLAGS} -S -o $@ $<
${WORKDIR}/%.o: ${LIBSRCDIR}/%.c
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
