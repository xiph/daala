lib_LTLIBRARIES = libdaalabase.la libdaaladec.la libdaalaenc.la

noinst_HEADERS = \
	adapt.h \
	block_size.h \
	dct.h \
	decint.h \
	encint.h \
	entcode.h \
	entdec.h \
	entenc.h \
	filter.h \
	generic_code.h \
	internal.h \
	intra.h \
        logging.h \
	mc.h \
	odintrin.h \
	pvq.h \
	pvq_code.h \
	state.h \
	tf.h

libdaalabase_la_CFLAGS = $(OGG_CFLAGS)
libdaalabase_la_LIBADD = $(OGG_LIBS)
if DUMP_IMAGES
  libdaalabase_la_LIBADD += $(PNG_LIBS)
endif
libdaalabase_la_LDFLAGS = -no-undefined \
 -version-info @OP_LT_CURRENT@:@OP_LT_REVISION@:@OP_LT_AGE@
libdaalabase_la_SOURCES = \
	adapt.c \
	block_size.c \
	entcode.c \
	entdec.c \
	entenc.c \
	filter.c \
	generic_code.c \
	generic_decoder.c \
	generic_encoder.c \
	info.c \
	internal.c \
	intra.c \
	intradata.c \
	laplace_tables.c \
	logging.c \
	mc.c \
	newdct.c \
	pvq.c \
	pvq_decoder.c \
	pvq_encoder.c \
	state.c \
	tf.c
if ENABLE_X86ASM
libdaalabase_la_SOURCES += \
	x86/cpu.c \
	x86/cpu.h \
	x86/x86int.h \
	x86/sse2mc.c \
	x86/x86state.c
endif

libdaaladec_la_CFLAGS = $(OGG_CFLAGS)
libdaaladec_la_LIBADD = libdaalabase.la $(OGG_LIBS)
if DUMP_IMAGES
  libdaaladec_la_LIBADD += $(PNG_LIBS)
endif
libdaaladec_la_LDFLAGS = -no-undefined \
 -version-info @OP_LT_CURRENT@:@OP_LT_REVISION@:@OP_LT_AGE@
libdaaladec_la_SOURCES = \
	decode.c

libdaalaenc_la_CFLAGS = $(OGG_CFLAGS)
libdaalaenc_la_LIBADD = libdaalabase.la $(OGG_LIBS)
if DUMP_IMAGES
  libdaalaenc_la_LIBADD += $(PNG_LIBS)
endif
libdaalaenc_la_LDFLAGS = -no-undefined \
 -version-info @OP_LT_CURRENT@:@OP_LT_REVISION@:@OP_LT_AGE@
libdaalaenc_la_SOURCES = \
	encode.c \
	decode.c \
	infoenc.c \
	mcenc.c
