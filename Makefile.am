# Libraries

# Provide the full test output for failed tests when using the parallel
# test suite (which is enabled by default with automake 1.13+).
export VERBOSE = yes

lib_LTLIBRARIES = src/libdaalabase.la src/libdaaladec.la src/libdaalaenc.la

AM_CPPFLAGS = -I$(srcdir)/include

noinst_HEADERS = \
	src/adapt.h \
	src/block_size.h \
	src/block_size_enc.h \
	src/dct.h \
	src/decint.h \
	src/encint.h \
	src/entcode.h \
	src/entdec.h \
	src/entenc.h \
	src/filter.h \
	src/generic_code.h \
	src/internal.h \
	src/intra.h \
	src/laplace_code.h \
	src/logging.h \
	src/mc.h \
	src/mcenc.h \
	src/odintrin.h \
	src/partition.h \
	src/pvq.h \
	src/pvq_code.h \
	src/quantizer.h \
	src/state.h \
	src/tf.h \
	src/zigzag.h \
	src/accounting.h \
	src/x86/cpu.h \
	src/x86/x86enc.h \
	src/x86/x86dct.h \
	src/x86/x86int.h \
	src/arm/cpu.h \
	src/arm/armint.h

src_dct_SOURCES = src/dct.c src/internal.c
if ENABLE_X86ASM
src_dct_SOURCES += src/x86/cpu.c
if ENABLE_SSE2_INTRINSICS
src_dct_SOURCES += src/x86/sse2dct.c
%sse2dct.o %sse2dct.lo: CFLAGS += -msse2
endif
if ENABLE_SSE41_INTRINSICS
src_dct_SOURCES += src/x86/sse41dct.c
%sse41dct.o %sse41dct.lo: CFLAGS += -msse4.1
endif
if ENABLE_AVX2_INTRINSICS
src_dct_SOURCES += src/x86/avx2dct.c
%avx2dct.o %avx2dct.lo: CFLAGS += -mavx2
endif
endif
if ENABLE_ARMASM
src_dct_SOURCES += src/arm/cpu.c
endif


src_libdaalabase_la_CFLAGS = $(OGG_CFLAGS)
src_libdaalabase_la_LIBADD = $(OGG_LIBS) $(LIBM)
if DUMP_IMAGES
  src_libdaalabase_la_LIBADD += $(PNG_LIBS)
endif
src_libdaalabase_la_LDFLAGS = -no-undefined \
 -version-info @OD_LT_CURRENT@:@OD_LT_REVISION@:@OD_LT_AGE@
src_libdaalabase_la_LDFLAGS += 
src_libdaalabase_la_SOURCES = \
	src/accounting.c \
	src/adapt.c \
	src/entcode.c \
	src/entdec.c \
	src/filter.c \
	src/generic_code.c \
	src/generic_decoder.c \
	src/info.c \
	src/intra.c \
	src/laplace_decoder.c \
	src/laplace_tables.c \
	src/logging.c \
	src/mc.c \
	src/partition.c \
	src/pvq.c \
	src/pvq_code.c \
	src/pvq_decoder.c \
	src/quantizer.c \
	src/state.c \
	src/switch_table.c \
	src/tf.c \
	src/zigzag4.c \
	src/zigzag8.c \
	src/zigzag16.c \
	src/zigzag32.c \
	src/zigzag64.c \
	$(src_dct_SOURCES)
if ENABLE_X86ASM
src_libdaalabase_la_SOURCES += \
	src/x86/sse2mc.c \
	src/x86/x86state.c
endif

src_libdaaladec_la_CFLAGS = $(OGG_CFLAGS)
src_libdaaladec_la_LIBADD = src/libdaalabase.la $(OGG_LIBS) $(LIBM)
src_libdaaladec_la_LDFLAGS = -no-undefined \
 -version-info @OD_LT_CURRENT@:@OD_LT_REVISION@:@OD_LT_AGE@
src_libdaaladec_la_SOURCES = \
	src/decode.c \
	src/infodec.c

src_libdaalaenc_la_CFLAGS = $(OGG_CFLAGS)
src_libdaalaenc_la_LIBADD = src/libdaalabase.la $(OGG_LIBS) $(LIBM)
if ENCODER_CHECK
  src_libdaalaenc_la_LIBADD += src/libdaaladec.la
endif
src_libdaalaenc_la_LDFLAGS = -no-undefined \
 -version-info @OD_LT_CURRENT@:@OD_LT_REVISION@:@OD_LT_AGE@
src_libdaalaenc_la_SOURCES = \
	src/block_size_enc.c \
	src/encode.c \
	src/entenc.c \
	src/generic_encoder.c \
	src/infoenc.c \
	src/laplace_encoder.c \
	src/mcenc.c \
	src/pvq_encoder.c
if ENABLE_X86ASM
src_libdaalaenc_la_SOURCES += \
        src/x86/x86enc.c \
        src/x86/x86mcenc.c
endif

# Example programs

noinst_PROGRAMS = \
	examples/dump_video \
	examples/encoder_example

if ENABLE_PLAYER_EXAMPLE
  noinst_PROGRAMS += examples/player_example
endif

if ENABLE_ANALYZER
  noinst_PROGRAMS += examples/analyzer
endif

examples_dump_video_SOURCES = examples/dump_video.c
examples_dump_video_CFLAGS = $(OGG_CFLAGS)
examples_dump_video_LDADD = src/libdaalabase.la src/libdaaladec.la  $(OGG_LIBS)

examples_encoder_example_SOURCES = examples/encoder_example.c
examples_encoder_example_CFLAGS = $(OGG_CFLAGS)
examples_encoder_example_LDADD = src/libdaalabase.la src/libdaalaenc.la $(OGG_LIBS) $(LIBM)

if ENABLE_PLAYER_EXAMPLE
examples_player_example_SOURCES = examples/player_example.c
examples_player_example_CFLAGS = $(OGG_CFLAGS) $(SDL_CFLAGS)
examples_player_example_LDADD = src/libdaalabase.la src/libdaaladec.la $(OGG_LIBS) $(SDL_LIBS)
endif

if ENABLE_ANALYZER
examples_analyzer_SOURCES = examples/analyzer.cc
examples_analyzer_CXXFLAGS = $(OGG_CFLAGS) $(WX_CPPFLAGS)
examples_analyzer_LDADD = src/libdaalabase.la src/libdaaladec.la $(OGG_LIBS) $(WX_LIBS)
endif

# Includes

daalaincludedir = $(includedir)/daala
daalainclude_HEADERS = \
	include/daala/codec.h \
	include/daala/daaladec.h \
	include/daala/daalaenc.h \
	include/daala/daala_integer.h

# Tools

CLEANFILES =
EXTRA_PROGRAMS =

if ENABLE_TOOLS

tools_TARGETS = \
	tools/png2y4m \
	tools/y4m2png \
	tools/dump_psnrhvs \
	tools/block_size_analysis \
	tools/to_monochrome \
	tools/downsample \
	tools/upsample \
	tools/divu_const \
	tools/trans \
	tools/trans2d \
	tools/trans_gain \
	tools/gen_cdf \
	tools/compute_basis \
	tools/gen_laplace_tables \
	tools/daalainfo \
	tools/dump_ssim \
	tools/dump_fastssim \
	tools/bjontegaard \
	tools/yuvjpeg \
	tools/jpegyuv \
	tools/yuv2yuv4mpeg \
	tools/dump_psnr \
	tools/vq_train \
	tools/draw_zigzags

noinst_HEADERS += \
	tools/cholesky.h \
	tools/int_search.h \
	tools/kiss99.h \
	tools/matidx.h \
	tools/od_defs.h \
	tools/od_filter.h \
	tools/pythag.h \
	tools/svd.h \
	tools/trans_tools.h \
	tools/vidinput.h

EXTRA_PROGRAMS += $(tools_TARGETS)
CLEANFILES += $(tools_TARGETS)

tools_daalainfo_SOURCES = tools/daalainfo.c
tools_daalainfo_CFLAGS = $(OGG_CFLAGS)
tools_daalainfo_LDADD = $(OGG_LIBS) src/libdaalabase.la src/libdaaladec.la

# png2y4m
tools_png2y4m_SOURCES = \
	tools/kiss99.c \
	tools/png2y4m.c
tools_png2y4m_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_png2y4m_LDADD = $(OGG_LIBS) $(PNG_LIBS)

# y4m2png
tools_y4m2png_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/y4m2png.c
tools_y4m2png_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_y4m2png_LDADD = $(OGG_LIBS) $(PNG_LIBS)

# dump_psnrhvs
tools_dump_psnrhvs_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	$(src_dct_SOURCES) \
	tools/dump_psnrhvs.c
tools_dump_psnrhvs_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_dump_psnrhvs_LDADD = $(OGG_LIBS) $(LIBM)

# dump_ssim
tools_dump_ssim_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/dump_ssim.c
tools_dump_ssim_CFLAGS = $(OGG_CFLAGS)
tools_dump_ssim_LDADD = $(OGG_LIBS) $(LIBM)

# dump_fastssim
tools_dump_fastssim_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/dump_fastssim.c
tools_dump_fastssim_CFLAGS = $(OGG_CFLAGS)
tools_dump_fastssim_LDADD = $(OGG_LIBS) $(LIBM)

# dump_psnr
tools_dump_psnr_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/dump_psnr.c
tools_dump_psnr_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_dump_psnr_LDADD = $(OGG_LIBS) $(LIBM)

# block_size_analysis
tools_block_size_analysis_SOURCES = \
	tools/block_size_analysis.c \
	src/block_size_enc.c \
	tools/vidinput.c \
	tools/y4m_input.c \
	$(src_dct_SOURCES) \
	src/generic_encoder.c \
	src/generic_code.c \
	src/laplace_encoder.c \
	src/logging.c \
	src/laplace_tables.c \
	src/entcode.c \
	src/entenc.c \
	src/filter.c \
	src/switch_table.c
tools_block_size_analysis_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_block_size_analysis_LDADD = $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

# to_monochrome
tools_to_monochrome_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/to_monochrome.c
tools_to_monochrome_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_to_monochrome_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS)

# downsample
tools_downsample_SOURCES = \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/downsample.c
tools_downsample_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_downsample_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS)

# upsample
tools_upsample_SOURCES = \
	$(src_dct_SOURCES) \
	src/generic_code.c \
	src/switch_table.c \
	src/logging.c \
	src/info.c \
	src/mc.c \
	src/pvq.c \
	src/state.c \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/upsample.c
if ENABLE_X86ASM
tools_upsample_SOURCES += \
	src/x86/sse2mc.c \
	src/x86/x86state.c
endif
tools_upsample_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_upsample_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

# divu_const
tools_divu_const_SOURCES = \
	tools/divu_const.c
tools_divu_const_CFLAGS = $(OGG_CFLAGS)
tools_divu_const_LDADD = $(OGG_LIBS) $(LIBM)

# trans
tools_trans_SOURCES = \
	tools/trans.c \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/od_filter.c \
	tools/trans_tools.c \
	tools/int_search.c \
	tools/trans_data.c \
	tools/kiss99.c \
	tools/svd.c \
	tools/cholesky.c \
	src/filter.c \
	$(src_dct_SOURCES)
tools_trans_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS) $(OPENMP_CFLAGS)
tools_trans_LDADD = $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

# trans_gain
tools_trans_gain_SOURCES = \
	tools/trans_gain.c \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/od_filter.c \
	tools/trans_tools.c \
	src/filter.c \
	$(src_dct_SOURCES)
tools_trans_gain_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS) $(OPENMP_CFLAGS)
tools_trans_gain_LDADD = $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

#trans2d
tools_trans2d_SOURCES = \
	tools/trans2d.c \
	tools/vidinput.c \
	tools/y4m_input.c \
	tools/od_filter.c \
	tools/trans_tools.c \
	tools/int_search.c \
	tools/trans_data.c \
	tools/kiss99.c \
	tools/svd.c \
	tools/cholesky.c \
	src/filter.c \
	$(src_dct_SOURCES)
tools_trans2d_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS) $(OPENMP_CFLAGS)
tools_trans2d_LDADD = $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

# gen_cdf
tools_gen_cdf_SOURCES = \
	tools/gen_cdf.c
tools_gen_cdf_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_gen_cdf_LDADD = $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

# compute_basis
tools_compute_basis_SOURCES = \
	tools/compute_basis.c \
	src/dct.c \
	src/filter.c \
	src/internal.c
tools_compute_basis_CFLAGS = $(OGG_CFLAGS)
tools_compute_basis_LDADD = $(OGG_LIBS) $(LIBM)

# gen_laplace_tables
tools_gen_laplace_tables_SOURCES = \
	tools/gen_laplace_tables.c
tools_gen_laplace_tables_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
tools_gen_laplace_tables_LDADD = $(OGG_LIBS) $(PNG_LIBS) $(LIBM)

# bjontegaard
tools_bjontegaard_SOURCES = \
	tools/bjontegaard.c \
	tools/cholesky.c \
	tools/qr.c \
	tools/svd.c
tools_bjontegaard_CFLAGS =
tools_bjontegaard_LDADD = $(LIBM)

# yuvjpeg
tools_yuvjpeg_SOURCES = \
	tools/yuvjpeg.c
tools_yuvjpeg_CFLAGS =
tools_yuvjpeg_LDADD = $(JPEG_LIBS)

# jpegyuv
tools_jpegyuv_SOURCES = \
	tools/jpegyuv.c
tools_jpegyuv_CFLAGS =
tools_jpegyuv_LDADD = $(JPEG_LIBS)

# yuv2yuv4mpeg
tools_yuv2yuv4mpeg_SOURCES = \
	tools/yuv2yuv4mpeg.c
tools_yuv2yuv4mpeg_CFLAGS =
tools_yuv2yuv4mpeg_LDADD =

# vq_train
tools_vq_train_SOURCES = \
	tools/vq_train.c
tools_vq_train_CFLAGS = $(OPENMP_CFLAGS)
tools_vq_train_LDADD = $(LIBM)

# draw_zigzags
tools_draw_zigzags_SOURCES = \
	tools/draw_zigzags.c \
	src/zigzag4.c \
	src/zigzag8.c \
	src/zigzag16.c \
	src/zigzag32.c \
	src/zigzag64.c
tools_draw_zigzags_CFLAGS =
tools_draw_zigzags_LDADD =


endif # ENABLE_TOOLS

# Tests

if ENABLE_UNIT_TESTS
noinst_PROGRAMS += \
	src/tests/dcttest \
	src/tests/test_32x32 \
	src/tests/ectest \
	src/tests/test_coef_coder \
	src/tests/logging_test \
	src/tests/test_divu_small \
	src/tests/check_tests

TESTS = \
	src/tests/dcttest \
	src/tests/test_32x32 \
	src/tests/ectest \
	src/tests/test_coef_coder \
	src/tests/logging_test \
	src/tests/test_divu_small \
	src/tests/check_tests

src_tests_dcttest_SOURCES = $(src_dct_SOURCES) src/filter.c
src_tests_dcttest_CFLAGS = $(OGG_CFLAGS) \
 -DOD_DCT_CHECK_OVERFLOW -DOD_DCT_TEST -DOD_ENABLE_ASSERTIONS
src_tests_dcttest_LDADD = $(LIBM)

src_tests_test_32x32_SOURCES = src/tests/test_32x32.c
src_tests_test_32x32_CFLAGS = $(OGG_CFLAGS) -D_GNU_SOURCE
src_tests_test_32x32_LDADD = $(LIBM)

src_tests_ectest_SOURCES = src/tests/ectest.c
src_tests_ectest_CFLAGS = $(OGG_CFLAGS)
src_tests_ectest_LDADD = \
 src/libdaalaenc.la \
 src/libdaalabase.la \
 $(LIBM)

src_tests_test_coef_coder_SOURCES = src/tests/test_coef_coder.c
src_tests_test_coef_coder_CFLAGS = $(OGG_CFLAGS)
src_tests_test_coef_coder_LDADD = \
 src/libdaalabase.la \
 src/libdaaladec.la \
 src/libdaalaenc.la \
 $(OGG_LIBS) \
 $(LIBM)

src_tests_logging_test_SOURCES = src/tests/logging_test.c
src_tests_logging_test_CFLAGS = $(OGG_CFLAGS)
src_tests_logging_test_LDADD = \
 src/libdaalabase.la \
 src/libdaalaenc.la \
 $(OGG_LIBS)

src_tests_test_divu_small_SOURCES = src/tests/test_divu_small.c
src_tests_test_divu_small_CFLAGS = $(OGG_CFLAGS)
src_tests_test_divu_small_LDADD = \
 src/libdaalabase.la \
 $(OGG_LIBS)

src_tests_check_tests_SOURCES = \
 src/tests/check_main.c \
 src/tests/headerencode_test.c
src_tests_check_tests_CFLAGS = $(OGG_CFLAGS) $(CHECK_CFLAGS) -Wno-variadic-macros
src_tests_check_tests_LDADD = \
 src/libdaalabase.la \
 src/libdaaladec.la \
 src/libdaalaenc.la \
 $(OGG_LIBS) \
 $(CHECK_LIBS) \
 $(LIBM)

endif

SUBDIRS=doc

ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS = -I$(top_srcdir)/include $(DEPS_CFLAGS)

dist_doc_DATA = COPYING AUTHORS

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = daalaenc.pc daaladec.pc

debug: DEBUG_OPTIONS = -DOD_ENABLE_ASSERTIONS -DOD_LOGGING_ENABLED
debug:
	$(MAKE) CFLAGS="$(CFLAGS) -O0 -ggdb3 $(DEBUG_OPTIONS)" all

EXTRA_DIST = \
	daalaenc.pc.in \
	daaladec.pc.in \
	daalaenc-uninstalled.pc.in \
	daaladec-uninstalled.pc.in \
	tools/unix/Makefile \
	unix/Makefile

# Targets to build and install just the library without the docs
daala install-daala: NO_DOXYGEN = 1

daala: all
install-daala: install

if ENABLE_TOOLS
tools: $(tools_TARGETS)
else
tools:
	@echo "error: Tools support not enabled. Try running 'configure --enable-tools'."
endif

# Or just the docs
docs: doc/doxygen-build.stamp

install-docs:
	@if [ -z "$(NO_DOXYGEN)" ]; then \
	  ( cd doc && \
	    echo "Installing documentation in $(DESTDIR)$(docdir)"; \
	    $(INSTALL) -d $(DESTDIR)$(docdir)/html/search; \
	    for f in `find html -type f \! -name "installdox"` ; do \
	      $(INSTALL_DATA) $$f $(DESTDIR)$(docdir)/$$f;     \
	    done ) \
	fi

doc/doxygen-build.stamp: doc/Doxyfile $(daalainclude_HEADERS)
	@[ -n "$(NO_DOXYGEN)" ] || ( cd doc && doxygen && touch $(@F) )

if HAVE_DOXYGEN

# Or everything (by default)
all-local: docs

install-data-local: install-docs

clean-local:
	$(RM) -r doc/html
	$(RM) -r doc/latex
	$(RM) -r doc/man
	$(RM) doc/doxygen-build.stamp
	$(RM) doc/doxygen_sqlite3.db

uninstall-local:
	$(RM) -r $(DESTDIR)$(docdir)/html

endif

# We check this every time make is run, with configure.ac being touched to
# trigger an update of the build system files if update_version changes the
# current PACKAGE_VERSION (or if package_version was modified manually by a
# user with either AUTO_UPDATE=no or no update_version script present - the
# latter being the normal case for tarball releases).
#
# We can't just add the package_version file to CONFIGURE_DEPENDENCIES since
# simply running autoconf will not actually regenerate configure for us when
# the content of that file changes (due to autoconf dependency checking not
# knowing about that without us creating yet another file for it to include).
#
# The MAKECMDGOALS check is a gnu-make'ism, but will degrade 'gracefully' for
# makes that don't support it.  The only loss of functionality is not forcing
# an update of package_version for `make dist` if AUTO_UPDATE=no, but that is
# unlikely to be a real problem for any real user.
$(top_srcdir)/configure.ac: force
	@case "$(MAKECMDGOALS)" in \
	    dist-hook)                             exit 0       ;; \
	    dist-* | dist | distcheck | distclean) _arg=release ;; \
	esac; \
	if ! $(top_srcdir)/update_version $$_arg 2> /dev/null; then \
	    if [ ! -e $(top_srcdir)/package_version ]; then \
		echo 'PACKAGE_VERSION="unknown"' > $(top_srcdir)/package_version; \
	    fi; \
	    . $(top_srcdir)/package_version || exit 1; \
	    [ "$(PACKAGE_VERSION)" != "$$PACKAGE_VERSION" ] || exit 0; \
	fi; \
	touch $@

force:

# Create a minimal package_version file when make dist is run.
dist-hook:
	echo 'PACKAGE_VERSION="$(PACKAGE_VERSION)"' > $(top_distdir)/package_version

.PHONY: daala install-daala docs install-docs figures tools
