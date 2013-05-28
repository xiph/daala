#Daala video codec
#Copyright (c) 2002-2012 Daala project contributors.  All rights reserved.
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



AM_CFLAGS=-DOD_ENABLE_ASSERTIONS -Werror -fopenmp


noinst_PROGRAMS = png2y4m y4m2png dump_psnrhvs block_size_analysis plot_intra_maps \
		  init_intra_maps intra_stats intra_pred intra_trace trans \
		  trans2d init_intra_xform gen_cdf gen_laplace_tables

noinst_HEADERS = \
	cholesky.h \
	image.h \
	image_tools.h \
	int_search.h \
	intra_fit_tools.h \
	kiss99.h \
	matidx.h \
	od_covmat.h \
	od_defs.h \
	od_filter.h \
	od_intra.h \
	pythag.h \
	stats_tools.h \
	svd.h \
	trans_tools.h \
	vidinput.h

# png2y4m
png2y4m_SOURCES = \
	kiss99.c \
	png2y4m.c
png2y4m_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
png2y4m_LDADD = $(OGG_LIBS) $(PNG_LIBS)


# y4m2png

y4m2png_SOURCES = \
	vidinput.c \
	y4m2png.c
y4m2png_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
y4m2png_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS)


# dump_psnrhvs

dump_psnrhvs_SOURCES = \
	vidinput.c \
	$(top_srcdir)/src/newdct.c \
	dump_psnrhvs.c
dump_psnrhvs_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
dump_psnrhvs_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm


# block_size_analysis

block_size_analysis_SOURCES = \
	block_size_analysis.c \
	$(top_srcdir)/src/block_size.c \
	vidinput.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/internal.c \
	$(top_srcdir)/src/entcode.c \
	$(top_srcdir)/src/entenc.c \
	$(top_srcdir)/src/filter.c
block_size_analysis_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
block_size_analysis_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm


# plot_intra_maps
plot_intra_maps_SOURCES = \
	image.c \
	plot_intra_maps.c \
	intra_fit_tools.c \
	vidinput.c
plot_intra_maps_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
plot_intra_maps_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm



# init_intra_maps
init_intra_maps_SOURCES = \
	init_intra_maps.c \
	intra_fit_tools.c \
	vidinput.c
init_intra_maps_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
init_intra_maps_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm



# intra_stats
intra_stats_SOURCES = \
	intra_stats.c \
	intra_fit_tools.c \
	od_covmat.c \
	od_filter.c \
	od_intra.c \
	od_intra_data.c \
	image_tools.c \
	stats_tools.c \
	vidinput.c \
	image.c \
	$(top_srcdir)/src/filter.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/intra.c \
	$(top_srcdir)/src/tf.c \
	$(top_srcdir)/src/internal.c \
	$(top_srcdir)/src/intradata.c
intra_stats_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS) -fopenmp
intra_stats_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm


# intra_pred
intra_pred_SOURCES = \
	intra_pred.c \
	intra_fit_tools.c \
	od_covmat.c \
	od_filter.c \
	od_intra.c \
	od_intra_data.c \
	image_tools.c \
	stats_tools.c \
	vidinput.c \
	image.c \
	svd.c \
	cholesky.c \
	$(top_srcdir)/src/filter.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/intra.c \
	$(top_srcdir)/src/intradata.c \
	$(top_srcdir)/src/internal.c \
	$(top_srcdir)/src/tf.c
intra_pred_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS) -fopenmp
intra_pred_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm


# intra_trace
intra_trace_SOURCES = \
	intra_trace.c \
	intra_fit_tools.c \
	od_covmat.c \
	od_filter.c \
	od_intra.c \
	od_intra_data.c \
	image_tools.c \
	stats_tools.c \
	trans_tools.c \
	vidinput.c \
	image.c \
	$(top_srcdir)/src/filter.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/intra.c \
	$(top_srcdir)/src/tf.c \
	$(top_srcdir)/src/internal.c \
	$(top_srcdir)/src/intradata.c
intra_trace_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS) -fopenmp
intra_trace_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm




# trans
trans_SOURCES = \
	trans.c \
	$(top_srcdir)/src/internal.c \
	intra_fit_tools.c \
	vidinput.c \
	od_covmat.c \
	od_filter.c \
	od_intra.c \
	od_intra_data.c \
	stats_tools.c \
	trans_tools.c \
	int_search.c \
	trans_data.c \
	kiss99.c \
	image.c \
	svd.c \
	cholesky.c \
	$(top_srcdir)/src/filter.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/intradata.c
trans_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS) -fopenmp
trans_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm



#trans2d
trans2d_SOURCES = \
	trans2d.c \
	$(top_srcdir)/src/internal.c \
	intra_fit_tools.c \
	vidinput.c \
	od_covmat.c \
	od_filter.c \
	od_intra.c \
	od_intra_data.c \
	stats_tools.c \
	trans_tools.c \
	int_search.c \
	trans_data.c \
	kiss99.c \
	image.c \
	svd.c \
	cholesky.c \
	$(top_srcdir)/src/filter.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/intradata.c
trans2d_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS) -fopenmp
trans2d_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm



# init_intra_xform
init_intra_xform_SOURCES = \
	cholesky.c \
	init_intra_xform.c \
	intra_fit_tools.c \
	svd.c \
	vidinput.c \
	$(top_srcdir)/src/filter.c \
	$(top_srcdir)/src/newdct.c \
	$(top_srcdir)/src/intra.c \
	$(top_srcdir)/src/tf.c \
	$(top_srcdir)/src/internal.c \
	$(top_srcdir)/src/intradata.c
init_intra_xform_CFLAGS = $(THEORA_CFLAGS) $(OGG_CFLAGS) $(PNG_CFLAGS)
init_intra_xform_LDADD = $(THEORA_LIBS) $(OGG_LIBS) $(PNG_LIBS) -lm



# gen_cdf
gen_cdf_SOURCES = \
	gen_cdf.c
gen_cdf_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
gen_cdf_LDADD = $(OGG_LIBS) $(PNG_LIBS) -lm


# gen_laplace_tables
gen_laplace_tables_SOURCES = \
	gen_laplace_tables.c
gen_laplace_tables_CFLAGS = $(OGG_CFLAGS) $(PNG_CFLAGS)
gen_laplace_tables_LDADD = $(OGG_LIBS) $(PNG_LIBS) -lm



