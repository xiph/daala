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

MAKEDEPEND = $(CC) -MM
CINCLUDE = `pkg-config ogg libpng --cflags`
CFLAGS = -Wall -Wno-overlength-strings -Wno-long-long -Wshadow
CFLAGS += -std=c89 -pedantic -fopenmp
CFLAGS += -O3
CFLAGS += -g
CFLAGS += -DOD_ENABLE_ASSERTIONS

# You could instead set this to HAVE_GETTIMEOFDAY if you don't have ftime.
CFLAGS += -DHAVE_FTIME

CINCLUDE := -I../../include ${CINCLUDE}
LIBSRCDIR=../../src
LOCALSRCDIR=..
WORKDIR=build

#png2y4m
P2Y_LOCAL_CSOURCES = \
kiss99.c \
png2y4m.c

P2Y_LIB_CSOURCES =

P2Y_LDFLAGS = `pkg-config ogg libpng --libs`

P2Y_OBJS:=${P2Y_LOCAL_CSOURCES} ${P2Y_LIB_CSOURCES}
P2Y_OBJS:=${P2Y_OBJS:%.c=${WORKDIR}/%.o}
P2Y_CSOURCES:=${P2Y_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${P2Y_LIB_CSOURCES:%=${LIBSRCDIR}/%}
P2Y_TARGET:=png2y4m

#y4m2png
Y2P_LOCAL_CSOURCES = \
vidinput.c \
y4m_input.c \
y4m2png.c

Y2P_LIB_CSOURCES =

Y2P_LDFLAGS = `pkg-config ogg libpng --libs`

Y2P_OBJS:=${Y2P_LOCAL_CSOURCES} ${Y2P_LIB_CSOURCES}
Y2P_OBJS:=${Y2P_OBJS:%.c=${WORKDIR}/%.o}
Y2P_CSOURCES:=${Y2P_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${Y2P_LIB_CSOURCES:%=${LIBSRCDIR}/%}
Y2P_TARGET:=y4m2png

#dump_psnr
PSNR_LOCAL_CSOURCES = \
vidinput.c \
y4m_input.c \
dump_psnr.c

PSNR_LIB_CSOURCES =

PSNR_LDFLAGS = `pkg-config ogg libpng --libs` -lm

PSNR_OBJS:=${PSNR_LOCAL_CSOURCES} ${PSNR_LIB_CSOURCES}
PSNR_OBJS:=${PSNR_OBJS:%.c=${WORKDIR}/%.o}
PSNR_CSOURCES:=${PSNR_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${PSNR_LIB_CSOURCES:%=${LIBSRCDIR}/%}
PSNR_TARGET:=dump_psnr

#dump_psnrhvs
PSNRHVS_LOCAL_CSOURCES = \
vidinput.c \
y4m_input.c \
dct.c \
dump_psnrhvs.c \
tf.c \
internal.c

PSNRHVS_LIB_CSOURCES =

PSNRHVS_LDFLAGS = `pkg-config ogg libpng --libs` -lm

PSNRHVS_OBJS:=${PSNRHVS_LOCAL_CSOURCES} ${PSNRHVS_LIB_CSOURCES}
PSNRHVS_OBJS:=${PSNRHVS_OBJS:%.c=${WORKDIR}/%.o}
PSNRHVS_CSOURCES:=${PSNRHVS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${PSNRHVS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
PSNRHVS_TARGET:=dump_psnrhvs

#block_size_analysis
BSIZE_LOCAL_CSOURCES = \
block_size_analysis.c \
block_size_enc.c \
generic_encoder.c \
generic_code.c \
laplace_encoder.c \
logging.c \
laplace_tables.c \
vidinput.c \
y4m_input.c \
dct.c \
internal.c \
entcode.c \
entenc.c \
filter.c \
switch_table.c \
tf.c

BSIZE_LIB_CSOURCES =

BSIZE_LDFLAGS = `pkg-config ogg libpng --libs` -lm

BSIZE_OBJS:=${BSIZE_LOCAL_CSOURCES} ${BSIZE_LIB_CSOURCES}
BSIZE_OBJS:=${BSIZE_OBJS:%.c=${WORKDIR}/%.o}
BSIZE_CSOURCES:=${BSIZE_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${BSIZE_LIB_CSOURCES:%=${LIBSRCDIR}/%}
BSIZE_TARGET:=block_size_analysis

#to_monochrome
MONO_LOCAL_CSOURCES = \
vidinput.c \
y4m_input.c \
to_monochrome.c

MONO_LIB_CSOURCES =

MONO_LDFLAGS = `pkg-config ogg libpng --libs`

MONO_OBJS:=${MONO_LOCAL_CSOURCES} ${MONO_LIB_CSOURCES}
MONO_OBJS:=${MONO_OBJS:%.c=${WORKDIR}/%.o}
MONO_CSOURCES:=${MONO_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${MONO_LIB_CSOURCES:%=${LIBSRCDIR}/%}
MONO_TARGET:=to_monochrome


#downsample
DS_LOCAL_CSOURCES = \
vidinput.c \
y4m_input.c \
downsample.c

DS_LIB_CSOURCES =

DS_LDFLAGS = `pkg-config ogg libpng --libs`

DS_OBJS:=${DS_LOCAL_CSOURCES} ${DS_LIB_CSOURCES}
DS_OBJS:=${DS_OBJS:%.c=${WORKDIR}/%.o}
DS_CSOURCES:=${DS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${DS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
DS_TARGET:=downsample


#upsample
US_LOCAL_CSOURCES = \
vidinput.c \
y4m_input.c \
upsample.c \
filter.c \
util.c

US_LIB_CSOURCES = \
dct.c \
internal.c \
generic_code.c \
switch_table.c \
logging.c \
info.c \
mc.c \
pvq.c \
state.c

US_LDFLAGS = `pkg-config ogg libpng --libs` -lm

US_OBJS:=${US_LOCAL_CSOURCES} ${US_LIB_CSOURCES}
US_OBJS:=${US_OBJS:%.c=${WORKDIR}/%.o}
US_CSOURCES:=${US_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${US_LIB_CSOURCES:%=${LIBSRCDIR}/%}
US_TARGET:=upsample

#trans
T_LOCAL_CSOURCES = \
trans.c \
internal.c \
vidinput.c \
y4m_input.c \
od_filter.c \
trans_tools.c \
int_search.c \
trans_data.c \
kiss99.c \
svd.c \
cholesky.c \
tf.c

T_LIB_CSOURCES = \
filter.c \
dct.c

T_LDFLAGS = `pkg-config ogg libpng --libs` -lm

T_OBJS:=${T_LOCAL_CSOURCES} ${T_LIB_CSOURCES}
T_OBJS:=${T_OBJS:%.c=${WORKDIR}/%.o}
T_CSOURCES:=${T_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${T_LIB_CSOURCES:%=${LIBSRCDIR}/%}
T_TARGET:=trans

#trans_gain
TG_LOCAL_CSOURCES = \
trans_gain.c \
internal.c \
vidinput.c \
y4m_input.c \
od_filter.c \
trans_tools.c \
tf.c

TG_LIB_CSOURCES = \
filter.c \
dct.c

TG_LDFLAGS = `pkg-config ogg libpng --libs` -lm

TG_OBJS:=${TG_LOCAL_CSOURCES} ${TG_LIB_CSOURCES}
TG_OBJS:=${TG_OBJS:%.c=${WORKDIR}/%.o}
TG_CSOURCES:=${TG_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${TG_LIB_CSOURCES:%=${LIBSRCDIR}/%}
TG_TARGET:=trans_gain

#compute_basis
CBASIS_LOCAL_CSOURCES = \
compute_basis.c

CBASIS_LIB_CSOURCES = \
dct.c \
filter.c \
internal.c

CBASIS_LDFLAGS = `pkg-config ogg libpng --libs` -lm

CBASIS_OBJS:=${CBASIS_LOCAL_CSOURCES} ${CBASIS_LIB_CSOURCES}
CBASIS_OBJS:=${CBASIS_OBJS:%.c=${WORKDIR}/%.o}
CBASIS_CSOURCES:=${CBASIS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${CBASIS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
CBASIS_TARGET:=compute_basis

#compute_haar_basis
CHBASIS_LOCAL_CSOURCES = \
compute_haar_basis.c

CHBASIS_LIB_CSOURCES = \
dct.c \
filter.c \
internal.c

CHBASIS_LDFLAGS = `pkg-config ogg libpng --libs` -lm

CHBASIS_OBJS:=${CHBASIS_LOCAL_CSOURCES} ${CHBASIS_LIB_CSOURCES}
CHBASIS_OBJS:=${CHBASIS_OBJS:%.c=${WORKDIR}/%.o}
CHBASIS_CSOURCES:=${CHBASIS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${CHBASIS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
CHBASIS_TARGET:=compute_haar_basis

#trans2d
T2D_LOCAL_CSOURCES = \
trans2d.c \
internal.c \
vidinput.c \
y4m_input.c \
od_filter.c \
trans_tools.c \
int_search.c \
trans_data.c \
kiss99.c \
svd.c \
cholesky.c \
tf.c

T2D_LIB_CSOURCES = \
filter.c \
dct.c

T2D_LDFLAGS = `pkg-config ogg libpng --libs` -lm

T2D_OBJS:=${T2D_LOCAL_CSOURCES} ${T2D_LIB_CSOURCES}
T2D_OBJS:=${T2D_OBJS:%.c=${WORKDIR}/%.o}
T2D_CSOURCES:=${T2D_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${T2D_LIB_CSOURCES:%=${LIBSRCDIR}/%}
T2D_TARGET:=trans2d

#gen_cdf
GEN_CDF_LOCAL_CSOURCES = \
gen_cdf.c

GEN_CDF_LIB_CSOURCES =

GEN_CDF_LDFLAGS = -lm

GEN_CDF_OBJS:=${GEN_CDF_LOCAL_CSOURCES} ${GEN_CDF_LIB_CSOURCES}
GEN_CDF_OBJS:=${GEN_CDF_OBJS:%.c=${WORKDIR}/%.o}
GEN_CDF_CSOURCES:=${GEN_CDF_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} \
 ${GEN_CDF_LIB_CSOURCES:%=${LIBSRCDIR}/%}
GEN_CDF_TARGET:=gen_cdf

#gen_laplace_tables
GEN_LAPLACE_TABLES_LOCAL_CSOURCES = \
gen_laplace_tables.c

GEN_LAPLACE_TABLES_LIB_CSOURCES =

GEN_LAPLACE_TABLES_LDFLAGS = -lm

GEN_LAPLACE_TABLES_OBJS:=${GEN_LAPLACE_TABLES_LOCAL_CSOURCES} ${GEN_LAPLACE_TABLES_LIB_CSOURCES}
GEN_LAPLACE_TABLES_OBJS:=${GEN_LAPLACE_TABLES_OBJS:%.c=${WORKDIR}/%.o}
GEN_LAPLACE_TABLES_CSOURCES:= \
 ${GEN_LAPLACE_TABLES_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} \
 ${GEN_LAPLACE_TABLES_LIB_CSOURCES:%=${LIBSRCDIR}/%}
GEN_LAPLACE_TABLES_TARGET:=gen_laplace_tables

#bjontegaard
BD_LOCAL_CSOURCES = \
bjontegaard.c \
cholesky.c \
svd.c \
qr.c

BD_LIB_CSOURCES =

BD_LDFLAGS = -lm

BD_OBJS:=${BD_LOCAL_CSOURCES} ${BD_LIB_CSOURCES}
BD_OBJS:=${BD_OBJS:%.c=${WORKDIR}/%.o}
BD_CSOURCES:=${BD_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} \
 ${BD_LIB_CSOURCES:%=${LIBSRCDIR}/%}
BD_TARGET:=bjontegaard

#yuvjpeg
YJ_LOCAL_CSOURCES = \
yuvjpeg.c

YJ_LIB_CSOURCES =

YJ_LDFLAGS = -ljpeg

YJ_OBJS:=${YJ_LOCAL_CSOURCES} ${YJ_LIB_CSOURCES}
YJ_OBJS:=${YJ_OBJS:%.c=${WORKDIR}/%.o}
YJ_CSOURCES:=${YJ_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} \
 ${YJ_LIB_CSOURCES:%=${LIBSRCDIR}/%}
YJ_TARGET:=yuvjpeg

#jpegyuv
JY_LOCAL_CSOURCES = \
jpegyuv.c

JY_LIB_CSOURCES =

JY_LDFLAGS = -ljpeg

JY_OBJS:=${JY_LOCAL_CSOURCES} ${JY_LIB_CSOURCES}
JY_OBJS:=${JY_OBJS:%.c=${WORKDIR}/%.o}
JY_CSOURCES:=${JY_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} \
 ${JY_LIB_CSOURCES:%=${LIBSRCDIR}/%}
JY_TARGET:=jpegyuv

#draw_zigzags
DZZ_LOCAL_CSOURCES = \
draw_zigzags.c

DZZ_LIB_CSOURCES = \
zigzag4.c \
zigzag8.c \
zigzag16.c \
zigzag32.c \
zigzag64.c

DZZ_LDFLAGS =

DZZ_OBJS:=${DZZ_LOCAL_CSOURCES} ${DZZ_LIB_CSOURCES}
DZZ_OBJS:=${DZZ_OBJS:%.c=${WORKDIR}/%.o}
DZZ_CSOURCES:=${DZZ_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} \
 ${DZZ_LIB_CSOURCES:%=${LIBSRCDIR}/%}
DZZ_TARGET:=draw_zigzags

ALL_OBJS:=${P2Y_OBJS} ${Y2P_OBJS} ${PSNR_OBJS} ${PSNRHVS_OBJS} \
 ${GEN_CDF_OBJS} ${GEN_LAPLACE_TABLES_OBJS} ${T_OBJS} \
 ${TG_OBJS} ${CBASIS_OBJS} ${CHBASIS_OBJS} ${T2D_OBJS} ${BSIZE_OBJS} \
 ${BD_OBJS} ${YJ_OBJS} ${JY_OBJS} ${DZZ_OBJS} ${MONO_OBJS} ${DS_OBJS} \
 ${US_OBJS}

ALL_ASMS:=${ALL_OBJS:%.o=%.s}

ALL_DEPS:=${ALL_OBJS:%.o=%.d}

ALL_TARGETS:=${P2Y_TARGET} ${Y2P_TARGET} ${PSNR_TARGET} ${PSNRHVS_TARGET} \
 ${GEN_CDF_TARGET} ${GEN_LAPLACE_TABLES_TARGET} ${T_TARGET} \
 ${TG_TARGET} ${CBASIS_TARGET} ${CHBASIS_TARGET} ${T2D_TARGET} \
 ${BSIZE_TARGET} ${BD_TARGET} ${YJ_TARGET} ${JY_TARGET} ${DZZ_TARGET} \
 ${MONO_TARGET} ${DS_TARGET} ${US_TARGET}

all: ${ALL_TARGETS}

asm: ${ALL_ASM}

check: all

clean:
	${RM} ${ALL_ASM} ${ALL_OBJS} ${ALL_DEPS}
	${RM} ${ALL_TARGETS}
	-rmdir ${WORKDIR}

${BSIZE_TARGET}: ${BSIZE_OBJS}
	${CC} ${CFLAGS} ${BSIZE_OBJS} ${BSIZE_LDFLAGS} -o $@

${P2Y_TARGET}: ${P2Y_OBJS}
	${CC} ${CFLAGS} ${P2Y_OBJS} ${P2Y_LDFLAGS} -o $@

${Y2P_TARGET}: ${Y2P_OBJS}
	${CC} ${CFLAGS} ${Y2P_OBJS} ${Y2P_LDFLAGS} -o $@

${PSNR_TARGET}: ${PSNR_OBJS}
	${CC} ${CFLAGS} ${PSNR_OBJS} ${PSNR_LDFLAGS} -o $@

${PSNRHVS_TARGET}: ${PSNRHVS_OBJS}
	${CC} ${CFLAGS} ${PSNRHVS_OBJS} ${PSNRHVS_LDFLAGS} -o $@

${T_TARGET}: ${T_OBJS}
	${CC} ${CFLAGS} ${T_OBJS} ${T_LDFLAGS} -o $@

${TG_TARGET}: ${TG_OBJS}
	${CC} ${CFLAGS} ${TG_OBJS} ${TG_LDFLAGS} -o $@

${CBASIS_TARGET}: ${CBASIS_OBJS}
	${CC} ${CFLAGS} ${CBASIS_OBJS} ${CBASIS_LDFLAGS} -o $@

${CHBASIS_TARGET}: ${CHBASIS_OBJS}
	${CC} ${CFLAGS} ${CHBASIS_OBJS} ${CHBASIS_LDFLAGS} -o $@

${T2D_TARGET}: ${T2D_OBJS}
	${CC} ${CFLAGS} ${T2D_OBJS} ${T2D_LDFLAGS} -o $@

${GEN_CDF_TARGET}: ${GEN_CDF_OBJS}
	${CC} ${CFLAGS} ${GEN_CDF_OBJS} ${GEN_CDF_LDFLAGS} -o $@

${GEN_LAPLACE_TABLES_TARGET}: ${GEN_LAPLACE_TABLES_OBJS}
	${CC} ${CFLAGS} ${GEN_LAPLACE_TABLES_OBJS} ${GEN_LAPLACE_TABLES_LDFLAGS} -o $@

${BD_TARGET}: ${BD_OBJS}
	${CC} ${CFLAGS} ${BD_OBJS} ${BD_LDFLAGS} -o $@

${YJ_TARGET}: ${YJ_OBJS}
	${CC} ${CFLAGS} ${YJ_OBJS} ${YJ_LDFLAGS} -o $@

${JY_TARGET}: ${JY_OBJS}
	${CC} ${CFLAGS} ${JY_OBJS} ${JY_LDFLAGS} -o $@

${DZZ_TARGET}: ${DZZ_OBJS}
	${CC} ${CFLAGS} ${DZZ_OBJS} ${DZZ_LDFLAGS} -o $@

${MONO_TARGET}: ${MONO_OBJS}
	${CC} ${CFLAGS} ${MONO_OBJS} ${MONO_LDFLAGS} -o $@

${DS_TARGET}: ${DS_OBJS}
	${CC} ${CFLAGS} ${DS_OBJS} ${DS_LDFLAGS} -o $@

${US_TARGET}: ${US_OBJS}
	${CC} ${CFLAGS} ${US_OBJS} ${US_LDFLAGS} -o $@

${ALL_OBJS} ${ALL_ASMS} ${ALL_DEPS} ${ALL_TARGETS}: Makefile

.PHONY: all asm clean check

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

${WORKDIR}/%.d: ${LOCALSRCDIR}/%.c
	mkdir -p ${dir $@}
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.o} > $@
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.s} >> $@
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT $@ >> $@
${WORKDIR}/%.s: ${LOCALSRCDIR}/%.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -S -o $@ $<
${WORKDIR}/%.o: ${LOCALSRCDIR}/%.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -c -o $@ $<

# Include header file dependencies
include ${ALL_DEPS}
