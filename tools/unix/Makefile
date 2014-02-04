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
CFLAGS = -Wall -Wno-parentheses -Wno-overlength-strings -Wno-long-long
CFLAGS += -std=c89 -pedantic -fopenmp
CFLAGS += -O3
CFLAGS += -g
CFLAGS += -DOD_ENABLE_ASSERTIONS

# You could instead set this to HAVE_GETTIMEOFDAY if you don't have ftime.
CFLAGS += -DHAVE_FTIME

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
newdct.c \
dump_psnrhvs.c \

PSNRHVS_LIB_CSOURCES =

PSNRHVS_LDFLAGS = `pkg-config ogg libpng --libs` -lm

PSNRHVS_OBJS:=${PSNRHVS_LOCAL_CSOURCES} ${PSNRHVS_LIB_CSOURCES}
PSNRHVS_OBJS:=${PSNRHVS_OBJS:%.c=${WORKDIR}/%.o}
PSNRHVS_CSOURCES:=${PSNRHVS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${PSNRHVS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
PSNRHVS_TARGET:=dump_psnrhvs

#block_size_analysis
BSIZE_LOCAL_CSOURCES = \
block_size_analysis.c \
block_size.c \
block_size_enc.c \
vidinput.c \
y4m_input.c \
newdct.c \
internal.c \
entcode.c \
entenc.c \
filter.c \
switch_table.c


BSIZE_LIB_CSOURCES =

BSIZE_LDFLAGS = `pkg-config ogg libpng --libs` -lm

BSIZE_OBJS:=${BSIZE_LOCAL_CSOURCES} ${BSIZE_LIB_CSOURCES}
BSIZE_OBJS:=${BSIZE_OBJS:%.c=${WORKDIR}/%.o}
BSIZE_CSOURCES:=${BSIZE_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${BSIZE_LIB_CSOURCES:%=${LIBSRCDIR}/%}
BSIZE_TARGET:=block_size_analysis


#plot_intra_maps
PIM_LOCAL_CSOURCES = \
image.c \
plot_intra_maps.c \
intra_fit_tools.c \
vidinput.c \
y4m_input.c \

PIM_LIB_CSOURCES =

PIM_LDFLAGS = `pkg-config ogg libpng --libs` -lm

PIM_OBJS:=${PIM_LOCAL_CSOURCES} ${PIM_LIB_CSOURCES}
PIM_OBJS:=${PIM_OBJS:%.c=${WORKDIR}/%.o}
PIM_CSOURCES:=${PIM_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${PIM_LIB_CSOURCES:%=${LIBSRCDIR}/%}
PIM_TARGET:=plot_intra_maps


#init_intra_maps
IIM_LOCAL_CSOURCES = \
init_intra_maps.c \
intra_fit_tools.c \
vidinput.c \
y4m_input.c \

IIM_LIB_CSOURCES =

IIM_LDFLAGS = `pkg-config ogg libpng --libs` -lm

IIM_OBJS:=${IIM_LOCAL_CSOURCES} ${IIM_LIB_CSOURCES}
IIM_OBJS:=${IIM_OBJS:%.c=${WORKDIR}/%.o}
IIM_CSOURCES:=${IIM_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${IIM_LIB_CSOURCES:%=${LIBSRCDIR}/%}
IIM_TARGET:=init_intra_maps


#intra_stats
IS_LOCAL_CSOURCES = \
intra_stats.c \
intra_fit_tools.c \
od_covmat.c \
od_filter.c \
od_intra.c \
od_intra_data.c \
image_tools.c \
stats_tools.c \
vidinput.c \
y4m_input.c \
image.c

IS_LIB_CSOURCES = \
block_size.c \
block_size_enc.c \
entenc.c \
entcode.c \
filter.c \
newdct.c \
intra.c \
switch_table.c \
tf.c \
internal.c \
intradata.c

IS_LDFLAGS = `pkg-config ogg libpng --libs` -lm

IS_OBJS:=${IS_LOCAL_CSOURCES} ${IS_LIB_CSOURCES}
IS_OBJS:=${IS_OBJS:%.c=${WORKDIR}/%.o}
IS_CSOURCES:=${IS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${IS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
IS_TARGET:=intra_stats


# cfl_stats
CS_LOCAL_CSOURCES = \
cfl_stats.c \
intra_fit_tools.c \
od_covmat.c \
od_filter.c \
od_intra.c \
od_intra_data.c \
image_tools.c \
stats_tools.c \
vidinput.c \
y4m_input.c \
image.c

CS_LIB_CSOURCES = \
block_size.c \
block_size_enc.c \
entenc.c \
entcode.c \
filter.c \
newdct.c \
intra.c \
switch_table.c \
tf.c \
internal.c \
intradata.c

CS_LDFLAGS = `pkg-config ogg libpng --libs` -lm

CS_OBJS:=${CS_LOCAL_CSOURCES} ${CS_LIB_CSOURCES}
CS_OBJS:=${CS_OBJS:%.c=${WORKDIR}/%.o}
CC_CSOURCES:=${CS_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${CS_LIB_CSOURCES:%=${LIBSRCDIR}/%}
CS_TARGET:=cfl_stats


#intra_pred
IP_LOCAL_CSOURCES = \
intra_pred.c \
intra_fit_tools.c \
od_covmat.c \
od_filter.c \
od_intra.c \
od_intra_data.c \
image_tools.c \
stats_tools.c \
vidinput.c \
y4m_input.c \
image.c \
svd.c \
cholesky.c

IP_LIB_CSOURCES = \
block_size.c \
block_size_enc.c \
entenc.c \
entcode.c \
filter.c \
newdct.c \
intra.c \
intradata.c \
internal.c \
switch_table.c \
tf.c

IP_LDFLAGS = `pkg-config ogg libpng --libs` -lm

IP_OBJS:=${IP_LOCAL_CSOURCES} ${IP_LIB_CSOURCES}
IP_OBJS:=${IP_OBJS:%.c=${WORKDIR}/%.o}
IP_CSOURCES:=${IP_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${IP_LIB_CSOURCES:%=${LIBSRCDIR}/%}
IP_TARGET:=intra_pred


#intra_trace
IT_LOCAL_CSOURCES = \
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
y4m_input.c \
image.c

IT_LIB_CSOURCES = \
block_size.c \
block_size_enc.c \
entenc.c \
entcode.c \
filter.c \
newdct.c \
intra.c \
switch_table.c \
tf.c \
internal.c \
intradata.c

IT_LDFLAGS = `pkg-config ogg libpng --libs` -lm

IT_OBJS:=${IT_LOCAL_CSOURCES} ${IT_LIB_CSOURCES}
IT_OBJS:=${IT_OBJS:%.c=${WORKDIR}/%.o}
IT_CSOURCES:=${IT_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${IT_LIB_CSOURCES:%=${LIBSRCDIR}/%}
IT_TARGET:=intra_trace


#trans
T_LOCAL_CSOURCES = \
trans.c \
internal.c \
intra_fit_tools.c \
vidinput.c \
y4m_input.c \
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
cholesky.c

T_LIB_CSOURCES = \
filter.c \
newdct.c \
intradata.c

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
trans_tools.c 

TG_LIB_CSOURCES = \
filter.c \
newdct.c \
intradata.c

TG_LDFLAGS = `pkg-config ogg libpng --libs` -lm

TG_OBJS:=${TG_LOCAL_CSOURCES} ${TG_LIB_CSOURCES}
TG_OBJS:=${TG_OBJS:%.c=${WORKDIR}/%.o}
TG_CSOURCES:=${TG_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${TG_LIB_CSOURCES:%=${LIBSRCDIR}/%}
TG_TARGET:=trans_gain


#trans2d
T2D_LOCAL_CSOURCES = \
trans2d.c \
internal.c \
intra_fit_tools.c \
vidinput.c \
y4m_input.c \
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
cholesky.c

T2D_LIB_CSOURCES = \
filter.c \
newdct.c \
intradata.c

T2D_LDFLAGS = `pkg-config ogg libpng --libs` -lm

T2D_OBJS:=${T2D_LOCAL_CSOURCES} ${T2D_LIB_CSOURCES}
T2D_OBJS:=${T2D_OBJS:%.c=${WORKDIR}/%.o}
T2D_CSOURCES:=${T2D_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${T2D_LIB_CSOURCES:%=${LIBSRCDIR}/%}
T2D_TARGET:=trans2d


#init_intra_xform
IIX_LOCAL_CSOURCES = \
cholesky.c \
init_intra_xform.c \
intra_fit_tools.c \
svd.c \
vidinput.c \
y4m_input.c \

IIX_LIB_CSOURCES = \
filter.c \
newdct.c \
intra.c \
tf.c \
internal.c \
intradata.c

IIX_LDFLAGS = `pkg-config ogg libpng --libs` -lm

IIX_OBJS:=${IIX_LOCAL_CSOURCES} ${IIX_LIB_CSOURCES}
IIX_OBJS:=${IIX_OBJS:%.c=${WORKDIR}/%.o}
IIX_CSOURCES:=${IIX_LOCAL_CSOURCES:%=${LOCALSRCDIR}/%} ${IIX_LIB_CSOURCES:%=${LIBSRCDIR}/%}
IIX_TARGET:=init_intra_xform


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

ALL_OBJS:=${P2Y_OBJS} ${Y2P_OBJS} ${PSNR_OBJS} ${PSNRHVS_OBJS} ${PIM_OBJS} ${IIM_OBJS} ${IIX_OBJS} \
 ${GEN_CDF_OBJS} ${GEN_LAPLACE_TABLES_OBJS} ${IS_OBJS} ${IP_OBJS} ${IT_OBJS} ${T_OBJS} \
 ${TG_OBJS} ${T2D_OBJS} ${BSIZE_OBJS} ${BD_OBJS} ${YJ_OBJS} ${JY_OBJS}

ALL_ASMS:=${ALL_OBJS:%.o=%.s}

ALL_DEPS:=${ALL_OBJS:%.o=%.d}

ALL_TARGETS:=${P2Y_TARGET} ${Y2P_TARGET} ${PSNR_TARGET} ${PSNRHVS_TARGET} \
 ${PIM_TARGET} ${IIM_TARGET} ${IIX_TARGET} ${GEN_CDF_TARGET} \
 ${GEN_LAPLACE_TABLES_TARGET} ${IS_TARGET} ${IP_TARGET} ${IT_TARGET} \
 ${T_TARGET} ${TG_TARGET} ${T2D_TARGET} ${BSIZE_TARGET} ${BD_TARGET} \
 ${YJ_TARGET} ${JY_TARGET}

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

${PIM_TARGET}: ${PIM_OBJS}
	${CC} ${CFLAGS} ${PIM_OBJS} ${PIM_LDFLAGS} -o $@

${IIM_TARGET}: ${IIM_OBJS}
	${CC} ${CFLAGS} ${IIM_OBJS} ${IIM_LDFLAGS} -o $@

${IS_TARGET}: ${IS_OBJS}
	${CC} ${CFLAGS} ${IS_OBJS} ${IS_LDFLAGS} -o $@

${IP_TARGET}: ${IP_OBJS}
	${CC} ${CFLAGS} ${IP_OBJS} ${IP_LDFLAGS} -o $@

${IT_TARGET}: ${IT_OBJS}
	${CC} ${CFLAGS} ${IT_OBJS} ${IT_LDFLAGS} -o $@

${T_TARGET}: ${T_OBJS}
	${CC} ${CFLAGS} ${T_OBJS} ${T_LDFLAGS} -o $@

${TG_TARGET}: ${TG_OBJS}
	${CC} ${CFLAGS} ${TG_OBJS} ${TG_LDFLAGS} -o $@

${T2D_TARGET}: ${T2D_OBJS}
	${CC} ${CFLAGS} ${T2D_OBJS} ${T2D_LDFLAGS} -o $@

${CS_TARGET}: ${CS_OBJS}
	${CC} ${CFLAGS} ${CS_OBJS} ${CS_LDFLAGS} -o $@

${IIX_TARGET}: ${IIX_OBJS}
	${CC} ${CFLAGS} ${IIX_OBJS} ${IIX_LDFLAGS} -o $@

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
