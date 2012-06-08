CC = gcc
MAKEDEPEND = gcc -MM
CFLAGS = -Wall -Wno-parentheses -Wno-overlength-strings -Wno-long-long
CFLAGS += -std=c89 -pedantic
CFLAGS += -O3
#CFLAGS += -g

LIBSRCDIR=../src
WORKDIR=build


#png2y4m
P2Y_LOCAL_CSOURCES = \
kiss99.c \
png2y4m.c

P2Y_LIB_CSOURCES =

P2Y_LDFLAGS = -ltheoradec -logg -lpng -lz

P2Y_OBJS:=${P2Y_LOCAL_CSOURCES} ${P2Y_LIB_CSOURCES}
P2Y_OBJS:=${P2Y_OBJS:%.c=${WORKDIR}/%.o}
P2Y_CSOURCES:=${P2Y_LOCAL_CSOURCES} ${P2Y_LOCAL_CSOURCES:%=${LIBSRCDIR}/%}
P2Y_TARGET:=png2y4m

#y4m2png
Y2P_LOCAL_CSOURCES = \
vidinput.c \
y4m2png.c

Y2P_LIB_CSOURCES =

Y2P_LDFLAGS = -ltheoradec -logg -lpng -lz

Y2P_OBJS:=${Y2P_LOCAL_CSOURCES} ${Y2P_LIB_CSOURCES}
Y2P_OBJS:=${Y2P_OBJS:%.c=${WORKDIR}/%.o}
Y2P_CSOURCES:=${Y2P_LOCAL_CSOURCES} ${Y2P_LOCAL_CSOURCES:%=${LIBSRCDIR}/%}
Y2P_TARGET:=y4m2png


#plot_intra_maps
PIM_LOCAL_CSOURCES = \
image.c \
plot_intra_maps.c \
intra_fit_tools.c \
vidinput.c

PIM_LIB_CSOURCES =

PIM_LDFLAGS = -ltheoradec -logg -lpng -lz -lm

PIM_OBJS:=${PIM_LOCAL_CSOURCES} ${PIM_LIB_CSOURCES}
PIM_OBJS:=${PIM_OBJS:%.c=${WORKDIR}/%.o}
PIM_CSOURCES:=${PIM_LOCAL_CSOURCES} ${PIM_LOCAL_CSOURCES:%=${LIBSRCDIR}/%}
PIM_TARGET:=plot_intra_maps


#init_intra_maps
IIM_LOCAL_CSOURCES = \
init_intra_maps.c \
intra_fit_tools.c \
svd.c \
vidinput.c

IIM_LIB_CSOURCES =

IIM_LDFLAGS = -ltheoradec -logg -lpng -lz -lm

IIM_OBJS:=${IIM_LOCAL_CSOURCES} ${IIM_LIB_CSOURCES}
IIM_OBJS:=${IIM_OBJS:%.c=${WORKDIR}/%.o}
IIM_CSOURCES:=${IIM_LOCAL_CSOURCES} ${IIM_LOCAL_CSOURCES:%=${LIBSRCDIR}/%}
IIM_TARGET:=init_intra_maps


#init_intra_xform
IIX_LOCAL_CSOURCES = \
init_intra_xform.c \
intra_fit_tools.c \
svd.c \
vidinput.c

IIX_LIB_CSOURCES = \
filter.c \
newdct.c

IIX_LDFLAGS = -ltheoradec -logg -lpng -lz -lm

IIX_OBJS:=${IIX_LOCAL_CSOURCES} ${IIX_LIB_CSOURCES}
IIX_OBJS:=${IIX_OBJS:%.c=${WORKDIR}/%.o}
IIX_CSOURCES:=${IIX_LOCAL_CSOURCES} ${IIX_LOCAL_CSOURCES:%=${LIBSRCDIR}/%}
IIX_TARGET:=init_intra_xform


ALL_OBJS:=${P2Y_OBJS} ${Y2P_OBJS} ${PIM_OBJS} ${IIM_OBJS} ${IIX_OBJS}

ALL_ASMS:=${ALL_OBJS:%.o=%.s}

ALL_DEPS:=${ALL_OBJS:%.o=%.d}

ALL_TARGETS:=${P2Y_TARGET} ${Y2P_TARGET} \
 ${PIM_TARGET} ${IIM_TARGET} ${IIX_TARGET}

all: ${ALL_TARGETS}

asm: ${ALL_ASM}

check: all

clean:
	${RM} ${ALL_ASM} ${ALL_OBJS} ${ALL_DEPS}
	${RM} ${ALL_TARGETS}
	-rmdir ${WORKDIR}

${P2Y_TARGET}: ${P2Y_OBJS}
	${CC} ${CFLAGS} ${P2Y_OBJS} ${P2Y_LDFLAGS} -o $@

${Y2P_TARGET}: ${Y2P_OBJS}
	${CC} ${CFLAGS} ${Y2P_OBJS} ${Y2P_LDFLAGS} -o $@

${PIM_TARGET}: ${PIM_OBJS}
	${CC} ${CFLAGS} ${PIM_OBJS} ${PIM_LDFLAGS} -o $@

${IIM_TARGET}: ${IIM_OBJS}
	${CC} ${CFLAGS} ${IIM_OBJS} ${IIM_LDFLAGS} -o $@

${IIX_TARGET}: ${IIX_OBJS}
	${CC} ${CFLAGS} ${IIX_OBJS} ${IIX_LDFLAGS} -o $@

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

${WORKDIR}/%.d: %.c
	mkdir -p ${dir $@}
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.o} > $@
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT ${@:%.d=%.s} >> $@
	${MAKEDEPEND} ${CINCLUDE} ${CFLAGS} $< -MT $@ >> $@
${WORKDIR}/%.s: %.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -S -o $@ $<
${WORKDIR}/%.o: %.c
	mkdir -p ${dir $@}
	${CC} ${CINCLUDE} ${CFLAGS} -c -o $@ $<
