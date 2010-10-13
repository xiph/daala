# Makefile to generate the paper from its sources.
# Requires transfig and pdflatex.

# Alias some shell commands, so we might have a chance of running on non-Unix
#  platforms some day.
RM = rm -f

PAPER_SRCS = mc.tex daala.bib

FIG_SRCS = 4-8sub.fig mvpred.fig
FIG_TEX_SRCS = blockidx.fig vbunsplit.fig

FIG_AUXS =

FIG_OBJS = $(FIG_SRCS:.fig=.pdf) \
 $(FIG_TEX_SRCS:.fig=.pdftex) $(FIG_TEX_SRCS:.fig=.pdftex_t)

EXTRA_FIGS =

all: mc.pdf

mc.pdf : $(PAPER_SRCS) $(FIG_OBJS) $(EXTRA_FIGS) mc.bbl
	pdflatex mc
	pdflatex mc

mc.bbl : $(PAPER_SRCS)
	if [ mc.tex -nt mc.aux ] ; then pdflatex mc ; fi
	bibtex mc

figures : $(FIG_OBJS)

# rule to generate latex versions of the xfig figures
%.tex : %.fig
	fig2dev -L latex $< $@

%.pdf : %.fig
	fig2dev -L pdf -p 0 $< $@

%.pdftex : %.fig
	fig2dev -L pdftex -p 0 $< $@

%.pdftex_t : %.fig
	fig2dev -L pdftex_t -p $(<:.fig=.pdftex) $< $@

.PHONY: all figures clean distclean maintainer-clean

# clean targets
clean:
	-$(RM) mc.pdf
	-$(RM) $(FIG_OBJS)
	-$(RM) $(FIG_AUXS)
	-$(RM) mc.aux
	-$(RM) mc.log
	-$(RM) mc.lof
	-$(RM) mc.lot
	-$(RM) mc.out
	-$(RM) mc.bbl
	-$(RM) mc.blg
	-$(RM) mc.toc

distclean: clean

maintainer-clean: distclean

maintainerclean: maintainer-clean
