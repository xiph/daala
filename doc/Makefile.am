# Makefile to generate the paper from its sources.
# Requires transfig and pdflatex.


PAPER_SRCS = mc.tex daala.bib

FIG_SRCS = 4-8sub.fig mvpred.fig
FIG_TEX_SRCS = blockidx.fig vbunsplit.fig

FIG_AUXS =

FIG_OBJS = $(FIG_SRCS:.fig=.pdf) \
 $(FIG_TEX_SRCS:.fig=.pdftex) $(FIG_TEX_SRCS:.fig=.pdftex_t)

EXTRA_FIGS =

noinst_DATA = mc.pdf

CLEANFILES = \
	mc.pdf \
	mc.out mc.log mc.blg mc.bbl mc.aux \
	blockidx.pdftex blockidx.pdftex_t \
	vbunsplit.pdftex vbunsplit.pdftex_t \
	mvpred.pdf 4-8sub.pdf

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

.PHONY: figures

