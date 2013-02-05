
######################
######################
# makefile generated 
# dude makefile 
# created by s.e.pav 
# $Id: Makefile 89 2006-01-26 20:59:08Z spav $
######################
######################

############### FLAGS ###############

# these don't work at times
# in that case set these by hand?

LATEX       := $(shell which latex)
BIBTEX      := $(shell which bibtex)
PDFLATEX    := $(shell which pdflatex)
HTLATEX     := $(shell which htlatex.sh)
MAKEINDEX   := $(shell which makeindex)
GHOSTVIEW   := $(shell which gv)
PDF_VIEWER  := $(shell which okular)
TODAY       := $(shell date +'%Y%b%d%H%m')
TIME        := $(shell perl -e 'print time();')
PAGER   		:= $(shell which less)
TAR         := $(shell which tar)
ASPELL  		:= $(shell which aspell)
MAXIMA      := $(shell which maxima)

RLIB         = /usr/lib64/R

#TEXINPADD    = .:$(HOME)/sys/etc/tex:$(HOME)/sys/etc/tex/SEPtex:$(HOME)/work/math/TEX:$(RLIB)/share/texmf/tex/latex
TEXINPADD    = .:$(RLIB)/share/texmf/tex/latex

PRETEX       = TEXINPUTS=$(TEXINPADD):$$TEXINPUTS
PREBIB       = BSTINPUTS=$(TEXINPADD):$$BSTINPUTS \
               BIBINPUTS=$(TEXINPADD):$$BIBINPUTS 

PREIDX       = INDEXSTYLE=$(TEXINPADD):$$INDEXSTYLE

#undoes psfrag for pdf
UNPSFRAG		 = perl $(HOME)/sys/bin/unpsfrag.pl
#unroll commands
DETEXIFY		 = perl $(HOME)/sys/perl/detexify.pl

SCREEN_SIZE  = normal
include	$(HOME)/sys/etc/.Makefile.local

PROJECT      = dude
TEX_SOURCE   = $(PROJECT).tex
BIB_SOURCE   = $(PROJECT).bib
DVI_TARGET   = $(PROJECT).dvi
PS_TARGET    = $(PROJECT).ps 
PDF_TARGET   = $(PROJECT).pdf
BBLS         = $(PROJECT).bbl

#SAVE
# tracked projects
PROJECTS     = $(PROJECT) 
#UNSAVE
# add on dependencies (subchapters of dude)
R_DEPS 			 = defs.R
TEX_EXTRAS   = sharpe_shortcuts.tex mac_overfit.tex skew_study_show.tex
TEX_EXTRAS  += $(R_DEPS)
# nonlocal dependencies
STY_FILES    = 

# ACK!
DIR_STUB     = work
DOC_DIR      = $(DIR_STUB)/$(PROJECT)/DOC

# more worry
OWNER        = $(shell whoami)
THIS_HOST    = $(shell hostname)
CWD          = $(shell pwd)

# push versions of the pdf to math
SHOW_TODAY   = dude_show_$(TODAY).pdf

# tar junk
BACKUP       = $(PROJECT)_backup_$(TODAY)
TAR_FILE     = $(BACKUP).tar
TARBZ_FILE   = $(TAR_FILE).bz2
TAR_MAJORS   = Makefile $(TEX_SOURCE) $(BIB_SOURCE) $(TEX_EXTRAS) *.Rnw
TAR_MINORS   = *.ps *.eps *.fig *.pro *.sty *.R 

STUB         = $(HOME)/sys/STUB
EMAIL        = spav@ucsd.edu

# X geometry
X_GEOM_TINY  = -geometry 975x760+15+10
X_GEOM_NORM  = -geometry 1066x920+15+10
X_GEOM_WIDE  = -geometry 1230x920+15+10
X_GEOM_HUGE  = -geometry 1800x1450+15+10
X_GEOM_HUGE  = -geometry 1514x1165+1920+0
X_GEOM_HUGE  = -geometry 1834x1107+1920+0

# xdvi
XDVI_NORM    = -s 5
XDVI_HUGE    = -s 3

# now select
ifeq ($(SCREEN_SIZE),big)
X_GEOM_PORT  = $(X_GEOM_HUGE)
X_GEOM_LAND  = $(X_GEOM_WIDE)
XDVI_SCALE   = $(XDVI_HUGE)
else
X_GEOM_PORT  = $(X_GEOM_NORM)
X_GEOM_LAND  = $(X_GEOM_WIDE)
XDVI_SCALE   = $(XDVI_NORM)
endif

# dvips
#DVIPS_FLAGS  = -R -t letter
#add the -z to get hyperlinks
DVIPS_FLAGS  = -R -t letter -z
DVIEPS_FLAGS = -R -t letter -Pcmz -Pamz
DVIPDF_FLAGS = -R -t letter -Pcmz -Pamz

# xdvi
XDVI_FLAGS   = -safer $(XDVI_SCALE) -expert -hush $(X_GEOM_PORT) -paper us \
								-keep -postscript 0

# pdfviewer
XPDF_FLAGS   = $(X_GEOM_PORT) -z width -bg grey
OKULAR_FLAGS = 
PDFSEE_FLAGS = -z width -q -bg black -fullscreen
PDFSEE_FLAGS = -z width -q
PDFSEE_FLAGS = $(XPDF_FLAGS)
PDFSEE_FLAGS = $(OKULAR_FLAGS)

#aspell
ASPELL_FLAGS = 

#dot
DOT_FLAGS    = 

SHOW_PNG     = $(wildcard dude_show*.png)
SHOW_EPS     = $(patsubst %.png,%.eps,$(SHOW_PNG))


############## DEFAULT ##############

default : all

############## MARKERS ##############

.PHONY   : 
.SUFFIXES: .tex .bib .dvi .ps .pdf .eps
.PRECIOUS: %.dvi %.ps %.pdf %.jpg %.gif %.eps

############ BUILD RULES ############

## the old way using Sweave
#%.tex : %.Rnw $(R_DEPS)
		#R CMD Sweave $<

#%.R : %.Rnw
		#R CMD Stangle $<

%.tex : %.Rnw $(R_DEPS)
		Rscript -e 'require(knitr);knit("$<")'

%.R : %.Rnw
		Rscript -e 'require(knitr);knit("$<",tangle=TRUE)'

%.pdf : %.tex
		$(PRETEX) $(PDFLATEX) $<
		if grep Citation $*.log > /dev/null; then $(PREBIB) $(BIBTEX) $*; $(PRETEX) $(PDFLATEX) $*; fi

# compile and convert
%.dvi : %.tex $(STY_FILES) 
		$(PRETEX) $(LATEX) $<
		if grep Citation $*.log > /dev/null; then $(PREBIB) $(BIBTEX) $*; $(PRETEX) $(LATEX) $*; fi
		if grep Rerun $*.log > /dev/null; then $(PRETEX) $(LATEX) $*; fi
		if grep Citation $*.log > /dev/null; then $(PREBIB) $(BIBTEX) $*; $(PRETEX) $(LATEX) $*; fi
		if grep Rerun $*.log > /dev/null; then $(PRETEX) $(LATEX) $*; fi

%.txt : %.mx
		$(MAXIMA) -b $< > $@

%.dep.dot : %.tex %.dvi $(STY_FILES)
		echo digraph G \{ > $@;
		$(PRETEX) $(LATEX) $< | grep 'DEP:' | perl -pe 's/DEP://;s/\\hbox {}//g;' >> $@
		echo \} >> $@;

%.ps : %.dvi
		dvips $(DVIPS_FLAGS) -o $@ $< 

#%.pdf : %.dvi
		#dvips $(DVIPDF_FLAGS) -o $*.ps $< 
		#ps2pdf $*.ps
		#-rm -f $*.ps

# converting to pdf the hard way

%.pdftex : %.tex
		$(UNPSFRAG) $<
		cp _$< $@

%.pdf : %.pdftex
		$(PDFLATEX) $<

%.mergepdf :
		gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@ $(*)*.pdf

# roll commands from sty files

%.detex : %.tex
		@-if ! [ -r $@ ]; then cp $< $@; fi
		$(DETEXIFY) $@ >> $@;

# crappy ps from pdf.
%.craps : %.pdf
		a2ps -X'ascii' -M'Letter' -R --columns=1 -o $*.ps $<

# crappy ps from edp
%.craps : %.edp
		a2ps -X'ascii' -M'Letter' -R --columns=2 -o $*.ps $<

# reverse?
%.pdf : %.eps
		epstopdf $<

# crappy ascii anyone?
%.txt : %.ps
		ps2ascii $< > $@

# crappy html anyone?
%.html : %.tex
		$(HTLATEX) $< 
		-rm -f tex4ht.fls $*.4ct $*.4tc $*.lg $*.idv $*.xref $*.tmp

# tex extras
%.bbl : %.bib
		$(PREBIB) $(BIBTEX) $*

%.bbl : %.aux
		$(PREBIB) $(BIBTEX) $*

%.ind : %.idx
		$(PREIDX) $(MAKEINDEX) $*

# convert from xfig
%.eps : %.fig
		/bin/rm -f $@ 
		fig2dev -L ps $< > $@

# see http://electron.mit.edu/~gsteele/pdf/ 
%.eps : %.png
		convert $< $@

# make a graph
%.dot : %.flc
		perl $(HOME)/sys/perl/flowcharter.pl $< > $@

%.ps : %.dot
		dot $(DOT_FLAGS) -Tps $< -o $@

%.fig : %.dot
		dot $(DOT_FLAGS) -Tfig $< -o $@

# viewing targets
%.xsee : %.dvi
		xdvi $(XDVI_FLAGS) $*.dvi
		
%.pdfsee : %.pdf
		$(PDF_VIEWER) $(PDFSEE_FLAGS) $*.pdf

# check a document
%.chk : %.dup %.spell

# check spelling
%.spell : %.tex
		$(ASPELL) $(ASPELL_FLAGS) --dont-tex-check-comments -t -l < $< | sort | uniq | $(PAGER)

# check duplicate words
%.dup : %.tex
		perl -an -F/\\s+/ -e 'BEGIN { $$last = q[]; $$line = 0; $$prevline = q[];}\
		$$line++;$$first = 1;\
		foreach $$word (@F) {\
		if ($$word eq $$last) {\
		if ($$first) { print qq[duplicate $$word, lines ],($$line-1),qq[-$$line:\n$$prevline$$_]; }\
		else { print qq[duplicate $$word, line $$line:\n$$_]; } }\
		$$last = $$word; $$first = 0; } \
		$$prevline = $$_;' < $< | $(PAGER)

# put into a tar file
%.tar : %.manif
		@ $(TAR) -cvf $@ -T $< $< Makefile

%.tgz : %.manif
		@ $(TAR) -czvf $@ -T $< $< Makefile


############# CLEAN UP ##############

# clean up
%.clean : 
		-rm -f $*.aux $*.log $*.dvi $*.bbl $*.blg $*.toc $*.ilg $*.ind
		-rm -f $*.out $*.idx $*.lot $*.lof $*.brf $*.nav $*.snm
%.realclean : %.clean
		-rm -f $*.ps $*.pdf
		-rm -f $*-[0-9][0-9][0-9]*.eps $*-[0-9][0-9][0-9]*.pdf


############### RULES ###############

# an easy target
all: $(DVI_TARGET) $(BBLS) figure/Curve_fitting.eps 

ratarb: ratarb.dvi $(BBLS)


dvis: $(patsubst %,%.dvi,$(PROJECTS))
pss: $(patsubst %,%.ps,$(PROJECTS))
pdfs: $(patsubst %,%.pdf,$(PROJECTS))

# targets
$(DVI_TARGET): $(TEX_SOURCE) $(TEX_EXTRAS) $(STY_FILES)
#$(PDF_TARGET): $(PS_TARGET)

hints :
	@-grep slides Makefile
	@-grep today Makefile

# see them 
presee: $(PROJECT).xsee
xsee: $(PROJECT).xsee
see: $(PROJECT).psee
pdfsee: $(PROJECT).pdfsee

# check it

spell: $(PROJECT).spell 

# clean up
clean: $(patsubst %,%.clean,$(PROJECTS))
realclean: $(patsubst %,%.realclean,$(PROJECTS))
		-rm -f Rplots.pdf

cleancache: 
		echo "killing knitr cache! ack!"
		-rm -rf cache

superclean: realclean cleancache

safe: realclean
		@$(TAR) -cvf $(TAR_FILE) $(TAR_MAJORS) $(SUB_FIGS)
		@ -$(TAR) -uvf $(TAR_FILE) $(DVI_TARGET)
		@ -$(TAR) -uvf $(TAR_FILE) $(TAR_MINORS) 
		@ -$(TAR) -tvf $(TAR_FILE)
		rm -f $(TARBZ_FILE)
		bzip2 $(TAR_FILE)

# see http://www.latex-community.org/forum/viewtopic.php?f=5&t=2070 
# sigh
show_extras: $(SHOW_EPS)

.tags :
	nice -n 18 ctags -f .tmp_tags --recurse --language-force=R --fields=+i `find . -regextype posix-egrep -regex '.*.R(nw)?'`;
	mv .tmp_tags $@

############# ONE OFFS ##############

figure/Curve_fitting.jpg :
	#wget -O $@ "http://www.automatedtrader.net/Files/images/Mags/Issue_04_Jan_07/Strategies/Curve%20fitting.jpg"
	wget -O $@ "www.imgur.com/5oOkkSR.jpg"

figure/Curve_fitting.eps : figure/Curve_fitting.jpg
	convert $^ $@

######################
######################

######################
######################
# remake the makefile from the stub directory
makefile:
		cp Makefile .oldMakefile;
		perl -pe 's/\s+PROJ_NAME/ $(PROJECT)/;s{\s+CWD_NAME}{ $(DIR_STUB)};' < $(STUB)/Makefile > Makefile;
		-diff Makefile .oldMakefile;
		-rm -i .oldMakefile;

#for vim modeline: (do not edit)
# vim:ts=2:sw=2:tw=79:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:tags=tags;:syn=make:ft=make:ai:si:cin:nu:fo=croqt:cino=p0t0c5(0:
