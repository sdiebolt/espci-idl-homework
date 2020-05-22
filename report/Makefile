# Makefile to create PDF documents from LaTeX files
# Author: Samuel Diebolt <samuel@diebolt.io>
# No license

DOC=main

COMPILER=/Library/TeX/texbin/pdflatex
ARGS=--shell-escape
BIBTEX=/Library/TeX/texbin/bibtex

RM=/bin/rm

$(DOC).pdf: $(DOC).bbl $(DOC).tex
	$(COMPILER) $(ARGS) $(DOC).tex -draftmode
	$(COMPILER) $(ARGS) $(DOC).tex 

$(DOC).bbl: $(DOC).aux
	$(BIBTEX) $(DOC).aux

$(DOC).aux: $(DOC).bib
	$(COMPILER) $(ARGS) $(DOC).tex -draftmode
	$(COMPILER) $(ARGS) $(DOC).tex -draftmode

.PHONY: all clean cleanall

all: $(DOC).pdf

clean:
	$(RM) -f *.tox *.aux *.log *.out *.glg *.gls *.glo *.ist *.log *.toc *.acn *.acr *.alg *.bbl *.blg

cleanall: clean
	$(RM) -f *.pdf
