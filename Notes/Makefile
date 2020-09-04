T=1

PAPER=filename

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex
	 latexmk -pvc -pvctimeout -pvctimeoutmins=$(T) -pdf $(PAPER).tex
	 latexmk -c -bibtex
