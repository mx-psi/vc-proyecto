## PANDOC VARIABLES
PANDOC:=pandoc # pandoc executable
#FILTERS:= filters/env.hs filters/pandoc-crossref pandoc-citeproc # in order
FILTERS:= filters/env.hs pandoc-citeproc # in order
PFLAGS:= $(foreach filter,$(FILTERS),-F $(filter))  # Pandoc flags


## OUTPUT-INDEPENDENT VARIABLES
SRCS:=memoria.md # Source files in order
DEPS:=src/citations.bib src/style.csl src/header.md src/footer.md

IMGS:= $(wildcard Implementacion/imagenes/*)

## PDF-SPECIFIC VARIABLES
OUTPDF:=memoria.pdf
DEPSPDF:=$(DEPS) src/latex.sty
PDFFLAGS:=-H src/latex.sty

.PHONY: all clean

all: proyecto.zip

$(OUTPDF): $(SRCS) $(DEPSPDF)
	$(PANDOC) $(PFLAGS) $(PDFFLAGS) src/header.md $(SRCS) src/footer.md -o $@

proyecto.zip: $(OUTPDF) Implementacion/main.py Implementacion/iterativo.py Implementacion/auxiliar.py $(IMGS)
	zip -9 $@ $^

clean:
	rm $(OUTPDF) proyecto.zip
