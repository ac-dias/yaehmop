# Build the developer notes: make (needs pdflatex + latexmk); make clean
DOC = advances_in_yaehmop

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex ../new_tests/MoS2_SOC/out/bands_soc.png ../new_tests/Fe_spin/out/bands_spin.png
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(DOC).tex

clean:
	latexmk -c $(DOC).tex

.PHONY: all clean
