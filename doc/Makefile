
all: main.pdf.open

preprint: preprint.pdf.open

%.open: %
	open $<

figure1.pdf: gp120tree.pdf
	cp $< $@

figure2.pdf: sim/results.pdf
	cp $< $@

%.pdf: %.tex paper.tex references.bib figure1.pdf figure2.pdf
	pdflatex $<
	bibtex $*
	pdflatex $<
	pdflatex $<
