

open: main.pdf
	open $<

main.pdf: main.tex references.bib
	pdflatex $<
	bibtex main
	pdflatex $<
	pdflatex $<
