

open: main.pdf
	open $<

main.pdf: main.tex
	pdflatex $<
	bibtex main
	pdflatex $<
	pdflatex $<
