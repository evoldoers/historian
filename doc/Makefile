

open: main.pdf
	open $<

figure1.pdf: gp120tree.pdf
	cp $< $@

figure2.pdf: sim/results.pdf
	cp $< $@

main.pdf: main.tex references.bib figure1.pdf figure2.pdf
	pdflatex $<
	bibtex main
	pdflatex $<
	pdflatex $<
