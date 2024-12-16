LATEX=pdflatex
FLAGS=-interaction=nonstopmode -halt-on-error -file-line-error
TARGET=main

.PHONY: ${TARGET}.pdf
${TARGET}.pdf:
	${LATEX} ${FLAGS} ${TARGET}.tex
	bibtex ${TARGET}
	${LATEX} ${FLAGS} ${TARGET}.tex
	${LATEX} ${FLAGS} ${TARGET}.tex

.PHONY: ${TARGET}.docx
${TARGET}.docx:
	pandoc -o $@ --verbose --citeproc --bibliography=ref.bib ${TARGET}.tex

clean:
	rm -f ${TARGET}.pdf
	rm -f *.dvi
	rm -f *.aux
	rm -f *.tex.backup
	rm -f *.loc
	rm -f *.soc
	rm -f *.lof
	rm -f *.log
	rm -f *.lot
	rm -f *.toc
	rm -f *.blg
	rm -f *.bbl
	rm -f *.bak
	rm -f *.out
	rm -f *.tex~
	rm -f *.cut
