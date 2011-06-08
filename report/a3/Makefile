TEX=pdflatex
OBJ=document.pdf
TRASH=$(shell ls *.aux *.idx *.log *.toc *.ilg *.ind *.backup *.nav *.out *.snm 2>/dev/null)

.SUFFIXES : .tex .pdf

.tex.pdf: 
	$(TEX) $<
	$(TEX) $<

all:${OBJ}


view:all
	evince  ${OBJ}&

clean:
	rm -f $(OBJ)
	rm -f $(TRASH)
test:view
