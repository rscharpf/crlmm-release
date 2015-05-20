All: vignettes clean

## all pdfs must be generated from the first target of the Makefile
vignettes: AffyGW.tex
	cp ../inst/scripts/AffyGW.pdf .
	cp ../inst/scripts/IlluminaPreprocessCN.pdf .
	cp ../inst/scripts/Infrastructure.pdf .
	cp ../inst/scripts/gtypeDownstream.pdf .
	cp ../inst/scripts/crlmmIllumina.pdf .
	"$(R_HOME)/bin/R" CMD texi2dvi --pdf genotyping.tex
	"$(R_HOME)/bin/R" CMD texi2dvi --pdf CopyNumberOverview.tex

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex
