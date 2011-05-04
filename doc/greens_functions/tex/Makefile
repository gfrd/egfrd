# Todo. Use that Makefile.am magic. And why is it recompiling even though 
# nothing changed?

# make should produce:
#   * 1 pdf for the derivation of the Greens function of each diffusion 
#   problem. (Todo)
#   * 1 pdf for each problem giving the results: Greens function, survival 
#   probability and propensity function.
#   * greensFunctions.pdf, all results combined.


PDFLATEX=pdflatex -output-directory=../build

greensFunctions.pdf: greensFunctions.tex greensFunction1DAbsAbs.pdf greensFunction1DAbsInf.pdf greensFunction1DRadAbs.pdf greensFunction1DRadInf.pdf greensFunction2DRadAbs.pdf
	$(PDFLATEX) greensFunctions.tex

greensFunction1DAbsAbs.pdf: greensFunction1DAbsAbs.tex include/greensFunction1DAbsAbs.include.tex 
	$(PDFLATEX) greensFunction1DAbsAbs.tex 

greensFunction1DAbsInf.pdf: greensFunction1DAbsInf.tex include/greensFunction1DAbsInf.include.tex
	$(PDFLATEX) greensFunction1DAbsInf.tex

greensFunction1DRadAbs.pdf: greensFunction1DRadAbs.tex include/greensFunction1DRadAbs.include.tex
	$(PDFLATEX) greensFunction1DRadAbs.tex

greensFunction1DRadInf.pdf: greensFunction1DRadInf.tex include/greensFunction1DRadInf.include.tex
	$(PDFLATEX) greensFunction1DRadInf.tex

greensFunction2DRadAbs.pdf: greensFunction2DRadAbs.tex include/greensFunction2DRadAbs.include.tex
	$(PDFLATEX) greensFunction2DRadAbs.tex

