####################################################################
#  SHHEIG usersguide makefile                                      #
#  Makefile for creating/updating the SHHEIG users's guide on Unix # 
#  machines.                                                       #
#  SHHEIG, Release 1.0                      ./shheig/src/makefile  #
#  Matthias Voigt, MPI Magdeburg                                   #
#  July 23, 2014.                                                  #
####################################################################
#
#  This is the makefile to create/update the users' guide for SHHEIG.
#
#  The command
#       make
#  without any arguments creates or updates a PDF file called
#       userManual.pdf.
#
#  To remove auxiliary files after the PDF is created, enter
#       make clean.
#
#######################################################################

usersguide:
	pdflatex userManual.tex
	bibtex userManual.aux
	pdflatex userManual.tex
	pdflatex userManual.tex

clean:
	rm -f userManual.aux
	rm -f userManual.bbl
	rm -f userManual.blg
	rm -f userManual.log
	rm -f userManual.out
