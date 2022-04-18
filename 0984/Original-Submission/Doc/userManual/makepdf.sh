#!/bin/bash

pdflatex ADiGatorUserGuide.tex
bibtex ADiGatorUserGuide
pdflatex ADiGatorUserGuide.tex
pdflatex ADiGatorUserGuide.tex