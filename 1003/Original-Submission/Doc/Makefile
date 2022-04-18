userManual.pdf: userManual.tex references.bib \
    Figures/BrotherlyMatching.eps \
    Figures/CommunityMatching.eps \
    Figures/MultilevelOverview.eps \
    title-info.tex Makefile
	xelatex userManual
	xelatex userManual
	bibtex userManual
	xelatex userManual
	xelatex userManual
