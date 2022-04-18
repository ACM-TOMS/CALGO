ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR := .
endif

DOCDIR     := $(DIR)
EXTRACTDOC := awk -f $(TOOLS)/extractdoc.awk

L1BLAS     := rotg rotmg rot rotm swap scal copy axpy dot dsdot nrm2 asum
L2BLAS     := gemv gbmv symv sbmv spmv trmv tbmv tpmv trsv tbsv tpsv \
              ger syr spr syr2 spr2
L3BLAS     := gemm symm syrk syr2k trmm trsm
RMDSCA     := scal axpy gemv gbmv symv sbmv spmv ger syr spr syr2 spr2 \
              gemm symm syrk syr2k trmm trsm
OTHER      := potrf

SCABLAS    := $(RMDSCA:%=%-scalars)

DOCL1      := $(L1BLAS:%=%.txt)
DOCL2      := $(L2BLAS:%=%.txt)
DOCL3      := $(L3BLAS:%=%.txt)
DOCSC      := $(SCABLAS:%=%.txt)
DOCOTHER   := $(OTHER:%=%.txt)

ALLTMP     := $(addprefix $(DIR)/,$(DOCL1) $(DOCL2) $(DOCL3) $(DOCOTHER))
ALLTMPS    := $(addprefix $(DIR)/,$(DOCSC))

REFDOC     := $(DIR)/refdoc.tex
REFTEX     := $(DIR)/reference-manual.tex
USERTEX    := $(DIR)/user-manual.tex
REFPDF     := $(REFTEX:%.tex=%.pdf)
USERPDF    := $(USERTEX:%.tex=%.pdf)

doc: $(REFPDF) $(USERPDF)

$(ALLTMP): $(DIR)/%.txt: $(SRC)/s%_rmd.f90
	$(EXTRACTDOC) $< > $@
$(ALLTMPS): $(DIR)/%-scalars.txt: $(SRC)/s%_rmds.f90
	$(EXTRACTDOC) $< > $@

$(REFDOC): $(ALLTMP) $(ALLTMPS)
	echo $(L1BLAS) | $(TOOLS)/build-ref-man.sh 1 > $(REFDOC)
	echo $(L2BLAS) | $(TOOLS)/build-ref-man.sh 2 >> $(REFDOC)
	echo $(L3BLAS) | $(TOOLS)/build-ref-man.sh 3 >> $(REFDOC)
	echo $(SCABLAS)| $(TOOLS)/build-ref-man.sh 4 >> $(REFDOC)
	echo $(OTHER)  | $(TOOLS)/build-ref-man.sh other >> $(REFDOC)

$(REFPDF): $(REFTEX) $(REFDOC)
	cd $(DOCDIR) && latexmk -pdf $(notdir $<)

$(USERPDF): $(USERTEX) $(USERDOC)
	cd $(DOCDIR) && latexmk -pdf $(notdir $<)

clean::
	rm -f $(ALLTMP) $(ALLTMPS) $(DOCDIR)/{refdoc.tex,*.aux,*.log,*.out,*.fls,*.toc}
	rm -f $(DOCDIR)/{*.dvi,*.fdb_latexmk,*.synctex.gz}

clean-all:: clean
	rm -f $(REFPDF) $(USERPDF)

.PHONY: doc clean clean-all
