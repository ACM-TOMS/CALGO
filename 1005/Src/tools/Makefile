ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR:=.
endif

STOOLS     := $(wildcard $(DIR)/rmd_s*.f90)
DTOOLS     := $(subst rmd_s,rmd_d,$(STOOLS))
TOOLSRC    := $(STOOLS) $(DTOOLS) $(DIR)/printutil.f90
TOOLOBJ    := $(TOOLSRC:.f90=$(OEXT))
TOOLLIB    := $(DIR)/libtools.a
TESTFOBJ   := $(DIR)/rmd_stesttools$(OEXT) $(DIR)/rmd_dtesttools$(OEXT)

tools: $(TOOLLIB)

$(DTOOLS): $(DIR)/rmd_d%: $(DIR)/rmd_s%
	python $(TOOLS)/rmd-convert-type.py $< -o $@

$(TOOLOBJ): %$(OEXT): %.f90
	$(FC) $(FOPTS) -I $(TOOLS) -c $< -o $@ $(MOD) $(TOOLS)

$(TOOLLIB): $(TOOLOBJ)
	$(AR) -cr $@ $^
	ranlib $@

clean::
	rm -f $(DTOOLS) $(TOOLS)/*$(OEXT) $(TOOLS)/*.mod $(TOOLLIB)

.PHONY: tools clean
