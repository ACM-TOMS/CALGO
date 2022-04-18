ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR := .
endif
TESTDIR := $(DIR)
MORE    := $(DIR)/more
#DISPMODULE := $(TOOLS)/dispmodule.o

TESTBLAS := srotmg.f90
ALLRMD   := $(notdir $(wildcard $(SRC)/s*_rmd.f90))
ALLTEST  := $(ALLRMD) $(TESTBLAS)
ALLMORE  := $(filter-out srot srotg srotm srotmg sdsdot,$(ALLRMD:%_rmd.f90=%))
STEST    := $(addprefix $(DIR)/test_,$(ALLTEST))
DTEST    := $(subst test_s,test_d,$(STEST))
DTEST    := $(filter-out %/test_ddsdot_rmd.f90,$(DTEST))
TEST     := $(STEST) $(DTEST)
MORETEST := $(addprefix $(MORE)/test_more_,$(ALLMORE))
MORETEST := $(MORETEST:%=%.f90)
TESTEXE     := $(TEST:%.f90=%$(EEXT))
STESTEXE    := $(STEST:%.f90=%$(EEXT))
DTESTEXE    := $(DTEST:%.f90=%$(EEXT))
MORETESTEXE := $(MORETEST:%.f90=%$(EEXT))
TESTLIBDEP  := $(TOOLS)/libtools.a $(LIB)/libblasrmd.a $(BLASDEP)
TESTLIBDIR  := -L$(TOOLS) -L$(LIB)
TESTLIBS    := -ltools -lblasrmd $(BLASLIB)

BUILDMORE   := $(TOOLS)/build-more-tests.awk
MOREPAR     := $(TOOLS)/test-more-param.txt
MORETEMPL   := $(TOOLS)/test-more-template.txt

test: $(TESTEXE) $(MORETESTEXE) # exe files for all tests incl. more
stest: $(STESTEXE) # base single precision tests
dtest: $(DTESTEXE) # base double precision tests
moretest: $(MORETESTEXE)

$(DTEST): $(DIR)/test_d%: $(DIR)/test_s% # .f90 files with dp tests
	python $(TOOLS)/rmd-convert-type.py $< -o $@

$(MORETEST): $(BUILDMORE) $(MOREPAR) $(MORETEMPL) # .f90 files with more tests
$(MORETEST): $(MORE)/test_more_%.f90: $(SRC)/%_rmd.f90
	@mkdir -p $(MORE)
	awk -f $(BUILDMORE) -v ofile=$@ -v idir=$(TOOLS)

$(TESTEXE): %$(EEXT): %.f90 $(TESTLIBDEP)
	$(FC) $(FOPTSTEST) $(TESTLIBDIR) -I$(TOOLS) $< -o $@ $(TESTLIBS)

$(MORETESTEXE): %$(EEXT): %.f90 $(TESTLIBDEP)
	$(FC) $(FOPTSTEST) $(TESTLIBDIR) -I$(TOOLS) $< -o $@ $(TESTLIBS)

basetestrun: stestrun dtestrun
testrun: basetestrun moretestrun

stestrun: $(STESTEXE)
	for t in $(STESTEXE); do \
	   echo $$t; \
	   ./$$t || exit 1;\
	done

dtestrun: $(DTESTEXE)
	for t in $(DTESTEXE); do \
	   echo $$t; \
	   ./$$t || exit 1;\
	done
moretestrun: $(MORETESTEXE)
	for t in $(MORETESTEXE); do \
	   echo $$t; \
	   ./$$t || exit 1;\
	done

clean::
	rm -f $(TESTEXE) $(DTEST)
	rm -f $(TESTDIR)/*.mod
	rm -rf $(MORE)
	rm -rf $(TESTDIR)/*.dSYM

.PHONY: test clean testrun testmore
