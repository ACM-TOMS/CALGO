ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR := .
endif

RMDLIB        := $(LIB)/libblasrmd.a
SHARED_RMDLIB := $(LIB)/libblasrmd.so

SSRC    := $(wildcard $(DIR)/s*_rmd.f90)
SSRC    += $(wildcard $(DIR)/s*_rmds.f90)
DSRC    := $(SSRC:$(DIR)/s%=$(DIR)/d%)
DSRC    := $(filter-out %/ddsdot_rmd.f90,$(DSRC))
RMDSRC  := $(SSRC) $(DSRC)
RMDOBJ  := $(RMDSRC:%.f90=%$(OEXT))

src: rmdlib # sharedlib

rmdlib: $(RMDLIB)

sharedlib: $(SHARED_RMDLIB)

$(DSRC): $(DIR)/d%.f90: $(DIR)/s%.f90
	python $(TOOLS)/rmd-convert-type.py $< -o $@

$(RMDLIB): $(RMDOBJ)
	mkdir -p $(LIB)
	$(AR) -cr $@ $^
	ranlib $@

$(RMDOBJ): %$(OEXT): %.f90
	$(FC) -c $(FOPTS) -I$(TOOLS) $< -o $@

$(SHARED_RMDLIB): $(RMDSRC)
	$(FC) -fpic -shared $(FOPTS) $(RMDSRC) -o $@ $(BLASDIR) $(BLASLIB)

clean::
	rm -f $(RMDOBJ) $(DSRC)
	rm -rf $(LIB)
	rm -rf $(DIR)/*.dSYM

.PHONY: src rmdlib sharedlib clean
