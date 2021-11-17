ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR := .
endif

REFBLASSRC := $(wildcard $(DIR)/*.f)
REFBLASOBJ := $(REFBLASSRC:.f=$(OEXT))
REFBLASLIB := $(LIB)/librefblas.a

refblas: $(REFBLASLIB)

$(REFBLASLIB): $(REFBLASOBJ)
	mkdir -p $(LIB)
	$(AR) -cr $@ $^
	ranlib $@

$(REFBLASOBJ): %$(OEXT): %.f
	$(FC) $(FOPTSREFBLAS) -c $< -o $@

clean::
	rm -f $(REFBLASOBJ) $(REFBLASLIB)

.PHONY: refblas clean
