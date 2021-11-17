ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR:=.
endif

ALLSRC  := $(wildcard $(DIR)/*.f90)
ALLOBJ  := $(ALLSRC:%.f90=%$(OEXT))
DEMODEP := $(LIB)/libblasrmd.a $(TOOLS)/libtools.a $(BLASDEP)
DEMODIR := $(DIR)

DEMOUTIL := $(TOOLS)/printutil$(OEXT)
DEMO1SRC := $(DIR)/demo1_sumsq.f90
DEMO1OBJ := $(DEMO1SRC:%.f90=%$(OEXT))
TEST1OBJ := $(DIR)/test_demo1$(OEXT)
TEST1EXE := $(DIR)/test_demo1$(EEXT)
DRIV1OBJ := $(DIR)/demo1_driver$(OEXT)
DRIV1EXE := $(DIR)/demo1_driver$(EEXT)

TEST2OBJ  := $(DIR)/test_demo2$(OEXT)
TEST2EXE  := $(DIR)/test_demo2$(EEXT)
DRIV2OBJ  := $(DIR)/demo2_driver$(OEXT)
DRIV2EXE  := $(DIR)/demo2_driver$(EEXT)

DEMO2SRC  := $(DIR)/demo2_likeli.f90
DEMO2OBJ  := $(DEMO2SRC:%.f90=%$(OEXT))

ifneq ($(slicot),)
  SB03MD   := $(DIR)/sb03md.f
  DLYAP    := $(DIR)/dlyap_slicot.f90
  DEMO2OBJ += $(SB03MD:%.f=%$(OEXT))
else  
  DLYAP := $(DIR)/dlyap.f90
endif

DEMO2OBJ += $(DLYAP:%.f90=%$(OEXT))
TESTOBJ := $(TEST1OBJ) $(TEST2OBJ)
TESTEXE := $(TEST1EXE) $(TEST2EXE) $(TESTLYAPEXE)
ALLEXE  := $(TESTEXE) $(DRIV1EXE) $(DRIV2EXE)

#cat slicot/{sb03md,mb01rd,sb03mx,sb03my,sb03mv,sb03mw,sb04px,select}.f > sb03md.f

demo: $(ALLEXE)

$(ALLEXE): $(DEMODEP)

$(TEST1OBJ) $(TEST2OBJ): $(TOOLS)/rmd_dtesttools$(OEXT)

$(ALLOBJ): %$(OEXT): %.f90 $(DEMOUTIL)
	$(FC) -c $(FOPTS) -I$(TOOLS) $< -o $@

$(TEST1EXE): $(TEST1OBJ) $(DEMO1OBJ)
	$(FC) $(LOPTS) -L$(TOOLS) -L$(LIB) $^ -o $@ -ltools -lblasrmd $(BLASLIB)

$(TEST2EXE): $(TEST2OBJ) $(DEMO2OBJ)
	$(FC) $(LOPTS) -L$(TOOLS) -L$(LIB) $^ -o $@ -ltools -lblasrmd $(BLASLIB)

$(TESTLYAPEXE): $(DIR)/test_dlyap$(OEXT) $(DIR)/dlyap$(OEXT)
	$(FC) $(LOPTS) -L$(TOOLS) $^ -o $@ -ltools $(BLASLIB)

$(DRIV1EXE): $(DRIV1OBJ) $(DEMO1OBJ)
	$(FC) $(LOPTS) -L$(TOOLS) -L$(LIB) $^ -o $@ -ltools -lblasrmd $(BLASLIB)

$(DRIV2EXE): $(DRIV2OBJ) $(DEMO2OBJ)
	$(FC) $(LOPTS) -L$(TOOLS) -L$(LIB) $^ -o $@ -ltools -lblasrmd $(BLASLIB)

demorun: $(DRIV1EXE) $(DRIV2EXE)
	$(DRIV1EXE)
	$(DRIV2EXE) $(DEMODIR)/var_2_3.dat

testdemo: testlyap testdemo1 testdemo2

testdemo1: $(TEST1EXE)
	@printf '%s' 'running test_demo1...'
	@$(TEST1EXE)
	@echo OK

testdemo2: $(TEST2EXE)
	@printf '%s' 'running test_demo2...'
	@$(TEST2EXE) $(DEMODIR)/var_1_4.dat
	@$(TEST2EXE) $(DEMODIR)/var_2_3.dat
	@$(TEST2EXE) $(DEMODIR)/var_2_4.dat
	@$(TEST2EXE) $(DEMODIR)/var_2_40.dat
	@echo OK

testlyap: $(TESTLYAPEXE)
	@printf '%s' 'running test_dlyap...'
	@$(TESTLYAPEXE)
	@echo OK

clean::
	rm -f $(DEMODIR)/*$(OEXT) $(ALLEXE)

.PHONY:  demo demorun testdemo testdemo1 testdemo2 testlyap clean
