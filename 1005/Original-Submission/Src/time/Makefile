ifneq ($(TOPLEVEL),true)
  include ../make.inc
  DIR := .
endif

TIME := timermd.f90 time_auxil.f90 $(NUMTHRD).f90
TIMESRC := $(addprefix $(DIR)/,$(TIME))
TIMEOBJ := $(TIMESRC:.f90=$(OEXT))
TIMEEXE := $(DIR)/timermd$(EEXT)
TIMEDIR := $(DIR)/

time: $(TIMEEXE)

timerun: time
	@./$(TIMEEXE) -n100 -s2

$(DIR)/time_auxil$(OEXT): $(DIR)/timermd$(OEXT)
$(DIR)/timermd$(OEXT): $(TOOLS)/printutil$(OEXT)

$(TIMEEXE): $(LIB)/libblasrmd.a $(TOOLS)/libtools.a $(BLASDEP)

$(TIMEOBJ) : %$(OEXT) : %.f90
	$(FC) $(FOPTS) -I$(TOOLS) -I$(TIMEDIR) -c $< -o $@ $(MOD) $(TIMEDIR)

$(TIMEEXE) : $(TIMEOBJ)
	$(FC) $(LOPTS) -L$(LIB) -L$(TOOLS) -o $@ $^ -lblasrmd -ltools $(BLASLIB)

clean::
	@rm -f $(TIMEEXE) $(TIMEDIR)/*$(OEXT) $(TIMEDIR)/*.mod

.PHONY: time clean
