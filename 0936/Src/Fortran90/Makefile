#Time-stamp:  <2016-07-21 18:02:31 m>
# # # This make file by Fred T. Krogh, benefiting from lots of testing help from
# # # Philip W. Sharp and Richard J. Hanson.

# Variables are defined in a parent makefile, ..Makefile
all: CHECK
.PHONY : all
$(LIBS) :
ifndef FC
	@echo; echo "ERROR ##### You need to define the compiler (FC)";\
  echo;exit -100
else
	mkdir -p $(R);\
 mkdir -p $(LIBS);\
 $(ENVOUT)
endif

# Some notes on make syntax:
#  for loops have the form
# for <var> in <something>
#   do
#     <stuff(var)>
#   done
# In bash the <var> is preceded with a single $ when inside <stuff(var)>, but
# here it must be preceded by $$.

# make1, make2, and make3 are empty files used to tell whether something
# needs recompiling.  By using these any combination of s, d, and q for PLETS
# should compile and test only what is specified.

# Substitutions such as: $(MESS:_x=_$$x), replace every _x in the make variable
# MESS (defined as messy_x) with _$$x where x is the variable in the "do" that
# is part of a for loop, with the result that finally the $$x is replaced by the
# x which in our case is a letter in PLETS.

MESS=messy_x
SRC_FILES=$(SRC)/messy_gen_m.F90 $(SRC)/messy_bod.F90 \
 $(SRC)/messy_s_m.F90 $(SRC)/messy_d_m.F90 $(SRC)/messy_q_m.F90
$(LIBS)/precision_m.o : $(SRC)/precision_m.F90 $(LIBS) Makefile
	$(FC) $(FFLAGS) -c -o $(LIBS)/precision_m.o $(SRC)/precision_m.F90 \
  $(MODU)
$(LIBS)/libmathf.a : $(LIBS)/precision_m.o $(SRC_FILES)
	for x in $(PLETS);\
  do\
    $(FC) $(FFLAGS) -c -o $(LIBS)/$(MESS:_x=_$$x)"_m.o" \
    $(SRC)/$(MESS:_x=_$$x)"_m.F90" $(MODU) $(INCL);\
  done;\
  $(FC) $(FFLAGS) -c -o $(LIBS)/messy_gen_m.o $(SRC)/messy_gen_m.F90\
   $(MODU) $(INCL);\
    ar -r $(LIBS)/libmathf.a $(LIBS)/*.o;\
    rm $(LIBS)/*.o

SAMP=sample_x
$(DRV)/sample_gen_m.o : $(LIBS)/libmathf.a $(DRV)/sample_gen_m.F90 \
  $(DRV)/sample_s_m.F90 $(DRV)/sample_d_m.F90 $(DRV)/sample_q_m.F90 \
  $(DRV)/sample_bod.F90
	for x in $(PLETS);\
 do\
   $(FC) $(FFLAGS) -c -o $(DRV)/$(SAMP:_x=_$$x)"_m.o" \
   $(DRV)/$(SAMP:_x=_$$x)"_m.F90" $(INCL) $(INCLD) $(MODUD);\
 done;\
 $(FC) $(FFLAGS) -c -o $(DRV)/sample_gen_m.o $(DRV)/sample_gen_m.F90 $(INCL)\
  $(INCLD) $(MODUD)



$(R)/result.g : $(DRV)/sample_gen_m.o $(DRV)/tmessy.F90\
   $(DRV)/tmessyg.F90
	for x in $(PLETS);\
 do\
   $(FC) $(FFLAGS) -Dmessy_=messy_$$x"_m" -Dsample_=sample_$$x"_m" -o $(DRV)/xtest\
    $(DRV)/tmessy.F90 $(DRV)/sample_$$x"_m.o" $(INCL) $(INCLD) $(LIBF);\
   $(RUN) $(DRV)/xtest >$(R)/result.$$x;\
   $(FC) $(FFLAGS) -Dplet_=\"$$x\" -Dmessy_=messy_$$x"_m" -o $(DRV)/xtest\
     $(DRV)/thrdtmessy.F90 $(INCL) $(INCLD) $(LIBF);\
     $(DRV)/xtest;\
 done;\
 $(FC) $(FFLAGS) -o $(DRV)/xtest $(DRV)/tmessyg.F90 $(INCL) $(LIBF);\
 $(DRV)/xtest  >$(R)/result.g

DIFF2=
ifneq ($(NUMT),1)
  DIFF2= for x in $(RO)/*1.$$y;\
   do\
     diff -C0 -w -p $$x $(R)/`basename $$x`;\
   done;
endif

DIFF=for y in $(PLETS);\
 do\
   diff -C0 -p $(RO)/result.$$y $(R)/result.$$y;\
   for x in $(RO)/*0.$$y;\
   do\
     diff -C0 -p $$x $(R)/`basename $$x`;\
   done;\
   $(DIFF2)\
 done;\
 diff -C0 -p $(RO)/result.g $(R)/result.g


CHECK: $(R)/result.g
	@echo "###### Differences for Fortran Codes are computed here:#####";\
  $(DIFF); echo; echo\
  "If differences indicated are insignifiant, then Fortran tests Passed.";\
echo "Significant differences indicate that a test did NOT PASS"; echo;

clean: # Used prior to zipping or to check after major changes.
	@rm -rfv $(LIBS) $(DRV)/NewResults $(R) $(DRV)/*.o $(DRV)/*.mod \
   *.mod $(DRV)/x* $(DOC)/auto; \
 for x in aux bbl blg dvi log out toc;do rm -vf  $(DOC)/*.$$x;\
 done;\
 find ../ -name "*~" -delete # -delete can be replaced with -exec "{}" \;

replace_orig: $(DRV)/NewResults # Used if different results are to be expected
	rm -fv $(DRV)/Results/*; cp -fv $(DRV)/NewResults/* $(DRV)/Results
