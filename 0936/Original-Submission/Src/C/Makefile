# Time-stamp:  <2016-07-21 18:03:03 m>
# # # This make file by Fred T. Krogh, benefiting from lots of testing help from
# # # Philip W. Sharp and Richard J. Hanson.

all: CHECK
.PHONY : all
$(R) :
ifndef CC
	@echo; echo "You need to define the compiler (CC)";\
 echo;exit 100
else
	mkdir -p $(R);\
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

# Substitutions such as: $(MESS:_x=_$$x), replace every _x in the make variable
# MESS (defined as messy_x) with _$$x where x is the variable in the "do" that
# is part of a for loop, with the result that finally the $$x is replaced by the
# x which in our case is a letter in CPLETS.

$(LIBS)/libmathc.a : $(SRC)/cmessy.h $(R) Makefile $(SRC)/cmessycall_bod.F90\
  $(SRC)/smessycall_m.F90 $(SRC)/dmessycall_m.F90 $(SRC)/qmessycall_m.F90
	for x in $(CPLETS);\
 do\
   $(FC) $(FFLAGS) -DL_=d -c -o $(LIBS)/$$x"messycall_m.o"\
  $(SRC)/$$x"messycall_m.F90" $(MODU) $(INCL);\
  $(CC) $(CFLAGS) -Dplet_=\'$$x\' -Dqplet_=\"$$x\" -c\
  -o $(LIBS)/$$x"messy.o" $(SRC)/$$x"messy.c";\
 done;\
 cp -f $(SRC)/cmessy.h $(LIBS)/cmessy.h;\
 ar -r $(LIBS)/libmathc.a $(LIBS)/?messy*.o;\
 rm -f $(LIBS)?messy*.o

$(DRV)/dsample.o : $(DRV)/csample_bod.c $(DRV)/ssample.c $(DRV)/dsample.c\
  $(DRV)/qsample.c $(LIBS)/libmathc.a $(DRV)/csample.h
	for x in $(CPLETS);\
 do\
   $(CC) $(CFLAGS) -Dplet_=\'$$x\' -Dqplet_=\"$$x\" -c -o $(DRV)/$$x"sample.o"\
    $(DRV)/$$x"sample.c" $(INCL);\
 done

ifdef SPTEST
$(R)/result.g: $(DRV)/ctmessy.c $(DRV)/dsample.o $(LIBS)/cmessy.h
	$(CC) $(CFLAGS) -Dplet_=\'d\' -Dqplet_=\"d\" -c -o $(DRV)/ctmessy.o\
 $(DRV)/ctmessy.c $(INCL) $(INCLD);\
 $(FC) $(FFLAGS) -Dplet_=\'d\' -Dqplet_=\"d\" -o $(DRV)/ctmessy $(DRV)/ctmessy.o\
   $(DRV)/dsample.o $(LIBC) $(LIBF);\
 $(DRV)/ctmessy
# After this is should be possible to run "ddd Drivers/ctmessy" from the C
# directory
else
$(R)/result.g: $(DRV)/ctmessy.c $(DRV)/dsample.o $(LIBS)/cmessy.h
	for x in $(CPLETS);\
  do\
    $(CC) $(CFLAGS) -Dplet_=\'$$x\' -Dqplet_=\"$$x\" -c -o $(DRV)/xtest.o\
 $(DRV)/ctmessy.c $(INCL) $(INCLD);\
 $(FC) $(FFLAGS) -Dplet_=\'$$x\' -Dqplet_=\"$$x\" -o $(DRV)/xtest $(DRV)/xtest.o\
 $(DRV)/$$x"sample.o" $(LIBC) $(LIBF);\
 $(RUN) $(DRV)/xtest >$(R)/result.$$x;\
  done;\
  touch $(R)/result.g

DIFF=for y in $(CPLETS);\
 do\
   diff -C0 -p $(RO)/result.$$y $(R)/result.$$y;\
 done
endif

CHECK: $(R)/result.g
	@echo "###### Differences for C Codes are computed here:#####";\
  $(DIFF);echo; echo\
  "If differences indicated are insignifiant, then C tests Passed.";\
echo "Significant differences indicate that a test did NOT PASS"; echo;


$(DOC)/messy_doc.pdf: $(DOC)/messy_doc.tex $(DOC)/messy.bib $(DOC)/remark_936.tex
	pdflatex $(DOC)/messy_doc.tex; bibtex $(DOC)/messy_doc;\
  pdflatex $(DOC)/messy_doc.tex; pdflatex $(DOC)/remark_936.tex;\
  bibtex $(DOC)/remark_936; pdflatex $(DOC)/remark_936.tex; mv *.pdf $(DOC)

clean: # Used prior to zipping or to check after major changes.
	rm -rf  $(DRV)/NewResults $(DRV)/*.o $(LIBS)/libmathc.a $(R) $(DRV)/x*\
  $(SRC)/*.mod;\
  find ../ -name "*~" -delete;  find ../ -name "\#*" -delete
# -delete above can be replaced with -exec "{}" \;

replace_orig: $(DRV)/NewResults # Used if different results are to be expected
	rm -fv $(RO)/*; cp -fv $(R)/* $(RO)
