# Suffix for Fortran 90 programs
.SUFFIXES:
.SUFFIXES : .f90 $(SUFFIXES)

#-----------------------------------------------------------------------
# The following are platform dependent

# Name of Fortran 90  (or 95) compiler
FC = f95

#  Select your compiler flags
CFLAGS = -ieee=full -nan -C=all -gline
#CFLAGS = -nan -gline
#CFLAGS = -fast

# Suffix for Modules (created by the Fortran compiler)
MSUFF = mod

AR= ar

#-----------------------------------------------------------------------

.f90.o:
	$(FC) -c $(CFLAGS) $<

#-----------------------------------------------------------------------

OBJS1= internal_types.o ds_routines.o divide.o \
	rule_tn.o rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_cn.o  \
	rule_general.o rule_1.o region_processor.o volume.o \
	check.o global_all.o error_handling.o
OBJS2=buckley.o cui.o

MOD1= internal_types.$(MSUFF) ds_routines.$(MSUFF) divide.$(MSUFF) \
	rule_t3.$(MSUFF) rule_t2.$(MSUFF) rule_c2.$(MSUFF) rule_c3.$(MSUFF)  \
	rule_general.$(MSUFF) rule_1.$(MSUFF) region_processor.$(MSUFF) volume.$(MSUFF) \
	check.$(MSUFF) global_all.$(MSUFF) error_handling.$(MSUFF)

SRC=cui.f90 check.f90 global_all.f90 error_handling.f90 \
	buckley.f90 internal_types.f90 ds_routines.f90 divide.f90 \
	rule_tn.f90 rule_t3.f90 rule_t2.f90 rule_c2.f90 rule_c3.f90 rule_cn.f90 \
	rule_general.f90 rule_1.f90 region_processor.f90 volume.f90

#-----------------------------------------------------------------------

all: libcubpack.a

#-----------------------------------------------------------------------

# dependencies
rule_1.o: buckley.o
rule_c2.o: buckley.o
rule_t2.o: buckley.o
rule_c3.o: buckley.o
rule_t3.o: buckley.o
rule_tn.o: buckley.o
rule_cn.o: buckley.o
volume.o: buckley.o
divide.o: buckley.o internal_types.o
rule_general.o: rule_t3.o rule_t2.o rule_c2.o rule_c3.o rule_1.o \
	        internal_types.o
region_processor.o : divide.o rule_general.o
global_all.o: region_processor.o volume.o ds_routines.o internal_types.o
cui.o: global_all.o check.o error_handling.o ds_routines.o internal_types.o
internal_types.o: buckley.o

#-----------------------------------------------------------------------

libcubpack.a: $(OBJS1) $(OBJS2)
	$(AR) rvu libcubpack.a $(OBJS1) $(OBJS2)

veryclean:
	/bin/rm -f $(OBJS1) $(OBJS2) *.$(MSUFF) libcubpack.a 

clean:
	/bin/rm -f $(OBJS1)
