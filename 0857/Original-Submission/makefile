# Normally this file need not be edited. Architecture-dependent
# make commands are in makefile.in.

include makefile.in

.SUFFIXES: .f .f90 .o

OBJECTS		= polsys_glp.o lapack_glp.o main_template.o
TEST_OBJECTS    = polsys_glp.o lapack_glp.o test_install.o

# -----------------------------------------------------------------------------
# Targets.
# -----------------------------------------------------------------------------
default: polsys_glp

polsys_glp:	$(OBJECTS)
		$(LINKER) -o $@ $(OBJECTS)

test_install:   $(TEST_OBJECTS)
		$(LINKER) -o $@ $(TEST_OBJECTS)

.f90.o:
	$(F90_COMPILER) -o $@ $(F_FLAGS) $(F90_FLAGS) $<
.f.o:
	$(F77_COMPILER) -o $@ $(F_FLAGS) $<


# -----------------------------------------------------------------------------
# Clean-up.
# -----------------------------------------------------------------------------
clean:
	rm -f *.o *.mod polsys_glp test_install

