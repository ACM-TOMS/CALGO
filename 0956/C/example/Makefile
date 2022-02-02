################################################################################
# Macro definitions
################################################################################

TOPLEVEL = ..
include ${TOPLEVEL}/Make.config  # Parse for generic macros
include ./driver.config

CSRC = $(wildcard *.c)
OBJS = $(notdir $(CSRC:.c=.o))

.PHONY: all run trees clean
################################################################################
# Targets
################################################################################
all: driver

run: driver
	${MPIEXEC} -n $(NPROCS) $< ${PARAMETERFILE}

driver: $(OBJS)
	${MPICC} $(LDPATHS) -o $@ $(OBJS) $(LDLIBS)

trees:
	cd ${DATADIR}; make trees

clean:
	${RM} -f *~ *.o driver
	cd ${DATADIR}; make clean
