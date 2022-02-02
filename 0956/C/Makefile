################################################################################
# Macro definitions
################################################################################

TOPLEVEL = .
include ${TOPLEVEL}/Make.config  # Parse for generic macros
################################################################################
# Targets
################################################################################
.PHONY: all lib example clean remake

all:
	${MAKE} lib
	${MAKE} example

lib:
	cd src; ${MAKE} lib

example:
	cd example; ${MAKE} all

clean:
	rm -f *~
	cd src; ${MAKE} clean
	cd example; ${MAKE} clean

remake:
	${MAKE} clean
	${MAKE} all
################################################################################
