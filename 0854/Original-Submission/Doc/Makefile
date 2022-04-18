include make.inc

all: lib example

clean: cleanlib cleanexample

lib:
	( cd src; $(MAKE) )

example: 
	( cd examples; $(MAKE) )

cleanlib:
	( cd src; $(MAKE) clean )

cleanexample:
	( cd examples; $(MAKE) clean )

