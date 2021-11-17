### GEMM-Based Level 3 BLAS ############################################

all: gblib

clean: cleangblib

gblib:
	( cd SGBL3B; $(MAKE) )
	( cd DGBL3B; $(MAKE) )
	( cd CGBL3B; $(MAKE) )
	( cd ZGBL3B; $(MAKE) )

cleangblib:
	( cd SGBL3B; $(MAKE) clean )
	( cd DGBL3B; $(MAKE) clean )
	( cd CGBL3B; $(MAKE) clean )
	( cd ZGBL3B; $(MAKE) clean )
