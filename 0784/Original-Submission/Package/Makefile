#
#  Top level Makefile for the GEMM-based level 3 BLAS.
#  March, 1997.
#

include make.inc

all: install $(CRETMG) gblib $(CREULIB) testing timing

clean: cleaninstall cleangblib cleanunderlib cleantesting cleantiming

install:
	( cd INSTALL; $(MAKE); \
	  cp lsame.f ../GBL3B/SGBL3B; cp lsame.f ../GBL3B/DGBL3B; \
	  cp lsame.f ../GBL3B/CGBL3B; cp lsame.f ../GBL3B/ZGBL3B; \
	  cp lsame.f ../GBBEN/SBENCH; cp lsame.f ../GBBEN/DBENCH; \
	  cp lsame.f ../GBBEN/CBENCH; cp lsame.f ../GBBEN/ZBENCH; \
	  cp lsame.f ../UNDERLIB; )

$(CRETMG):
	( cd INSTALL; \
	  cp second.f ../GBBEN/TMGLIB; cp dsecnd.f ../GBBEN/TMGLIB )

gblib:
	( cd GBL3B; $(MAKE) )

$(CREULIB):
	( cd UNDERLIB; $(MAKE) )

testing:
	( cd TESTING; $(MAKE) )

timing:
	( cd GBBEN; $(MAKE) )

cleaninstall:
	( cd INSTALL; $(MAKE) clean )

cleangblib:
	( cd GBL3B; $(MAKE) clean )

cleanunderlib:
	( cd UNDERLIB; $(MAKE) clean )

cleantesting:
	( cd TESTING; $(MAKE) clean )

cleantiming:
	( cd GBBEN; $(MAKE) clean )
