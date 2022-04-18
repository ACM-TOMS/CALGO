### GEMM-Based Level 3 BLAS Benchmark ##################################

all: tmg timing

clean: cleantmg cleantiming

tmg:
	( cd TMGLIB; $(MAKE) )

timing:
	( cd SBENCH; $(MAKE) )
	( cd DBENCH; $(MAKE) )
	( cd CBENCH; $(MAKE) )
	( cd ZBENCH; $(MAKE) )

cleantmg:
	( cd TMGLIB; $(MAKE) clean )

cleantiming:
	( cd SBENCH; $(MAKE) clean )
	( cd DBENCH; $(MAKE) clean )
	( cd CBENCH; $(MAKE) clean )
	( cd ZBENCH; $(MAKE) clean )
