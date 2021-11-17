#-------------------------------------------------------------------------------
# Makefile for CHOLMOD_ACM_TOMS packages
#-------------------------------------------------------------------------------

include UFconfig/UFconfig.mk

# Compile the default rules for each package
default:
	( cd UFconfig/xerbla ; $(MAKE) )
	( cd AMD ; $(MAKE) )
	( cd CAMD ; $(MAKE) )
	( cd COLAMD ; $(MAKE) )
	( cd CCOLAMD ; $(MAKE) )
	( cd CHOLMOD ; $(MAKE) )

library: default

# Compile the MATLAB mexFunctions
mex:
	( cd AMD ; $(MAKE) mex )
	( cd CAMD ; $(MAKE) mex )
	( cd COLAMD ; $(MAKE) mex )
	( cd CCOLAMD ; $(MAKE) mex )
	( cd CHOLMOD ; $(MAKE) mex )

# Remove all files not in the original distribution
purge:
	( cd UFconfig/xerbla ; $(MAKE) purge )
	( cd AMD ; $(MAKE) purge )
	( cd CAMD ; $(MAKE) purge )
	( cd COLAMD ; $(MAKE) purge )
	( cd CCOLAMD ; $(MAKE) purge )
	( cd CHOLMOD ; $(MAKE) purge )

# Remove all files not in the original distribution, but keep the libraries
clean:
	( cd UFconfig/xerbla ; $(MAKE) clean )
	( cd AMD ; $(MAKE) clean )
	( cd CAMD ; $(MAKE) clean )
	( cd COLAMD ; $(MAKE) clean )
	( cd CCOLAMD ; $(MAKE) clean )
	( cd CHOLMOD ; $(MAKE) clean )

distclean: purge

# statement coverage (Linux only); this requires a lot of time.
cov:
	( cd CHOLMOD ; $(MAKE) cov )
