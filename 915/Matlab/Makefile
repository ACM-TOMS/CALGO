#-------------------------------------------------------------------------------
# Makefile for all UF sparse matrix packages (SuiteSparseQR subset)
#-------------------------------------------------------------------------------

include UFconfig/UFconfig.mk

# Compile the default rules for each package
default:
	( cd UFconfig/xerbla ; $(MAKE) )
	( cd UFconfig ; $(MAKE) )
#	( cd metis-4.0 ; $(MAKE) )
	( cd AMD ; $(MAKE) )
	( cd CAMD ; $(MAKE) )
	( cd COLAMD ; $(MAKE) )
	( cd CCOLAMD ; $(MAKE) )
	( cd CHOLMOD ; $(MAKE) )
	( cd SPQR ; $(MAKE) )

# install all packages in /usr/local/lib and /usr/local/include
install:
	( cd UFconfig ; $(MAKE) install )
	( cd AMD ; $(MAKE) install )
	( cd CAMD ; $(MAKE) install )
	( cd COLAMD ; $(MAKE) install )
	( cd CCOLAMD ; $(MAKE) install )
	( cd CHOLMOD ; $(MAKE) install )
	( cd SPQR ; $(MAKE) install )

# uninstall all packages
uninstall:
	( cd UFconfig ; $(MAKE) uninstall )
	( cd AMD ; $(MAKE) uninstall )
	( cd CAMD ; $(MAKE) uninstall )
	( cd COLAMD ; $(MAKE) uninstall )
	( cd CCOLAMD ; $(MAKE) uninstall )
	( cd CHOLMOD ; $(MAKE) uninstall )
	( cd SPQR ; $(MAKE) uninstall )

library:
	( cd UFconfig/xerbla ; $(MAKE) )
	( cd UFconfig ; $(MAKE) )
#	( cd metis-4.0 ; $(MAKE) )
	( cd AMD ; $(MAKE) library )
	( cd CAMD ; $(MAKE) library )
	( cd CCOLAMD ; $(MAKE) library )
	( cd COLAMD ; $(MAKE) library )
	( cd CHOLMOD ; $(MAKE) library )
	( cd SPQR ; $(MAKE) library )

# Remove all files not in the original distribution
purge:
	( cd UFconfig/xerbla ; $(MAKE) purge )
	( cd UFconfig ; $(MAKE) purge )
#	( cd metis-4.0 ; $(MAKE) realclean )
	( cd AMD ; $(MAKE) purge )
	( cd CAMD ; $(MAKE) purge )
	( cd COLAMD ; $(MAKE) purge )
	( cd CCOLAMD ; $(MAKE) purge )
	( cd CHOLMOD ; $(MAKE) purge )
	( cd SPQR ; $(MAKE) purge )

# Remove all files not in the original distribution, but keep the libraries
clean:
	( cd UFconfig/xerbla ; $(MAKE) clean )
	( cd UFconfig ; $(MAKE) clean )
#	( cd metis-4.0 ; $(MAKE) clean )
	( cd AMD ; $(MAKE) clean )
	( cd CAMD ; $(MAKE) clean )
	( cd COLAMD ; $(MAKE) clean )
	( cd CCOLAMD ; $(MAKE) clean )
	( cd CHOLMOD ; $(MAKE) clean )
	( cd SPQR ; $(MAKE) clean )

distclean: purge

# statement coverage (Linux only); this requires a lot of time.
# The umfpack tcov requires a lot of disk space
cov:
	( cd CHOLMOD ; $(MAKE) cov )
	( cd SPQR ; $(MAKE) cov )
