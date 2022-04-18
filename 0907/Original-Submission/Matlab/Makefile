#-------------------------------------------------------------------------------
# Makefile for KLU, BTF, AMD, and COLAMD
#-------------------------------------------------------------------------------

include UFconfig/UFconfig.mk

# Compile the default rules for each package
default:
	( cd AMD ; $(MAKE) )
	( cd COLAMD ; $(MAKE) )
	( cd BTF ; $(MAKE) )
	( cd KLU ; $(MAKE) )

# install all packages in /usr/local/lib and /usr/local/include
install:
	( cd AMD ; $(MAKE) install )
	( cd COLAMD ; $(MAKE) install )
	( cd BTF ; $(MAKE) install )
	( cd KLU ; $(MAKE) install )

# uninstall all packages
uninstall:
	( cd AMD ; $(MAKE) uninstall )
	( cd COLAMD ; $(MAKE) uninstall )
	( cd BTF ; $(MAKE) uninstall )
	( cd KLU ; $(MAKE) uninstall )

library:
	( cd AMD ; $(MAKE) library )
	( cd BTF ; $(MAKE) library )
	( cd COLAMD ; $(MAKE) library )
	( cd KLU ; $(MAKE) library )

# Remove all files not in the original distribution
purge:
	( cd UFconfig/xerbla ; $(MAKE) purge )
	( cd UFconfig ; $(MAKE) purge )
	( cd AMD ; $(MAKE) purge )
	( cd COLAMD ; $(MAKE) purge )
	( cd BTF ; $(MAKE) purge )
	( cd KLU ; $(MAKE) purge )

# Remove all files not in the original distribution, but keep the libraries
clean:
	( cd UFconfig/xerbla ; $(MAKE) clean )
	( cd UFconfig ; $(MAKE) clean )
	( cd AMD ; $(MAKE) clean )
	( cd COLAMD ; $(MAKE) clean )
	( cd BTF ; $(MAKE) clean )
	( cd KLU ; $(MAKE) clean )

distclean: purge

# statement coverage (Linux only); this requires a lot of time.
# The umfpack tcov requires a lot of disk space
cov:
	( cd KLU ; $(MAKE) cov )
