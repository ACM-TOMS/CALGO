###########################################################################
# RealPaver v. 0.4                                                        #
###########################################################################
# Author: Laurent Granvilliers                                            #
# Copyright (c) 1999-2003 Institut de Recherche en Informatique de Nantes #
# Copyright (c) 2004      Laboratoire d'Informatique de Nantes Atlantique #
###########################################################################
# RealPaver is distributed WITHOUT ANY WARRANTY. Read the associated      #
# COPYRIGHT file for more details.                                        #
###########################################################################
# Makefile.am                                                             #
###########################################################################


SUBDIRS=. src
DIST_SUBDIRS=src

# Rule to build tar-bzipped distribution package
$(PACKAGE)-$(VERSION).tar.gz: dist

# Rule to build RPM distribution package
rpm: $(PACKAGE)-$(VERSION).tar.gz
	rpm --tarball a $(PACKAGE)-$(VERSION).tar.gz
