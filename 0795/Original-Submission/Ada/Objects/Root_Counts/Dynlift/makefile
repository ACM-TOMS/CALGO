all: drivdynl babldmvc climp

# possible options for the compiler (normal,optim,suppress):
#  normal : without optimizations
#  optim : optimize as far as possible
#  suppress : same as optim, but with the suppression of all checks

normal= 
optim= -O 
suppress= -O -S 
options=$(suppress)
makeoptions=$(options) -v -f

# compiling command :
#  ada    : compile the file
#  a.make : verifies dependencies before compiling

compile=ada $(options)
# compile=a.make $(makeoptions)

# Making an Ada library :

ada.lib:
	@-../../makelib
	a.path -a ../../System
	a.path -a ../../Math_Lib/Numbers
	a.path -a ../../Math_Lib/Matrices
	a.path -a ../../Math_Lib/Polynomials
	a.path -a ../../Math_Lib/Supports
	a.path -a ../../Homotopy
	a.path -a ../../Continuation
	a.path -a ../Product
	a.path -a ../Implift
	a.path -a ../Stalift

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Root_Counts/Dynlift

links:
	@-ln -s $(linkrc)/simplex.a      simplex.a
	@-ln -s $(linkrc)/simplexB.a     simplexB.a
	@-ln -s $(linkrc)/simplex_io.a   simplex_io.a
	@-ln -s $(linkrc)/simplex_ioB.a  simplex_ioB.a
	@-ln -s $(linkrc)/triangle.a     triangle.a
	@-ln -s $(linkrc)/triangleB.a    triangleB.a
	@-ln -s $(linkrc)/triangle_io.a  triangle_io.a
	@-ln -s $(linkrc)/triangle_ioB.a triangle_ioB.a
	@-ln -s $(linkrc)/glodyntri.a    glodyntri.a
	@-ln -s $(linkrc)/glodyntriB.a   glodyntriB.a
	@-ln -s $(linkrc)/dynlift.a      dynlift.a
	@-ln -s $(linkrc)/dynliftB.a     dynliftB.a
	@-ln -s $(linkrc)/dyntri.a       dyntri.a
	@-ln -s $(linkrc)/dyntriB.a      dyntriB.a
	@-ln -s $(linkrc)/cayemb.a       cayemb.a
	@-ln -s $(linkrc)/cayembB.a      cayembB.a
	@-ln -s $(linkrc)/cayley.a       cayley.a
	@-ln -s $(linkrc)/cayleyB.a      cayleyB.a
	@-ln -s $(linkrc)/minkpoly.a     minkpoly.a
	@-ln -s $(linkrc)/minkpolyB.a    minkpolyB.a
	@-ln -s $(linkrc)/drivmink.a     drivmink.a
	@-ln -s $(linkrc)/drivminkB.a    drivminkB.a
	@-ln -s $(linkrc)/commfaces.a    commfaces.a
	@-ln -s $(linkrc)/commfacesB.a   commfacesB.a
	@-ln -s $(linkrc)/enumfaces.a    enumfaces.a
	@-ln -s $(linkrc)/enumfacesB.a   enumfacesB.a
	@-ln -s $(linkrc)/freqgraph.a    freqgraph.a
	@-ln -s $(linkrc)/freqgraphB.a   freqgraphB.a
	@-ln -s $(linkrc)/initmice.a     initmice.a
	@-ln -s $(linkrc)/initmiceB.a    initmiceB.a
	@-ln -s $(linkrc)/flatmisu.a     flatmisu.a
	@-ln -s $(linkrc)/flatmisuB.a    flatmisuB.a
	@-ln -s $(linkrc)/unfolding.a    unfolding.a
	@-ln -s $(linkrc)/unfoldingB.a   unfoldingB.a
	@-ln -s $(linkrc)/triamisu.a     triamisu.a
	@-ln -s $(linkrc)/triamisuB.a    triamisuB.a
	@-ln -s $(linkrc)/dymisudi.a     dymisudi.a
	@-ln -s $(linkrc)/dymisudiB.a    dymisudiB.a
	@-ln -s $(linkrc)/dynpolco.a     dynpolco.a
	@-ln -s $(linkrc)/dynpolcoB.a    dynpolcoB.a
	@-ln -s $(linkrc)/drivdynl.a     drivdynl.a
	@-ln -s $(linkrc)/drivdynlB.a    drivdynlB.a
	@-ln -s $(linkrc)/dbkkcomp.a     dbkkcomp.a
	@-ln -s $(linkrc)/dbkkcompB.a    dbkkcompB.a
	@-ln -s $(linkrc)/blackmvc.a     blackmvc.a
	@-ln -s $(linkrc)/blackmvcB.a    blackmvcB.a
	@-ln -s $(linkrc)/babldmvc.a     babldmvc.a
	@-ln -s $(linkrc)/babldmvcB.a    babldmvcB.a

# simplices :

simplex: ada.lib simplex.a simplexB.a simplex_io.a simplex_ioB.a
	$(compile) simplex.a simplexB.a simplex_io.a simplex_ioB.a

# triangulations :

triangle: simplex triangle.a triangleB.a triangle_io.a triangle_ioB.a
	$(compile) triangle.a triangleB.a triangle_io.a triangle_ioB.a

# utilities and dynamic lifting 

dyntri: triangle glodyntri.a glodyntriB.a dynlift.a dynliftB.a dyntri.a dyntriB.a
	$(compile) glodyntri.a glodyntriB.a dynlift.a dynliftB.a
	$(compile) dyntri.a dyntriB.a

# The Cayley trick :

flatmisu: ada.lib flatmisu.a flatmisuB.a
	$(compile) flatmisu.a flatmisuB.a

cayley: dyntri flatmisu cayemb.a cayembB.a cayley.a cayleyB.a
	$(compile) cayemb.a cayembB.a cayley.a cayleyB.a

minkpoly: cayley minkpoly.a minkpolyB.a drivmink.a drivminkB.a
	$(compile) minkpoly.a minkpolyB.a drivmink.a drivminkB.a

# Dynamic mixed subdivisions :

fac_util: ada.lib commfaces.a commfacesB.a enumfaces.a enumfacesB.a
	$(compile) commfaces.a commfacesB.a enumfaces.a enumfacesB.a

initmice: ada.lib freqgraph.a freqgraphB.a initmice.a initmiceB.a
	$(compile) freqgraph.a freqgraphB.a initmice.a initmiceB.a

unfolding: flatmisu unfolding.a unfoldingB.a triamisu.a triamisuB.a
	$(compile) unfolding.a unfoldingB.a triamisu.a triamisuB.a

dymisudi: dyntri fac_util initmice unfolding dymisudi.a dymisudiB.a
	$(compile) dymisudi.a dymisudiB.a

dynpolco: dymisudi dynpolco.a dynpolcoB.a
	$(compile) dynpolco.a dynpolcoB.a

# The drivers :

drivdynl: minkpoly dynpolco drivdynl.a drivdynlB.a
	$(compile) drivdynl.a drivdynlB.a

blackmvc: dynpolco dbkkcomp.a dbkkcompB.a blackmvc.a blackmvcB.a
	ada dbkkcomp.a dbkkcompB.a
	$(compile) blackmvc.a blackmvcB.a

babldmvc: blackmvc babldmvc.a babldmvcB.a
	$(compile) babldmvc.a babldmvcB.a
