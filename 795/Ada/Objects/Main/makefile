all: drivroco drivrore climp

# possible options for the compiler (normal,optim,suppress):
#  normal : without optimizations
#  optim : optimize as far as possible
#  suppress : same as optim, but with the suppression of all checks

normal=
optim= -O 
suppress= -O -S
options=$(suppress)

# Making an Ada library :

ada.lib:
	@-../makelib
	a.path -a ../System
	a.path -a ../Math_Lib/Numbers
	a.path -a ../Math_Lib/Matrices
	a.path -a ../Math_Lib/Supports
	a.path -a ../Math_Lib/Polynomials
	a.path -a ../Homotopy
	a.path -a ../Continuation
	a.path -a ../Root_Counts/Product
	a.path -a ../Root_Counts/Implift
	a.path -a ../Root_Counts/Stalift
	a.path -a ../Root_Counts/Dynlift
	a.path -a ../Root_Counts/Symmetry

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning imports :

climp:
	@-rm -f .imports/*

# Cleaning superfluous instantiations :

cleaninst:
	@-a.cleaninst

# Installing everything :

install:
	@-./makeall

# Cleaning everything :

cleanall:
	@-./makeclean

# Establishing the links :

linkrc = ../../Main

links:
	@-ln -s $(linkrc)/drivowst.a    drivowst.a
	@-ln -s $(linkrc)/drivowstB.a   drivowstB.a
	@-ln -s $(linkrc)/drivroco.a    drivroco.a
	@-ln -s $(linkrc)/drivrocoB.a   drivrocoB.a
	@-ln -s $(linkrc)/blackroco.a   blackroco.a
	@-ln -s $(linkrc)/blackrocoB.a  blackrocoB.a
	@-ln -s $(linkrc)/bablroco.a    bablroco.a
	@-ln -s $(linkrc)/bablrocoB.a   bablrocoB.a
	@-ln -s $(linkrc)/mainroco.a    mainroco.a
	@-ln -s $(linkrc)/mainrocoB.a   mainrocoB.a
	@-ln -s $(linkrc)/drivrore.a    drivrore.a
	@-ln -s $(linkrc)/drivroreB.a   drivroreB.a
	@-ln -s $(linkrc)/mainvali.a    mainvali.a
	@-ln -s $(linkrc)/mainvaliB.a   mainvaliB.a
	@-ln -s $(linkrc)/bablvali.a    bablvali.a
	@-ln -s $(linkrc)/bablvaliB.a   bablvaliB.a
	@-ln -s $(linkrc)/bablphc.a     bablphc.a
	@-ln -s $(linkrc)/bablphcB.a    bablphcB.a
	@-ln -s $(linkrc)/mainphc.a     mainphc.a
	@-ln -s $(linkrc)/mainphcB.a    mainphcB.a
	@-ln -s $(linkrc)/dispatch.a    dispatch.a
	@-ln -s $(linkrc)/dispatchB.a   dispatchB.a
	@-ln -s $(linkrc)/phcpack.a     phcpack.a
	@-ln -s $(linkrc)/phcpackB.a    phcpackB.a
	@-ln -s $(linkrc)/use_phc.a     use_phc.a

# The interactive driver for the preprocessing :

drivprep: ada.lib drivprep.a drivprepB.a
	ada $(options) drivprep.a drivprepB.a

# The interactive driver for the root counting methods :

drivowst: ada.lib drivowst.a drivowstB.a
	ada $(options) drivowst.a drivowstB.a

drivroco: drivowst drivroco.a drivrocoB.a mainroco.a mainrocoB.a
	ada $(options) drivroco.a drivrocoB.a mainroco.a mainrocoB.a

roco: drivroco blackroco.a blackrocoB.a bablroco.a bablrocoB.a
	ada $(options) blackroco.a bablroco.a bablrocoB.a
	ada blackrocoB.a

# The driver for root refining :

validate: ada.lib mainvali.a mainvaliB.a bablvali.a bablvaliB.a
	ada $(options) mainvali.a mainvaliB.a bablvali.a bablvaliB.a

drivrore: validate drivrore.a drivroreB.a
	ada $(options) drivrore.a drivroreB.a

# The main interactive dispatcher:

bablphc: roco drivrore mainphc.a mainphcB.a bablphc.a bablphcB.a
	ada $(options) mainphc.a mainphcB.a bablphc.a bablphcB.a

phc: bablphc dispatch.a dispatchB.a
	ada $(options) dispatch.a dispatchB.a
	@-make climp
	a.ld -o ../../../bin/phc dispatch
