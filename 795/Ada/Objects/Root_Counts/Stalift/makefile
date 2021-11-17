all: bkkcomp drivstal climp

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

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Root_Counts/Stalift

links:
	@-ln -s $(linkrc)/intlift.a      intlift.a
	@-ln -s $(linkrc)/intliftB.a     intliftB.a
	@-ln -s $(linkrc)/intlift_ut.a   intlift_ut.a
	@-ln -s $(linkrc)/intlift_utB.a  intlift_utB.a
	@-ln -s $(linkrc)/intmisu.a      intmisu.a
	@-ln -s $(linkrc)/intmisuB.a     intmisuB.a
	@-ln -s $(linkrc)/intmisu_io.a   intmisu_io.a
	@-ln -s $(linkrc)/intmisu_ioB.a  intmisu_ioB.a
	@-ln -s $(linkrc)/intprune.a     intprune.a
	@-ln -s $(linkrc)/intpruneB.a    intpruneB.a
	@-ln -s $(linkrc)/comisudi.a     comisudi.a
	@-ln -s $(linkrc)/comisudiB.a    comisudiB.a
	@-ln -s $(linkrc)/mivoco.a       mivoco.a
	@-ln -s $(linkrc)/mivocoB.a      mivocoB.a
	@-ln -s $(linkrc)/intpolco.a     intpolco.a
	@-ln -s $(linkrc)/intpolcoB.a    intpolcoB.a
	@-ln -s $(linkrc)/bkkcomp.a      bkkcomp.a
	@-ln -s $(linkrc)/bkkcompB.a     bkkcompB.a
	@-ln -s $(linkrc)/fltmisu.a      fltmisu.a
	@-ln -s $(linkrc)/fltmisuB.a     fltmisuB.a
	@-ln -s $(linkrc)/fltmisu_io.a   fltmisu_io.a
	@-ln -s $(linkrc)/fltmisu_ioB.a  fltmisu_ioB.a
	@-ln -s $(linkrc)/fltpolco.a     fltpolco.a
	@-ln -s $(linkrc)/fltpolcoB.a    fltpolcoB.a
	@-ln -s $(linkrc)/fltlift.a      fltlift.a
	@-ln -s $(linkrc)/fltliftB.a     fltliftB.a
	@-ln -s $(linkrc)/fltlift_ut.a   fltlift_ut.a
	@-ln -s $(linkrc)/fltlift_utB.a  fltlift_utB.a
	@-ln -s $(linkrc)/fltprune.a     fltprune.a
	@-ln -s $(linkrc)/fltpruneB.a    fltpruneB.a
	@-ln -s $(linkrc)/cofltint.a     cofltint.a
	@-ln -s $(linkrc)/cofltintB.a    cofltintB.a
	@-ln -s $(linkrc)/inncones.a     inncones.a
	@-ln -s $(linkrc)/innconesB.a    innconesB.a
	@-ln -s $(linkrc)/noconint.a     noconint.a
	@-ln -s $(linkrc)/noconintB.a    noconintB.a
	@-ln -s $(linkrc)/contrimv.a     contrimv.a
	@-ln -s $(linkrc)/contrimvB.a    contrimvB.a
	@-ln -s $(linkrc)/drivcomv.a     drivcomv.a
	@-ln -s $(linkrc)/drivcomvB.a    drivcomvB.a
	@-ln -s $(linkrc)/drivcrit.a     drivcrit.a
	@-ln -s $(linkrc)/drivcritB.a    drivcritB.a
	@-ln -s $(linkrc)/drivlift.a     drivlift.a
	@-ln -s $(linkrc)/drivliftB.a    drivliftB.a
	@-ln -s $(linkrc)/prunstat.a     prunstat.a
	@-ln -s $(linkrc)/prunstatB.a    prunstatB.a
	@-ln -s $(linkrc)/drivstal.a     drivstal.a
	@-ln -s $(linkrc)/drivstalB.a    drivstalB.a

# Integer regular mixed subdivisions :

lift: ada.lib intlift.a intliftB.a
	$(compile) intlift.a intliftB.a

misudi: ada.lib intmisu.a mivoco.a intmisuB.a intmisu_io.a intmisu_ioB.a
	$(compile) intmisu.a mivoco.a intmisuB.a
	$(compile) intmisu_io.a intmisu_ioB.a

intlift_ut: misudi lift intlift_ut.a intlift_utB.a
	$(compile) intlift_ut.a intlift_utB.a

comisudi: intlift_ut intprune.a intpruneB.a comisudi.a comisudiB.a
	$(compile) intprune.a intpruneB.a comisudi.a comisudiB.a

# Floating-point regular mixed subdivisions :

fltmisu: ada.lib fltmisu.a fltmisuB.a
	$(compile) fltmisu.a fltmisuB.a

fltpolco: fltmisu fltpolco.a fltpolcoB.a
	$(compile) fltpolco.a fltpolcoB.a

fltprune: fltpolco fltprune.a fltpruneB.a fltlift.a fltliftB.a
	$(compile) fltprune.a fltpruneB.a fltlift.a fltliftB.a

cofltint: fltmisu cofltint.a cofltintB.a fltlift_ut.a fltlift_utB.a
	$(compile) cofltint.a cofltintB.a fltlift_ut.a fltlift_utB.a

fltmisu_io: cofltint fltmisu_io.a fltmisu_ioB.a
	$(compile) fltmisu_io.a fltmisu_ioB.a

drivlift: fltmisu_io fltprune drivlift.a drivliftB.a
	$(compile) drivlift.a drivliftB.a

# Computing the mixed volume, based on a coherent mixed subdivision :

mixvol: comisudi mixvol.a
	$(compile) mixvol.a
	@-make climp
	a.ld -o /tmp/mixvol mixed_volume

mivoco: comisudi misudi mivocoB.a
	$(compile) mivocoB.a

# Polyhedral homotopy continuation :

intpolco: misudi bkkcomp.a intpolco.a intpolcoB.a
	$(compile) bkkcomp.a
	$(compile) intpolco.a intpolcoB.a

# Black box routines for computing BKK bound :

bkkcomp: mivoco intpolco bkkcompB.a
	ada bkkcompB.a

# note that optimization for bkkcompB.a causes an assertion_error

# Computing zero contributions and essential sets :

norcones: ada.lib inncones.a innconesB.a noconint.a noconintB.a
	$(compile) inncones.a innconesB.a noconint.a noconintB.a

drivcomv: norcones contrimv.a contrimvB.a drivcomv.a drivcomvB.a
	$(compile) contrimv.a contrimvB.a drivcomv.a drivcomvB.a

drivcrit: drivcomv drivcrit.a drivcritB.a
	$(compile) drivcrit.a drivcritB.a

# Drivers for constructing a regular mixed subdivision :

prunstat: ada.lib prunstat.a prunstatB.a
	$(compile) prunstat.a prunstatB.a

drivstal: intlift_ut bkkcomp drivcrit prunstat drivlift drivstal.a drivstalB.a
	$(compile) drivstal.a drivstalB.a
