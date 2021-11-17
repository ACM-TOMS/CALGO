all: drivsss drivsyml climp

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

# NOTE :
# drivgrp_io raises exceptions when optimized

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
	a.path -a ../Dynlift

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Root_Counts/Symmetry

links:
	@-ln -s $(linkrc)/perms.a        perms.a
	@-ln -s $(linkrc)/permsB.a       permsB.a
	@-ln -s $(linkrc)/permops.a      permops.a
	@-ln -s $(linkrc)/permopsB.a     permopsB.a
	@-ln -s $(linkrc)/templates.a    templates.a
	@-ln -s $(linkrc)/templatesB.a   templatesB.a
	@-ln -s $(linkrc)/symgrp.a       symgrp.a
	@-ln -s $(linkrc)/symgrpB.a      symgrpB.a
	@-ln -s $(linkrc)/symgrp_io.a    symgrp_io.a
	@-ln -s $(linkrc)/symgrp_ioB.a   symgrp_ioB.a
	@-ln -s $(linkrc)/sbsymgrp_io.a  sbsymgrp_io.a
	@-ln -s $(linkrc)/sbsymgrp_ioB.a sbsymgrp_ioB.a
	@-ln -s $(linkrc)/orbits.a       orbits.a
	@-ln -s $(linkrc)/orbitsB.a      orbitsB.a
	@-ln -s $(linkrc)/orbits_io.a    orbits_io.a
	@-ln -s $(linkrc)/orbits_ioB.a   orbits_ioB.a
	@-ln -s $(linkrc)/equpsys.a      equpsys.a
	@-ln -s $(linkrc)/equpsysB.a     equpsysB.a
	@-ln -s $(linkrc)/sym_ss.a       sym_ss.a
	@-ln -s $(linkrc)/sym_ssB.a      sym_ssB.a
	@-ln -s $(linkrc)/lsymred.a      lsymred.a
	@-ln -s $(linkrc)/lsymredB.a     lsymredB.a
	@-ln -s $(linkrc)/drivgrp_io.a   drivgrp_io.a
	@-ln -s $(linkrc)/drivgrp_ioB.a  drivgrp_ioB.a
	@-ln -s $(linkrc)/drivorbi.a     drivorbi.a
	@-ln -s $(linkrc)/drivorbiB.a    drivorbiB.a
	@-ln -s $(linkrc)/drivsss.a      drivsss.a
	@-ln -s $(linkrc)/drivsssB.a     drivsssB.a
	@-ln -s $(linkrc)/symlift.a      symlift.a
	@-ln -s $(linkrc)/symliftB.a     symliftB.a
	@-ln -s $(linkrc)/symrand.a      symrand.a
	@-ln -s $(linkrc)/symrandB.a     symrandB.a
	@-ln -s $(linkrc)/symbkk.a       symbkk.a
	@-ln -s $(linkrc)/symbkkB.a      symbkkB.a
	@-ln -s $(linkrc)/sympolco.a     sympolco.a
	@-ln -s $(linkrc)/sympolcoB.a    sympolcoB.a
	@-ln -s $(linkrc)/faceperm.a     faceperm.a
	@-ln -s $(linkrc)/facepermB.a    facepermB.a
	@-ln -s $(linkrc)/facesypo.a     facesypo.a
	@-ln -s $(linkrc)/facesypoB.a    facesypoB.a
	@-ln -s $(linkrc)/gencells.a     gencells.a
	@-ln -s $(linkrc)/gencellsB.a    gencellsB.a
	@-ln -s $(linkrc)/drivsyml.a     drivsyml.a
	@-ln -s $(linkrc)/drivsymlB.a    drivsymlB.a
	@-ln -s $(linkrc)/mainsmvc.a     mainsmvc.a
	@-ln -s $(linkrc)/mainsmvcB.a    mainsmvcB.a

# Working with permutation groups :

perms: ada.lib perms.a permsB.a
	$(compile) perms.a permsB.a

symgrp: perms symgrp.a symgrpB.a symgrp_io.a symgrp_ioB.a
	$(compile) symgrp.a symgrpB.a symgrp_io.a symgrp_ioB.a

sbsymgrp_io: symgrp sbsymgrp_io.a sbsymgrp_ioB.a
	ada sbsymgrp_io.a sbsymgrp_ioB.a

drivgrp_io: sbsymgrp_io drivgrp_io.a drivgrp_ioB.a
	ada drivgrp_io.a drivgrp_ioB.a

# Applying the symmetry :

permops: symgrp permops.a permopsB.a
	$(compile) permops.a permopsB.a

orbits: permops orbits.a orbitsB.a orbits_io.a orbits_ioB.a
	$(compile) orbits.a orbitsB.a orbits_io.a orbits_ioB.a

drivorbi: orbits drivorbi.a drivorbiB.a
	$(compile) drivorbi.a drivorbiB.a
 
equpsys: permops equpsys.a equpsysB.a
	$(compile) equpsys.a equpsysB.a

# Symmetric set structure analysis :

sym_ss: permops templates.a templatesB.a sym_ss.a sym_ssB.a
	$(compile) templates.a templatesB.a
	$(compile) sym_ss.a sym_ssB.a

lsymred: permops lsymred.a lsymredB.a
	$(compile) lsymred.a lsymredB.a

# The interactive driver needed in the main program :

drivsss: drivgrp_io sym_ss lsymred drivorbi equpsys drivsss.a drivsssB.a
	$(compile) drivsss.a drivsssB.a

# Symmetric lifting algorithm :

symlift: ada.lib symlift.a symliftB.a
	$(compile) symlift.a symliftB.a

sympolco: ada.lib symbkk.a sympolco.a sympolcoB.a
	$(compile) symbkk.a sympolco.a sympolcoB.a

symbkk: symlift sympolco symbkk.a symbkkB.a
	ada symbkkB.a

# Symmetric randomize :

symrand: ada.lib symrand.a symrandB.a
	$(compile) symrand.a symrandB.a

# Computing the mixed volume, based on a symmetric subdivision :

facesypo: ada.lib faceperm.a facepermB.a facesypo.a facesypoB.a
	$(compile) faceperm.a facepermB.a facesypo.a facesypoB.a

gencells: facesypo gencells.a gencellsB.a
	$(compile) gencells.a gencellsB.a

drivsyml: symrand gencells symbkk drivsyml.a drivsymlB.a mainsmvc.a mainsmvcB.a
	$(compile) drivsyml.a drivsymlB.a mainsmvc.a mainsmvcB.a
