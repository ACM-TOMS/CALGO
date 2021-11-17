all: startsys drivmhom drivmuho drivss climp

# possible options for the compiler (normal,optim,suppress):
#  normal : without optimizations
#  optim : optimize as far as possible
#  suppress : same as optim, but with the suppression of all checks

# WARNING : Optimizer produces wrong results with enumeration strategy
#           in the package m_Homogeneous_Bezout_Numbers!

normal= 
optim= -O 
suppress= -O -S
options=$(normal)
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

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Root_Counts/Product

links:
	@-ln -s $(linkrc)/startsys.a      startsys.a
	@-ln -s $(linkrc)/startsysB.a     startsysB.a
	@-ln -s $(linkrc)/sets.a          sets.a
	@-ln -s $(linkrc)/setsB.a         setsB.a
	@-ln -s $(linkrc)/sets_io.a       sets_io.a
	@-ln -s $(linkrc)/sets_ioB.a      sets_ioB.a
	@-ln -s $(linkrc)/partsets.a      partsets.a
	@-ln -s $(linkrc)/partsetsB.a     partsetsB.a
	@-ln -s $(linkrc)/partsets_io.a   partsets_io.a
	@-ln -s $(linkrc)/partsets_ioB.a  partsets_ioB.a
	@-ln -s $(linkrc)/degsets.a       degsets.a
	@-ln -s $(linkrc)/degsetsB.a      degsetsB.a
	@-ln -s $(linkrc)/mhombez.a       mhombez.a
	@-ln -s $(linkrc)/mhombezB.a      mhombezB.a
	@-ln -s $(linkrc)/mhomstart.a     mhomstart.a
	@-ln -s $(linkrc)/mhomstartB.a    mhomstartB.a
	@-ln -s $(linkrc)/generate.a      generate.a
	@-ln -s $(linkrc)/generateB.a     generateB.a
	@-ln -s $(linkrc)/intehom.a       intehom.a
	@-ln -s $(linkrc)/intehomB.a      intehomB.a
	@-ln -s $(linkrc)/drivinho.a      drivinho.a
	@-ln -s $(linkrc)/drivinhoB.a     drivinhoB.a
	@-ln -s $(linkrc)/drivmhom.a      drivmhom.a
	@-ln -s $(linkrc)/drivmhomB.a     drivmhomB.a
	@-ln -s $(linkrc)/rps.a           rps.a
	@-ln -s $(linkrc)/rpsB.a          rpsB.a
	@-ln -s $(linkrc)/rps_io.a        rps_io.a
	@-ln -s $(linkrc)/rps_ioB.a       rps_ioB.a
	@-ln -s $(linkrc)/rpss.a          rpss.a
	@-ln -s $(linkrc)/rpssB.a         rpssB.a
	@-ln -s $(linkrc)/ds.a            ds.a
	@-ln -s $(linkrc)/dsB.a           dsB.a
	@-ln -s $(linkrc)/drivmuho.a      drivmuho.a
	@-ln -s $(linkrc)/drivmuhoB.a     drivmuhoB.a
	@-ln -s $(linkrc)/ss.a            ss.a
	@-ln -s $(linkrc)/ssB.a           ssB.a
	@-ln -s $(linkrc)/ss_io.a         ss_io.a
	@-ln -s $(linkrc)/ss_ioB.a        ss_ioB.a
	@-ln -s $(linkrc)/crpss.a         crpss.a
	@-ln -s $(linkrc)/crpssB.a        crpssB.a
	@-ln -s $(linkrc)/persets.a       persets.a
	@-ln -s $(linkrc)/persetsB.a      persetsB.a
	@-ln -s $(linkrc)/drivss.a        drivss.a
	@-ln -s $(linkrc)/drivssB.a       drivssB.a

# Start systems based on the total degree :

startsys: ada.lib startsys.a startsysB.a
	$(compile) startsys.a startsysB.a

# Sets of unknowns :

sets: ada.lib sets.a setsB.a sets_io.a sets_ioB.a
	$(compile) sets.a setsB.a sets_io.a sets_ioB.a

partsets: sets partsets.a partsetsB.a partsets_io.a partsets_ioB.a
	$(compile) partsets.a partsetsB.a partsets_io.a partsets_ioB.a

degsets: partsets degsets.a degsetsB.a generate.a generateB.a
	$(compile) degsets.a degsetsB.a generate.a generateB.a

# Random product systems :

rps: ada.lib rps.a rpsB.a rps_io.a rps_ioB.a
	$(compile) rps.a rpsB.a rps_io.a rps_ioB.a

# m-homogenization (interpolating and linear-product start systems) :

mhom: degsets mhombez.a mhombezB.a rps mhomstart.a mhomstartB.a
	$(compile) mhombez.a mhombezB.a
	$(compile) mhomstart.a mhomstartB.a

intehom: degsets intehom.a intehomB.a drivinho.a drivinhoB.a
	$(compile) intehom.a intehomB.a drivinho.a drivinhoB.a

drivmhom: mhom intehom drivmhom.a drivmhomB.a
	$(compile) drivmhom.a drivmhomB.a

# Projective transformations:

mprojtrans: sets mprojtrans.a mprojtransB.a
	$(compile) mprojtrans.a mprojtransB.a

# multi-homogenization (the GBQ-Algorithm) :

muho: degsets ds.a dsB.a rps rpss.a rpssB.a
	$(compile) ds.a dsB.a
	$(compile) rpss.a rpssB.a

drivmuho: muho drivmuho.a drivmuhoB.a
	$(compile) drivmuho.a drivmuhoB.a

# set-structure analysis (Symbolic Homotopy Construction) :

ss: degsets ss.a ssB.a ss_io.a ss_ioB.a rps crpss.a crpssB.a
	$(compile) ss.a ssB.a ss_io.a ss_ioB.a crpss.a crpssB.a

drivss: ss persets.a persetsB.a drivss.a drivssB.a
	$(compile) persets.a persetsB.a drivss.a drivssB.a
