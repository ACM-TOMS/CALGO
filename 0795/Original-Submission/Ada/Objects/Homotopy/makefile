all: solutions projtrans drivhoco drivscal drivred climp

# possible options for the compiler (normal,optim,suppress):
#   normal : without optimizations
#   optim : optimize as far as possible
#   suppress : same as optim, but with the suppression of all checks

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
	@-../makelib
	a.path -a ../System
	a.path -a ../Math_Lib/Numbers
	a.path -a ../Math_Lib/Matrices
	a.path -a ../Math_Lib/Polynomials

# Cleaning unnecessary instantiations :

cleaninst:
	a.cleaninst

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../Homotopy

links:
	@-ln -s $(linkrc)/equals.a         equals.a
	@-ln -s $(linkrc)/equalsB.a        equalsB.a
	@-ln -s $(linkrc)/solutions.a      solutions.a
	@-ln -s $(linkrc)/solutionsB.a     solutionsB.a
	@-ln -s $(linkrc)/solutions_io.a   solutions_io.a
	@-ln -s $(linkrc)/solutions_ioB.a  solutions_ioB.a
	@-ln -s $(linkrc)/homotopy.a       homotopy.a
	@-ln -s $(linkrc)/homotopyB.a      homotopyB.a
	@-ln -s $(linkrc)/homzers.a        homzers.a
	@-ln -s $(linkrc)/homzersB.a       homzersB.a
	@-ln -s $(linkrc)/drivhoco.a       drivhoco.a
	@-ln -s $(linkrc)/drivhocoB.a      drivhocoB.a
	@-ln -s $(linkrc)/projtrans.a      projtrans.a
	@-ln -s $(linkrc)/projtransB.a     projtransB.a
	@-ln -s $(linkrc)/scaling.a        scaling.a
	@-ln -s $(linkrc)/scalingB.a       scalingB.a
	@-ln -s $(linkrc)/drivscals.a      drivscals.a
	@-ln -s $(linkrc)/drivscalsB.a     drivscalsB.a
	@-ln -s $(linkrc)/mainscal.a       mainscal.a
	@-ln -s $(linkrc)/mainscalB.a      mainscalB.a
	@-ln -s $(linkrc)/redupoly.a       redupoly.a
	@-ln -s $(linkrc)/redupolyB.a      redupolyB.a
	@-ln -s $(linkrc)/redufull.a       redufull.a
	@-ln -s $(linkrc)/redufullB.a      redufullB.a
	@-ln -s $(linkrc)/reduover.a       reduover.a
	@-ln -s $(linkrc)/reduoverB.a      reduoverB.a
	@-ln -s $(linkrc)/drivreds.a       drivreds.a
	@-ln -s $(linkrc)/drivredsB.a      drivredsB.a
	@-ln -s $(linkrc)/mainred.a        mainred.a
	@-ln -s $(linkrc)/mainredB.a       mainredB.a

# The solutions :
#   for IBM RS/6000 optimizer on solutions causes STORAGE_ERROR 
#   when clustering of solutions during continuation appears.

sols: ada.lib equals.a equalsB.a solutions.a solutionsB.a
	$(compile) equals.a equalsB.a solutions.a
	ada solutionsB.a

solutions: sols solutions_io.a solutions_ioB.a
	$(compile) solutions_io.a solutions_ioB.a

# The homotopy :

homotopy: ada.lib homotopy.a homotopyB.a
	$(compile) homotopy.a homotopyB.a

# parametric homotopies :

parahoms: ada.lib parahoms.a parahomsB.a parahoms_io.a parahoms_ioB.a
	ada parahoms.a parahomsB.a parahoms_io.a parahoms_ioB.a

# Adding homogeneous equations :

homzers: ada.lib homzers.a homzersB.a
	$(compile) homzers.a homzersB.a

projtrans: homzers projtrans.a projtransB.a
	$(compile) projtrans.a projtransB.a

# The interactive driver for the construction of the homotopy :

drivhoco: homotopy homzers drivhoco.a drivhocoB.a
	$(compile) drivhoco.a drivhocoB.a

# Scaling and reduction :

scaling: solutions scaling.a scalingB.a
	$(compile) scaling.a scalingB.a

drivscal: scaling drivscals.a drivscalsB.a mainscal.a mainscalB.a
	$(compile) drivscals.a drivscalsB.a mainscal.a mainscalB.a

redufull: ada.lib redupoly.a redupolyB.a redufull.a redufullB.a
	$(compile) redupoly.a redupolyB.a redufull.a redufullB.a

reduover: redufull reduover.a reduoverB.a
	$(compile) reduover.a reduoverB.a

drivred: reduover drivreds.a drivredsB.a mainred.a mainredB.a
	$(compile) drivreds.a drivredsB.a mainred.a mainredB.a
