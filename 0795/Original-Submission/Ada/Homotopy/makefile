ADAPATH=/cw/vads.dec/v6.1.0d

all: solutions sparsity projtrans drivhoco drivscal drivred climp

# possible options for the compiler (normal,optim,suppress):
#   normal : without optimizations
#   optim : optimize as far as possible
#   suppress : same as optim, but with the suppression of all checks

normal= 
optim= -O
suppress= -O -S
options=$(suppress)
makeoptions=$(options) -v -f

# compiling commands :
#  ada    : invoke the Ada compiler directly
#  a.make : check dependencies before compiling

compile=ada $(options)

# Making an Ada library :

ada.lib:
	a.mklib . $(ADAPATH)/verdixlib
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

# The solutions :

sols: ada.lib equals.a equalsB.a solutions.a solutionsB.a
	$(compile) equals.a equalsB.a solutions.a solutionsB.a

solutions: sols solutions_io.a solutions_ioB.a condsols.a condsolsB.a
	$(compile) solutions_io.a solutions_ioB.a
	$(compile) condsols.a condsolsB.a

# The homotopy :

homotopy: ada.lib homotopy.a homotopyB.a
	$(compile) homotopy.a homotopyB.a

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

drivred: reduover ratio.a drivreds.a drivredsB.a mainred.a mainredB.a
	$(compile) ratio.a drivreds.a drivredsB.a mainred.a mainredB.a

# Some tools

usepert: ada.lib perturb.a pert.a
	$(compile) perturb.a pert.a
	a.ld -o /tmp/pert Use_Perturb

evalfj: ada.lib evalfj.a
	$(compile) evalfj.a
	a.ld -o /tmp/evalfj evalfj

# Sampling polynomial systems :

samplers: ada.lib samplers.a samplersB.a
	$(compile) samplers.a samplersB.a

# Measuring the sparsity of a polynomial system :

sparsity: ada.lib sparsity.a sparsityB.a
	$(compile) sparsity.a sparsityB.a

avsp: sparsity avsp.a
	$(compile) avsp.a
	@-make climp
	a.ld -o /tmp/avsp avsp

# Test routines :

ts_drivred: drivred ts_drivred.a
	$(compile) ts_drivred.a
	@-make climp
	a.ld -o /tmp/ts_drivred ts_drivred

ts_drivscal: drivscal ts_drivscal.a
	$(compile) ts_drivscal.a
	@-make climp
	a.ld -o /tmp/ts_drivscal ts_drivscal

ts_nterms: sparsity ts_nterms.a
	$(compile) ts_nterms.a
	@-make climp
	a.ld -o /tmp/ts_nterms ts_nterms

ts_sam: samplers ts_samplers.a
	$(compile) ts_samplers.a
	@-make climp
	a.ld -o /tmp/ts_sam ts_samplers
