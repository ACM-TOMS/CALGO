all: bablpoco drivwind valipoco climp

# possible options for the compiler (normal,optim,suppress):
#  normal : without optimizations
#  optim : optimize as far as possible
#  suppress : same as optim, but with the suppression of all checks

normal= 
optim= -O
suppress= -O -S
options=$(suppress)
makeoptions=$(options) -v -f

# compiling commands :
#  ada    : invoke the Ada compiler directly
#  a.make : check dependencies before compiling

compile=ada $(options)
# compile=a.make $(makeoptions)

# Making an Ada library :

ada.lib:
	@-../makelib
	a.path -a ../System
	a.path -a ../Math_Lib/Numbers
	a.path -a ../Math_Lib/Matrices
	a.path -a ../Math_Lib/Polynomials
	a.path -a ../Homotopy

# Cleaning superfluous instantations :

cleaninst:
	@-a.cleaninst

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../Continuation

links:
	@-ln -s $(linkrc)/process_io.a      process_io.a
	@-ln -s $(linkrc)/process_ioB.a     process_ioB.a
	@-ln -s $(linkrc)/cont_data.a       cont_data.a
	@-ln -s $(linkrc)/cont_dataB.a      cont_dataB.a
	@-ln -s $(linkrc)/conpar.a          conpar.a
	@-ln -s $(linkrc)/conparB.a         conparB.a
	@-ln -s $(linkrc)/conpar_io.a       conpar_io.a
	@-ln -s $(linkrc)/conpar_ioB.a      conpar_ioB.a
	@-ln -s $(linkrc)/predictors.a      predictors.a
	@-ln -s $(linkrc)/predictorsB.a     predictorsB.a
	@-ln -s $(linkrc)/dispred.a         dispred.a
	@-ln -s $(linkrc)/dispredB.a        dispredB.a
	@-ln -s $(linkrc)/correctors.a      correctors.a
	@-ln -s $(linkrc)/correctorsB.a     correctorsB.a
	@-ln -s $(linkrc)/pathtrac.a        pathtrac.a
	@-ln -s $(linkrc)/pathtracB.a       pathtracB.a
	@-ln -s $(linkrc)/continue.a        continue.a
	@-ln -s $(linkrc)/continueB.a       continueB.a
	@-ln -s $(linkrc)/rootrefi.a        rootrefi.a
	@-ln -s $(linkrc)/rootrefiB.a       rootrefiB.a
	@-ln -s $(linkrc)/drivpoco.a        drivpoco.a
	@-ln -s $(linkrc)/drivpocoB.a       drivpocoB.a
	@-ln -s $(linkrc)/vlprstab.a        vlrpstab.a
	@-ln -s $(linkrc)/vlprstabB.a       vlrpstabB.a
	@-ln -s $(linkrc)/vlprsalg.a        vlrpsalg.a
	@-ln -s $(linkrc)/vlprsalgB.a       vlrpsalgB.a
	@-ln -s $(linkrc)/dirpaths.a        dirpaths.a
	@-ln -s $(linkrc)/dirpathsB.a       dirpathsB.a
	@-ln -s $(linkrc)/drivpadi.a        drivpadi.a
	@-ln -s $(linkrc)/drivpadiB.a       drivpadiB.a
	@-ln -s $(linkrc)/drivwind.a        drivwind.a
	@-ln -s $(linkrc)/drivwindB.a       drivwindB.a
	@-ln -s $(linkrc)/mainpoco.a        mainpoco.a
	@-ln -s $(linkrc)/mainpocoB.a       mainpocoB.a
	@-ln -s $(linkrc)/scanpoco.a        scanpoco.a
	@-ln -s $(linkrc)/scanpocoB.a       scanpocoB.a
	@-ln -s $(linkrc)/blackpoco.a       blackpoco.a
	@-ln -s $(linkrc)/blackpocoB.a      blackpocoB.a
	@-ln -s $(linkrc)/bablpoco.a        bablpoco.a
	@-ln -s $(linkrc)/bablpocoB.a       bablpocoB.a
	@-ln -s $(linkrc)/valipoco.a        valipoco.a
	@-ln -s $(linkrc)/valipocoB.a       valipocoB.a

# the package the manages the i/o during continuation :

proc_io: ada.lib process_io.a process_ioB.a
	$(compile) process_io.a process_ioB.a

# The data and parameters for the continuation :

condat: ada.lib cont_data.a cont_dataB.a
	$(compile) cont_data.a cont_dataB.a

conpar: condat conpar.a conparB.a conpar_io.a conpar_ioB.a
	$(compile) conpar.a conparB.a conpar_io.a conpar_ioB.a

# predictor-corrector based on increment-and-fix :

pc_fix: proc_io condat predictors.a predictorsB.a correctors.a correctorsB.a
	$(compile) predictors.a predictorsB.a correctors.a correctorsB.a

# The vLpRs-Algorithm for extrapolating path directions :

vlprs: ada.lib vlprstab.a vlprstabB.a vlprsalg.a vlprsalgB.a
	$(compile) vlprstab.a vlprstabB.a vlprsalg.a vlprsalgB.a

dirpaths: vlprs dirpaths.a dirpathsB.a
	$(compile) dirpaths.a dirpathsB.a

# the path trackers and the main continuation routines :
# ( Note : optimizing causes run-time error on IBM RS/6000...)

pathtrac: conpar pc_fix dirpaths dispred.a dispredB.a pathtrac.a pathtracB.a
	ada dispred.a dispredB.a pathtrac.a pathtracB.a

continue: pathtrac continue.a continueB.a
	ada continue.a continueB.a

# roots refinement :

roots: ada.lib rootrefi.a rootrefiB.a
	$(compile) rootrefi.a rootrefiB.a

# The interactive drivers :

drivpadi : continue drivpadi.a drivpadiB.a
	$(compile) drivpadi.a drivpadiB.a

drivpoco: drivpadi proc_io roots drivpoco.a drivpocoB.a mainpoco.a mainpocoB.a
	$(compile) drivpoco.a drivpocoB.a mainpoco.a mainpocoB.a

drivwind: drivpoco drivwind.a drivwindB.a
	$(compile) drivwind.a drivwindB.a

blackpoco: drivpoco scanpoco.a scanpocoB.a blackpoco.a blackpocoB.a
	$(compile) scanpoco.a scanpocoB.a blackpoco.a blackpocoB.a

bablpoco: blackpoco bablpoco.a bablpocoB.a
	$(compile) bablpoco.a bablpocoB.a

valipoco: ada.lib valipoco.a valipocoB.a
	$(compile) valipoco.a valipocoB.a
