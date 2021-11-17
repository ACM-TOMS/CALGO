all: drivimpl climp

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

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Root_Counts/Implift

links:
	@-ln -s $(linkrc)/transfo.a      transfo.a
	@-ln -s $(linkrc)/transfoB.a     transfoB.a
	@-ln -s $(linkrc)/transfo_io.a   transfo_io.a
	@-ln -s $(linkrc)/transfo_ioB.a  transfo_ioB.a
	@-ln -s $(linkrc)/intvec_ut.a    intvec_ut.a
	@-ln -s $(linkrc)/intvec_utB.a   intvec_utB.a
	@-ln -s $(linkrc)/lstivc_ut.a    lstivc_ut.a
	@-ln -s $(linkrc)/lstivc_utB.a   lstivc_utB.a
	@-ln -s $(linkrc)/arrlivc_ut.a   arrlivc_ut.a
	@-ln -s $(linkrc)/arrlivc_utB.a  arrlivc_utB.a
	@-ln -s $(linkrc)/transols.a     transols.a
	@-ln -s $(linkrc)/transolsB.a    transolsB.a
	@-ln -s $(linkrc)/tranlaur.a     tranlaur.a
	@-ln -s $(linkrc)/tranlaurB.a    tranlaurB.a
	@-ln -s $(linkrc)/tranlists.a    tranlists.a
	@-ln -s $(linkrc)/tranlistsB.a   tranlistsB.a
	@-ln -s $(linkrc)/binom.a        binom.a
	@-ln -s $(linkrc)/binomB.a       binomB.a
	@-ln -s $(linkrc)/fewnom.a       fewnom.a
	@-ln -s $(linkrc)/fewnomB.a      fewnomB.a
	@-ln -s $(linkrc)/tv.a           tv.a
	@-ln -s $(linkrc)/tvB.a          tvB.a
	@-ln -s $(linkrc)/tv_io.a        tv_io.a
	@-ln -s $(linkrc)/tv_ioB.a       tv_ioB.a
	@-ln -s $(linkrc)/powlis.a       powlis.a
	@-ln -s $(linkrc)/powlisB.a      powlisB.a
	@-ln -s $(linkrc)/durker.a       durker.a
	@-ln -s $(linkrc)/durkerB.a      durkerB.a
	@-ln -s $(linkrc)/mihoco.a       mihoco.a
	@-ln -s $(linkrc)/mihocoB.a      mihocoB.a
	@-ln -s $(linkrc)/facenusu.a     facenusu.a
	@-ln -s $(linkrc)/facenusuB.a    facenusuB.a
	@-ln -s $(linkrc)/volumes.a      volumes.a
	@-ln -s $(linkrc)/volumesB.a     volumesB.a
	@-ln -s $(linkrc)/vertices.a     vertices.a
	@-ln -s $(linkrc)/verticesB.a    verticesB.a
	@-ln -s $(linkrc)/setsvol.a      setsvol.a
	@-ln -s $(linkrc)/setsvolB.a     setsvolB.a
	@-ln -s $(linkrc)/genpos.a       genpos.a
	@-ln -s $(linkrc)/genposB.a      genposB.a
	@-ln -s $(linkrc)/drivpts.a      drivpts.a
	@-ln -s $(linkrc)/drivptsB.a     drivptsB.a
	@-ln -s $(linkrc)/drivplc.a      drivplc.a
	@-ln -s $(linkrc)/drivplcB.a     drivplcB.a
	@-ln -s $(linkrc)/drivimpl.a     drivimpl.a
	@-ln -s $(linkrc)/drivimplB.a    drivimplB.a

# Unimodular transformations :

transfo: ada.lib transfo.a transfoB.a transfo_io.a transfo_ioB.a
	$(compile) transfo.a transfoB.a transfo_io.a transfo_ioB.a

intvec_ut: transfo intvec_ut.a intvec_utB.a
	$(compile) intvec_ut.a intvec_utB.a

tranlists: intvec_ut tranlists.a tranlistsB.a
	$(compile) tranlists.a tranlistsB.a

transform: tranlists transols.a transolsB.a tranlaur.a tranlaurB.a
	$(compile) transols.a transolsB.a tranlaur.a tranlaurB.a

# Utilities for lists of integer vectors :

lstivc_ut: intvec_ut lstivc_ut.a lstivc_utB.a
	$(compile) lstivc_ut.a lstivc_utB.a

arrlivc_ut: transform lstivc_ut arrlivc_ut.a arrlivc_utB.a
	$(compile) arrlivc_ut.a arrlivc_utB.a

# solving binomial and fewnomial systems :

bi: ada.lib transform binom.a binomB.a
	$(compile) binom.a binomB.a

few: bi fewnom.a fewnomB.a
	$(compile) fewnom.a fewnomB.a

# auxiliairy data structure :

tv: ada.lib tv.a tvB.a tv_io.a tv_ioB.a
	$(compile) tv.a tvB.a tv_io.a tv_ioB.a

# computing volumes and mixed volumes :

vol: arrlivc_ut tv facenusu.a facenusuB.a volumes.a volumesB.a
	$(compile) facenusu.a facenusuB.a volumes.a volumesB.a

# constructing support sets

powlis: ada.lib powlis.a powlisB.a
	$(compile) powlis.a powlisB.a

# The method of Durand-Kerner :

durker: ada.lib durker.a durkerB.a
	$(compile) durker.a durkerB.a

# mixed homotopy continuation :

mihoco: few vol powlis durker mihoco.a mihocoB.a
	$(compile) mihoco.a mihocoB.a

# extracting the vertices :

vertices: ada.lib vertices.a verticesB.a drivpts.a drivptsB.a
	$(compile) vertices.a verticesB.a drivpts.a drivptsB.a

# checking generic position :

genpos: powlis tv genpos.a genposB.a
	$(compile) genpos.a genposB.a

# an interactive driver for computing the BKK bound and for
# constructing a random coefficient start system :

drivplc: ada.lib drivplc.a drivplcB.a
	$(compile) drivplc.a drivplcB.a

setsvol: mihoco vol setsvol.a setsvolB.a
	$(compile) setsvol.a setsvolB.a

drivimpl: setsvol vertices genpos drivplc drivimpl.a drivimplB.a
	$(compile) drivimpl.a drivimplB.a
