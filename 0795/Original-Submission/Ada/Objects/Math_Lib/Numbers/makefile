all: complex_io random nat_inst int_inst flt_inst cmp_inst strnum climp

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

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc=../../../Math_Lib/Numbers

links:
	@-ln -s $(linkrc)/float.a       float.a
	@-ln -s $(linkrc)/mathfun.a     mathfun.a
	@-ln -s $(linkrc)/mathfunB.a    mathfunB.a
	@-ln -s $(linkrc)/numbers_io.a  numbers_io.a
	@-ln -s $(linkrc)/numbers_ioB.a numbers_ioB.a
	@-ln -s $(linkrc)/complex.a     complex.a
	@-ln -s $(linkrc)/complexB.a    complexB.a
	@-ln -s $(linkrc)/complex_io.a  complex_io.a
	@-ln -s $(linkrc)/complex_ioB.a complex_ioB.a
	@-ln -s $(linkrc)/random.a      random.a
	@-ln -s $(linkrc)/randomB.a     randomB.a
	@-ln -s $(linkrc)/nat_inst.a    nat_inst.a
	@-ln -s $(linkrc)/nat_instB.a   nat_instB.a
	@-ln -s $(linkrc)/int_inst.a    int_inst.a
	@-ln -s $(linkrc)/int_instB.a   int_instB.a
	@-ln -s $(linkrc)/flt_inst.a    flt_inst.a
	@-ln -s $(linkrc)/flt_instB.a   flt_instB.a
	@-ln -s $(linkrc)/cmp_inst.a    cmp_inst.a
	@-ln -s $(linkrc)/cmp_instB.a   cmp_instB.a
	@-ln -s $(linkrc)/strnum.a      strnum.a
	@-ln -s $(linkrc)/strnumB.a     strnumB.a

# Floating point numbers and user friendly i/o of numbers :

float: ada.lib float.a numbers_io.a numbers_ioB.a mathfun.a mathfunB.a
	$(compile) float.a numbers_io.a numbers_ioB.a mathfun.a mathfunB.a

# Complex numbers :

complex: float complex.a complexB.a
	$(compile) complex.a complexB.a

complex_io: complex complex_io.a complex_ioB.a
	$(compile) complex_io.a complex_ioB.a

# Random number generators :

random: complex random.a randomB.a
	$(compile) random.a randomB.a

# instantiation packages :

nat_inst: ada.lib nat_inst.a nat_instB.a
	$(compile) nat_inst.a nat_instB.a

int_inst: ada.lib int_inst.a int_instB.a
	$(compile) int_inst.a int_instB.a

flt_inst: float flt_inst.a flt_instB.a
	$(compile) flt_inst.a flt_instB.a

cmp_inst: complex cmp_inst.a cmp_instB.a
	$(compile) cmp_inst.a cmp_instB.a

# string conversions :

strnum: ada.lib strnum.a strnumB.a
	$(compile) strnum.a strnumB.a
