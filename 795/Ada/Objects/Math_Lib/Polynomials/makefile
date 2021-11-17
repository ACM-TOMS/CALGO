all: cmp_polsys_io jacobi pollaco tools climp

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
	a.path -a ../Numbers
	a.path -a ../Matrices

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Math_Lib/Polynomials

links:
	@-ln -s $(linkrc)/nat_grad_lex.a   nat_grad_lex.a
	@-ln -s $(linkrc)/nat_grad_lexB.a  nat_grad_lexB.a
	@-ln -s $(linkrc)/int_grad_lex.a   int_grad_lex.a
	@-ln -s $(linkrc)/int_grad_lexB.a  int_grad_lexB.a
	@-ln -s $(linkrc)/lists.a          lists.a
	@-ln -s $(linkrc)/listsB.a         listsB.a
	@-ln -s $(linkrc)/expvec.a         expvec.a
	@-ln -s $(linkrc)/expvecB.a        expvecB.a
	@-ln -s $(linkrc)/mul_poly.a       mul_poly.a
	@-ln -s $(linkrc)/mul_polyB.a      mul_polyB.a
	@-ln -s $(linkrc)/mul_laur.a       mul_laur.a
	@-ln -s $(linkrc)/mul_laurB.a      mul_laurB.a
	@-ln -s $(linkrc)/cmp_mpoly.a      cmp_mpoly.a
	@-ln -s $(linkrc)/cmp_mlaur.a      cmp_mlaur.a
	@-ln -s $(linkrc)/symtab.a         symtab.a
	@-ln -s $(linkrc)/symtabB.a        symtabB.a
	@-ln -s $(linkrc)/symtab_io.a      symtab_io.a
	@-ln -s $(linkrc)/symtab_ioB.a     symtab_ioB.a
	@-ln -s $(linkrc)/cmp_mpoly_io.a   cmp_mpoly_io.a
	@-ln -s $(linkrc)/cmp_mpoly_ioB.a  cmp_mpoly_ioB.a
	@-ln -s $(linkrc)/cmp_polsys.a     cmp_polsys.a
	@-ln -s $(linkrc)/cmp_polsysB.a    cmp_polsysB.a
	@-ln -s $(linkrc)/cmp_laursys.a    cmp_laursys.a
	@-ln -s $(linkrc)/cmp_laursysB.a   cmp_laursysB.a
	@-ln -s $(linkrc)/pollaco.a        pollaco.a
	@-ln -s $(linkrc)/pollacoB.a       pollacoB.a
	@-ln -s $(linkrc)/lapolco.a        lapolco.a
	@-ln -s $(linkrc)/lapolcoB.a       lapolcoB.a
	@-ln -s $(linkrc)/cmp_polsys_io.a  cmp_polsys_io.a
	@-ln -s $(linkrc)/cmp_polsys_ioB.a cmp_polsys_ioB.a
	@-ln -s $(linkrc)/jacobi.a         jacobi.a
	@-ln -s $(linkrc)/jacobiB.a        jacobiB.a
	@-ln -s $(linkrc)/laurjaco.a       laurjaco.a
	@-ln -s $(linkrc)/laurjacoB.a      laurjacoB.a
	@-ln -s $(linkrc)/polrand.a        polrand.a
	@-ln -s $(linkrc)/polrandB.a       polrandB.a
	@-ln -s $(linkrc)/laurrand.a       laurrand.a
	@-ln -s $(linkrc)/laurrandB.a      laurrandB.a
	@-ln -s $(linkrc)/substits.a       substits.a
	@-ln -s $(linkrc)/substitsB.a      substitsB.a

# Term ordenings :

nat_grad_lex: ada.lib nat_grad_lex.a nat_grad_lexB.a
	$(compile) nat_grad_lex.a nat_grad_lexB.a

int_grad_lex: ada.lib int_grad_lex.a int_grad_lexB.a
	$(compile) int_grad_lex.a int_grad_lexB.a

# Generic lists :

lists: ada.lib lists.a listsB.a
	$(compile) lists.a listsB.a

# Generic multivariate polynomials :

mul_poly: lists mul_poly.a mul_polyB.a 
	$(compile) mul_poly.a mul_polyB.a

mul_laur: lists mul_laur.a mul_laurB.a
	$(compile) mul_laur.a mul_laurB.a

# Complex multivariate polynomials :

cmp_mpoly: mul_poly nat_grad_lex cmp_mpoly.a
	$(compile) cmp_mpoly.a

cmp_mlaur: mul_laur int_grad_lex cmp_mlaur.a
	$(compile) cmp_mlaur.a

# Input/Output operations of multivariate polynomials :

symtab: ada.lib symtab.a symtabB.a symtab_io.a symtab_ioB.a
	ada symtab.a symtabB.a symtab_io.a symtab_ioB.a

# Optimizing causes a run time error :

cmp_mpoly_io: symtab cmp_mpoly cmp_mpoly_io.a cmp_mpoly_ioB.a
	ada cmp_mpoly_io.a cmp_mpoly_ioB.a

# Complex polynomial systems :

cmp_polsys: cmp_mpoly cmp_polsys.a cmp_polsysB.a
	$(compile) cmp_polsys.a cmp_polsysB.a

cmp_laursys: cmp_mlaur cmp_laursys.a cmp_laursysB.a expvec.a expvecB.a
	$(compile) cmp_laursys.a cmp_laursysB.a expvec.a expvecB.a

# Converters between the polynomial systems :

pollaco: cmp_polsys cmp_laursys pollaco.a pollacoB.a lapolco.a lapolcoB.a
	$(compile) pollaco.a pollacoB.a lapolco.a lapolcoB.a

# Input/Output operations of complex polynomial systems :

cmp_polsys_io: cmp_polsys cmp_mpoly_io cmp_polsys_io.a cmp_polsys_ioB.a
	$(compile) cmp_polsys_io.a cmp_polsys_ioB.a

# The package Jacobi_Matrices :

jacobi: cmp_polsys cmp_laursys jacobi.a jacobiB.a laurjaco.a laurjacoB.a
	$(compile) jacobi.a jacobiB.a laurjaco.a laurjacoB.a

# Some tools: randomizers and substitutors :

tools: rand substits

rand: cmp_polsys cmp_laursys polrand.a polrandB.a laurrand.a laurrandB.a
	$(compile) polrand.a polrandB.a
	$(compile) laurrand.a laurrandB.a

substits: cmp_polsys substits.a substitsB.a
	$(compile) substits.a substitsB.a
