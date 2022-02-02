all: vec vec_io mat_io lss vvc vvc_io climp

vec: vectors nat_vec int_vec flt_vec cmp_vec cmp_norms
vec_io: vec nat_vec_io int_vec_io flt_vec_io cmp_vec_io
mat: vec nat_mat int_mat flt_mat cmp_mat
mat_io: mat int_mat_io flt_mat_io cmp_mat_io
lss: mat gcd int_lss int_lis flt_lss cmp_lss
vvc: vec nat_vvc int_vvc flt_vvc cmp_vvc
vvc_io: vvc nat_vvc_io int_vvc_io flt_vvc_io cmp_vvc_io

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
	a.path -a ../Numbers

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Math_Lib/Matrices

links:
	@-ln -s $(linkrc)/gcd.a          gcd.a
	@-ln -s $(linkrc)/gcdB.a         gcdB.a
	@-ln -s $(linkrc)/vectors.a      vectors.a
	@-ln -s $(linkrc)/vectorsB.a     vectorsB.a
	@-ln -s $(linkrc)/nat_vec.a      nat_vec.a
	@-ln -s $(linkrc)/nat_vec_io.a   nat_vec_io.a
	@-ln -s $(linkrc)/nat_vec_ioB.a  nat_vec_ioB.a
	@-ln -s $(linkrc)/int_vec.a      int_vec.a
	@-ln -s $(linkrc)/int_vec_io.a   int_vec_io.a
	@-ln -s $(linkrc)/int_vec_ioB.a  int_vec_ioB.a
	@-ln -s $(linkrc)/flt_vec.a      flt_vec.a
	@-ln -s $(linkrc)/flt_vec_io.a   flt_vec_io.a
	@-ln -s $(linkrc)/flt_vec_ioB.a  flt_vec_ioB.a
	@-ln -s $(linkrc)/cmp_vec.a      cmp_vec.a
	@-ln -s $(linkrc)/cmp_vec_io.a   cmp_vec_io.a
	@-ln -s $(linkrc)/cmp_vec_ioB.a  cmp_vec_ioB.a
	@-ln -s $(linkrc)/cmp_norms.a    cmp_norms.a
	@-ln -s $(linkrc)/cmp_normsB.a   cmp_normsB.a
	@-ln -s $(linkrc)/nat_vvc.a      nat_vvc.a
	@-ln -s $(linkrc)/nat_vvc_io.a   nat_vvc_io.a
	@-ln -s $(linkrc)/nat_vvc_ioB.a  nat_vvc_ioB.a
	@-ln -s $(linkrc)/nat_mat.a      nat_mat.a
	@-ln -s $(linkrc)/nat_matB.a     nat_matB.a
	@-ln -s $(linkrc)/nat_mat_io.a   nat_mat_io.a
	@-ln -s $(linkrc)/nat_mat_ioB.a  nat_mat_ioB.a
	@-ln -s $(linkrc)/int_mat.a      int_mat.a
	@-ln -s $(linkrc)/int_matB.a     int_matB.a
	@-ln -s $(linkrc)/int_mat_io.a   int_mat_io.a
	@-ln -s $(linkrc)/int_mat_ioB.a  int_mat_ioB.a
	@-ln -s $(linkrc)/int_lss.a      int_lss.a
	@-ln -s $(linkrc)/int_lssB.a     int_lssB.a
	@-ln -s $(linkrc)/int_lis.a      int_lis.a
	@-ln -s $(linkrc)/int_lisB.a     int_lisB.a
	@-ln -s $(linkrc)/int_vvc.a      int_vvc.a
	@-ln -s $(linkrc)/int_vvc_io.a   int_vvc_io.a
	@-ln -s $(linkrc)/int_vvc_ioB.a  int_vvc_ioB.a
	@-ln -s $(linkrc)/flt_mat.a      flt_mat.a
	@-ln -s $(linkrc)/flt_matB.a     flt_matB.a
	@-ln -s $(linkrc)/flt_mat_io.a   flt_mat_io.a
	@-ln -s $(linkrc)/flt_mat_ioB.a  flt_mat_ioB.a
	@-ln -s $(linkrc)/flt_lss.a      flt_lss.a
	@-ln -s $(linkrc)/flt_lssB.a     flt_lssB.a
	@-ln -s $(linkrc)/flt_vvc.a      flt_vvc.a
	@-ln -s $(linkrc)/flt_vvc_io.a   flt_vvc_io.a
	@-ln -s $(linkrc)/flt_vvc_ioB.a  flt_vvc_ioB.a
	@-ln -s $(linkrc)/cmp_mat.a      cmp_mat.a
	@-ln -s $(linkrc)/cmp_matB.a     cmp_matB.a
	@-ln -s $(linkrc)/cmp_mat_io.a   cmp_mat_io.a
	@-ln -s $(linkrc)/cmp_mat_ioB.a  cmp_mat_ioB.a
	@-ln -s $(linkrc)/cmp_lss.a      cmp_lss.a
	@-ln -s $(linkrc)/cmp_lssB.a     cmp_lssB.a
	@-ln -s $(linkrc)/cmp_vvc.a      cmp_vvc.a
	@-ln -s $(linkrc)/cmp_vvc_io.a   cmp_vvc_io.a
	@-ln -s $(linkrc)/cmp_vvc_ioB.a  cmp_vvc_ioB.a

# Greatest common divisor :

gcd: ada.lib gcd.a gcdB.a
	$(compile) gcd.a gcdB.a

# The vectors package :

vectors: ada.lib vectors.a vectorsB.a
	$(compile) vectors.a vectorsB.a

# Instantiation of natural, integer, float and complex vectors :

nat_vec: vectors nat_vec.a
	$(compile) nat_vec.a

int_vec: vectors int_vec.a
	$(compile) int_vec.a

flt_vec: vectors flt_vec.a
	$(compile) flt_vec.a

cmp_vec: vectors cmp_vec.a
	$(compile) cmp_vec.a

# Norms of complex vectors :

cmp_norms: cmp_vec cmp_norms.a cmp_normsB.a
	$(compile) cmp_norms.a cmp_normsB.a

# Input/Output operations of all the vector packages :

nat_vec_io: nat_vec nat_vec_io.a nat_vec_ioB.a
	$(compile) nat_vec_io.a nat_vec_ioB.a

int_vec_io: int_vec int_vec_io.a int_vec_ioB.a
	$(compile) int_vec_io.a int_vec_ioB.a

flt_vec_io: flt_vec flt_vec_io.a flt_vec_ioB.a
	$(compile) flt_vec_io.a flt_vec_ioB.a

cmp_vec_io: cmp_vec cmp_vec_io.a cmp_vec_ioB.a
	$(compile) cmp_vec_io.a cmp_vec_ioB.a

# Matrices with natural, integer, float and complex components :

nat_mat: nat_vec nat_mat.a nat_matB.a
	$(compile) nat_mat.a nat_matB.a

int_mat: int_vec int_mat.a int_matB.a
	$(compile) int_mat.a int_matB.a

flt_mat: flt_vec flt_mat.a flt_matB.a
	$(compile) flt_mat.a flt_matB.a

cmp_mat: cmp_vec cmp_mat.a cmp_matB.a
	$(compile) cmp_mat.a cmp_matB.a

# Input/Output operations of all the vector packages :

nat_mat_io: nat_mat nat_mat_io.a nat_mat_ioB.a
	$(compile) nat_mat_io.a nat_mat_ioB.a

int_mat_io: int_mat int_mat_io.a int_mat_ioB.a
	$(compile) int_mat_io.a int_mat_ioB.a

flt_mat_io: flt_mat flt_mat_io.a flt_mat_ioB.a
	$(compile) flt_mat_io.a flt_mat_ioB.a

cmp_mat_io: cmp_mat cmp_mat_io.a cmp_mat_ioB.a
	$(compile) cmp_mat_io.a cmp_mat_ioB.a

# Linear system solvers :

int_lss: int_mat gcd int_lss.a int_lssB.a
	ada int_lss.a int_lssB.a

flt_lss: flt_mat int_vec flt_lss.a flt_lssB.a
	$(compile) flt_lss.a flt_lssB.a

cmp_lss: cmp_mat int_vec cmp_lss.a cmp_lssB.a
	$(compile) cmp_lss.a cmp_lssB.a

# Linear inequality solvers :

int_lis: int_lss int_lis.a int_lisB.a
	$(compile) int_lis.a int_lisB.a

# Vectors of Vectors :

nat_vvc: nat_vec nat_vvc.a
	$(compile) nat_vvc.a

int_vvc: int_vec int_vvc.a
	$(compile) int_vvc.a

flt_vvc: flt_vec flt_vvc.a
	$(compile) flt_vvc.a

cmp_vvc: cmp_vec cmp_vvc.a
	$(compile) cmp_vvc.a

# Input/Output operations of all the vectors of vectors packages :

nat_vvc_io: nat_vvc nat_vvc_io.a nat_vvc_ioB.a
	$(compile) nat_vvc_io.a nat_vvc_ioB.a

int_vvc_io: int_vvc int_vvc_io.a int_vvc_ioB.a
	$(compile) int_vvc_io.a int_vvc_ioB.a

flt_vvc_io: flt_vvc flt_vvc_io.a flt_vvc_ioB.a
	$(compile) flt_vvc_io.a flt_vvc_ioB.a

cmp_vvc_io: cmp_vvc cmp_vvc_io.a cmp_vvc_ioB.a
	$(compile) cmp_vvc_io.a cmp_vvc_ioB.a

