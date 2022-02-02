all: lp facenum fltfanum data climp

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
	a.path -a ../Polynomials

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# Establishing the links :

linkrc = ../../../Math_Lib/Supports

links:
	@-ln -s $(linkrc)/dictios.a       dictios.a
	@-ln -s $(linkrc)/dictiosB.a      dictiosB.a
	@-ln -s $(linkrc)/lp.a            lp.a
	@-ln -s $(linkrc)/lpB.a           lpB.a
	@-ln -s $(linkrc)/givens.a        givens.a
	@-ln -s $(linkrc)/givensB.a       givensB.a
	@-ln -s $(linkrc)/farkas.a        farkas.a
	@-ln -s $(linkrc)/farkasB.a       farkasB.a
	@-ln -s $(linkrc)/int_farkas.a    int_farkas.a
	@-ln -s $(linkrc)/int_farkasB.a   int_farkasB.a
	@-ln -s $(linkrc)/facenu_ut.a     facenu_ut.a
	@-ln -s $(linkrc)/facenu_utB.a    facenu_utB.a
	@-ln -s $(linkrc)/facenum.a       facenum.a
	@-ln -s $(linkrc)/facenumB.a      facenumB.a
	@-ln -s $(linkrc)/flt_lis.a       flt_lis.a
	@-ln -s $(linkrc)/flt_lisB.a      flt_lisB.a
	@-ln -s $(linkrc)/fltintro.a      fltintro.a
	@-ln -s $(linkrc)/fltintroB.a     fltintroB.a
	@-ln -s $(linkrc)/fltfanum.a      fltfanum.a
	@-ln -s $(linkrc)/fltfanumB.a     fltfanumB.a
	@-ln -s $(linkrc)/lstivc.a        lstivc.a
	@-ln -s $(linkrc)/lstivcB.a       lstivcB.a
	@-ln -s $(linkrc)/lstivc_io.a     lstivc_io.a
	@-ln -s $(linkrc)/lstivc_ioB.a    lstivc_ioB.a
	@-ln -s $(linkrc)/arrlivc.a       arrlivc.a
	@-ln -s $(linkrc)/arrlivcB.a      arrlivcB.a
	@-ln -s $(linkrc)/arrlivc_io.a    arrlivc_io.a
	@-ln -s $(linkrc)/arrlivc_ioB.a   arrlivc_ioB.a
	@-ln -s $(linkrc)/intsupfu.a      intsupfu.a
	@-ln -s $(linkrc)/intsupfuB.a     intsupfuB.a
	@-ln -s $(linkrc)/intfaces.a      intfaces.a
	@-ln -s $(linkrc)/intfacesB.a     intfacesB.a
	@-ln -s $(linkrc)/intfaces_io.a   intfaces_io.a
	@-ln -s $(linkrc)/intfaces_ioB.a  intfaces_ioB.a
	@-ln -s $(linkrc)/lstfvc.a        lstfvc.a
	@-ln -s $(linkrc)/lstfvcB.a       lstfvcB.a
	@-ln -s $(linkrc)/lstfvc_io.a     lstfvc_io.a
	@-ln -s $(linkrc)/lstfvc_ioB.a    lstfvc_ioB.a
	@-ln -s $(linkrc)/arrlfvc.a       arrlfvc.a
	@-ln -s $(linkrc)/arrlfvcB.a      arrlfvcB.a
	@-ln -s $(linkrc)/arrlfvc_io.a    arrlfvc_io.a
	@-ln -s $(linkrc)/arrlfvc_ioB.a   arrlfvc_ioB.a
	@-ln -s $(linkrc)/fltsupfu.a      fltsupfu.a
	@-ln -s $(linkrc)/fltsupfuB.a     fltsupfuB.a
	@-ln -s $(linkrc)/fltfaces.a      fltfaces.a
	@-ln -s $(linkrc)/fltfacesB.a     fltfacesB.a

# Linear programming and the simplex method :

lp: ada.lib dictios.a dictiosB.a lp.a lpB.a
	$(compile) dictios.a dictiosB.a lp.a lpB.a

# Farkas lemma and complementary slackness :

farkas: ada.lib givens.a givensB.a farkas.a farkasB.a flt_lis.a flt_lisB.a
	ada givens.a givensB.a flt_lis.a flt_lisB.a
	ada farkas.a farkasB.a

int_farkas: farkas int_farkas.a int_farkasB.a fltintro.a fltintroB.a
	$(compile) int_farkas.a int_farkasB.a
	ada fltintro.a fltintroB.a

# Enumeration of faces :

facenum: int_farkas facenu_ut.a facenu_utB.a facenum.a facenumB.a
	$(compile) facenu_ut.a facenu_utB.a facenum.a facenumB.a

fltfanum: farkas fltfanum.a fltfanumB.a
	$(compile) fltfanum.a fltfanumB.a

# data structures:

lstivc: ada.lib lstivc.a lstivcB.a lstivc_io.a lstivc_ioB.a
	$(compile) lstivc.a lstivcB.a lstivc_io.a lstivc_ioB.a

arrlivc: lstivc arrlivc.a arrlivcB.a arrlivc_io.a arrlivc_ioB.a
	$(compile) arrlivc.a arrlivcB.a arrlivc_io.a arrlivc_ioB.a

intfaces: lstivc intfaces.a intfacesB.a intfaces_io.a intfaces_ioB.a
	$(compile) intfaces.a intfacesB.a intfaces_io.a intfaces_ioB.a

intsupfu: lstivc intsupfu.a intsupfuB.a
	$(compile) intsupfu.a intsupfuB.a

lstfvc: ada.lib lstfvc.a lstfvcB.a lstfvc_io.a lstfvc_ioB.a
	$(compile) lstfvc.a lstfvcB.a lstfvc_io.a lstfvc_ioB.a

arrlfvc: lstfvc arrlfvc.a arrlfvcB.a arrlfvc_io.a arrlfvc_ioB.a
	$(compile) arrlfvc.a arrlfvcB.a arrlfvc_io.a arrlfvc_ioB.a

fltfaces: lstfvc fltsupfu.a fltsupfuB.a fltfaces.a fltfacesB.a
	$(compile) fltsupfu.a fltsupfuB.a fltfaces.a fltfacesB.a

data: arrlivc intfaces intsupfu arrlfvc fltfaces
