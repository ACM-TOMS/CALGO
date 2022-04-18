all: cmdline timing byemess communication file_ops climp

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
	@-../makelib

# Cleaning the Ada library :

clean:
	a.rmlib -f

# Cleaning the .imports directory :

climp:
	@-rm -f .imports/*

# command line package :

cmdline: ada.lib ucmdline.a ucmdlineB.a
	$(compile) ucmdline.a ucmdlineB.a

# The timing package :

timing: ada.lib rusage.a rusageB.a timing.a timingB.a
	$(compile) rusage.a rusageB.a timing.a timingB.a

# System calls :

machines: ada.lib system_call.a system_callB.a machines.a machinesB.a
	$(compile) system_call.a system_callB.a
	$(compile) machines.a machinesB.a

byemess: machines byemess.a byemessB.a
	$(compile) byemess.a byemessB.a

# Communications with user :

communication: ada.lib commuser.a commuserB.a
	$(compile) commuser.a commuserB.a

# file operations :

file_ops: ada.lib file_ops.a file_opsB.a
	$(compile) file_ops.a file_opsB.a

# Making the links with the original sources :

linkrc = ../../System

links:
	@-ln -s $(linkrc)/byemess.a      byemess.a
	@-ln -s $(linkrc)/byemessB.a     byemessB.a
	@-ln -s $(linkrc)/commuser.a     commuser.a
	@-ln -s $(linkrc)/commuserB.a    commuserB.a
	@-ln -s $(linkrc)/rusage.a       rusage.a
	@-ln -s $(linkrc)/rusageB.a      rusageB.a
	@-ln -s $(linkrc)/timing.a       timing.a
	@-ln -s $(linkrc)/timingB.a      timingB.a
	@-ln -s $(linkrc)/ucmdline.a     ucmdline.a
	@-ln -s $(linkrc)/ucmdlineB.a    ucmdlineB.a
	@-ln -s $(linkrc)/system_call.a  system_call.a
	@-ln -s $(linkrc)/system_callB.a system_callB.a
	@-ln -s $(linkrc)/machines.a     machines.a
	@-ln -s $(linkrc)/machinesB.a    machinesB.a
	@-ln -s $(linkrc)/file_ops.a     file_ops.a
	@-ln -s $(linkrc)/file_opsB.a    file_opsB.a
