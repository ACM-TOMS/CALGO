# Use gmake instead of make on Sun Solaris
# The default compiler name is f95, platform is unix, and fflags is -c
#
# USAGE EXAMPLES:
#
#    make demo
#    make check
#    make clean
#    make check-all-kinds
#    make dispmodule.obj
#    make check compiler=nag                     # also = gfortran, g95, ifort, absoft or pathf95
#    make check-naninf                           # works with g95, gfortran, ifort and Nag f95
#    make check-naninf-ieee                      # works with Nag f95 and newer versions of g95
#    make check-quadprec                         # works at least with ifort
#    make demo platform=windows compiler=ifort   # assumes gnu make is accessible (e.g. via cygwin)
#    make check debug=yes                        # use if test crashes, for more information
#    gmake check compiler=gfortran optimize=yes  # use gmake and gfortran -c -O2
#    make clean compiler=ifort                   # needed to clean for ifort (because obj defaults to .o)
#
# (alternatively edit the makefile)

# Set default compiler, platform, don't optimize or debug, compile flags = -c
compiler   = f95
platform   = unix
optimize   = no  # compile with -O2 flags?
debug      = no  # compile with debug flags?
fflags     = -c -O0
lnkflgs    =
ieeeflags  =

# Default file extensions are .o for object files, nothing for executables,
# use -o name to specify executable name:
obj  = .o
exe  =
out  = -o # (followed by space)
oout = -o #

compile   = $(compiler)
link      = $(compiler)

# The definition of add_ons must be changed according to the number of bytes in
# a default integer (and perhaps also in other situations)
dispmodule_src = dispmodule.f90
dispmodule_obj = dispmodule$(obj)
add_ons = disp_i1mod$(obj) disp_i2mod$(obj) disp_i8mod$(obj) disp_l1mod$(obj)
test_ikinds = test_kind_irange2$(exe) test_kind_irange4$(exe) \
              test_kind_irange9$(exe) test_kind_irange18$(exe)
test_rlkinds = test_kind_snglprec$(exe) test_kind_dbleprec$(exe) \
               test_kind_log1byte$(exe)

# Set flags and variables for different compilers
dflags = -c -g
lnkflgs =
ifeq ($(compiler),f95)
  dflags = -c -g
  ppflag = -fpp
endif
ifeq ($(compiler),nag)
  dflags = -c -g -C=all -gline
  ieeeflags = -ieee=full -w=x95
  ppflag = -fpp
  compile = f95
  link = f95
endif
ifeq ($(compiler),absoft)
  compile = af95
  link = af95
  add_ons = disp_i1mod_absoft$(obj)
  fflags = -c -O0 -Dabsoft
  test_ikinds = test_kind_irange2$(exe) test_kind_irange9$(exe)
  test_rlkinds = test_kind_snglprec$(exe) test_kind_dbleprec$(exe)
  dispmodule_src = dispmodule_absoft.f90
  dispmodule_obj = dispmodule_absoft$(obj)
endif
ifeq ($(compiler),g95)
  dflags = -c -g -ftrace=full
  ppflag = -cpp
endif
ifeq ($(compiler),gfortran)
  dflags = -c -g
  ppflag = -x f95-cpp-input
endif
ifeq ($(compiler),ifort)
  fflags = -c -nologo -O0
  dflags = $(fflags) -debug -traceback
  ppflag = -fpp
  lnkflgs = -nologo
  # (Od for no optimization)
endif

# Set flags for compilation (debugging or optimization)
ifeq ($(debug),yes)
  fflags = $(dflags)
endif
ifeq ($(optimize),yes)
  fflags = -c -O2
endif

# Change executeable file extension on windows (and obj and out for ifort)
ifeq ($(platform),windows)
  exe = .exe
  ifeq ($(compiler),ifort)
    obj = .obj
    out = -exe:
    oout = -object:
    fflags = -c -nologo -Od -stand:f95
    dflags = $(fflags) -debug -traceback
  endif
endif

objects = $(dispmodule_obj) $(add_ons)
otherobjects = test_dispmodule$(obj) test_naninf_mod$(obj) test_naninf$(obj) \
                test_naninf_ieee$(obj) dispdemo$(obj)

all : test_dispmodule$(exe) dispdemo$(exe) test_naninf$(exe)

demo : dispdemo$(exe)
	./dispdemo$(exe)

check : test_dispmodule$(exe)
	./test_dispmodule$(exe)

check-naninf : test_naninf$(exe)
	./test_naninf$(exe)

check-naninf-ieee : test_naninf_ieee$(exe)
	./test_naninf_ieee$(exe)

clean :
	rm -f test_dispmodule$(exe) dispdemo$(exe) test_naninf$(exe) test_naninf_ieee$(exe) \
	test_kind_quadprec$(exe) \
	$(test_ikinds) \
	$(test_rlkinds) \
	$(objects) \
	$(otherobjects)

check-all-kinds : $(test_ikinds) $(test_rlkinds)
	@./test_kind_irange2
	@./test_kind_irange4
	@./test_kind_irange9
	@./test_kind_irange18
	@./test_kind_snglprec
	@./test_kind_dbleprec
	@./test_kind_log1byte

check-all-absoft-kinds : $(test_ikinds) $(test_rlkinds)
	@./test_kind_irange2
	@./test_kind_snglprec
	@./test_kind_dbleprec

check-quadprec : test_kind_quadprec$(exe)
	@./test_kind_quadprec

test_dispmodule$(exe) : test_dispmodule$(obj) $(dispmodule_obj)
	$(link) $(lnkflgs) $(out)test_dispmodule$(exe) $< $(dispmodule_obj)

dispdemo$(exe) : dispdemo$(obj) $(dispmodule_obj)
	$(link) $(lnkflgs) $(out)$@ $< $(dispmodule_obj)

test_naninf$(exe) : test_naninf$(obj) test_naninf_mod$(obj) $(dispmodule_obj)
	$(link) $(lnkflgs) $(ieeeflags) $(out)test_naninf$(exe) $< test_naninf_mod$(obj) $(dispmodule_obj)

test_naninf_ieee$(exe) : test_naninf_ieee$(obj) test_naninf_mod$(obj) $(dispmodule_obj)
	$(link) $(lnkflgs) $(ieeeflags) $(out)test_naninf_ieee$(exe) $< test_naninf_mod$(obj) $(dispmodule_obj)

test_kind_quadprec$(exe) : test_dispmodule_fpp.F90 $(objects) disp_r16mod$(obj)
	$(compile) $(fflags) $(oout)test_kind_quadprec$(obj) -Dquadprec $(ppflag) $<
	$(link) $(lnkflgs) $(out)$@ test_kind_quadprec$(obj) disp_r16mod$(obj) $(objects)

test_kind_irange%$(exe) : test_dispmodule_fpp.F90 $(objects)
	$(compile) $(fflags) $(oout)test_kind_irange$*$(obj) -Dirange=$* $(ppflag) $<
	$(link) $(lnkflgs) $(out)$@ test_kind_irange$*$(obj) $(objects)

test_kind_%$(exe) : test_dispmodule_fpp.F90 $(objects)
	$(compile) $(fflags) $(oout)test_kind_$*$(obj) -D$* $(ppflag) $<
	$(link) $(lnkflgs) $(out)$@ test_kind_$*$(obj) $(objects)

$(dispmodule_obj) : $(dispmodule_src)
	$(compile) $(fflags) $<

.PRECIOUS : disp_%mod$(obj)

disp_%mod$(obj) : disp_%mod.f90 $(dispmodule_obj) # matches disp_i1mod, disp_i2mod etc.
	$(compile) $(fflags) $<

disp_%mod_absoft$(obj) : disp_%mod_absoft.f90 $(dispmodule_obj) # matches disp_i1mod_absoft
	$(compile) $(fflags) $<

test_dispmodule$(obj) : test_dispmodule.f90 $(dispmodule_obj)
	$(compile) $(fflags) $<

test_naninf_mod$(obj) : test_naninf_mod.f90 $(dispmodule_obj)
	$(compile) $(fflags) $<

test_naninf$(obj) : test_naninf.f90 test_naninf_mod$(obj)
	$(compile) $(fflags) $(ieeeflags) $<

test_naninf_ieee$(obj) : test_naninf_ieee.f90 test_naninf_mod$(obj)
	$(compile) $(fflags) $(ieeeflags) $<

dispdemo$(obj) : dispdemo.f90 $(dispmodule_obj)
	$(compile) $(fflags) $<
