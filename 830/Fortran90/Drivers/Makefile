#
# Makefile for
#  1.  battery tests on the updated *rotmg routines
#              make remark
#  2.  the downdating test
#              make downdate
#  3.  all the battery tests for both the remark and the new versions
#              make all
#
.SUFFIXES: .f90

# F90 should be set to the name of the Fortran 90 compiler to be used
F90         = f95
F90         = /opt/NAGWare/bin/f95

# F90OPTS should be set to compiler flags to be used
F90OPTS     = -C -g -u -xcheck=stkovf -ansi -XlistE
F90OPTS     = -g90 -gline -dcfuns -nan -save -C

# F90LINK and F90LINKOPTS should be set to the relevant linker and 
# linker options. For most systems these are the same as for the compiler
F90LINK     = $(F90)
F90LINKOPTS = $(F90OPTS)

.f90.o:
	$(F90) $(F90OPTS) -c $*.f90
.f.o:
	$(F90) $(F90OPTS) -c $*.f
D_ROTOBJS= Import_Kinds.o Givens_Rotations.o d_rot_battery.o
DROTOBJS= Import_Kinds.o blas.o drot_battery.o
D_ROTMOBJS= Import_Kinds.o Givens_Rotations.o d_rotm_battery.o
DROTMOBJS= Import_Kinds.o drotmg.o blas.o drotm_battery.o
DGROTMOBJS= Import_Kinds.o Givens_Rotations.o rotm_battery.f90
DGROTOBJS= Import_Kinds.o Givens_Rotations.o rot_battery.f90

S_ROTOBJS= Import_Kinds.o Givens_Rotations.o s_rot_battery.o
SROTOBJS= Import_Kinds.o blas.o srot_battery.o
S_ROTMOBJS= Import_Kinds.o Givens_Rotations.o s_rotm_battery.o
SROTMOBJS= Import_Kinds.o srotmg.o blas.o srotm_battery.o
SGROTMOBJS= Import_Kinds.o Givens_Rotations.o gsrotm_battery.f90
SGROTOBJS= Import_Kinds.o Givens_Rotations.o gsrot_battery.f90

D_ROTDATA= data1
DROTDATA= ../../Remark/Drivers/data2
D_ROTMDATA= data3
DROTMDATA= ../../Remark/Drivers/data4
DGROTMDATA= data3
DGROTDATA= data1

DDATE_OBJS=Import_Kinds.o Givens_Rotations.o DXrot.o drotmg.o blas.o
SDATE_OBJS=Import_Kinds.o Givens_Rotations.o SXrot.o srotmg.o blas.o

all: res1 res2 res3 res4 res5 res6 res7 res8 res9 res10 res11 res12\
     sdate_res ddate_res
remark: res4 res5 res10 res11
downdate: sdate_res ddate_res

sdate_res: $(SDATE_OBJS)
	$(F90LINK) -o sdate $(F90LINKOPTS) $(SDATE_OBJS)
	./sdate > sdate_res
ddate_res: $(DDATE_OBJS)
	$(F90LINK) -o ddate $(F90LINKOPTS) $(DDATE_OBJS)
	./ddate > ddate_res

res1: $(DGROTMOBJS) $(DGROTMDATA)
	$(F90LINK) $(F90LINKOPTS) $(DGROTMOBJS)
	./a.out < $(DGROTMDATA) > res1

res2: $(DGROTOBJS) $(DGROTDATA)
	$(F90LINK) $(F90LINKOPTS) $(DGROTOBJS)
	./a.out < $(DGROTDATA) > res2

res3: $(D_ROTOBJS) $(D_ROTDATA)
	$(F90LINK) $(F90LINKOPTS) $(D_ROTOBJS)
	./a.out < $(D_ROTDATA) > res3

res4: $(DROTOBJS) $(DROTDATA)
	$(F90LINK) $(F90LINKOPTS) $(DROTOBJS)
	./a.out < $(DROTDATA) > res4

res5: $(DROTMOBJS) $(DROTMDATA)
	$(F90LINK) $(F90LINKOPTS) $(DROTMOBJS)
	./a.out < $(DROTMDATA) > res5

res6: $(D_ROTMOBJS) $(D_ROTMDATA)
	$(F90LINK) $(F90LINKOPTS) $(D_ROTMOBJS)
	./a.out < $(D_ROTMDATA) > res6


res7: $(SGROTMOBJS) $(DGROTMDATA)
	$(F90LINK) $(F90LINKOPTS) $(SGROTMOBJS)
	./a.out < $(DGROTMDATA) > res7

res8: $(SGROTOBJS) $(DGROTDATA)
	$(F90LINK) $(F90LINKOPTS) $(SGROTOBJS)
	./a.out < $(DGROTDATA) > res8

res9: $(S_ROTOBJS) $(D_ROTDATA)
	$(F90LINK) $(F90LINKOPTS) $(S_ROTOBJS)
	./a.out < $(D_ROTDATA) > res9

res10: $(SROTOBJS) $(DROTDATA)
	$(F90LINK) $(F90LINKOPTS) $(SROTOBJS)
	./a.out < $(DROTDATA) > res10

res11: $(SROTMOBJS) $(DROTMDATA)
	$(F90LINK) $(F90LINKOPTS) $(SROTMOBJS)
	./a.out < $(DROTMDATA) > res11

res12: $(S_ROTMOBJS) $(D_ROTMDATA)
	$(F90LINK) $(F90LINKOPTS) $(S_ROTMOBJS)
	./a.out < $(D_ROTMDATA) > res12

Import_Kinds.o: ../Src/Import_Kinds.f90
	$(F90) $(F90OPTS) -c ../Src/Import_Kinds.f90

Givens_Rotations.o: Import_Kinds.o ../Src/Givens_Rotations.f90
	$(F90) $(F90OPTS) -c ../Src/Givens_Rotations.f90

DXrot.o: DXrot.f90
	$(F90) $(F90OPTS) -c DXrot.f90

SXrot.o: SXrot.f90
	$(F90) $(F90OPTS) -c SXrot.f90

gsrotm_battery.o:  gsrotm_battery.f90
	$(F90) $(F90OPTS) -c gsrotm_battery.f90

rotm_battery.o:  rotm_battery.f90
	$(F90) $(F90OPTS) -c rotm_battery.f90

rot_battery.o:  rot_battery.f90
	$(F90) $(F90OPTS) -c rot_battery.f90

d_rot_battery.o:  d_rot_battery.f90
	$(F90) $(F90OPTS) -c d_rot_battery.f90

s_rot_battery.o:  s_rot_battery.f90
	$(F90) $(F90OPTS) -c s_rot_battery.f90

d_rotm_battery.o:  d_rotm_battery.f90
	$(F90) $(F90OPTS) -c d_rotm_battery.f90

s_rotm_battery.o:  s_rotm_battery.f90
	$(F90) $(F90OPTS) -c s_rotm_battery.f90

gsrot_battery.o:  gsrot_battery.f90
	$(F90) $(F90OPTS) -c gsrot_battery.f90

blas.o: ../../Remark/Drivers/blas.f
	$(F90) $(F90OPTS) -c ../../Remark/Drivers/blas.f

drot_battery.o:  ../../Remark/Drivers/drot_battery.f90
	$(F90) $(F90OPTS) -c ../../Remark/Drivers/drot_battery.f90

drotm_battery.o:  ../../Remark/Drivers/drotm_battery.f90
	$(F90) $(F90OPTS) -c ../../Remark/Drivers/drotm_battery.f90

srot_battery.o:  ../../Remark/Drivers/srot_battery.f90
	$(F90) $(F90OPTS) -c ../../Remark/Drivers/srot_battery.f90

srotm_battery.o:  ../../Remark/Drivers/srotm_battery.f90
	$(F90) $(F90OPTS) -c ../../Remark/Drivers/srotm_battery.f90

drotmg.o: ../../Remark/Src/drotmg.f
	$(F90) $(F90OPTS) -c ../../Remark/Src/drotmg.f

srotmg.o: ../../Remark/Src/srotmg.f
	$(F90) $(F90OPTS) -c ../../Remark/Src/srotmg.f

clean:
	# remove executables and possible debris
	rm -rf ./a.out core ./sdate ./ddate
	# remove various module files -- depends on compiler
	rm -rf *.o *.g90 *.d *.mod
	# Remove work files for epcf90 compiler
	rm -rf work.pc work.pcl
