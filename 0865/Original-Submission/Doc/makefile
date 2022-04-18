f90= frt -X9 -Am -Kauto
f90= frt -Am -Kfast,eval,gs,preex
f90= f90  -fast 
f90= f90  -fast 
f90= frt -Am -Kfast,eval,gs,preex
f90= f90  -fast -ldxml
f90= frt -Hasu -X9 -Kauto -Am -Kcover
f90= ifort -O3
f90= f95  -O4  
#f90= f95  -C=all -V -dcfuns -nan -gline  -O
#f90= lf95 -O2
#f90= g95  

DBLAS= $(blas)dscal.f $(blas)dcopy.f $(blas)dsyrk.f $(lapack)dpotrf.f\
       $(blas)dgemm.f $(blas)dtrsm.f  $(blas)dgemv.f $(blas)xerbla.f \
       $(blas)lsame.f $(blas)ilaenv.f  $(lapack)dpotf2.f\
       $(blas)ddot.f $(blas)dtpsv.f  $(blas)dtrsv.f $(blas)dtpmv.f\
       /numerical/num/LAPACK3/BLAS/SRC/dspmv.f

all: lower_ex upper_ex kernel_ex full_test

speed:ddeps.o  
	$(f90) block_hybrid_Cholesky.f90 block_hybrid_speed.f90  ddeps.o 
	a.out 


full_test:ddeps.o full_test.o 
	$(f90) block_hybrid_Cholesky.f90 full_test.o ddeps.o 
	a.out < full_test.data >temp
	diff temp full_test.output

cover:ddeps.o all_together.f90 ddeps.f
#	nag_coverage95 all_together.f90 -first ddeps.o 
	nag_coverage95 all_together.f90 ddeps.f -first  
	all_together.inst.exe < full_test.data >temp
	diff temp program_results
#	gedit --new-window all_together.out all_together.listing.f90

lower_ex: 	ddeps.o  
	$(f90) block_hybrid_Cholesky.f90 lower_hybrid_ex.f90 ddeps.o 
	a.out <lower_hybrid_ex.data >temp
	diff temp lower_hybrid_ex.output

upper_ex: 	ddeps.o  
	$(f90) block_hybrid_Cholesky.f90 upper_hybrid_ex.f90 $(lapack)dpptrs.f ddeps.o 
	a.out <upper_hybrid_ex.data > temp
	diff temp upper_hybrid_ex.output

kernel_ex: 	ddeps.o  
	$(f90) block_hybrid_Cholesky.f90 kernel_ex.f90  ddeps.o 
	a.out <kernel_ex.data >temp
	diff temp kernel_ex.output

full_test.o:	full_test.f90
	 $(f90) -c full_test.f90
ddeps.o:	ddeps.f
	 $(f90) -c ddeps.f
block_hybrid_Cholesky.o:	block_hybrid_Cholesky.f90
	 $(f90) -c block_hybrid_Cholesky.f90

all_together.f90: block_hybrid_Cholesky.f90 full_test.f90
	cat  block_hybrid_Cholesky.f90 full_test.f90 > all_together.f90
ddeps.f:
	cat  $(DBLAS)  > ddeps.f

clean:
	rm a.out *.o  temp  *deps.f *.mod

lapack= /numerical/num/LAPACK3/SRC/
linpack= /numerical/num/linpack/
blas= /numerical/num/hsl2004/packages/blas/
