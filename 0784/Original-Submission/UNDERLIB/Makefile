include ../make.inc

### Underlying BLAS ####################################################

ULIB =   $(CREULIB)

DGBS =	 saxpy.f scopy.f sgemm.f sgemv.f sger.f  sscal.f ssyr.f \
         strmv.f strsv.f daxpy.f dcopy.f dcabs1.f dgemm.f dgemv.f \
         dger.f dscal.f dsyr.f  dtrmv.f dtrsv.f caxpy.f ccopy.f \
         cgemm.f cgemv.f cher.f  cscal.f ctrmv.f ctrsv.f zaxpy.f \
         zcopy.f zgemm.f zgemv.f zher.f  zscal.f ztrmv.f ztrsv.f \
         lsame.f xerbla.f

DGB =	 saxpy.o scopy.o sgemm.o sgemv.o sger.o  sscal.o ssyr.o \
         strmv.o strsv.o daxpy.o dcopy.o dcabs1.o dgemm.o dgemv.o \
         dger.o dscal.o dsyr.o  dtrmv.o dtrsv.o caxpy.o ccopy.o \
         cgemm.o cgemv.o cher.o  cscal.o ctrmv.o ctrsv.o zaxpy.o \
         zcopy.o zgemm.o zgemv.o zher.o  zscal.o ztrmv.o ztrsv.o \
         lsame.o xerbla.o

########################################################################
all: $(CREULIB)

$(CREULIB):	$(DGB)
	$(ARCH) $(ARCHFLAGS) $(ULIB) $(DGB)
	$(RANLIB) $(ULIB)

$(DGB):		$(DGBS)
	$(FORTRAN) -c $(UNDEROPT) $(DGBS)

clean:
	rm -f *.o
