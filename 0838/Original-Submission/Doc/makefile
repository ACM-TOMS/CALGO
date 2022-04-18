ASOURCE = .
AFUNS = $(ASOURCE)/airy_functions.f90
FC = f95 
#FC = g95 
#FC = nag_f95 

all: test_install parameter_info sample_airy

airy_functions.o: $(ASOURCE)/airy_functions.f90 $(ASOURCE)/airy_real $(ASOURCE)/airy_complex $(ASOURCE)/airy_head $(ASOURCE)/airy_parameters
	$(FC) -c -o airy_functions.o $(ASOURCE)/airy_functions.f90

test_install: test_install.f90 airy_functions.o
	$(FC) -o test_install test_install.f90 airy_functions.o

clean:
	rm -f test_install airy_functions.o
	rm -f *.mod airyerr.log

%.o: %.f90
	$(FC) -c $< -o $@ 
