# Makefile for the FRB fortran package

FC = /opt/NAGWare/bin/f95
FOPTS = -gline -u -dcfuns -nan -save -C -C=undefined

all: driver_example0 driver_example1 driver_example2 driver_example3 driver_example4

driver_example0 : driver_example0.f frb_step.f
	$(FC) $(FOPTS) -o driver_example0 driver_example0.f frb_step.f

driver_example1 : driver_example1.f frb_step.f
	$(FC) $(FOPTS) -o driver_example1 driver_example1.f frb_step.f

driver_example2 : driver_example2.f quat_step.f
	$(FC) $(FOPTS) -o driver_example2 driver_example2.f quat_step.f

driver_example3 : driver_example3.f frb_step.f
	$(FC) $(FOPTS) -o driver_example3 driver_example3.f frb_step.f

driver_example4 : driver_example4.f quat_step.f
	$(FC) $(FOPTS) -o driver_example4 driver_example4.f quat_step.f

clean:
	rm -f driver_example0 driver_example1 driver_example2 driver_example3
	rm -f driver_example4
