File README.txt

FRB FORTRAN package for the solution of Free Rigid Body equations
    Pi' = Pi x Omega
    Q' = Q skew(Omega)
where Pi := angular momentum, Omega := Angular velocity = aInert^{-1} Pi, Q := attitude matrix. 
aInert is the inertia tensor of the body.

1. System requirements

The software requires you have installed a recent FORTRAN compiler,
like FORTRAN 95 (f95, g95) or a compatible version (like gfortran), on which
it has been tested. 
The compiler can be downloaded for free at
    http://hpc.sourceforge.net/
    http://www.g95.org/
The software is compatible with FORTRAN77.

To compile, run the make command.
The compilation and linking are made using the'g95' compilator.
If you wish to use another fortran compilator, change (inline command)
the Fortran compilator FC variable.


2. Usage of the FRB FORTRAN package:

This package contains the following main FORTRAN subroutines
     frb_step.f
     quat_step.f
that perform the Free Rigid body computations. The first file uses
matrix-based rotations, while the second file uses quaternion-based
rotations, and a slightly different algorithm.

Subroutine calls:
Typical calls within another program:
	call frb_step(Piout, Qout, h, aInert, Piin, Qin, Np) 
	call quat_step(Piout, qout, h, aInert, Piin, qin, Np) 

Parameters description:
Piout := output vector (3x1) of angular momentum (of unit norm)
Qout := output attitude matrix (3x3) 
qout := output attitude quaternion vector (4x1)
h := stepsize of integration (scalar)
aInert:= vector (3x1) of inertia moments of the rigid body, with components in increasing
order
Piin := input vector (3x1) of angular momentum (of unit norm)
Qin := input attitude matrix (3x3)
qin := input attidude quaternion (4x1)
Np := integration order for the elliptic integral of third kind.
   Use Np=0 for exact integration
   Use 1<=Np<= 10 for Gaussian quadrature of order 2*Np.
   

N.B. For input angular momenta Pi of arbitrary norm, use Piin = Pi/norm(Pi) and 
h = norm(Pi)*h. The correct output is obtained as Piout*norm(Pi).


3. Compilation of driver examples:

The main computational subroutine need to be linked to a driver. We
have included the following example drivers:
     drivers			subroutine to be linked to
--------------------------------------------------------------------
     driver_example0.f     |	     frb_step.f
     driver_example1.f     |	     frb_step.f
     driver_example2.f     |	     quat_step.f
     driver_example3.f     |	     frb_step.f
     driver_example4.f     |	     quat_step.f
--------------------------------------------------------------------


4. All files (alphabetical order):

driver_example0.f	(example of FRB integration in a single step from 0 to Tfin). 
				 To be linked with frb_step.f
driver_example1.f	(example of rigid body with small torque, from 0 to Tfin, uses a 
				 splitting method). To be linked with frb_step.f
driver_example2.f   (Like driver_example0.f, but uses quaternions in place for rotations). 
				 To be linked with quat_step.f
driver_example3.f	(example for FRB integration with many steps
				 of lenght h from 0 to Tfin, using
				 semi-exact methods) 
				 To be linked with frb_step.f
driver_example4.f	(Like driver_example3.f but uses quaternions) 
				 To be linked with quat_step.f


frb_step.f		(Main computational subroutine, see above)

init_example0.dat 	(Data file for example 0)
init_example1.dat	(Data file for example 1)
init_example2.dat	(Data file for example 2)
init_example3.dat	(Data file for example 3)
init_example4.dat	(Data file for example 3)


makefile 		      	makefile compiling all the drivers
				using g95. See also 1. System Requirements 

out_example0.dat	(Output file for example 0)
out_example1.dat	(Output file for example 1)
out_example2.dat	(Output file for example 2)
out_example3.dat	(Output file for example 3)
out_example4.dat	(Output file for example 4)


quat_step.f		(Main computational subroutine, see above)


