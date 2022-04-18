#! /bin/sh

#default values:
#	    tol= depends on the function        	required accuracy
#	    ntf= 1					number of functions to test RELIADIFF with
#	    nmax= 2000					maximum number of series coefficients
#	    pcoeff= 0					don't print coefficients
#           x= 1, 5, 10, 15       			evaluation points

#YOU COULD GIVE LESS THAN ALL ARGUMENTS: IF ONE OF THE FIRST THREE ARGUMENTS IS n THE PROGRAM USE THE DEFAULT VALUE; 
#IF THERE ARE LESS THAN 4 ARGUMENTS, THE PROGRAM USE THE DEFAULT EVALUATION POINTS. 

# THIS SCRIPT SHELL ALLOWS TO  TEST THE SOFTWARE ON 10 FUNCTIONS OF THE DATABASE 
# WITH A  TOLERANCE TOL=1.E-04. THE INVERSE FUNCTION IS COMPUTED 5 points directly given by user
# NMAX HAS THE DEFAULT VALUE: 2000

tol=1.0e-04
nmax=n		#default: 2000 coefficients max
pcoeff=n	#default: do not print Taylor coefficients
ntf=10		#10 function to give at runtime
range=0		#user will give evaluation points one by ones
dim=5		#5 evaluation points

echo Running test_on_database with third default value \(nmax\), with tol=$tol on $ntf func and at $dim given points \(you shall choose $ntf func from the database\)...

../test_on_database $tol $ntf $nmax $pcoeff $range $dim 2.0 4.5 6.3 8.5 10.0


mkdir output_demo1_10_tol-4
mv fz*files output_demo1_10_tol-4