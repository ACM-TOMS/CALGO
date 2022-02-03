#! /bin/sh

#default values:
#	    tol= depends on the function        	required accuracy
#	    ntf= 1					number of functions to test RELIADIFF with
#	    nmax= 2000					maximum number of series coefficients
#	    pcoeff= 0					don't print coefficients
#           x= 1, 5, 10, 15       			evaluation points

#YOU COULD GIVE LESS THAN ALL ARGUMENTS: IF ONE OF THE FIRST THREE ARGUMENTS IS n THE PROGRAM USE THE DEFAULT VALUE; 
#IF THERE ARE LESS THAN 4 ARGUMENTS, THE PROGRAM USE THE DEFAULT EVALUATION POINTS. 


# THIS SCRIPT SHELL ALLOWS TO  TEST THE SOFTWARE ON ALL FUNCTIONS OF THE DATABASE 
# WITH A  TOLERANCE TOL=1.E-04. THE INVERSE FUNCTION IS COMPUTED ON [0.5,15] WITH STEP SIZE=0.5
# NMAX HAS THE DEFAULT VALUE: 2000




tol=1.0e-04
nmax=n		#default: 2000 coefficients max
pcoeff=1	#print the involved Taylor coefficients
ntf=a 		#all database functions
range=1		#user will give points in a range
a=0.5
b=15
step=0.5

../test_on_database $tol $ntf $nmax $pcoeff $range $a $b $step >screen_output.txt

mkdir output_demo1_all_tol-4
mv  *.txt output_demo1_all_tol-4
mv  fz*files output_demo1_all_tol-4
