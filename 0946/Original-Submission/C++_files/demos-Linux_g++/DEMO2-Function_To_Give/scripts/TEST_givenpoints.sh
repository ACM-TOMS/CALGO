#! /bin/sh

# default values:
#	    tol= depends on the function        	required accuracy
#	    sigma0= 1					convergence abscissa for F
#	    nmax= 2000					maximum number of series coefficients
#	    pcoeff= 0					don't print coefficients
#	    szero=0					the Transform has not a singularity at zero
#           x= 1, 5, 10, 15       			evaluation points

#YOU COULD GIVE LESS THAN ALL ARGUMENTS: IF ONE OF THE FIRST THREE ARGUMENTS IS n THE PROGRAM USE THE DEFAULT VALUE; 
#IF THERE ARE LESS THAN 5 ARGUMENTS, THE PROGRAM USE THE DEFAULT EVALUATION POINTS. 

##################################################################################################################
# THIS SHELL SCRIPT ALLOWS TO RUN THE DEMO WITH FIRST 3 DEFAULT VALUE IN ARGUMENTS AND 4 GIVEN EVALUATION POINTS #
##################################################################################################################

echo Running demo with first for default values on four points given by user...

range=0
dim=3

../test_functiontogive n n n n n $range $dim 1 1.8 2.5 >screen_output.txt

mkdir output_demo2_givenpoints
mv *.txt output_demo2_givenpoints
