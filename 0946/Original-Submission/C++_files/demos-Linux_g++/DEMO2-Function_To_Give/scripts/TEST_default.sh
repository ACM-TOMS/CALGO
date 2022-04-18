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
# THIS SHELL SCRIPT ALLOWS TO RUN THE DEMO PROVIDING NO ARGUMENTS                                                #
##################################################################################################################


echo Running demo with all default values...

../test_functiontogive >screen_output.txt

mkdir output_demo2_default
mv *.txt output_demo2_default
