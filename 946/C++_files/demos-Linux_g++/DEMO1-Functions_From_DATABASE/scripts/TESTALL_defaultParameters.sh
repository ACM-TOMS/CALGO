#! /bin/sh

#default values:
#	    tol= depends on the function        	required accuracy
#	    ntf= 1					number of functions to test RELIADIFF with
#	    nmax= 2000					maximum number of series coefficients
#	    pcoeff= 0					don't print coefficients
#           x= 1, 5, 10, 15       			evaluation points

#YOU COULD GIVE LESS THAN ALL ARGUMENTS: IF ONE OF THE FIRST THREE ARGUMENTS IS n THE PROGRAM USE THE DEFAULT VALUE; 
#IF THERE ARE LESS THAN 4 ARGUMENTS, THE PROGRAM USE THE DEFAULT EVALUATION POINTS. 


# THIS SCRIPT ALLOWS TO  TEST THE SOFTWARE ON ALL FUNCTIONS OF THE DATABASE 
# WITHOUT PROVIDING ANY ARGUMENT OTHER THAN THE MAX NUMBER OF FUNCTIONS. THE INVERSE FUNCTIONS ARE COMPUTED ON 5 DEFAULT POINTS
# NMAX HAS THE DEFAULT VALUE: 2000


tol=n
ntf=a 		#all database functions
#pcoeff=n	
#nmax=n
#range=1
#a=1
#b=15
#step=0.5

echo Running test_on_database with all default values on all database functions

../test_on_database $tol $ntf >screen_output.txt

mkdir output_demo1_all_default
mv *.txt fz*files output_demo1_all_default