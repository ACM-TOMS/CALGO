#! /bin/sh

#default values:
#	    tol= depends on the function        	required accuracy
#	    ntf= 1					number of functions to test RELIADIFF with
#	    nmax= 2000					maximum number of series coefficients
#	    pcoeff= 0					don't print coefficients
#           x= 1, 5, 10, 15       			evaluation points

#YOU COULD GIVE LESS THAN ALL ARGUMENTS: IF ONE OF THE FIRST THREE ARGUMENTS IS n THE PROGRAM USE THE DEFAULT VALUE; 
#IF THERE ARE LESS THAN 4 ARGUMENTS, THE PROGRAM USE THE DEFAULT EVALUATION POINTS. 


echo Running test_driver with all default values \(you shall choose a func from the database\)...

../test_on_database

mkdir output_demo1_one_default
mv fz*files output_demo1_one_default