#! /bin/sh

#default values:
#	    tol= depends on the function        	required accuracy
#	    ntf= 1					number of functions to test RELIADIFF with
#	    nmax= 2000					maximum number of series coefficients
#	    pcoeff= 0					don't print coefficients
#           x= 1, 5, 10, 15       			evaluation points

#YOU COULD GIVE LESS THAN ALL ARGUMENTS: IF ONE OF THE FIRST THREE ARGUMENTS IS n THE PROGRAM USE THE DEFAULT VALUE; 
#IF THERE ARE LESS THAN 4 ARGUMENTS, THE PROGRAM USE THE DEFAULT EVALUATION POINTS. 

echo Running test_driver with first three default values on four points given by user \(you shall choose a func from the database\)...

range=0		#user will give points one by ones
dim=4

../test_on_database n n n n $range $dim 1 1.5 1.8 2.5

mkdir output_demo1_one_givenpoints
mv fz*files output_demo1_one_givenpoints