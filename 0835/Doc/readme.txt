
This Matlab package is developed by Zhonggang Zeng (zzeng@neiu.edu)
dated Jan. 6, 2004. This package is released through ACM Transaction 
on Mathematical Software. The author is not responsible for any damage
that may be related to using this package. 

1. General information on the package

There are three directories when the package is unzipped:

      documentation
      multroot
      testsuite

There are two papers in "documentation" directory:

  [1]. "Computing multiple roots of inexact polynomials", Z. Zeng,
       Math. Comp., to appear.

  [2]. "MultRoot - A Matlab package for computing polynomial roots
       and multiplicities", ACM Trans. Math. Software, to appear.

The paper [1] presents the overall theory and algorithm, while [2]
describes the software package and the test suite in detail. 

The directory "multroot" contains the Matlab modules of the MultRoot
package. 

The directory "testsuite" contains the Matlab modules that generate
test polynomials. 

2. Installation

Unzip and save all the files, open Matlab, set proper path to 
the directory containing the Matlab modules, the package is
now ready to run.  

3. Simple applications

To find roots and multiplicities of a polynomial p in Matlab format,
simply call

>> z = multroot(p)

The first column of z lists distinct roots, the 2nd column shows
the corresponding multiplicities.

At the current version, the code works for polynomials with 
moderate root multiplicities (such as 30).

Example:

>> f = poly([ones(1,20),2*ones(1,20),3*ones(1,10),4*ones(1,5)]);

>> z = multroot(f)

THE  STRUCTURE-PRESERVING CONDITION NUMBER:     	73.6821 
THE BACKWARD ERROR:                    	    8.63e-016 
THE ESTIMATED FORWARD ROOT ERROR:      	    1.27e-013 

        computed roots         	multiplicities

        4.000000000000002 	 	 	   5 
        2.999999999999999 	 	 	  10 
        2.000000000000000 	 	 	  20 
        1.000000000000000 	 	 	  20 


4. Sophisticated applications

There are five main programs:

1. MultRoot:	calculates roots and corresponding multiplicities,
		by calling the next two expert programs.


2. GcdRoot:   	calculates the multiplicity structure and initial 
               	estimates of the roots. In many cases, the estimates
               	are accurate enough.

		The minimum input is the polynomial coefficient vector

3. PejRoot:  	from given multiplicity structure and initial root estimates,
              	PejRoot calculate the root accurately.

                The minimum input is the polyn. coef. vector, initial 
                iterate and multiplicity structure.

4. mroot:       Same as multroot except that mroot returns a flag and 
                terminates when a nontrivial multiplicity structure is 
                not found. 

5. spcond:      calculates the structure-preserving condition number.


Expert users should read the papers attached in the package. At least
the paper [2] above. All modules contain "help" message that can be
accessed in Matlab. You may try

>> help multroot

Zhonggang Zeng
Associate Professor


