/********************************************************************************/
/*        GKLS-Generator of Classes of ND  (non-differentiable),                */
/*                                 D  (continuously differentiable), and        */
/*                                 D2 (twice continuously differentiable)       */
/*                     Test Functions for Global Optimization                   */
/*                                                                              */
/*   Authors:                                                                   */
/*                                                                              */
/*   M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev                           */
/*                                                                              */
/*   (C) 2002-2003                                                              */
/*                                                                              */
/*	 References:                                                            */
/*                                                                              */
/*   1. M.Gaviano, D.E.Kvasov, D.Lera, and Ya.D.Sergeyev (2003),                */
/*   Software for Generation of Classes of Test Functions                       */
/*   with Known Local and Global Minima for Global Optimization                 */
/*                                                                              */
/*   2. D.Knuth (1997), The Art of Computer Programming, Vol. 2:                */
/*   Seminumerical Algorithms (Third Edition). Reading, Massachusetts:          */
/*   Addison-Wesley.                                                            */
/*                                                                              */
/*   The software constructs a convex quadratic function (paraboloid) and then  */
/*   systematically distorts randomly selected parts of this function           */
/*   by polynomials in order to introduce local minima and to construct test    */
/*   functions which are non-differentiable (ND-type), continuously             */
/*   differentiable (D-type), and twice continuously differentiable (D2-type)   */
/*   at the feasible region.                                                    */
/*                                                                              */
/*   Each test class is defined by the following parameters:                    */
/*  (1) problem dimension                                                       */
/*  (2) number of local minima including the paraboloid min and the global min  */
/*  (3) global minimum value                                                    */
/*  (3) distance from the paraboloid vertex to the global minimizer             */
/*  (4) radius of the attraction region of the global minimizer                 */
/********************************************************************************/

This version is available also from the following WWW-address:
http://wwwinfo.deis.unical.it/~yaro

It includes: 
  \Package -- Directory containing the C code of the generator. 
              The package consists of the following files:
               gkls.c -- the main file;
	       gkls.h -- the header file that users should include in their 
	                 application projects in order to call subroutines from
	                 the file gkls.c;
	       rnd_gen.c -- the file containing the uniform random number
	                    generator proposed by D.Knuth (1997);
	       rnd_gen.h -- the header file for linkage to the file rnd_gen.c;
	       example.c -- an example of the GKLS-generator usage;
	       Makefile  -- an example of a UNIX makefile provided to UNIX users 
	                    for a simple compilation and linkage of separate files
	                    of the application project.

  UserManual.pdf -- A User Manual. 
  
  \Tex -- LaTeX-files of the User Manual.


