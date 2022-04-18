 Comments on the use of the DoubleIntegral class.
 -----------------------------------------------
 Written by Ian Robinson and Michael Hill
 Last updated: 23 April, 2001

 This file contains comments on, and examples of, the use of a C++ class
 called DoubleIntegral.  It includes a description of a C++ program called
 r2d2lri_sample.cpp for evaluating the given sample integrals.  Output from
 the program is written to a file called r2d2lri_sample.out.

 The main purpose of the DoubleIntegral class is to evaluate double integrals
 of the form:

          b   h(x)
        I   I     f(x,y) dy dx.
         a   g(x)

 The method used for the evaluation of DoubleIntegrals of this form is
 described in the paper written by Michael Hill and Ian Robinson: "r2d2lri:
 An algorithm for automatic two-dimensional cubature", (submitted for
 publication).  The definition of the class (along with a complete
 description of all the constructors and member functions) can be found in
 the file r2d2lri.h; its implementation is in the file r2d2lri.cpp.

 To evaluate a particular integral, one must first define the integral (in
 the form described below) and then use one of the "evaluate" functions
 contained in the file r2d2lri.cpp.  Including the directive

        #include "r2d2lri.h"

 in a program is sufficient to make the DoubleIntegral class available to
 that program.

 To use the class, integrals must be defined in the form:

   const double a = ...;
   const double b = ...;
   double g(double x) {return (some function of x);}
   double h(double x) {return (some function of x);}
   double f(double x, double y) {return (some function of x and y);}

 For instance, the integral

             1  x^2
       I1 = I  I     x.exp(y) dy dx
             0  0

 could be defined as follows:

   const double a1 = 0.0;
   const double b1 = 1.0;
   double g1(double x) {return 0.0;}
   double h1(double x) {return x*x;}
   double f1(double x, double y) {return x*exp(y);}

 Any of a, b, g(x) and h(x) may be infinite.  The value INFINITY = HUGE_VAL
 is defined in r2d2lri.h for this purpose.  (If use of HUGE_VAL causes problems
 with some compilers, the definition of INFINITY could be changed to a large
 value, such as 10^30.)  An example utilizing this value is the following
 integral over a semi-infinite domain:

               0   oo         dy dx
       I30 = I    I    ------------------- .
              -oo  0   sqrt(-xy).(y-x+1)^2

 This integral might be coded as follows:

   const double a30 = -INFINITY;
   const double b30 = 0.0;
   double g30(double x) {return 0.0;}
   double h30(double x) {return INFINITY;}
   double f30(double x, double y)
   {
      double z1 = y - x + 1.0;
      z1 = z1*z1;
      if (z1 == 0.0) return 0.0;
      double z2 = -x*y;
      return z2 <= 0.0 ? 0.0 : 1.0/(z1 * sqrt(z2));
   }

 (Integrals I1 and I30 are drawn from the Robinson and De Doncker testbed
 referred to in Hill, M. and Robinson, I., "d2lri: A non-adaptive algorithm
 for two-dimensional cubature", J.Comput.Appl.Math., 112, 1999, pp 121-145.)

 It must be remembered that access to the value INFINITY is available only to
 programs which contain the statement: #include "r2d2lri.h".

 There are a number of ways to now define the above sample integrals as
 members of the DoubleIntegral class in order to evaluate them.  Here, we
 describe two straightforward methods of doing so.  (Brief descriptions of
 all of the features of the class are given in the file r2d2lri.h.)

 1. The following two statements suffice to evaluate a stand-alone integral
    to, say, 10 significant figures:

          DoubleIntegral I(a1,b1,g1,h1,f1);
          value = I.evaluate(0.5E-10);

    Calls to I.rel_err_est() and I.error_flag() then provide an estimate of
    the relative error in the computed value of the integral and an indication
    of the likely success of the algorithm, respectively.  (A value other than
    zero for the error flag indicates that the requested accuracy may not have
    been achieved.)

    The same outcome to the above could be achieved with the statements:

          DoubleIntegral I;
          value = I.evaluate(a1,b1,g1,h1,f1,rel_err_est,error_flag,0.5E-10);

    In this case, the error estimate and error flag indicator are returned as
    rel_err_est and error_flag, respectively.

 2. To evaluate more than one integral to, say, 8-figure accuracy, one could
    proceed as follows:

          DoubleIntegral J(0.5E-8);
          J.set_new_integral(a1,b1,g1,h1,f1);
          J.evaluate();

           Output results for integral I1.

          J.set_new_integral(a30,b30,g30,h30,f30);
          J.evaluate();

           Output results for integral I30.

 A full set of functions (for example, set_outer_interval, set_rel_tol,
 get_max_evals, etc) is available for accessing private member values of
 an object of type DoubleIntegral and for altering these values (either
 individually or in sensible groupings).

 (CAUTION:

 In C++, the statement

        cout << I.evaluate() << I.rel_err_est() << endl;

 has a different outcome to that of the two statements

        cout << I.evaluate();
        cout << I.rel_err_est() << endl;

 In the first example, the member function rel_err_est() for the object I
 is called BEFORE the function evaluate() is called.  The opposite order of
 evaluation occurs in the second example.)

 The file r2d2lri_sample.cpp contains a program which provides several examples
 of the use of the DoubleIntegral class and its associated algorithm, r2d2lri.
 Note that the file r2d2lri.h contains #include directives for the math.h and
 stdlib.h library header files.



