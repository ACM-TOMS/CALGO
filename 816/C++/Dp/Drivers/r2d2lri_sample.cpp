//*****************************************************************************
// r2d2lri_sample.cpp
// Sample program demonstrating use of the DoubleIntegral class.
//
// Written by Ian Robinson and Michael Hill.
// Last updated 23 April, 2001.
//
// This program demonstrates use of various forms of the "evaluate" function
// and other public functions provided by the DoubleIntegral class to evaluate
// integrals 1, 21 and 30 from the Robinson and De Doncker testbed to various
// accuracies.  Refer to the README.txt file for comments on the set-up of this
// program.
//
// To run the program, ensure that the files r2d2lri.h and r2d2lri.cpp (or
// their corresponding object files) are in the same directory as this file.
// Then compile r2d2lri_sample.cpp and run the resulting executable file.
//
// Output from the program is directed to a file called r2d2lri_sample.out.
//
//*****************************************************************************
#include "r2d2lri.h"
#include <fstream.h>
#include <iomanip.h>

// Definitions of Integrals 1, 21 and 30 from the RD testbed.

// Integral 1.
const double a1 = 0.0;
const double b1 = 1.0;
double g1(double x) {return 0.0;}
double h1(double x) {return x*x;}
double f1(double x, double y) {return x*exp(y);}
const double exact1 = exp(1.0)/2.0 - 1.0;

// Integral 21.
const double a21 = 0.0;
const double b21 = 2.0;
double g21(double x) {return 0.0;}
double h21(double x)
{
   const double TWO_THIRDS = 2.0/3.0;
   double z = 1.0 - pow((x/2.0),1.5);
   return z <= 0.0 ? 0.0 : 3.0*pow(z,TWO_THIRDS);
};
double   f21(double x, double y)
{
   double z = x*y;
   return z == 0.0 ? 0.0 : pow(z,-0.1);
};
const double exact21 = 4.486951668283621;

// Integral 30.
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
const double exact30 = 4.0*atan(1.0);

int main()
{
  ofstream out;
  out.open("r2d2lri_sample.out");
  out << setiosflags(ios::scientific);

  out << "OUTPUT FILE FOR r2d2lri_sample.cpp" << endl << endl;

  // First, evaluate I30 to 2, 4 and 6 significant figures.

  out << "COMPUTATION OF I30 TO 2, 4, 6 AND 10 FIGURES USING r2d2lri"
      << endl << endl;

  out << "           0   oo         dy dx          " << endl
      << "   I30 = I    I    ------------------- . " << endl
      << "          -oo  0   sqrt(-xy).(y-x+1)^2   " << endl << endl;


  DoubleIntegral I(a30,b30,g30,h30,f30);
  out << "Requested accuracy = 2 significant figures" << endl;
  out << setprecision(14) << "The computed value of I30 is "
      << I.evaluate(0.5E-2) << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact30)/exact30 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations used was " << I.evals()
      << endl << endl;

  out << "Requested accuracy = 4 significant figures" << endl;
  out << setprecision(14) << "The computed value of I30 is "
      << I.evaluate(0.5E-4) << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact30)/exact30 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations was " << I.evals()
      << endl << endl;

  out << "Requested accuracy = 6 significant figures" << endl;
  out << setprecision(14) << "The computed value of I30 is "
      << I.evaluate(0.5E-6) << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact30)/exact30 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations was " << I.evals()
      << endl << endl;

  // Now try to evaluate I30 to 10 significant figures. (The algorithm is
  // unable to achieve this accuracy and sets the error flag accordingly.)

  out << "Requested accuracy = 10 significant figures" << endl;
  out << setprecision(14) << "The computed value of I30 is "
      << I.evaluate(0.5E-10) << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact30)/exact30 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations was " << I.evals()
      << endl << endl;

  // Now, evaluate I1 to 12 significant figures.

  DoubleIntegral I1;
  double value, error;
  int flag;

  out << endl << "COMPUTATION OF I1 TO 12 FIGURES" << endl << endl;

  out << "         1  x^2 " << endl
      << "   I1 = I  I    x.exp(y) dy dx " << endl
      << "         0  0 " << endl << endl;

  value = I1.evaluate(a1,b1,g1,h1,f1,error,flag,0.5E-12);
  out << "Requested accuracy = 12 significant figures" << endl;
  out << setprecision(14) << "The computed value of I1 is "
      << value << endl;
  out << setprecision(1) << "The estimated relative error is "
      << error << endl;
  out << "The actual relative error is "
      << fabs(value - exact1)/exact1 << endl;
  out << "The error flag has been set to " << flag << endl;
  out << "The number of function evaluations was " << I1.evals()
      << endl << endl;

  // Finally, evaluate each of the three integrals to 7 significant figures.

  out << endl << "COMPUTATION OF I30, I21 AND I1 TO 7 FIGURES" << endl << endl;

  out << "          2  h(x) " << endl
      << "   I21 = I  I    (xy)^(-0.1) dy dx " << endl
      << "          0  0 " << endl << endl
      << "   where " << endl
      << "      h(x) = 3(1 - (x/2)^1.5)^(2/3) " << endl << endl;

  I.set_rel_tol(0.5E-7);
  out << "Requested accuracy = 7 significant figures" << endl;
  out << setprecision(14) << "The computed value of I30 is "
      << I.evaluate() << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact30)/exact30 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations was " << I.evals()
      << endl << endl;

  I.set_new_integral(a21,b21,g21,h21,f21);
  out << "Requested accuracy = 7 significant figures" << endl;
  out << setprecision(14) << "The computed value of I21 is "
      << I.evaluate() << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact21)/exact21 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations was " << I.evals()
      << endl << endl;

  I = I1;
  out << "Requested accuracy = 7 significant figures" << endl;
  out << setprecision(14) << "The computed value of I1 is "
      << I.evaluate() << endl;
  out << setprecision(1) << "The estimated relative error is "
      << I.rel_err_est() << endl;
  out << "The actual relative error is "
      << fabs(I.value() - exact1)/exact1 << endl;
  out << "The error flag has been set to " << I.error_flag() << endl;
  out << "The number of function evaluations was " << I.evals() << endl;

  return 0;
}
