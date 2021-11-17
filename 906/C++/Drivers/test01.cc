#include 	<iostream>
#include 	<cmath>
#include 	<elrint3d.h>

using namespace std;

/*********************************************************************
  This example demonstrates the use of elrint3d for integration over 
  the unit cube. The integrand is a Gaussian function, defined using 
  a function pointer.  The requested relative error is 1.0e-12.
 *********************************************************************/

double gaussian(double x, double y, double z)
{
    // Normally, the following parameters would be
    // available as externally-defined variables

    double c1 = 0.2;
    double c2 = 0.5;
    double c3 = 0.8;
    double w1 = 0.4;
    double w2 = 0.7;
    double w3 = 0.25;

    return exp(- c1 * c1 * (x - w1)*(x - w1)
               - c2 * c2 * (y - w2)*(y - w2)
               - c3 * c3 * (z - w3)*(z - w3));
}

int main()
{
    Elrint3d I1(gaussian, 0, 1, 0, 1, 0, 1, 1e-12);
    cout.precision(17);
    cout << "Cubature        = " << I1.evaluate() << endl;
    cout << "Error Flag      = " << I1.errFlag() << endl;
    cout << "Estimated Error = " << I1.estErr() << endl;
    cout << "No of Fun Vals  = " << I1.evals() << endl;
    return 0;
}
