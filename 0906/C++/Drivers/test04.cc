#include 	<iostream>
#include 	<cmath>
#include 	<elrint3d.h>

using namespace std;

/*********************************************************************
  This example demonstrates the use of elrint3d for integration over 
  a variable domain.  The integrand is a Gaussian function, defined
  using a class.  The requested relative error is 1.0e-12.
 *********************************************************************/


class GaussianIntegrand : public Integrand<double>
{
private:
    double c1, c2, c3;
    double w1, w2, w3;

public:
    GaussianIntegrand(double C1, double C2, double C3, double W1, double W2, double W3) :
    c1(C1), c2(C2), c3(C3), w1(W1), w2(W2), w3(W3)
    {
    }
    
    double fun(const double x[]) const
    {
        return exp(- c1 * c1 * (x[0] - w1)*(x[0] - w1)
                   - c2 * c2 * (x[1] - w2)*(x[1] - w2)
                   - c3 * c3 * (x[2] - w3)*(x[2] - w3));
    }
};

//Define the lower boundary in the y direction

double yLower(double x)
{
    return 0;
}
//Define the upper boundary in the y direction

double yUpper(double x)
{
    return 1 - x;
}
//Define the lower boundary in the z direction

double zLower(double x, double y)
{
    return 0;
}
//Define the upper boundary in the z direction

double zUpper(double x, double y)
{
    return 1 - x - y;
}

int main()
{
    GaussianIntegrand f(0.2, 0.5, 0.8, 0.4, 0.7, 0.25);
    Elrint3d I2(f, 0, 1, yLower, yUpper, zLower, zUpper, 1e-12);
    cout.precision(17);
    cout << "Cubature        = " << I2.evaluate() << endl;
    cout << "Error Flag      = " << I2.errFlag() << endl;
    cout << "Estimated Error = " << I2.estErr() << endl;
    cout << "No of Fun Vals  = " << I2.evals() << endl;
    return 0;
}
