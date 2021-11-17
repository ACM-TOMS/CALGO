
// Example suggested by a referee, who called this `reasonable'.

// One can easily see, e.g. by substituting y=1, that this
// integral involves severe difficulties such as fast oscillations
// in the neighbourhood of the line x+y=2, and discontinuous derivatives.


#include <cubpack.h>
#include <iostream.h>

#define sqr(x) ((x)*(x))

static const real sqrt_PI = sqrt(M_PI),
                  p35_PI  = pow(M_PI,0.35e0);

real f(const Point& p)
 { real x=p.X() , y=p.Y() ;
     return (  sqrt( sqr(x-sqrt_PI) + sqr(y-1) + 0.00001)
              + pow(pow(y-p35_PI,4) + sqr(x-1) + 0.00002,1.0/3.0)
            ) * sin( (x-y-sqrt_PI+p35_PI)/(sqr(x+y-2) + 0.01) );
 }

int main ()
 {
   TRIANGLE Roof( Point(0,1), Point(2,1), Point(1,2));
   RECTANGLE Walls( Point(0,0), Point(0,1), Point(2,0));
   REGION_COLLECTION House = Walls + Roof;
   EvaluationCounter count;

   count.Start();
   cout <<"The integral is " << Integrate(f,House,0,0.5e-1,1000000);
   cout <<" with absolute error " << House.AbsoluteError() << endl;
   count.Stop();
   cout << count.Read() << " function evaluations were used." << endl;

   return 0;
 }
