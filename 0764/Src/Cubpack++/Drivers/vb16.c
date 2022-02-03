
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
   TRIANGLE T1( Point(0,0), Point(1,0), Point(0,1)),
            T2( Point(1,0), Point(2,0), Point(1.5,0.5)),
            T3( Point(2,0), Point(2,1), Point(1.5,0.5));
   RECTANGLE R1( Point(1,0), Point(0,1), Point(2,1));
   REGION_COLLECTION House1 = T1 + T2 + T3 + R1;

   RECTANGLE R3( Point(1,0), Point(0,1), Point(2,1)),
             R2( Point(1.8,0), Point(1.4,0.4), Point(2,0.2));
   TRIANGLE T4( Point(1,0), Point(1.8,0), Point(1.4,0.4)),
            T7( Point(0,0), Point(1,0), Point(0,1)),
            T5( Point(1.8,0), Point(2,0), Point(2,0.2)),
            T6( Point(2,0.2), Point(2,1), Point(1.6,0.6));
   REGION_COLLECTION House2 = T7 + R3 + R2 + T4 + T5 + T6;
   EvaluationCounter count;

   count.Start();
   cout <<"The integral is " << Integrate(f,House1,0,0.5e-1,1000000);
   cout <<" with absolute error " << House1.AbsoluteError() << endl;
   count.Stop();
   cout << count.Read() << " function evaluations were used." << endl;
   cout << "Contribution of the different parts:" << endl;
   cout << "R1: integral and error are " << R1.Integral() << " and "
                                        << R1.AbsoluteError() <<endl;
   cout << "T1: integral and error are " << T1.Integral() << " and "
                                        << T1.AbsoluteError() <<endl;
   cout << "T2: integral and error are " << T2.Integral() << " and "
                                        << T2.AbsoluteError() <<endl;
   cout << "T3: integral and error are " << T3.Integral() << " and "
                                        << T3.AbsoluteError() <<endl;
   cout << endl;
   count.Start();
   cout <<"The integral is " << Integrate(f,House2,0,0.5e-1,1000000);
   cout <<" with absolute error " << House2.AbsoluteError() << endl;
   count.Stop();
   cout << count.Read() << " function evaluations were used." << endl;

   return 0;
 }
