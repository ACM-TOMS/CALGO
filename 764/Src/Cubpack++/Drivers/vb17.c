// Example 23 from ditamo
// The function is rescaled such that the solution is 1.

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 {
   real x=p.X() , y=p.Y();
   return 0.25/sqrt(x*(2.0-x));
 }

real Height(const Point& p)
  { return sqrt(4-2*p.X()); }

int main ()
 {
   real Integral_est, Error_est;
   Boolean Success;
   Point A(0,0), B(2,0);
   GENERALIZED_RECTANGLE parabola(Height,A,B);

   do
     {
      Integrate(f,parabola,Integral_est,Error_est,Success,0,0.5e-1,10000);
      cout <<"The integral is " << Integral_est <<endl;
      cout <<" with estimated absolute error "
           << Error_est << endl;
      cout << "-------------------------------------------------"<<endl;
     }
   while ( ! Success );

   return 0;
 }
