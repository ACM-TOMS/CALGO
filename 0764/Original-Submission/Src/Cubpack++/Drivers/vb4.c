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
   Point A(0,0), B(2,0);
   GENERALIZED_RECTANGLE parabola(Height,A,B);
   Chrono TikTak;

   TikTak.Start();
   cout <<"The integral is " << Integrate(f,parabola,0,0.5e-4);
   TikTak.Stop();
   cout <<" with absolute error " << parabola.AbsoluteError() << endl;
   cout <<"Elapsed Time: " << TikTak.Read()/1000<<" seconds" << endl;

   return 0;
 }
