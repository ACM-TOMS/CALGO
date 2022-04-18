// Example 34 from DITAMO
// The function is rescaled such that the solution is 1.

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { real x=p.X() , y=p.Y(),z , a=1.0/3.0 , b=2.0/3.0 ;
   if (x*y < 1e-15) return 0.0;
   z = x + y + 1/(x*y);
   if ((z > 100) || (z < -100))
     { return 0.0;}
   else
     { return
       exp(-z) / (pow(x,b)*pow(y,a)*0.1806075059054343159);}
 }

int main ()
 {
// Using the first constructor:
      Point origin(0,0),p1(1,0),p2(0,1);
      PLANE_SECTOR quadrant(origin,p1,p2);
// Using the second constructor:
//    Point origin(0,0);
//    PLANE_SECTOR quadrant(origin,0.0, 0.0, M_PI/2);

   EvaluationCounter count;

   count.Start();
   cout <<"The integral is " << Integrate(f,quadrant,0,0.5e-4);
   cout <<" with absolute error " << quadrant.AbsoluteError() << endl;
   count.Stop();
   cout << count.Read() << " function evaluations were used." << endl;

   return 0;
 }
