// Example 20 from ditamo
// exact value = 4/15

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { real x=p.X() , y=p.Y() ;
   return
       sqrt(1-x-y);
 }

int main ()
 {
   Point origin(0,0),p1(1,0),p2(0,1);
   TRIANGLE ditamo20(origin,p1,p2);

   EvaluationCounter count;

   count.Start();
   cout <<"The integral is " << Integrate(f,ditamo20,0,0.5e-4,10000);
   cout <<" with absolute error " << ditamo20.AbsoluteError() << endl;
   count.Stop();
   cout << count.Read() << " function evaluations were used." << endl;

   return 0;
 }
