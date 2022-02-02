#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { real r=p.Length();
   return 1.0/(1.0+r*r*r*r); }

int main ()
 { Point origin(0,0);
   real radius=1;
   CIRCLE cir(origin,radius);

   cout <<"The integral is " << Integrate(f, cir, 0, 1.0e-6, 10000);
   cout <<" with absolute error " << cir.AbsoluteError() << endl;

   return 0;
 }
