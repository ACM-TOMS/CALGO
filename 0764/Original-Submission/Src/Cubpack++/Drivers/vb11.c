// This example assumes that math.h contains the function erfc.
// The erfc function is only used to compute the exact value.

// Problem posted in sci.math.num-analysis on 9 Feb 1995
// by N. Sukumar (E-mail: n-sukumar@nwu.edu),
// Theoretical & Applied Mechanics, Northwestern U, Evanston IL 60208.
// Exact value = 1/2 - (1/2) * erf (h/sqrt(2))
//          where h = b / sqrt (1 + a^2)

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { real x=p.X() , y=p.Y();
       return
       exp(-(x*x+y*y)/2.0)/(2*M_PI);
 }

int main ()
 {
      real a=1,b=1,h=b/sqrt(1.0+a*a),t=h/sqrt(2.0);
      Point A(0,b),B(1,a+b),D=2*A-B;
      PLANE_SECTOR halfplane(A,B,D);

   EvaluationCounter count;

   count.Start();
   cout <<"The estimated integral is " << Integrate(f,halfplane,0,0.5e-6)<<endl;
   cout <<"The exact value is        "<< erfc(t)/2.0 <<endl;
   count.Stop();
   cout << count.Read() << " function evaluations were used." << endl;

   return 0;
 }
