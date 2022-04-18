// This example comes from Problem 94-7 by D.E. Loper
// Siam Review, Vol. 36 (1994) page 277

#include <iostream.h>
#include <cubpack.h>
#include <math.h>

static real q;

real A(const Point& p)
 { real m=p.X() , f=p.Y();
   real t=(1-m*m)*(1-m*m), s = sin(f);
   s = s*s;
   return
          (t*s)/(t*s*s+q*q*m*m);
 }

real B(const Point& p)
 { real m=p.X() , f=p.Y();
   real t=(1-m*m), s = sin(f);
   s = s*s;
   return
          (t*s*(t*s+m*m))/(t*t*s*s+q*q*m*m);
 }

int main ()
 { Point Origin(0,0), mu(1,0),phi(0,M_PI/2);

   cout<<"q	A(q)		B(q)		A(q)/B(q)"<<endl;
   cout.precision(7);
   cout.setf(ios::fixed,ios::adjustfield);
   for ( q=1; q<20 ; q++)
   {
     RECTANGLE Loper1(Origin,mu,phi);
     RECTANGLE Loper2(Origin,mu,phi);
     real int1,int2;
     int1 = Integrate(A,Loper1, 0.0, 0.5e-7, 100000);
     int2 = Integrate(B,Loper2, 0.0, 0.5e-7, 100000);
     cout << q<<"	";
     cout << int1 << "	";
     cout << int2 << "	";
     cout << int1/int2 << endl;
   }

   return 0;
 }
