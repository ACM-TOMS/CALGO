// Example

#include <cubpack.h>
#include <iostream.h>
#include <math.h>

real f(const Point& p)
 { 
   return exp(-p.Y()*p.Y());
 }

int main ()
 { Point A(0,0), B(1,0);
   INFINITE_STRIP langelat(A,B);

   EvaluationCounter TikTak; TikTak.Start();
 //  cout.setf(ios::left,ios::adjustfield);
   cout.setf(ios::scientific,ios::floatfield);

   cout <<"The integral is " << Integrate(f,langelat,0,0.5e-3,10000);
   cout <<" with estimated absolute error "
        << langelat.AbsoluteError() << endl;
   cout <<"The real error is "
        << langelat.Integral() - sqrt(M_PI) <<endl;

   TikTak.Stop(); cout<<"Number of evaluations = "<<TikTak.Read()<<endl;

   return 0;
 }
