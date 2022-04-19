// Example 29 from ditamo

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { 
   return fabs(p.Y())*exp(-p.Y()*p.Y()/2);
 }

int main ()
 { Point A(0,0), B(1.107148717794090503,0);
   SEMI_INFINITE_STRIP langelat(A,B);

   EvaluationCounter TikTak; TikTak.Start();
 //  cout.setf(ios::left,ios::adjustfield);
   cout.setf(ios::scientific,ios::floatfield);

   cout <<"The integral is " << Integrate(f,langelat,0,0.5e-4,10000);
   cout <<" with estimated absolute error "
        << langelat.AbsoluteError() << endl;

   TikTak.Stop(); cout<<"Number of evaluations = "<<TikTak.Read()<<endl;

   return 0;
 }
