// Example 28 from ditamo

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { real y=p.Y();
   return 1/(y*y);
 }

int main ()
 { Point A(2,-1), B(1,-1);
   SEMI_INFINITE_STRIP langelat(A,B);

   EvaluationCounter TikTak; TikTak.Start();
 //  cout.setf(ios::left,ios::adjustfield);
   cout.setf(ios::scientific,ios::floatfield);

   cout <<"The integral is " << Integrate(f,langelat,0,0.5e-7,20000);
   cout <<" with estimated absolute error "
        << langelat.AbsoluteError() << endl;
   cout <<"The real error is "
        << langelat.Integral() - 1 <<endl;

   TikTak.Stop(); cout<<"Number of evaluations = "<<TikTak.Read()<<endl;

   return 0;
 }
