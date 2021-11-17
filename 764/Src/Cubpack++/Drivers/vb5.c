// Example 5 from ditamo

#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
 { real x=p.X() , y=p.Y();
   return p.Y();
 }

real Height(const Point& p)
  { real x=p.X();
    return sqrt(4*cos(x)*cos(x)/9 - 16*sin(x)*sin(x)/25);
  }

int main ()
 { Point A(0,0), B(atan(5.0/6.0),0);
   GENERALIZED_RECTANGLE lemniscate(Height,A,B);

   EvaluationCounter TikTak; TikTak.Start();
 //  cout.setf(ios::left,ios::adjustfield);
   cout.setf(ios::scientific,ios::floatfield);

   cout <<"The integral is " << Integrate(f,lemniscate,0,0.5e-7,10000);
   cout <<" with estimated absolute error "
        << lemniscate.AbsoluteError() << endl;
   cout <<"The real error is "
        << lemniscate.Integral() - (2-11*atan(5.0/6.0)/15)/15 <<endl;

   TikTak.Stop(); cout<<"Number of evaluations = "<<TikTak.Read()<<endl;

   return 0;
 }
