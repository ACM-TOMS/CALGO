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
   real Integral_est, Error_est;
   Boolean Success;

   EvaluationCounter TikTak; TikTak.Start();
   cout.setf(ios::scientific,ios::floatfield);

   do
     {
      Integrate(f,langelat,Integral_est,Error_est,Success,0,0.5e-7,10000);
      cout <<"The integral is " << Integral_est;
      cout <<" with estimated absolute error "
           << Error_est << endl;
      cout <<"The real error is "
           << Integral_est - 1 <<endl;
      cout << "-------------------------------------------------"<<endl;
     }
   while ( ! Success );

   TikTak.Stop(); cout<<"Total number of evaluations = "
                      <<TikTak.Read()<<endl;

   return 0;
 }
