// Example 33 from ditamo
// Exact value = arctan(2)
// Note that the "exact" value in the original paper is wrong.

#include <cubpack.h>
#include <iostream.h>
#include <iomanip.h>

real f(const Point& p)
 {
   real x=p.X() , y=p.Y();
   return exp(0.5*(-x*x -y*y));
 }

int main ()
 {
   Point origin(0,0);
   real innerradius=0, alfa=0, beta=atan(2.0);
   PLANE_SECTOR wedge(origin,innerradius,alfa,beta);
   EvaluationCounter count;

   cout.setf(ios::scientific,ios::floatfield);
   cout<<"req. rel. error    est integral    est error   abs error  evaluations" <<endl;

   count.Start();
   for (real req_err=0.05; req_err>1e-12; req_err/=10)
   {
     cout << setprecision(1) << "   <" << req_err <<"      " 
          << setprecision(10) << Integrate(f,wedge,0,req_err) << "   ";
     cout << setprecision(1) << wedge.AbsoluteError() << "     "
          << fabs(wedge.Integral() - atan(2.0)) << "     "
          << count.Read() << endl;
   }
   count.Stop();

   return 0;
 }
