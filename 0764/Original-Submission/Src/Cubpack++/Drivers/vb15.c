// The `House' example with on each subregion a different integrand.

#include <cubpack.h>
#include <iostream.h>

real f1(const Point& p)
 { 
     return ( p.X()*p.X() );
 }

real f2(const Point& p)
 {
     return ( p.Y()*p.Y() );
 }

int main ()
 {
   TRIANGLE Roof( Point(0,1), Point(2,1), Point(1,2));
   RECTANGLE Walls( Point(0,0), Point(0,1), Point(2,0));
   REGION_COLLECTION House = Walls + Roof;
   Roof.LocalIntegrand(f1);
   Walls.LocalIntegrand(f2);
   EvaluationCounter count;

   count.Start();
   cout <<"The integral is " << Integrate(House);
   cout <<" with estimated absolute error " << House.AbsoluteError()
        << endl;
   cout <<"The contributions of the subregions are:"<<endl;
   cout <<"  Walls: integral = "<<Walls.Integral()
        <<", error = " << Walls.AbsoluteError() << endl;
   cout <<"  Roof:  integral = "<<Roof.Integral()
        <<" , error = " << Roof.AbsoluteError() << endl;
   count.Stop();
   cout <<"The exact value is 2/3 + 7/6 = 11/6" << endl;
   cout << count.Read() << " function evaluations were used." << endl;

   return 0;
 }
