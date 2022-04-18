#include <iostream.h>
#include <cubpack.h>

extern "C" {void zbesh_ (real&, real&, real&, int&, int&, int&, 
  real[], real[], int&, int&);}

real AbsHankel1 ( const Point& z)
  { real x=z.X(), y=z.Y(), cr[10], ci[10], fnu=0;
    int kode=1, m=1, n=1, nz, ierr;
    zbesh_(x,y,fnu,kode,m,n,cr,ci,nz,ierr);
    return sqrt(cr[0]*cr[0]+ci[0]*ci[0]);
  }

real AbsHankel2 ( const Point& z)
  { real x=(z.X()*0.8+z.Y()*0.6), y=(0.8*z.Y()-z.X()*0.6), cr[10], ci[10], fnu=0;
    int kode=1, m=1, n=1, nz, ierr;
    zbesh_(x,y,fnu,kode,m,n,cr,ci,nz,ierr);
    return sqrt(cr[0]*cr[0]+ci[0]*ci[0]);
  }

int main()
  { Point O(0,0), A(-1,0),B(-0.8,-0.6);
    EvaluationCounter count;
    cout << "   Omega = 1"<<endl;
    count.Start();
    cout << "Integral over circle: " 
         << Integrate(AbsHankel1, CIRCLE(O,A));
    count.Stop();  
    cout << "(" << count.Read() << " function evaluations used)" << endl;
    count.Reset(); count.Start();
    cout << "Integral over cut circle: "
         << Integrate(AbsHankel1, POLAR_RECTANGLE(O,A,A));
    count.Stop();
    cout << "(" << count.Read() << " function evaluations used)" << endl;
    cout << endl<< "   Omega = (4+3i)/5"<<endl;
    count.Start();
    cout << "Integral over circle: " 
         << Integrate(AbsHankel2, CIRCLE(O,B));
    count.Stop();  
    cout << "(" << count.Read() << " function evaluations used)" << endl;
    count.Reset(); count.Start();
    cout << "Integral over cut circle: "
         << Integrate(AbsHankel2, POLAR_RECTANGLE(O,B,B));
    count.Stop();
    cout << "(" << count.Read() << " function evaluations used)" << endl;
  }
