#include <cubpack.h>
#include <iostream.h>

real f(const Point& p)
  { real x=p.X();
    return x*x;
  }

int main()
  { Point p1(0,0), p2(1,1), p3(2,0);
    TRIANGLE T(p1,p2,p3);
    cout << "The integral is " << Integrate(f,T) << endl;

    return 0;
  }
