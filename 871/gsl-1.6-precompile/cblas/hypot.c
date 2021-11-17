#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "ArithmosIO.hh"

#include <math.h>

static MpIeee xhypot(const MpIeee x, const MpIeee y);

static MpIeee xhypot(const MpIeee x, const MpIeee y)
{
  MpIeee xabs=  fabs(x) ;
  MpIeee yabs=  fabs(y) ;
  MpIeee min;MpIeee  max;

  if (xabs < yabs) {
    min = xabs ;
    max = yabs ;
  } else {
    min = yabs ;
    max = xabs ;
  }

  if (min == 0) 
    {
      return max ;
    }

  {
    MpIeee u=  min / max ;
    return max * sqrt (1 + u * u) ;
  }
}
