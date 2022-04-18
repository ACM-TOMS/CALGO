#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

#define DBL_MEMCPY(dest,src,n) memcpy((dest),(src),(n)*sizeof(double))
#define DBL_ZERO_MEMSET(dest,n) memset((dest),0,(n)*sizeof(double))
