#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

static inline void
reset_nrmax (gsl_integration_workspace * workspace);

static inline void
reset_nrmax (gsl_integration_workspace * workspace)
{
  workspace->nrmax = 0;
  workspace->i = workspace->order[0] ;
}
