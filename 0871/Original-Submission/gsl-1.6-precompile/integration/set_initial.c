#include <iostream>
#include <iomanip>
using namespace std;

#include "MpIeee.hh"
#include "CMpIeee.hh"
#include "ArithmosIO.hh"

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
                         MpIeee result, MpIeee error);

static inline
void set_initial_result (gsl_integration_workspace * workspace, 
                         MpIeee result, MpIeee error)
{
  workspace->size = 1;
  workspace->rlist[0] = result;
  workspace->elist[0] = error;
}
