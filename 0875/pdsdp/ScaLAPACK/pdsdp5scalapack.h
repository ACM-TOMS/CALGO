#include "dsdp5.h"

#include "mpi.h"

int PDSDPUseSCALAPACKLinearSolver(DSDP);
void DSDPSetRank(int);

/*
To use PDSDP, 
1. Get a problem working using DSDP.
2. Add the routines:
   -- MPI_Initialize()
   -- MPI_Finalize()
   -- PDSDPUseSCALAPACKLinearSolver()
*/
