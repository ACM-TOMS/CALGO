#include <mpi.h>
#include "pampac.h"
/**********************************************************************/
/* Cleanup routine to terminate MPI slave processes.                  */
/**********************************************************************/
void
stop_slaves (int N_p) {
  double temp = 0.0;
  int pid;
  for (pid = 1; pid < N_p; pid++)
    MPI_Send (&temp, 1, MPI_DOUBLE, pid, EXIT_TAG, MPI_COMM_WORLD);
  return;
}
