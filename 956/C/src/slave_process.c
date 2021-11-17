#include <gsl_cblas.h>
#include <mpi.h>
#include "pampac.h"
/**********************************************************************/
/* Slave processors wait to receive the current corrector iterate and */
/* tangent vector to compute the next corrector iterate. Two calls to */
/* MPI_Recv are needed for the current iterate and tangent vector.    */
/**********************************************************************/
void
slave_process (int N_dim) {
  int master = 0, istatus;
  double z[2*N_dim], *T, res_norm;
  MPI_Status status;

  int pid;
  MPI_Comm_rank(MPI_COMM_WORLD,&pid);

  T = z+N_dim;
  while (true) {
    /* Corresponding MPI_Send calls in compute_corrector_steps */
    MPI_Recv (z, N_dim, MPI_DOUBLE, master, MPI_ANY_TAG, MPI_COMM_WORLD,
              &status);
    if (status.MPI_TAG == EXIT_TAG) {
      break;
    }
    MPI_Recv (T, N_dim, MPI_DOUBLE, master, MPI_ANY_TAG, MPI_COMM_WORLD,
              &status);

    istatus = single_corrector_step (N_dim, z, T);
    if (istatus!=0) {
		fprintf (stderr, "%s: Failed corrector step on PID %i\n",
		         __func__, pid);
		fprintf (stderr, "%s: error status = %i\n", __func__, istatus);
	    break;
    }
    /* residual vector stored in T after call to compute residual */
    istatus = compute_residual (N_dim, z, T);
    if (istatus!=0) {
		fprintf (stderr, "%s: Error computing residual, status = %i\n",
		         __func__, istatus);
		fprintf (stderr, "%s: error status = %i\n", __func__, istatus);
		break;
    }
    res_norm = cblas_dnrm2 (N_dim-1, T, 1);
    /* Corresponding MPI_Recv calls in assess_residuals */
    MPI_Send (z, N_dim, MPI_DOUBLE, master, CONTINUE_TAG, MPI_COMM_WORLD);
    MPI_Send (&res_norm, 1, MPI_DOUBLE, master, CONTINUE_TAG, MPI_COMM_WORLD);
  }
  return;
}
