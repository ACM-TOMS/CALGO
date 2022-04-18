#include <mpi.h>
#include "pampac.h"
/**********************************************************************/
/* Traverse through tree sending vectors to distinct processors for   */
/* computing corrector steps concurrently.                            */
/**********************************************************************/
void
compute_corrector_steps (PTnode * alpha, int N_p) {
  int k;
  bool needs_correction_step, has_data;

  needs_correction_step = (alpha->state==PROGRESSING || alpha->state==CONVERGING);
  has_data = (alpha->z!=NULL) && (alpha->T_init!=NULL);
  if (has_data && needs_correction_step) {
    /* Corresponding MPI_Recv calls in slave_process */
    MPI_Send (alpha->z, alpha->N_dim, MPI_DOUBLE, alpha->pid,
              CONTINUE_TAG, MPI_COMM_WORLD);
    MPI_Send (alpha->T_init, alpha->N_dim, MPI_DOUBLE, alpha->pid,
              CONTINUE_TAG, MPI_COMM_WORLD);
  }
  for (k = 0; k < alpha->max_children; k++)
    if (alpha->child[k] != NULL) /* Recursive call */
      compute_corrector_steps (alpha->child[k], N_p);
}
