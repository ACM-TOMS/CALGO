#include "mpi.h"
#include "pampac.h"
/**********************************************************************/
/* Blanks fields of PTnodes with default values (to make it easier in */
/* in principle to detect errors in which fields have not been        */
/* assigned). As such, pointers are initialised with NULL, integer-   */
/* valued fields are initialised with -1 and real-valued fields are   */
/* initialised with NAN.                                              */
/* NOTE: initialize_PTnode allocates memory for the PTnode itself and */
/* for the field child. The remaining fields to allocate (namely z,   */
/* z_init, and T_init) must be explicitly allocated and deallocated   */
/* elsewhere.                                                         */
/**********************************************************************/
PTnode *
initialize_PTnode (int max_children) {
  PTnode *alpha;
  int k;
  alpha = malloc (sizeof (*alpha));
  alpha->N_dim = -1;
  /* label is initialised with negative value */
  alpha->label = -1;
  alpha->pid = MPI_PROC_NULL;
  alpha->state = PROGRESSING;
  alpha->nu = -1;
  alpha->nu_init = -1;
  alpha->nu_valid = -1;
  alpha->nu_viable = -1;
  alpha->h = NAN;
  alpha->res_norm = NAN;	/* Assume large initially */
  alpha->valid_path_length = NAN;
  alpha->viable_path_length = NAN;
  alpha->valid_index = -1;
  alpha->viable_index = -1;
  alpha->depth = -1;
  /* Initialise all pointers to NULL */
  alpha->z = NULL;
  alpha->T_init = NULL;
  alpha->z_init = NULL;
  alpha->max_children = max_children;
  alpha->child = NULL;
  alpha->child = malloc (max_children * sizeof(PTnode*));
  for (k=0; k<max_children; k++)
    alpha->child[k] = NULL;
  return (alpha);
}
