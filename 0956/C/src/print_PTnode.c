#include "pampac.h"
/**********************************************************************/
/* Given a PTnode, print its contents to stdout.                      */
/**********************************************************************/
void
print_PTnode (PTnode *alpha) {
  printf ("\nN_dim = %d\n", alpha->N_dim);
  printf ("depth = %d\n", alpha->depth);
  printf ("max_children = %d\n", alpha->max_children);
  printf ("label = %d\n", alpha->label);
  printf ("pid = %d\n", alpha->pid);
  printf ("state = ");
  print_state (alpha,stdout);
  printf("\n");
  printf ("nu = %d\n", alpha->nu);
  printf ("nu_init = %d\n", alpha->nu_init);
  printf ("nu_valid = %d\n", alpha->nu_valid);
  printf ("nu_viable = %d\n", alpha->nu_viable);
  printf ("h_init = %12.5g\n", alpha->h_init);
  printf ("h = %12.5g\n", alpha->h);
  printf ("res_norm = %12.5g\n", alpha->res_norm);
  printf ("valid_path_length = %12.5g\n", alpha->valid_path_length);
  printf ("valid_index = %d\n", alpha->valid_index);
  printf ("viable_path_length = %12.5g\n", alpha->viable_path_length);
  printf ("viable_index = %d\n", alpha->viable_index);

  /* For sufficiently small problems, print vectors in full */
  int k, nt=4;
  if (alpha->z==NULL)
    printf ("alpha->z=NULL\n");
  else if (alpha->N_dim>=nt) {
    printf("z = (");
    for (k=0; k<nt; k++)
      printf(" %10.3e,", alpha->z[k]);
    printf("... )\n");
  } else 
    printf ("alpha->z=%p\n",alpha->z);

  if (alpha->z_init==NULL)
    printf ("alpha->z_init=NULL\n");
  else if (alpha->N_dim>=nt) {
    printf("z_init = (");
    for (k=0; k<nt; k++)
      printf(" %10.3e,", alpha->z_init[k]);
    printf("... )\n");
  } else
    printf("alpha->z_init=%p\n",alpha->z_init);

  if (alpha->T_init==NULL)
    printf ("alpha->T_init=NULL\n");
  else if (alpha->N_dim>=nt) {
    printf("T_init = (");
    for (k=0; k<nt; k++)
      printf(" %10.3e,", alpha->T_init[k]);
    printf("... )\n");
  } else
    printf("alpha->T_init=%p\n",alpha->T_init);

  if (alpha->child==NULL)
    printf ("alpha->child=NULL\n");
  else {
    printf("alpha->child=%p\n",alpha->child);
    for (k=0; k<alpha->max_children; k++)
      printf("alpha->child[%i]=%p\n", k, alpha->child[k]);
  }
}
