#include "pampac.h"

/* want to check that copy_PTnode make a shallow copy and to find out how copies/originals can be deallocated carfeully */

void main(int argc, char* argv[]) {
  PTnode *orig=NULL, *clone=NULL;
  int N = 6;
  orig = new_PTnode();
  clone = new_PTnode();
  orig->N_dim = N;
  orig->z = malloc(N*sizeof(double));
  orig->z[0]=1.7;
  orig->z_init = malloc(N*sizeof(double));
  orig->z_init[1]=-3.2;
  orig->T_init = malloc(N*sizeof(double));
  orig->T_init[2]=4.5;
  printf("***\nPrinting ORIGINAL node\n***\n");
  print_PTnode(orig);
  copy_PTnode(orig,clone);
  printf("NOW: assigning original node's allocated fields to NULL");
  orig->z = NULL;
  orig->z_init = NULL;
  orig->T_init = NULL;
  printf("***\nPrinting CLONED node\n***\n");
  print_PTnode(clone);
  printf("***\nPrinting ORIGINAL node AGAIN\n***\n");
  print_PTnode(orig);
  //  printf("Deallocating memory for original node explecitly.\n");
  // free(orig);
  printf("Deleting original using delete_tree\n");
  delete_tree(orig);
  printf("***\nPrinting CLONE node AGAIN\n***\n");
  print_PTnode(clone);
  printf("Deleting clone using delete_tree\n");
  delete_tree(clone);
}
