#include "mpi.h"
#include "pampac.h"
/**********************************************************************/
/* Breadth-first traversal of tree to assign processor IDs ("pid"s)   */
/* to nodes that have not converged or diverged (i.e., that still     */
/* need processors). Assign pid=MPI_PROC_NU for CONVERGED & FAILED nodes;  */
/* other nodes are assigned pids such that 0 < pid <= Np-1. As each   */
/* node is traversed, its pid, its label, and its depth is recorded.  */
/**********************************************************************/
void
assign_processes (PTnode *root, options_struct *opts, int N_proc) {
  Queue q;
  PTnode *node;
  int k, label = 0, pid = 1;

  /* Note: pid starts at 1. Process 0 is the master process. */
  root->depth = 0;
  root->state = CONVERGED;
  if (N_proc==1) /* Terminates if only one processor */
    return;

  /* Breadth-first traversal of tree */
  init_queue (&q);
  enqueue (&q, root);
  while (!empty_queue (&q)) {
    node = front_of_queue (&q);
    dequeue (&q); /* Mark node as visited by dequeuing */
    /* Assign processor ID=MPI_PROC_NULL to nodes:
       * that have previously converged (i.e., CONVERGED)
       * that have diverged (i.e., FAILED)
       * when there are no more processors available
    */
    node->label = label++;
    node->pid = MPI_PROC_NULL;
    if (node->state==PROGRESSING || node->state==CONVERGING)
      if (pid <= N_proc-1)
        node->pid = pid++;
    /* Append children of node just dequeued (if any) to end of queue */
    if (node->child != NULL) {
      for (k=0; k<node->max_children; k++) {
        if (node->child[k] != NULL) {
          enqueue (&q, node->child[k]);
        }
      }
    }
  }
}
