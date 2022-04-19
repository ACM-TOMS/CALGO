#include "mpi.h"
#include "string.h"
#include "KS_example.h"
/**********************************************************************/
/* Generic MPI main routine with input options communicated through a */
/* custom data structure opts.                                        */
/**********************************************************************/
int
main (int argc, char *argv[]) {
  int p_id;         /* rank of process */
  int N_p;          /* number of processes */
  int N_dim;        /* number of variables */
  bool success;     /* flag to check for successful setup */
  
  MPI_Init (NULL, NULL);                 /* start up MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &p_id); /* find process ID */
  MPI_Comm_size (MPI_COMM_WORLD, &N_p);  /* find # of processes */

  if (p_id == 0) {
	options_struct opts;        /* Problem/algorithm parameters */
    initialize_options (&opts); /* Set defaults before reading input */
    success = parse_options (argc, argv, &opts);
    if (success) {
      N_dim = opts.N_dim;
      /* Broadcast N_dim from master to other processes */
      MPI_Bcast (&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
      success = setup_globals (N_dim); /* Initialisation */
      if (success) {
        master_process (N_p, &opts);
        teardown_globals ();   /* Clean-up */
        printf ("PID %d: Cleaning up master process.\n", p_id);
        fflush (stdout);
      } else {
		fprintf (stderr, "%s: Error setting up globals.\n",
		         __func__);
	  }
    } else {
		/* Error parsing parameter file, so broadcast N_dim==-1 to
		 * other processes so all processes can exit gracefully. */
		N_dim = -1; 
        MPI_Bcast (&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
        fprintf (stderr, "%s: Error parsing parameters file.\n",
		         __func__);
    }
    MPI_Barrier (MPI_COMM_WORLD); /* Synchronise before shutting down */
  } else {
    /* N_dim is read from disk by master and broadcast to slaves */
    MPI_Bcast (&N_dim,1,MPI_INT,0,MPI_COMM_WORLD);
    success = (N_dim>0) && setup_globals (N_dim); /* Initialisation */
    if (success) {
      slave_process (N_dim);
      teardown_globals ();   /* Clean-up */
      printf ("PID %d: Cleaning up slave process.\n", p_id);
      fflush (stdout);
    }
    MPI_Barrier (MPI_COMM_WORLD); /* Synchronise before shutting down */
  }
  MPI_Finalize ();      /* Shut down MPI */
  return 0;
}
