// $Id: main.cpp,v 1.1.2.1 2005/06/28 23:42:49 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/example-custom-executor/Attic/main.cpp,v $ 

/*
  This is the main file. It starts MPI, appropriately calls the master
  and worker tasks, and terminates MPI when everything is finished.
*/

#include <iostream>		// <-- for cout and cerr
#include "mpi.h"                // <-- Provides MPI

// forward declarations for functions defined in master.cpp and
// worker.cpp, respectively.
void master(int argc, char* argv[], int nprocs);
void worker();

int main(int argc, char* argv[])
{
  // Start MPI
  MPI_Init(&argc, &argv);

  // Get my rank 
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Get number of processors
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // Check that there is at least one worker 
  if (nprocs < 2)
  {
    std::cerr << "Error: This program requires at least 2 processes." << std::endl;
    MPI_Finalize();
    return 1;
  }
  
  // Check the number of input arguments 
  if (argc < 2)
  {
    if (rank == 0)
      std::cout << "Usage: " << argv[0] << " <problem size>" << std::endl;
    MPI_Finalize();
    return 1;
  }

  // Run main program
  if (rank == 0)	
  {
    master(argc, argv, nprocs);
  }
  else			
  {
    worker();
  } 

  // Exit MPI
  MPI_Finalize();

} // main

