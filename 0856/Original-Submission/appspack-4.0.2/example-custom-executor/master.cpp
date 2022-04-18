// $Id: master.cpp,v 1.1.2.1 2005/06/28 23:42:49 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/example-custom-executor/Attic/master.cpp,v $ 

/* 
   This file provides the routine for the master. The master is
   responsible for the following: 

   o Parsing the command line to get the number of variables
   o Sending the number of variables to each worker (see note 1)
   o Creating the parameter list (see note 2)
   o Instantiating the executor
   o Instantiating the bounds
   o Creating the solver
   o Running the solver (see note 3)
   o Sending a termination message to each worker

   Note 1: It is possible for the workers to also read the number of
   variables from the command line; however, we do it this way to show
   an example of sending initialization information to the workers.

   Note 2: The parameter list can be read from a file or hard-coded or
   some combination. APPSPACK provides tools for reading the parameter
   list form a file; see the function APPSPACK::parseTextInputFile.

   Note 3: Once the solver has been run, the best point can be
   extracted. See APPSPACK::Solver for more details.
*/

#include "mpi.h"                            // <-- Provides MPI
#include "msgtags.hpp"		            // <-- Provides the message tags
#include "executor.hpp"             // <-- Provides user-defined custom executor
#include "APPSPACK_Solver.hpp"              // <-- Provides APPSPACK::Solver
#include "APPSPACK_Vector.hpp"              // <-- Provides APPSPACK::Vector 
#include "APPSPACK_Constraints_Bounds.hpp"  // <-- Provides APPSPACK::Constraints::Bounds

void master(int argc, char* argv[], int nprocs)
{
  // *** Parse the command line ***

  int n;			// Problem size
  sscanf(argv[1],"%d",&n);
  
  // *** Intialize the workers by sending the problem size to each ***
  // Note: Must have a corresponding recv on the worker

  for (int i = 1; i < nprocs; i ++)
    MPI_Send(&n, 1, MPI_INT, i, SIZE, MPI_COMM_WORLD);
  
  // *** Create and fill the parameter list ***
  
  APPSPACK::Parameter::List params;

  APPSPACK::Vector initialx(n, 1.0/n);
  params.sublist("Solver").setParameter("Initial X", initialx);

  // Be sure to do casts on any ambiguous constants, such as in this case

  params.sublist("Solver").setParameter("Debug", static_cast<int>(3));
  params.sublist("Solver").setParameter("Step Tolerance", static_cast<double>(0.0001));

  APPSPACK::Vector lowerbounds(n, -1.0);
  params.sublist("Bounds").setParameter("Lower", lowerbounds);

  APPSPACK::Vector upperbounds(n, 3.0);
  params.sublist("Bounds").setParameter("Upper", upperbounds);
  
  // *** Instantiate the custom executor ***

  CustomExecutor executor(n,nprocs);
  
  // *** Instantiate the bounds ***

  APPSPACK::Constraints::Bounds bounds(params.sublist("Bounds"));
  
  // *** Create the solver ***

  APPSPACK::Solver solver(params.sublist("Solver"), executor, bounds);
  
  // *** Run the solver ***

  APPSPACK::Solver::State state = solver.solve();
  
  // *** Send a termination command to every worker ***

  // Note: The worker needs to keep checking for the quit message

  for (int i = 1; i < nprocs; i ++)
    MPI_Send(&n, 1, MPI_INT, i, QUIT, MPI_COMM_WORLD);

}
