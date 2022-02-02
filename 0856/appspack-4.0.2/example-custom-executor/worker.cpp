// $Id: worker.cpp,v 1.1.2.1 2005/06/28 23:42:49 tgkolda Exp $ 
// $Source: /space/CVS-Acro/acro/packages/appspack/appspack/example-custom-executor/Attic/worker.cpp,v $ 

/* 
   This file provides the routine for the worker. The worker is
   responsible for the following: 

   o Receive set-up information
   o Evaluate points as they arrive
   o Return upon termination message

*/

#include "mpi.h"		// <-- Provides MPI
#include "msgtags.hpp"		// <-- For MPI tags

// Return true if the constraint is violated
bool constraint(int n, double* x)
{
  double tmp = 0;

  for (int i = 0; i < n; i ++)
    tmp += x[i] * x[i];

  return (tmp <= 1.0);
} 

// Function evaluation
double feval(int n, double* x)
{
  double f = 0;
  
  for (int i = 0; i < n; i ++)
    f += (i + 1) * x[i] * x[i];

  return(f);
} 

void worker()
{
  // **** Receive initial message with problem size ***
  int n;			// Problem size
  MPI_Status status;

  // NOTE: This matches a send by the master.
  MPI_Recv(&n, 1, MPI_INT, 0, SIZE, MPI_COMM_WORLD, &status);
  
  // *** Continuously receive and process incoming messages ***

  // Create other variables
  double x[n];			// Vector to be evaluated
  int msgtag;			// Message tag
  int tag;			// Vector tag
  int code;			// Solution code (0 = constraint violation, 1 = success)
  double f;			// Functin value
  
  while (1)
  {
    // Blocking probe for the next message
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    
    // Extract the message tag
    msgtag = status.MPI_TAG;
    
    // Check for termination tag and, if so, quit
    if (msgtag == QUIT)
      break;
    
    // Receive message in two parts: tag and x 
    // NOTE: This must match the spawn function in the evaluator
    MPI_Recv(&tag, 1, MPI_INT, 0, XTAG,  MPI_COMM_WORLD, &status);
    MPI_Recv(x, n, MPI_DOUBLE, 0, XVEC, MPI_COMM_WORLD, &status);
    
    // Evaluate the function at x if the constraint is not violated
    if (constraint(n,x))
    {
      code = 0;
    }
    else
    {
      code = 1;
      f = feval(n,x);
    }
    
    // Return the answer in three parts: tag, code, f. These are
    // interpreted by the recv function in the evaluator.  
    // NOTE: This must match the recv function in the evaluator
    MPI_Send(&tag,  1, MPI_INT, 0, XTAG, MPI_COMM_WORLD);
    MPI_Send(&code, 1, MPI_INT, 0, CODE, MPI_COMM_WORLD);
    MPI_Send(&f, 1, MPI_DOUBLE, 0, FVAL, MPI_COMM_WORLD);
  }

}

