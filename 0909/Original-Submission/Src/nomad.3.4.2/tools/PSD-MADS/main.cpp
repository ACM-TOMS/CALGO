/*-------------------------------------------------------------*/
/*                            PSD-MADS                         */
/*-------------------------------------------------------------*/
/*                                                             */
/*  usage:                                                     */
/*                                                             */
/*    "mpirun -np p psdmads param_file bbe ns"                 */
/*    with p > 2, bbe > 0, and 1 <= ns <= number of variables  */
/*    . ns is the number of free variables for each process    */
/*    . bbe is the max number of evaluations for each process  */
/*                                                             */
/*-------------------------------------------------------------*/
/*                                                             */
/*  processes:                                                 */
/*                                                             */
/*            0: master                                        */
/*            1: pollster slave (1 direction)                  */
/*     2 to p-2: regular slaves (2ns directions)               */
/*          p-1: cache server                                  */
/*                                                             */
/*-------------------------------------------------------------*/
/*  See the user guide for other details and the description   */
/*  of the algorithm                                           */
/*-------------------------------------------------------------*/

/*-----------------------------------------------------------*/
#include "Master_Slaves.hpp"
using namespace std;
using namespace NOMAD;

const bool DEBUG = false;

/*-----------------------------------*/
/*           main function           */
/*-----------------------------------*/
int main ( int argc , char ** argv ) {

  // MPI initialization:
  MPI_Init ( &argc , &argv );
  int rank , np;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &np   );

  // check the arguments and the number of processes:
  if ( np <= 2 || argc != 4 ) {
    if ( rank == 0 )
      cerr << "usage: mpirun -np p " << argv[0]
	   << " param_file bbe ns, with p>2,"
	   << " bbe>1, and 1<=ns<=n."
	   << endl;
    MPI_Finalize();
    return 1;
  }

  // display:
  Display out ( cout );
  out.precision ( 16 );

  // parameters:
  NOMAD::Parameters p ( out );
  int bbe = atoi ( argv[2] );
  int ns  = atoi ( argv[3] );

  try {
    
    // read the parameters file:
    p.read ( argv[1] );

    // check the parameters:
    p.check();

    if ( ns < 1 || ns > p.get_dimension() )
      throw Exception ( __FILE__ , __LINE__ ,
      "Bad value for ns the number of free variables for each process" );

    if ( p.get_nb_obj() > 1 )
      throw Exception ( __FILE__ , __LINE__ ,
      "PSD-MADS is not designed for multi-objective optimization" );
  }
  catch ( exception & e ) {
    if ( rank == 0 )
      cerr << "error with parameters" << endl;
    MPI_Finalize();
    return 1;
  }

  // start the master:
  Master_Slaves master_slaves ( rank , np , bbe , ns , p , DEBUG );
  master_slaves.start();

  // cache server:
  Cache_Server cache ( out                  ,
		       rank                 ,
		       np                   ,
		       p.get_h_min()        ,
		       p.get_max_bb_eval()  ,
		       false                ,  // ALLOW_MULTIPLE_EVALS
		       DEBUG                  );

  // start the cache server:
  if ( rank == np-1 ) {
    if ( !DEBUG )
      out << endl << "TIME\tBBE\tOBJ" << endl << endl;
    cache.start();
  }

  // slaves: algorithm creation and execution:
  if ( rank != 0 && rank != np-1 ) {

    // MADS run:
    master_slaves.mads_run ( cache );

    // stop the master:
    master_slaves.stop();
  }

  // stop the cache server:
  cache.stop();

  // display the final solution:
  if ( !DEBUG && rank == np-1 )
    cache.display_best_points ( out );

  // MPI finalization:
  MPI_Finalize();
  return 0;
}
