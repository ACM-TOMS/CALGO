#include "FLAME.h"
#include "MyFLA_Parallel_normal_random.h"


// ============================================================================
// Declaration of local variables.

static unsigned int  global_seed        = 3;
static int           global_num_threads = 4;
static double        rev2to15           = 1.0 / 32768.0;

// ============================================================================
// Declaration of local prototypes.

static double MyFLA_Normal_random_number( double mu, double sigma,
                  int * alternate, double * b1, double * b2, 
                  unsigned int * seedp );

static double myrand( unsigned int * seedp );


// ============================================================================
void MyFLA_Init_parallel_normal_random_generator( int num_threads, int seed ) {
  MyFLA_Reset_parallel_normal_random_generator( num_threads, seed );
}

// ============================================================================
void MyFLA_Reset_parallel_normal_random_generator( int num_threads, int seed ) {
  global_seed        = ( unsigned int ) ( seed > 0 ? seed : - seed );
  global_num_threads = num_threads;
}

// ============================================================================
void MyFLA_Close_parallel_normal_random_generator() {
}

// ============================================================================
FLA_Error MyFLA_Parallel_normal_random_matrix( FLA_Obj A ) {
// Set random matrix with normal distribution.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      double  * buff_A;
      int     m_A, n_A, ldim_A, i, j;

      int          alternate = 0;
      double       b1        = 0.0;
      double       b2        = 0.0;
      unsigned int seed      = global_seed; // + omp_get_thread_num();

      //// printf( "Thread: %d  Global seed: %d\n", 
      ////         omp_get_thread_num(), global_seed );

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width ( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      ldim_A  = FLA_Obj_col_stride( A );

      // Main loop.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i + j * ldim_A ] = MyFLA_Normal_random_number( 0.0, 1.0,
                                         & alternate, & b1, & b2, & seed );
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Parallel_normal_random_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  // Change global_seed to obtain different in the next call.
  #pragma omp atomic
    global_seed += global_num_threads;


  // Show the random block.
#if 0
  #pragma omp critical ( print_parallel_normal_random_matrix ) 
  {
    printf( "MyFLA_Parallel_normal_random_matrix:\n" );
    FLA_Obj_show( "  Arandom = [ ", A, "%18.12le", " ];" );
  }
#endif

  return FLA_SUCCESS;
}

// ============================================================================
static double MyFLA_Normal_random_number( double mu, double sigma,
                  int * alternate, double * b1, double * b2, 
                  unsigned int * seedp ) {
//
// It computes and returns a normal random number.
//
  double  c1, c2, a, factor;

  // Quick return.
  if( ( * alternate ) == 1 ) {
    * alternate = ! ( * alternate );
    return( mu + sigma * ( * b2 ) );
  }
  // Main loop.
  do {
    c1 = -1.0 + 2.0 * myrand( seedp );
    c2 = -1.0 + 2.0 * myrand( seedp );
    a  = c1 * c1 + c2 * c2;
  } while ( ( a == 0.0 )||( a >= 1.0 ) );
  factor                = sqrt( ( -2.0 * log( a ) ) / a );
  * b1                  = c1 * factor;
  * b2                  = c2 * factor;
  * alternate           = ! ( * alternate );
  return( mu + sigma * ( * b1 ) );
}


// ============================================================================
static double myrand( unsigned int * next ) {
  *next = *next * 1103515245 + 12345;
  //// unsigned int uival = (unsigned int) ( *next / 65536 ) % 32768;
  //// double       dval  = ( (double) uival ) * rev2to15;
  //// printf( "  uival: %u  dval: %lf\n", uival, dval );
  //// return dval;
  return ( (double) ( (unsigned int) ( *next / 65536 ) % 32768 ) ) * rev2to15;
}


