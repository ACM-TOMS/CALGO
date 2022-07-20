#include "MyFLA_Normal_random.h"


// ============================================================================
// Declaration of local variables.

static int     alternate = 0;
static double  b1, b2;


// ============================================================================
// Declaration of local prototypes.

static double MyFLA_Normal_random_number( double mu, double sigma );


// ============================================================================
FLA_Error MyFLA_Normal_random_matrix( FLA_Obj A ) {
// Set random matrix with normal distribution.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      double  * buff_A;
      dim_t     m_A, n_A, ldim_A, i, j;

      // Some initializations.
      m_A     = FLA_Obj_length( A ); //FLASH_Obj_scalar_length( A );
      n_A     = FLA_Obj_width( A ); //FLASH_Obj_scalar_width ( A );
      buff_A  = ( double * ) FLASH_Obj_extract_buffer( A );
      ldim_A  = FLA_Obj_col_stride( A );

      // Added critical section to speed up performances since "rand()" is
      // not thread-safe and a false sharing could be generated.
      // Surprisingly, a new parallel normal random generator does not seem to
      // faster than this code with the critical section.
      #pragma omp critical( normal_random )
      {
        // Main loop.
        for ( j = 0; j < n_A; j++ ) {
          for ( i = 0; i < m_A; i++ ) {
	    //printf("Iter %d, %d (ldim %d) -> buff_A[i+j*ldim_A] = %f\n", i, j, ldim_A, buff_A[i+j*ldim_A]);
            buff_A[ i + j * ldim_A ] = MyFLA_Normal_random_number( 0.0, 1.0 );
            //buff_A[ i + j * ldim_A ] = 77.01;
          }
        }
      }
    }

    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Normal_random_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
void MyFLA_Reset_normal_random_generator( int seed ) {
//
// It resets the normal random generator.
//
  srand( seed );
  alternate = 0;
}

// ============================================================================
static double MyFLA_Normal_random_number( double mu, double sigma ) {
//
// It computes and returns a normal random number.
//
  double  c1, c2, a, factor, denom;

  // Quick return.
  if( alternate == 1 ) {
    alternate = ! alternate;
    return( mu + sigma * b2 );
  }
  // Main loop.
  denom = ( ( double ) RAND_MAX ) + 1.0;
  do {
    c1 = -1.0 + 2.0 * rand() / denom;
    c2 = -1.0 + 2.0 * rand() / denom;
    a  = c1 * c1 + c2 * c2;
  } while ( ( a == 0.0 )||( a >= 1.0 ) );
  factor    = sqrt( ( -2.0 * log( a ) ) / a );
  b1        = c1 * factor;
  b2        = c2 * factor;
  alternate = ! alternate;
  return( mu + sigma * b1 );
}

