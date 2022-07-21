#include <stdlib.h>
#include <omp.h>
#include "flush_cache.h"


// ============================================================================
void flush_cache( void )
{
  int nsize = 6 * 1024 * 1024 + 10240;
  int ndim  = nsize / sizeof( double );
  int nrep  = 50;
  int irep, i;

  double * a;

  // Create object.
  a = (double *) malloc( nsize );

  // Initialize data.
  for( i = 0; i < ndim; i++ ) {
    a[ i ] = 0.0;
  }

  for( irep = 0; irep < nrep; irep++ ) {
    // Write data.
    for( i = 0; i < ndim; i++ ) {
      a[ i ] += 1.0;
    }
  }

  // Free object.
  free( a );
}


// ============================================================================
void flush_cache_parallel( void )
{
  int nsize = 6 * 1024 * 1024 + 10240;
  int ndim  = nsize / sizeof( double );
  int nrep  = 50;
  int irep, i;

  double * a;

  #pragma omp parallel private( a, i, irep )
  {
    // Create object.
    a = (double *) malloc( nsize );

    // Initialize data.
    for( i = 0; i < ndim; i++ ) {
      a[ i ] = 0.0;
    }

    for( irep = 0; irep < nrep; irep++ ) {
      // Write data.
      for( i = 0; i < ndim; i++ ) {
        a[ i ] += 1.0;
      }
    }

    // Free object.
    free( a );
  }
}

