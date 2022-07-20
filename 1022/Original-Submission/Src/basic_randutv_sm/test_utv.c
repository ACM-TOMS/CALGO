#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "FLAME.h"
#include "matrix_generate.h"
#include "rank_sv.h"
#include "time_utv.h"


#define fla_utv_var49h_q0_no
#define fla_utv_var49h_q1_no
#define fla_utv_var49h_q2_no

#define fla_utv_var49h_q0_or
#define fla_utv_var49h_q1_or
#define fla_utv_var49h_q2_or

// ============================================================================
// Declaration of local prototypes.

static int test_read_info( int * ndim, int ** mvalues, int ** nvalues,
    int * nnb,  int ** nbvalues,
    int * num_threads, int * num_execs, int * matrix_type,
    int * print_data, int * check_result,
    int * seeds );

static void print_result( char * name, int test, int m, int n, int nb,
    double time, double gflops, double resid );


// ============================================================================
int main( int argc, char *argv[] ) {
  FLA_Obj  A, Acopy, sv;
  int      ndim, idim, * mvalues, * nvalues, m, n,
           nnb,  inb,  * nbvalues, nb,
           num_threads, info, matrix_type, num_execs, print_data, check_result,
           seeds[4],
           variant, itest, check_result_mt,
           build_ort, q;
  double   dtime, gflops, resid, max_resid;
  double   rthresh = 1.0e-12;

  // Initialize FLAME.
  FLA_Init();

  // Check the number of arguments.
  if ( argc != 1 ) {
    printf( "Usage:  %s\n\n", argv[0] );
    exit( -1 );
  }

  // Read the input arguments.
  info = test_read_info( & ndim, & mvalues, & nvalues,
                         & nnb,  & nbvalues,
                         & num_threads, & num_execs, & matrix_type,
                         & print_data, & check_result, seeds );
  if( info != 0 ) {
    printf( "Error in test_read_info\n" );
    exit( -1 );
  }

  // Initialize max_resid.
  max_resid = -1.0;

  itest = 1;
  for ( idim = 0; idim < ndim; idim++ ) {
    m = mvalues[ idim ];
    n = nvalues[ idim ];

    // Allocate space for matrices.
    FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, & A );
    FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, & Acopy );
    FLA_Obj_create( FLA_DOUBLE, min( m, n ), 1, 0, 0, & sv );

    printf( "\n" );
    printf( "%% ===================================================\n" );
    printf( "%% Testing matrix type: %d\n", matrix_type );
    printf( "%% ===================================================\n" );

    // ========================================================================
    // Generate data: m-by-n matrix A.
    // ========================================================================

    check_result_mt = check_result;
    matrix_generate( matrix_type, & check_result_mt, seeds, rthresh, A, sv );
    FLA_Copy( A, Acopy );

    // Compute and print rank of matrix.
    if( matrix_type != 100 ) {
      printf( "%% Rank computed with singular values: %d\n",
              rank_sv( sv, rthresh ) );
    }

    printf( "%%\n" );
    printf( "%%                                m     " );
    printf( "n     nb     time     gflops   resid\n" );
    printf( "%% =====================================" );
    printf( "========================================\n" );
    
    // ========================================================================
    // Test:  UTV AB UT.  var49h. q = 0. Do not build U,V.
    // ========================================================================

#ifdef fla_utv_var49h_q0_no
    variant = 5049;

    for ( inb = 0; inb < nnb; inb++ ) {
      nb        = nbvalues[ inb ];
      build_ort = 0;
      q         = 0;

      time_utv( variant, num_threads, nb, num_execs, print_data,
          check_result_mt, A, Acopy, sv, build_ort, q,
          & dtime, & gflops, & resid );
      if( resid > max_resid ) max_resid = resid;

      print_result( "FLA_UTV_v49h_q0_no", itest, m, n, nb,
                    dtime, gflops, resid );
      fflush( stdout );
    }
    printf( "\n" );
    fflush( stdout );
#endif

    // ========================================================================
    // Test:  UTV AB UT.  var49h. q = 1. Do not build U,V.
    // ========================================================================

#ifdef fla_utv_var49h_q1_no
    variant = 5049;

    for ( inb = 0; inb < nnb; inb++ ) {
      nb        = nbvalues[ inb ];
      build_ort = 0;
      q         = 1;

      time_utv( variant, num_threads, nb, num_execs, print_data,
          check_result_mt, A, Acopy, sv, build_ort, q,
          & dtime, & gflops, & resid );
      if( resid > max_resid ) max_resid = resid;

      print_result( "FLA_UTV_v49h_q1_no", itest, m, n, nb,
                    dtime, gflops, resid );
      fflush( stdout );
    }
    printf( "\n" );
    fflush( stdout );
#endif

    // ========================================================================
    // Test:  UTV AB UT.  var49h. q = 2. Do not build U,V.
    // ========================================================================

#ifdef fla_utv_var49h_q2_no
    variant = 5049;

    for ( inb = 0; inb < nnb; inb++ ) {
      nb        = nbvalues[ inb ];
      build_ort = 0;
      q         = 2;

      time_utv( variant, num_threads, nb, num_execs, print_data,
          check_result_mt, A, Acopy, sv, build_ort, q,
          & dtime, & gflops, & resid );
      if( resid > max_resid ) max_resid = resid;

      print_result( "FLA_UTV_v49h_q2_no", itest, m, n, nb,
                    dtime, gflops, resid );
      fflush( stdout );
    }
    printf( "\n" );
    fflush( stdout );
#endif


    printf( "%% -------------------------------------" );
    printf( "----------------------------------------\n" );

    // ========================================================================
    // Test:  UTV AB UT.  var49h. q = 0. Build U,V.
    // ========================================================================

#ifdef fla_utv_var49h_q0_or
    variant = 5049;

    for ( inb = 0; inb < nnb; inb++ ) {
      nb        = nbvalues[ inb ];
      build_ort = 1;
      q         = 0;

      time_utv( variant, num_threads, nb, num_execs, print_data,
          check_result_mt, A, Acopy, sv, build_ort, q,
          & dtime, & gflops, & resid );
      if( resid > max_resid ) max_resid = resid;

      print_result( "FLA_UTV_v49h_q0_or", itest, m, n, nb,
                    dtime, gflops, resid );
      fflush( stdout );
    }
    printf( "\n" );
    fflush( stdout );
#endif

    // ========================================================================
    // Test:  UTV AB UT.  var49h. q = 1. Build U,V.
    // ========================================================================

#ifdef fla_utv_var49h_q1_or
    variant = 5049;

    for ( inb = 0; inb < nnb; inb++ ) {
      nb        = nbvalues[ inb ];
      build_ort = 1;
      q         = 1;

      time_utv( variant, num_threads, nb, num_execs, print_data,
          check_result_mt, A, Acopy, sv, build_ort, q,
          & dtime, & gflops, & resid );
      if( resid > max_resid ) max_resid = resid;

      print_result( "FLA_UTV_v49h_q1_or", itest, m, n, nb,
                    dtime, gflops, resid );
      fflush( stdout );
    }
    printf( "\n" );
    fflush( stdout );
#endif

    // ========================================================================
    // Test:  UTV AB UT.  var49h. q = 2. Build U,V.
    // ========================================================================

#ifdef fla_utv_var49h_q2_or
    variant = 5049;

    for ( inb = 0; inb < nnb; inb++ ) {
      nb        = nbvalues[ inb ];
      build_ort = 1;
      q         = 2;

      time_utv( variant, num_threads, nb, num_execs, print_data,
          check_result_mt, A, Acopy, sv, build_ort, q,
          & dtime, & gflops, & resid );
      if( resid > max_resid ) max_resid = resid;

      print_result( "FLA_UTV_v49h_q2_or", itest, m, n, nb,
                    dtime, gflops, resid );
      fflush( stdout );
    }
    printf( "\n" );
    fflush( stdout );
#endif

    itest++;

    // Remove matrices.
    FLA_Obj_free( & A );
    FLA_Obj_free( & Acopy );
    FLA_Obj_free( & sv );
  }

  // Remove dynamic vectors.
  free( mvalues );
  free( nvalues );
  free( nbvalues );

  // Finalize FLAME.
  printf( "\n" );
  printf( "%% Maximum residual: %10.2le\n", max_resid );
  printf( "%% End of Program\n" );
  FLA_Finalize();

  return 0;
}

// ============================================================================
static int test_read_info( int * ndim, int ** mvalues, int ** nvalues,
    int * nnb,  int ** nbvalues,
    int * num_threads, int * num_execs, int * matrix_type,
    int * print_data, int * check_result,
    int * seeds ) {

  FILE  * fp;
  int   MAX_LEN_LINE = 2048;
  char  myLine[ MAX_LEN_LINE ];
  int   i, rv;
  char  * rv_pc;

  // Open the file.
  if ( ( fp = fopen( "test.in", "r" ) ) == NULL ) {
    return -1;
  }

  // Read the number of dimensions to test.
  rv = fscanf( fp, "%d", ndim );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  if( * ndim == 0 ) {
    printf( "ERROR: Number of dimensions is zero.\n" );
    return (-1);
  }

  // Create the vector to hold the dimensions "m".
  *mvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    rv = fscanf( fp, "%d", (*mvalues+i) );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );

  // Create the vector to hold the dimensions "n".
  *nvalues = (int *) malloc( *ndim * sizeof( int ) );

  // Read the dimensions.
  for( i = 0; i < * ndim; i++ )
    rv = fscanf( fp, "%d", (*nvalues+i) );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );

  // Write the dimensions.
  printf( "%% Test %d dimensions:\n", * ndim );
  printf( "%%   m: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*mvalues+i) );
  printf( "\n" );

  printf( "%%   n: " );
  for( i = 0; i < * ndim; i++ )
    printf( "%4d ", *(*nvalues+i) );
  printf( "\n" );

  // Read the number of block sizes to test.
  rv = fscanf( fp, "%d", nnb );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  if( * nnb == 0 ) {
    printf( "ERROR: Number of block sizes is zero.\n" );
    return (-1);
  }

  // Create the vector to hold the block sizes.
  *nbvalues = (int *) malloc( *nnb * sizeof( int ) );

  // Read the block sizes.
  for( i = 0; i < * nnb; i++ )
    rv = fscanf( fp, "%d", (*nbvalues+i) );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );

  // Write the block sizes.
  printf( "%% Test %d block sizes:\n", * nnb );
  printf( "%%   nb: " );
  for( i = 0; i < * nnb; i++ )
    printf( "%4d ", *(*nbvalues+i) );
  printf( "\n" );

  // Read the number of threads to use.
  rv = fscanf (fp, "%d", num_threads );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  printf( "%% Number of threads:                %d\n", * num_threads );

  // Read the number of executions.
  rv = fscanf (fp, "%d", num_execs );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  printf( "%% Number of executions:             %d\n", * num_execs );

  // Read the matrix_type.
  rv = fscanf (fp, "%d", matrix_type );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  printf( "%% Matrix type:                      %d\n", * matrix_type );

  // Read whether print data.
  rv = fscanf( fp, "%d", print_data );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  printf( "%% Print_data (0=no;1=yes):          %d\n", * print_data );

  // Read whether check result.
  rv = fscanf( fp, "%d", check_result );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  printf( "%% Check_result (0=no;1=yes):        %d\n", * check_result );

  // Read seeds. Last value must be odd.
  rv = fscanf( fp, "%d", seeds );
  rv = fscanf( fp, "%d", seeds+1 );
  rv = fscanf( fp, "%d", seeds+2 );
  rv = fscanf( fp, "%d", seeds+3 );
  rv_pc = fgets( myLine, MAX_LEN_LINE, fp );
  if( *( seeds+3 ) % 2 == 0 ) {
    printf( "ERROR: Last value of seeds must be odd.\n" );
    return (-1);
  }
  printf( "%% Seeds:                            %d %d %d %d\n",
          *seeds, *(seeds+1), *(seeds+2), *(seeds+3) );

  // Close the file.
  fclose( fp );

  return 0;
}

// ============================================================================
static void print_result( char * name, int test, int m, int n, int nb,
    double time, double gflops, double resid ) {

  printf( "%s( %d, 1:6 )=[ %5d %5d %4d  %7.4le  %6.2lf  %8.1le ];\n",
          name, test, m, n, nb, time, gflops, resid );
}

