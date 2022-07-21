#include <math.h>
#include "FLAME.h"
#include "flush_cache.h"
#include "MyFLA_Utils.h"
#include "check_residuals.h"
#include "compute_svd.h"
#include "FLA_QR_Apply_QT.h"
#include "FLA_UTV_AB_UT_var49h.h"
#include "time_utv.h"


// ============================================================================
// Declaration of local macros and prototypes.

#define max(x,y) ( (x) > (y) ? (x) : (y) )

static double flops_qr( int m, int n );


// ============================================================================
void time_utv( int variant, int num_threads, dim_t nb, int num_execs,
    int print_data, int check_result,
    FLA_Obj A, FLA_Obj Acopy, FLA_Obj sv,
    int build_ort_matrices, int q,
    double * dtime, double * gflops, double * res ) {
//
//
  int      m_A, n_A, datatype, irep;
  double   t1, t2, res_qrf, res_ort, res_ort1, res_ort2, res_svd, res_dif;
  FLA_Obj  mtau, QT, resid;
  int      m_mtau, n_mtau;
  FLA_Obj  Ah, Uh, Vh;  // Matrices with hierarchical storage.
  FLA_Obj  U, V, Ucopy, Vcopy;
  int      build_u, build_v;

  // Some initializations.
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width ( A );
  datatype = FLA_Obj_datatype( A );

  build_u  = build_ort_matrices;
  build_v  = build_ort_matrices;

  // Creating and initializing objects Ucopy and Vcopy.
  FLA_Obj_create( datatype, m_A, m_A, 0, 0, & Ucopy );
  FLA_Obj_create( datatype, n_A, n_A, 0, 0, & Vcopy );
  MyFLA_Set_to_identity( Ucopy );
  MyFLA_Set_to_identity( Vcopy );

  // Creating and initializing objects U and V.
  FLA_Obj_create( datatype, m_A, m_A, 0, 0, & U );
  FLA_Obj_create( datatype, n_A, n_A, 0, 0, & V );

  // Create object for storing tau values.
  // Create matrix "mtau" for storing tau values.
  m_mtau = ( ( FLA_Obj_length( A ) + nb - 1 ) / nb )*nb;
  n_mtau = ( FLA_Obj_width( A ) + nb - 1 ) / nb ;
  FLA_Obj_create( FLA_Obj_datatype( A ), m_mtau, n_mtau, 0, 0, & mtau );
  //// printf( "Created object mtau. Dims: %d x %d \n", mtau.m, mtau.n );
  MyFLA_Obj_set_to_one( mtau );

  // Create hierarchical matrices.
  if( variant / 1000 == 2 ) {
    FLASH_Obj_create( FLA_DOUBLE, m_A, n_A, 1, &nb, & Ah );
  } else if( variant / 1000 == 5 ) {
    FLASH_Obj_create( FLA_DOUBLE, m_A, n_A, 1, &nb, & Ah );
    FLASH_Obj_create( FLA_DOUBLE, m_A, m_A, 1, &nb, & Uh );
    FLASH_Obj_create( FLA_DOUBLE, n_A, n_A, 1, &nb, & Vh );
  }

  // Perform "num_execs" executions.
  for ( irep = 0 ; irep < num_execs; irep++ ) {

    // Set number of threads.
    omp_set_num_threads( num_threads );

    // Copy initial matrix Acopy into working matrix (A or Ah).
    MyFLA_Obj_set_to_one( A );
    if( variant / 1000 == 2 ) {
      // Hierarchical-storage codes.
      FLASH_Obj_hierarchify( Acopy, Ah );
    } else if( variant / 1000 == 5 ) {
      // Hierarchical-storage codes.
      FLASH_Obj_hierarchify( Acopy, Ah );
      FLASH_Obj_hierarchify( Ucopy, Uh );
      FLASH_Obj_hierarchify( Vcopy, Vh );
    } else {
      // Usual-storage codes.
      FLA_Copy( Acopy, A );
      FLA_Copy( Ucopy, U );
      FLA_Copy( Vcopy, V );
    }

    // Print initial data.
    if ( print_data == 1 ) {
      FLA_Obj_show( " ai = [ ", Acopy, "%22.15le", " ];" );
    }

    // Flush cache memory.
    flush_cache_parallel();

    t1 = FLA_Clock();

    // Families of variants:
    //  1* : QR.
    //  2* : QR. Hierarchical storage.
    //  3* : UTV. WY transformations.
    //  4* : UTV. UT transformations (only one included in this driver).
    //  5* : UTV. WY transformations. Hierarchical storage.

    // Execute routine.
    switch( variant ){

    case 5049:
      // Multithreaded UTV: Wave Dataflow. var49a.
      FLA_UTV_AB_UT_var49h( m_A, n_A, Ah, build_u, Uh, build_v, Vh, nb, q, num_threads );
      break;

    default:
      printf( " ERROR in switch of time_utv: unknown case.!!!\n" );
    }

    t2 = FLA_Clock();
    if ( irep == 0 ) {
      *dtime = ( t2 > t1 ? t2 - t1 : 0.0 );
    } else {
      *dtime += ( t2 > t1 ? t2 - t1 : 0.0 );
    }
  }

  // Unpack data from hierarchical to plain storage, for hierarchical codes.
  if( variant / 1000 == 2 ) {
    MyFLA_Obj_set_to_one( A );
    FLASH_Obj_flatten( Ah, A );
  } else if( variant / 1000 == 5 ) {
    MyFLA_Obj_set_to_one( A );
    FLASH_Obj_flatten( Ah, A );
    MyFLA_Obj_set_to_one( U );
    FLASH_Obj_flatten( Uh, U );
    MyFLA_Obj_set_to_one( V );
    FLASH_Obj_flatten( Vh, V );
  }

  // Print final data.
  if( print_data == 1 ) {
    FLA_Obj_show( " af = [ ", A, "%22.15le", " ];" );
    if( ( variant / 1000 == 3 )||
        ( variant / 1000 == 4 )||
        ( variant / 1000 == 5 ) ) {
      FLA_Obj_show( " uf = [ ", U, "%22.15le", " ];" );
      FLA_Obj_show( " vf = [ ", V, "%22.15le", " ];" );
    }
  }

  // Set number of threads again after the execution of the factorizations.
  omp_set_num_threads( num_threads );

  //
  // Check if residuals must be computed.
  //
  if( check_result == 0 ) {
    *res = -1.0;
  } else if( variant < 3000 ) {
    //
    // Slab and incremental QR factorizations.
    //

    FLA_Obj_create( FLA_Obj_datatype( A ), m_A, m_A, 0, 0, & QT );
    FLA_Obj_create( FLA_Obj_datatype( A ),   1,   1, 0, 0, & resid );

    MyFLA_Set_to_identity( QT );

    FLA_QR_Apply_QT( A, mtau, QT, nb );

    check_qr_with_qt( Acopy, QT, A, resid );

    res_qrf = *((double *) FLA_Obj_buffer_at_view( resid ) );

    check_ort( QT, resid );

    res_ort = *((double *) FLA_Obj_buffer_at_view( resid ) );

    check_sv_of_upper_triang_matrix( A, sv, resid );

    res_svd = *((double *) FLA_Obj_buffer_at_view( resid ) );

    FLA_Obj_free( &QT );
    FLA_Obj_free( &resid );

    *res = max( res_qrf, max( res_ort, res_svd ) );

  } else {
    //
    // UTV factorization.
    //
    if( ( build_u == 1 )&&( build_v == 1 ) ) {
      // Full residuals are computed only if both matrices U and V are built.

      FLA_Obj_create( FLA_Obj_datatype( A ), 1,   1,   0, 0, & resid );

      check_utv( Acopy, U, A, V, resid );
      res_dif = *( ( double * ) FLA_Obj_buffer_at_view( resid ) );
      //// FLA_Obj_show( " resid_of_check_utf = [ ", resid, "%15.5e", " ];" );

      check_ort( U, resid );
      res_ort1 = *( ( double * ) FLA_Obj_buffer_at_view( resid ) );
      //// FLA_Obj_show( " resid_of_check_ort_1 = [ ", resid, "%15.5e", " ];" );

      check_ort( V, resid );
      res_ort2 = *( ( double * ) FLA_Obj_buffer_at_view( resid ) );
      //// FLA_Obj_show( " resid_of_check_ort_2 = [ ", resid, "%15.5e", " ];" );

      res_ort = max( res_ort1, res_ort2 );

      check_sv_of_dense_matrix( A, sv, resid );
      res_svd = *( ( double * ) FLA_Obj_buffer_at_view( resid ) );
      //// FLA_Obj_show( " resid_of_check_svd = [ ", resid, "%15.5e", " ];" );

      FLA_Obj_free( & resid );

      * res = max( res_dif, max( res_ort, res_svd ) );

    } else {
      // Only residuals of SVD are computed since matrices U,V have not been 
      // built.

      FLA_Obj_create( FLA_Obj_datatype( A ), 1,   1,   0, 0, & resid );

      check_sv_of_dense_matrix( A, sv, resid );
      res_svd = *( ( double * ) FLA_Obj_buffer_at_view( resid ) );
      //// FLA_Obj_show( " resid_of_check_svd = [ ", resid, "%15.5e", " ];" );

      FLA_Obj_free( & resid );

      * res = res_svd;

    }
  }

  // Free hierarchical matrices.
  if( variant / 1000 == 2 ) {
    FLASH_Obj_free( & Ah );
  } else if( variant / 1000 == 5 ) {
    FLASH_Obj_free( & Ah );
    FLASH_Obj_free( & Uh );
    FLASH_Obj_free( & Vh );
  }

  // Compute and return time and gigaflops per sec.
  * dtime = * dtime / num_execs;
  if( * dtime == 0.0 ) {
    * gflops = -1.0;
  } else {
    * gflops = flops_qr( m_A, n_A  ) / ( *dtime * 1.0e+9 );
  }

  // Remove object for storing tau values.
  FLA_Obj_free( & mtau );

  // Removing objects U and V.
  FLA_Obj_free( & U );
  FLA_Obj_free( & V );

  // Removing objects Ucopy and Vcopy.
  FLA_Obj_free( & Ucopy );
  FLA_Obj_free( & Vcopy );
}

// ============================================================================
static double flops_qr( int m, int n ) {
  double  num_flops, d_m, d_n, d_i;
  int     i;

  d_m = ( double ) m;
  d_n = ( double ) n;
  num_flops = 0.0;
  for( i = 0; i < min( m, n ); i++ ) {
    d_i = ( double ) i;
    num_flops += 4.0 * ( d_m - d_i ) * ( d_n - d_i );
  }
  return num_flops;
}


