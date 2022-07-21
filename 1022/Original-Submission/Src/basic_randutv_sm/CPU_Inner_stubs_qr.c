#include <stdlib.h>
#include <omp.h>
#include "FLAME.h"
#include "MyFLA_Tools_QR_WY_var101.h"
#include "CPU_Inner_stubs_qr.h"


#undef PRINT_MESSAGES


// ============================================================================
// Declaration of prototypes.

#ifdef PRINT_COMPUTING_TIMINGS
static double flops_comp_dense_qr( int m, int n );
static double flops_comp_td_qr( int m, int n );
static double flops_apply_left_Qt_of_dense_QR( int m, int n, int k );
static double flops_apply_left_Qt_of_td_QR( int m, int n, int k );
#endif


// ============================================================================
FLA_Error CPU_Compute_dense_QR_WY_inner_qr( FLA_Obj A, FLA_Obj t, FLA_Obj S,
    int nb_alg ) {
// Compute QR of dense matrix A.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_A, n_A;
#endif

#ifdef PRINT_MESSAGES
  printf( "  Starting CPU_Compute_dense_QR_WY_inner_qr\n" );
  FLA_Obj_show( " Ai = [ ", A, "%le", " ];" );
  FLA_Obj_show( " ti = [ ", t, "%le", " ];" );
  FLA_Obj_show( " Si = [ ", S, "%le", " ];" );
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  //// printf( "    A:      %5d x %5d \n", A.m, A.n );
  //// printf( "    t:      %5d x %5d \n", t.m, t.n );
  //// printf( "    S:      %5d x %5d \n", S.m, S.n );
  //// printf( "    nb_alg: %5d \n", nb_alg );

  // Compute operation.
  FLA_Compute_dense_QR_WY_var101( A, t, S, nb_alg );

#ifdef PRINT_MESSAGES
  printf( "  End of CPU_Compute_dense_QR_WY_inner_qr\n" );
  FLA_Obj_show( " Af = [ ", A, "%le", " ];" );
  FLA_Obj_show( " tf = [ ", t, "%le", " ];" );
  FLA_Obj_show( " Sf = [ ", S, "%le", " ];" );
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_A = FLA_Obj_length( A );
  n_A = FLA_Obj_width ( A );
  if( tt != 0.0 ) {
    gf = ( flops_comp_dense_qr( m_A, n_A ) / ( tt * 1.0e+9 ) );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Comp_De   %5d x %5d           Time: %le   Speed: %6.1lf \n",
          m_A, n_A, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Compute_td_QR_WY_inner_qr( FLA_Obj U, FLA_Obj D, FLA_Obj t,
    FLA_Obj S, int nb_alg ) {
// Compute:
//   [ U; D ] = QR( [ U; D ] );
// where U is upper triangular.
//
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_D, n_D;
#endif

#ifdef PRINT_MESSAGES
  printf( "CPU_Compute_td_QR_WY_inner_qr\n" );
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Compute operation.
  FLA_Compute_td_QR_WY_var101( U, D, t, S, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_D = FLA_Obj_length( D );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = ( flops_comp_td_qr( m_D, n_D ) / ( tt * 1.0e+9 ) );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Comp_TD   %5d x %5d           Time: %le   Speed: %6.1lf \n",
          m_D, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_left_Qt_of_dense_QR_WY_inner_qr( FLA_Obj U, FLA_Obj S, 
    FLA_Obj C, int nb_alg ) {
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from left side:
//   for i = 1, n_U/nb_alg
//     C := ( I - U_i * S_i * U_i' ) * C.
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_C, n_C, k_U;
#endif

#ifdef PRINT_MESSAGES
  printf( "  Starting CPU_Apply_left_Qt_of_dense_QR_WY_inner_qr\n" );
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Compute operation.
  FLA_Apply_left_Qt_of_dense_QR_WY_var101( U, S, C, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  k_U = FLA_Obj_width ( U );
  if( tt != 0.0 ) {
    gf = ( flops_apply_left_Qt_of_dense_QR( m_C, n_C, k_U ) / 
         ( tt * 1.0e+9 ) );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_De   %5d x %5d x %5d   Time: %le   Speed: %6.1lf\n",
          m_C, n_C, k_U, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_left_Qt_of_td_QR_WY_inner_qr( FLA_Obj D, FLA_Obj S, 
    FLA_Obj F, FLA_Obj G, int nb_alg ) {
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F; G ] from left side.
//   for i = 1, n_U/nb_alg
//     [ F; G ] := ( I - [ I; D_i ] * S_i * [ I; D_i ]' ) * [ F; G ].
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_G, n_G, n_D;
#endif

#ifdef PRINT_MESSAGES
  printf( "CPU_Apply_left_Qt_of_td_QR_WY_inner_qr\n" );
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_left_Qt_of_td_QR_WY_var101( D, S, F, G, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_G = FLA_Obj_length( G );
  n_G = FLA_Obj_width ( G );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = ( flops_apply_left_Qt_of_td_QR( m_G, n_G, n_D ) / ( tt * 1.0e+9 ) );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_TD   %5d x %5d x %5d   Time: %le   Speed: %6.1lf\n",
          m_G, n_G, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
// Copy the contents of A into B.
// This is a special case of general routine FLA_Copy.
// In this routine, A.m can be < B.m
FLA_Error CPU_Mycopy_inner_qr( FLA_Obj A, FLA_Obj B ) {
  FLA_Obj  dst, none;

  //// printf( "  MyFLA_Mycopy\n" );
  if( ( A.m == B.m )&&( A.n == B.n ) ) {
    // Same sizes.
    FLA_Copy( A, B );
  } else {
    // A has less rows than B. Copy only those rows.
    FLA_Part_2x1( B,  & dst,
                      & none,   A.m, FLA_TOP );
    FLA_Copy( A, dst );
  }

  return FLA_SUCCESS;
}

#ifdef PRINT_COMPUTING_TIMINGS
// ============================================================================
static double flops_comp_dense_qr( int m, int n ) {
  double  num_flops, d_m, d_n, d_i;
  int    i;

  d_m = ( double ) m;
  d_n = ( double ) n;
  num_flops = 0.0;
  for( i = 0; i < min( m, n ); i++ ) {
    d_i = ( double ) i;
    num_flops += 4.0f * ( d_m - d_i ) * ( d_n - d_i );
  }
  return num_flops;
}


// ============================================================================
static double flops_comp_td_qr( int m, int n ) {
  double  num_flops, d_m, d_n, d_i;
  int    i;

  d_m = ( double ) m;
  d_n = ( double ) n;
  num_flops = 0.0;
  for( i = 0; i < min( m, n ); i++ ) {
    d_i = ( double ) i;
    num_flops += 4.0f * ( d_m ) * ( d_n - d_i );
  }
  return num_flops;
}


// ============================================================================
static double flops_apply_left_Qt_of_dense_QR( int m, int n, int k ) {
  double  num_flops, d_m, d_n, d_i;
  int    i;

  d_m = ( double ) m;
  d_n = ( double ) n;
  num_flops = 0.0;
  for( i = 0; i < k; i++ ) {
    d_i = ( double ) i;
    num_flops += 4.0f * ( d_m - d_i ) * ( d_n );
  }
  return num_flops;
}


// ============================================================================
static double flops_apply_left_Qt_of_td_QR( int m, int n, int k ) {
  double  num_flops, d_m, d_n;
  int    i;

  d_m = ( double ) m;
  d_n = ( double ) n;
  num_flops = 0.0;
  for( i = 0; i < k; i++ ) {
    num_flops += 4.0f * d_m * d_n;
  }
  return num_flops;
}
#endif

