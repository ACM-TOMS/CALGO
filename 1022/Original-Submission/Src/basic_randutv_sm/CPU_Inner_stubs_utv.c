#include "MyFLA_Utils.h"
#include "MyFLA_Tools_QR_UT_var102.h"
#include "CPU_Inner_stubs_utv.h"


// ============================================================================
// Declaration of local prototypes.

static FLA_Error MyFLA_Copy_vector_into_diagonal( FLA_Obj v, FLA_Obj A );

#ifdef PRINT_COMPUTING_TIMINGS
static double flops_comp_dense_qr( int m, int n );
static double flops_comp_td_qr( int m, int n );
static double flops_apply_left_Qt_of_dense_QR( int m, int n, int k );
static double flops_apply_left_Qt_of_td_QR( int m, int n, int k );
static double flops_gemm( int m, int n, int k );
#endif

// ============================================================================
FLA_Error CPU_Gemm_abta_inner_utv( FLA_Obj A, FLA_Obj B ) {
// Compute:  A := B' * A.
//
  FLA_Obj  Acopy, Btl, none2, none3, none4;
  dim_t      m_A, n_A;
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     n_Btl;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  //// printf( "FLA_Gemm_abta_inner_utv.\n" );
  //// FLA_Obj_show( " ai  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " bi  = [ ", B, "%le", " ];" );

  // Some initializations.
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  // Object Btl is created so that the code works fine on rectangular matrices.
  FLA_Part_2x2( B,  & Btl,   & none2,
                    & none3, & none4,   m_A, m_A, FLA_TL );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( A ), m_A, n_A, 0, 0, & Acopy );

  // Compute A := B' * A.
  FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, Btl, A, FLA_ZERO, Acopy );
  FLA_Copy( Acopy, A );

  //// FLA_Obj_show( " acopyf  = [ ", Acopy, "%le", " ];" );
  //// FLA_Obj_show( " af  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " bf  = [ ", B, "%le", " ];" );

  // Remove auxiliary object.
  FLA_Obj_free( & Acopy );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  n_Btl = FLA_Obj_width( Btl );
  if( tt != 0.0 ) {
    gf = flops_gemm( m_A, n_A, n_Btl ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Gemm_abta %5d x %5d x %5d     Time: %le   Speed: %6.1lf \n",
          m_A, n_A, n_Btl, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Gemm_aabt_inner_utv( FLA_Obj A, FLA_Obj B ) {
// Compute:  A := A * B'.
//
  FLA_Obj  Acopy, Btl, none, none2, none3;
  dim_t      m_A, n_A;

#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_Btl;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  //// printf( "FLA_Gemm_aabt_inner_utv.\n" );
  //// FLA_Obj_show( " ai  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " bi  = [ ", B, "%le", " ];" );

  // Some initializations.
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  // Extract the top left part of B.
  FLA_Part_2x2( B,   &Btl,   &none,
                     &none2, &none3,   n_A, n_A, FLA_TL );

  //// FLA_Obj_show( " btl  = [ ", Btl, "%le", " ];" );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( A ), m_A, n_A, 0, 0, & Acopy );

  // Compute A := A * B'.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
            FLA_ONE, A, Btl, FLA_ZERO, Acopy );
  FLA_Copy( Acopy, A );

  //// FLA_Obj_show( " acopyf  = [ ", Acopy, "%le", " ];" );
  //// FLA_Obj_show( " af  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " bf  = [ ", B, "%le", " ];" );

  // Remove auxiliary object.
  FLA_Obj_free( & Acopy );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_Btl = FLA_Obj_length( Btl );
  if( tt != 0.0 ) {
    gf = flops_gemm( m_A, n_A, m_Btl ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Gemm_aabt %5d x %5d x %5d     Time: %le   Speed: %6.1lf \n",
          m_A, n_A, m_Btl, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Gemm_aab_inner_utv( FLA_Obj A, FLA_Obj B ) {
// Compute:  A := A * B.
//
  FLA_Obj  Acopy, Btl, none, none2, none3;
  dim_t      m_A, n_A;
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     n_Btl;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  //// printf( "FLA_Gemm_aab_inner_utv.\n" );
  //// FLA_Obj_show( " ai  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " bi  = [ ", B, "%le", " ];" );

  // Some initializations.
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  // Extract the top left part of B.
  FLA_Part_2x2( B,   &Btl,   &none,
                     &none2, &none3,   n_A, n_A, FLA_TL );

  //// FLA_Obj_show( " btl  = [ ", Btl, "%le", " ];" );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( A ), m_A, n_A, 0, 0, & Acopy );

  // Compute A := A * B.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, A, Btl, FLA_ZERO, Acopy );
  FLA_Copy( Acopy, A );

  //// FLA_Obj_show( " acopyf  = [ ", Acopy, "%le", " ];" );
  //// FLA_Obj_show( " af  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " bf  = [ ", B, "%le", " ];" );

  // Remove auxiliary object.
  FLA_Obj_free( & Acopy );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  n_Btl = FLA_Obj_width( Btl );
  if( tt != 0.0 ) {
    gf = flops_gemm( m_A, n_A, n_Btl ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Gemm_aab  %5d x %5d x %5d     Time: %le   Speed: %6.1lf \n",
          m_A, n_A, n_Btl, tt, gf );
#endif

  return FLA_SUCCESS;
}


// ============================================================================
// Copy the contents of A into B.
// This is a special case of general routine FLA_Copy.
// In this routine, A.m can be < B.m
FLA_Error CPU_Mycopy_inner_utv( FLA_Obj A, FLA_Obj B ) {
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

// ============================================================================
FLA_Error CPU_Compute_dense_QR_UT_inner_utv( FLA_Obj A, FLA_Obj S, dim_t nb_alg ) {
// Compute QR of dense matrix A.
  FLA_Obj  vt;
  dim_t      n_A;
#ifdef PRINT_COMPUTING_TIMINGS
  double   tt1, tt2, tt, gf;
  dim_t      mm_A, nn_A;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Some initializations.
  n_A = FLA_Obj_width ( A );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( A ), n_A, 1, 0, 0, & vt );

  // Perform operation.
  FLA_Compute_dense_QR_UT_var102( A, vt, S, nb_alg );

  // Remove auxiliary object.
  FLA_Obj_free( & vt );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  mm_A = FLA_Obj_length( A );
  nn_A = FLA_Obj_width ( A );
  if( tt != 0.0 ) {
    gf = flops_comp_dense_qr( mm_A, nn_A ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Comp_De   %5d x %5d             Time: %le   Speed: %6.1lf \n",
          mm_A, nn_A, tt, gf );
#endif

  return FLA_SUCCESS;
}


// ============================================================================
FLA_Error CPU_Apply_left_Qt_of_td_QR_UT_inner_utv( FLA_Obj D, FLA_Obj S,
   FLA_Obj F, FLA_Obj G, dim_t nb_alg ) {
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F; G ] from left side.
//   for i = 1, n_U/nb_alg
//     [ F; G ] := ( I - [ I; D_i ] * S_i * [ I; D_i ]' ) * [ F; G ].
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_G, n_G, n_D;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_left_Qt_of_td_QR_UT_var102( D, S, F, G, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_G = FLA_Obj_length( G );
  n_G = FLA_Obj_width ( G );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_td_QR( m_G, n_G, n_D ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_l_TD %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          m_G, n_G, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}


// ============================================================================
FLA_Error CPU_Apply_left_Qt_of_dense_QR_UT_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, dim_t nb_alg ) {
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
  dim_t     m_C, n_C, k_B;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_left_Qt_of_dense_QR_UT_var102( U, S, C, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  k_B = FLA_Obj_width ( U );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_dense_QR( m_C, n_C, k_B ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_l_De %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          m_C, n_C, k_B, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Compute_svd_inner_utv( FLA_Obj U, FLA_Obj A, FLA_Obj VT ) {
// Compute:  SVD of A.
// Return singular values in the diagonal of A,
// and orthogonal matrices in U and VT.
//
  FLA_Obj  sv, Workspace;
  double   * buff_A, * buff_U, * buff_sv, * buff_VT, * buff_Workspace, swork;
  int      nb_alg, info, m_A, n_A, min_mn_A, ldim_A, ldim_U, ldim_VT, lwork;
  char     all = 'A';
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  //// printf( "FLA_Compute_svd_inner_utv.\n" );
  //// FLA_Obj_show( " ai  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " ui  = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " vti = [ ", VT, "%le", " ];" );

  // Some initializations.
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  min_mn_A = min( m_A, n_A );
  buff_A   = ( double * ) FLA_Obj_buffer_at_view( A );
  ldim_A   = FLA_Obj_col_stride( A );
  buff_U   = ( double * ) FLA_Obj_buffer_at_view( U );
  ldim_U   = FLA_Obj_col_stride( U );
  buff_VT  = ( double * ) FLA_Obj_buffer_at_view( VT );
  ldim_VT  = FLA_Obj_col_stride( VT );

  // Create vector sv for storing singular values.
  FLA_Obj_create( FLA_Obj_datatype( A ), min_mn_A, 1, 0, 0, & sv );
  buff_sv  = ( double * ) FLA_Obj_buffer_at_view( sv );

  // Compute optimal workspace length.
  lwork = -1;
  dgesvd_( & all, & all, & m_A, & n_A,
           buff_A, & ldim_A, buff_sv,
           buff_U, & ldim_U, buff_VT, & ldim_VT,
           & swork, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_f: %d \n", info );
  }
  lwork = ( int ) swork;
  //// printf( "  Optimal lwork: %d\n", lwork );

  // Create object Workspace.
  FLA_Obj_create( FLA_Obj_datatype( A ), lwork, 1, 0, 0, & Workspace );
  buff_Workspace = ( double * ) FLA_Obj_buffer_at_view( Workspace );

  // Call to SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
  //                            WORK, LWORK, INFO )
  dgesvd_( & all, & all, & m_A, & n_A,
           buff_A, & ldim_A, buff_sv,
           buff_U, & ldim_U, buff_VT, & ldim_VT,
           buff_Workspace, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_f: %d \n", info );
  }

  // Set matrix A to zero.
  FLA_Set( FLA_ZERO, A );

  // Copy singular values into the diagonal of A.
  MyFLA_Copy_vector_into_diagonal( sv, A );

  //// FLA_Obj_show( " sv  = [ ", sv, "%le", " ];" );
  //// FLA_Obj_show( " af  = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " uf  = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " vtf = [ ", VT, "%le", " ];" );

  // Remove object sv.
  FLA_Obj_free( & sv );

  // Remove object Workspace.
  FLA_Obj_free( & Workspace );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  if( tt != 0.0 ) {
    gf = -1.0;
  } else {
    gf = -1.0;
  }
  printf( "%%%% Svd       %5d x %5d             Time: %le   Speed: %6.1lf \n",
          m_A, n_A, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Copy_vector_into_diagonal( FLA_Obj v, FLA_Obj A ) {
// Copy the contents of vector v into the diagonal of matrix A.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      dim_t    mn_A, rs_A, cs_A, i;
      double  * buff_v, * buff_A;

      buff_v  = ( double * ) FLA_Obj_buffer_at_view( v );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      mn_A    = FLA_Obj_min_dim( A );

      // Copy the vector into the diagonal.
      for ( i = 0; i < mn_A; i++ ) {
        buff_A[ i * rs_A + i * cs_A ] = buff_v[ i ];
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Copy_vector_into_diagonal:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}


// ============================================================================
FLA_Error CPU_Keep_upper_triang_inner_utv( FLA_Obj A ) {
// Keep upper triangular part, and zero strictly lower triangular part.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_A, n_A;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_A = FLA_Obj_length( A );
  n_A = FLA_Obj_width ( A );
  if( tt != 0.0 ) {
    gf = -1.0;
  } else {
    gf = -1.0;
  }
  printf( "%%%% Keep_upp  %5d x %5d             Time: %le   Speed: %6.1lf \n",
          m_A, n_A, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Set_to_zero_inner_utv( FLA_Obj A ) {
// Set to zero.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_A, n_A;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Set( FLA_ZERO, A );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_A = FLA_Obj_length( A );
  n_A = FLA_Obj_width ( A );
  if( tt != 0.0 ) {
    gf = -1.0;
  } else {
    gf = -1.0;
  }
  printf( "%%%% Zero      %5d x %5d             Time: %le   Speed: %6.1lf \n",
          m_A, n_A, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_right_Q_of_dense_QR_UT_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, dim_t nb_alg ) {
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from the right side:
//   for i = 1, n_U/nb_alg
//     C := C * ( I - U_i * S_i' * U_i' ).
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_C, n_C, k_B;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_right_Q_of_dense_QR_UT_var102( U, S, C, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  k_B = FLA_Obj_width ( U );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_dense_QR( n_C, m_C, k_B ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_r_De %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          n_C, m_C, k_B, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_right_Q_of_td_QR_UT_inner_utv( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, dim_t nb_alg ) {
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F, G ] from the right side.
//   for i = 1, n_U/nb_alg
//     [ F G ] := [ F G ] * ( I - [ I; D_i ] * S_i * [ I; D_i ]' ).
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_G, n_G, n_D;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_right_Q_of_td_QR_UT_var102( D, S, F, G, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_G = FLA_Obj_length( G );
  n_G = FLA_Obj_width ( G );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_td_QR( n_G, m_G, n_D ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  gf = -1.0;
  printf( "%%%% Appl_r_TD %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          n_G, m_G, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}


// ============================================================================
FLA_Error CPU_Compute_td_QR_UT_inner_utv( FLA_Obj U, FLA_Obj D, FLA_Obj S,
    dim_t nb_alg ) {
// Compute:
//   [ U; D ] = QR( [ U; D ] );
// where U is upper triangular.
//
  FLA_Obj  vt;
  dim_t      n_U;
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  dim_t     m_D, n_D;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Some initializations.
  n_U = FLA_Obj_width( U );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( U ), n_U, 1, 0, 0, & vt );

  // Perform operation.
  FLA_Compute_td_QR_UT_var102( U, D, vt, S, nb_alg );

  // Remove auxiliary object.
  FLA_Obj_free( & vt );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_D = FLA_Obj_length( D );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = flops_comp_td_qr( m_D, n_D ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Comp_TD   %5d x %5d             Time: %le   Speed: %6.1lf \n",
          m_D, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}


#if 0
// ============================================================================
FLA_Error CPU_Compute_dense_QR_WY_inner_utv( FLA_Obj A, FLA_Obj S, int nb_alg ) {
// Compute QR of dense matrix A.
  FLA_Obj  vt;
  int      n_A;
#ifdef PRINT_COMPUTING_TIMINGS
  double   tt1, tt2, tt, gf;
  int      mm_A, nn_A;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Some initializations.
  n_A = FLA_Obj_width ( A );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( A ), n_A, 1, 0, 0, & vt );

  // Perform operation.
  FLA_Compute_dense_QR_WY_var101( A, vt, S, nb_alg );

  // Remove auxiliary object.
  FLA_Obj_free( & vt );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  mm_A = FLA_Obj_length( A );
  nn_A = FLA_Obj_width ( A );
  if( tt != 0.0 ) {
    gf = flops_comp_dense_qr( mm_A, nn_A ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Comp_De   %5d x %5d             Time: %le   Speed: %6.1lf \n",
          mm_A, nn_A, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_left_Qt_of_dense_QR_WY_inner_utv( FLA_Obj U, FLA_Obj S,
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
  int     m_C, n_C, k_B;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_left_Qt_of_dense_QR_WY_var101( U, S, C, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  k_B = FLA_Obj_width ( U );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_dense_QR( m_C, n_C, k_B ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_l_De %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          m_C, n_C, k_B, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_right_Q_of_dense_QR_WY_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, int nb_alg ) {
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from the right side:
//   for i = 1, n_U/nb_alg
//     C := C * ( I - U_i * S_i' * U_i' ).
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_C, n_C, k_B;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_right_Q_of_dense_QR_WY_var101( U, S, C, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  k_B = FLA_Obj_width ( U );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_dense_QR( n_C, m_C, k_B ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_r_De %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          n_C, m_C, k_B, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Compute_td_QR_WY_inner_utv( FLA_Obj U, FLA_Obj D, FLA_Obj S,
    int nb_alg ) {
// Compute:
//   [ U; D ] = QR( [ U; D ] );
// where U is upper triangular.
//
  FLA_Obj  vt;
  int      n_U;
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_D, n_D;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Some initializations.
  n_U = FLA_Obj_width( U );

  // Create auxiliary object.
  FLA_Obj_create( FLA_Obj_datatype( U ), n_U, 1, 0, 0, & vt );

  // Perform operation.
  FLA_Compute_td_QR_WY_var101( U, D, vt, S, nb_alg );

  // Remove auxiliary object.
  FLA_Obj_free( & vt );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_D = FLA_Obj_length( D );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = flops_comp_td_qr( m_D, n_D ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Comp_TD   %5d x %5d             Time: %le   Speed: %6.1lf \n",
          m_D, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_left_Qt_of_td_QR_WY_inner_utv( FLA_Obj D, FLA_Obj S,
   FLA_Obj F, FLA_Obj G, int nb_alg ) {
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F; G ] from left side.
//   for i = 1, n_U/nb_alg
//     [ F; G ] := ( I - [ I; D_i ] * S_i * [ I; D_i ]' ) * [ F; G ].
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_G, n_G, n_D;
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
    gf = flops_apply_left_Qt_of_td_QR( m_G, n_G, n_D ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Appl_l_TD %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          m_G, n_G, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Apply_right_Q_of_td_QR_WY_inner_utv( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, int nb_alg ) {
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F, G ] from the right side.
//   for i = 1, n_U/nb_alg
//     [ F G ] := [ F G ] * ( I - [ I; D_i ] * S_i * [ I; D_i ]' ).
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_G, n_G, n_D;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Apply_right_Q_of_td_QR_WY_var101( D, S, F, G, nb_alg );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_G = FLA_Obj_length( G );
  n_G = FLA_Obj_width ( G );
  n_D = FLA_Obj_width ( D );
  if( tt != 0.0 ) {
    gf = flops_apply_left_Qt_of_td_QR( n_G, m_G, n_D ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  gf = -1.0;
  printf( "%%%% Appl_r_TD %5d x %5d x %5d     Time: %le   Speed: %6.1lf\n",
          n_G, m_G, n_D, tt, gf );
#endif

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error CPU_Gemm_tn_inner_utv( FLA_Obj alpha, FLA_Obj A, FLA_Obj B,
    FLA_Obj beta, FLA_Obj C ) {
// Compute matrix-matrix product.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_C, n_C, m_A;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  m_A = FLA_Obj_length( A );
  if( tt != 0.0 ) {
    gf = flops_gemm( m_C, n_C, m_A ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Gemm_tn   %5d x %5d x %5d     Time: %le   Speed: %6.1lf \n",
          m_C, n_C, m_A, tt, gf );
#endif

  return FLA_SUCCESS;

}

// ============================================================================
FLA_Error CPU_Gemm_nn_inner_utv( FLA_Obj alpha, FLA_Obj A, FLA_Obj B,
    FLA_Obj beta, FLA_Obj C ) {
// Compute matrix-matrix product.
#ifdef PRINT_COMPUTING_TIMINGS
  double  tt1, tt2, tt, gf;
  int     m_C, n_C, n_A;
#endif

#ifdef PRINT_COMPUTING_TIMINGS
  tt1 = FLA_Clock();
#endif

  // Perform operation.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, alpha, A, B, beta, C );

#ifdef PRINT_COMPUTING_TIMINGS
  tt2 = FLA_Clock();
  tt = tt2 - tt1;
  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width ( C );
  n_A = FLA_Obj_width ( A );
  if( tt != 0.0 ) {
    gf = flops_gemm( m_C, n_C, n_A ) / ( tt * 1.0e+9 );
  } else {
    gf = -1.0;
  }
  printf( "%%%% Gemm_nn   %5d x %5d x %5d     Time: %le   Speed: %6.1lf \n",
          m_C, n_C, n_A, tt, gf );
#endif

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

// ============================================================================
static double flops_gemm( int m, int n, int k ) {
  double  num_flops, d_m, d_n, d_k;

  d_m = ( double ) m;
  d_n = ( double ) n;
  d_k = ( double ) k;
  num_flops = 2.0f * d_m * d_n * d_k;
  return num_flops;
}
#endif
#endif

