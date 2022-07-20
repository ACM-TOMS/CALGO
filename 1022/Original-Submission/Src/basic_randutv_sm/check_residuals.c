// ============================================================================
// check_residuals:
//   Version:  0.03
//   Date:     2018-04-20
// ============================================================================
#include "MyFLA_Utils.h"
#include "compute_svd.h"
#include "check_residuals.h"


// ============================================================================
FLA_Error build_q_matrix( FLA_Obj U, FLA_Obj s, FLA_Obj C ) {
// Build matrix Q with LAPACK's dorgqr routine.
  FLA_Obj  Work, Ul, Ur;
  int      info, datatype, m_C, n_C, mn_C, ldim_C, lwork;
  double   * buff_Work, * buff_s, * buff_C, dwork;

  // Some initializations.
  datatype = FLA_Obj_datatype( C );
  m_C      = FLA_Obj_length( C );
  n_C      = FLA_Obj_width( C );
  mn_C     = FLA_Obj_min_dim( C );
  ldim_C   = FLA_Obj_col_stride( C );
  buff_s   = ( double * ) FLA_Obj_buffer_at_view( s );
  buff_C   = ( double * ) FLA_Obj_buffer_at_view( C );

  // Compute optimal length of workspace.
  //// lwork = n_C * nb_alg;
  lwork = -1;
  dorgqr_( & m_C, & n_C, & mn_C, buff_C, & ldim_C, buff_s,
           & dwork, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, "*** Error in dorgqr_. Info: %d\n", info );
  }
  lwork = ( int ) dwork;
  //// printf( "  Optimal lwork: %d\n", lwork );

  // Create workspace.
  FLA_Obj_create( datatype, lwork, 1, 0, 0, & Work );
  buff_Work = ( double * ) FLA_Obj_buffer_at_view( Work );

  // Copy left part of U into C.
  FLA_Part_1x2( U, & Ul, & Ur,   n_C, FLA_LEFT );
  FLA_Copy( Ul, C );
  //// printf( "Dims Ul: %ld x %ld \n", Ul.m, Ul.n );
  //// printf( "Dims C: %ld x %ld \n", m_C, n_C );
  //// printf( "mn_C: %ld \n", mn_C );

  // Call to dorgqr.
  dorgqr_( & m_C, & n_C, & mn_C, buff_C, & ldim_C, buff_s,
           buff_Work, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, "*** Error in dorgqr_. Info: %d\n", info );
  }
  //// FLA_Obj_show( " Q = [ ", C, "%le", " ];" );

  // Remove workspace.
  FLA_Obj_free( & Work );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_ort( FLA_Obj Q, FLA_Obj resid ) {
// Check orthogonality of a matrix.
// Compute: || In - Q' * Q || / || Q ||, where Q is m-by-n.
  FLA_Obj  Dif, nrm;
  double   dnorm_dif, dnorm_Q;
  int      n_Q;

  // Create objects Dif and nrm.
  n_Q = FLA_Obj_width( Q );
  FLA_Obj_create( FLA_Obj_datatype( Q ), n_Q, n_Q, 0, 0, & Dif );
  FLA_Obj_create( FLA_Obj_datatype( Q ),   1,   1, 0, 0, & nrm );

  // Compute the norm of the difference.
  MyFLA_Set_to_identity( Dif );
  FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, Q, Q, FLA_MINUS_ONE, Dif );
  FLA_Norm_frob( Dif, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );
  //// printf( "   My norm of difference: %le\n", dnorm_dif );

  // Compute the norm of Q.
  FLA_Norm_frob( Q, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_Q );
  //// printf( "   My norm of Q: %le\n", dnorm_Q );

  // Compute and return the quotient.
  if( dnorm_Q == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_Q;
  }

  // Remove objects.
  FLA_Obj_free( & Dif );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_qr_with_q( FLA_Obj A, FLA_Obj Q, FLA_Obj R, FLA_Obj resid ) {
// Check QR factorization, given Q.
// Compute: || A - Q * R || / || A ||, where QT is m-by-m.
  FLA_Obj  Dif, TriuR, nrm, TriuRt, None;
  double   dnorm_dif, dnorm_A;

  // Create objects Dif, TriuR, and nrm.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & Dif );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, & TriuR );
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );

  // Compute: Dif = Q * triu( R ) - A.
  FLA_Copy( A, Dif );
  FLA_Copy( R, TriuR );
  MyFLA_Zero_strict_lower_triangular( TriuR );
  FLA_Part_2x1( TriuR,  & TriuRt,
                        & None,   FLA_Obj_width( Q ), FLA_TOP );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, Q, TriuRt, FLA_MINUS_ONE, Dif );

  // Compute the norm of the difference.
  FLA_Norm_frob( Dif, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );
  //// printf( "Norma de dif: %le\n", dnorm_dif );

  // Compute the norm of A.
  FLA_Norm_frob( A, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_A );
  //// printf( "Norma de A  : %le\n", dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_A;
  }

  // Remove objects Dif, TriuR, and nrm.
  FLA_Obj_free( & Dif );
  FLA_Obj_free( & TriuR );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_qr_with_qt( FLA_Obj A, FLA_Obj QT, FLA_Obj R, FLA_Obj resid ) {
// Check QR factorization, given Q.
// Compute: || A - QT' * R || / || A ||, where QT is m-by-m.
  FLA_Obj  Dif, TriuR, nrm, TriuRt, None;
  double   dnorm_dif, dnorm_A;

  // Create objects Dif, TriuR, and nrm.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & Dif );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, & TriuR );
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );

  // Compute: Dif = Q * triu( R ) - A.
  FLA_Copy( A, Dif );
  FLA_Copy( R, TriuR );
  MyFLA_Zero_strict_lower_triangular( TriuR );
  FLA_Part_2x1( TriuR,  & TriuRt,
                        & None,   FLA_Obj_length( QT ), FLA_TOP );
  FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, QT, TriuRt, FLA_MINUS_ONE, Dif );

  // Compute the norm of the difference.
  FLA_Norm_frob( Dif, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );
  //// printf( "Norma de dif: %le\n", dnorm_dif );

  // Compute the norm of A.
  FLA_Norm_frob( A, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_A );
  //// printf( "Norma de A  : %le\n", dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_A;
  }

  // Remove objects Dif, TriuR, and nrm.
  FLA_Obj_free( & Dif );
  FLA_Obj_free( & TriuR );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_qrcp_with_q( FLA_Obj A, FLA_Obj Q, FLA_Obj R, FLA_Obj p,
    FLA_Obj resid ) {
// Check QR with column pivoting factorization, given Q.
// Compute: || A * P - Q * R || / || A ||.
  FLA_Obj  Dif, TriuR, nrm, TriuRt, None;
  double   dnorm_dif, dnorm_A;
  double   * buff_A, * buff_Dif;
  int      m_A, n_A, ldim_A, ldim_Dif, * buff_p, col, i, j;

  // Create objects Dif, TriuR, and nrm.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & Dif );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, & TriuR );
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );

  // Compute: Dif = Q * triu( R ).
  MyFLA_Obj_set_to_zero( Dif );
  FLA_Copy( R, TriuR );
  MyFLA_Zero_strict_lower_triangular( TriuR );
  FLA_Part_2x1( TriuR,  & TriuRt,
                        & None,   FLA_Obj_width( Q ), FLA_TOP );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, Q, TriuRt, FLA_ZERO, Dif );
  //// FLA_Obj_show( " Q*R = [ ", Dif, "%le", " ];" );

  // Compute: Dif = Dif - A * P.
  buff_A   = ( double * ) FLA_Obj_buffer_at_view( A );
  buff_Dif = ( double * ) FLA_Obj_buffer_at_view( Dif );
  ldim_A   = FLA_Obj_col_stride( A );
  ldim_Dif = FLA_Obj_col_stride( Dif );
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width ( A );
  buff_p   = ( int * ) FLA_Obj_buffer_at_view( p );

  for( j = 0; j < n_A; j++ ) {
    col = buff_p[ j ] - 1;
    for( i = 0; i < m_A; i++ ) {
      buff_Dif[ i + j * ldim_Dif ] -= buff_A[ i + col * ldim_A ];
    }
  }
  //// FLA_Obj_show( " Dif = [ ", Dif, "%le", " ];" );

  // Compute the norm of the difference.
  FLA_Norm_frob( Dif, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );
  //// printf( "Norma de dif: %le\n", dnorm_dif );

  // Compute the norm of A.
  FLA_Norm_frob( A, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_A );
  //// printf( "Norma de A  : %le\n", dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_A;
  }

  // Remove objects Dif, TriuR, and nrm.
  FLA_Obj_free( & Dif );
  FLA_Obj_free( & TriuR );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_sv_of_upper_triang_matrix( FLA_Obj R, FLA_Obj vs,
    FLA_Obj resid ) {
// Check singular values of upper triangular matrix.
// Compute: || singular values of R - vs ||_F,
// where R is upper triangular m-by-n.
  FLA_Obj  Rcopy, d, nrm;
  double   dnorm_vs, dnorm_dif;
  int      min_mn_R;

  // Create a copy of R into Rcopy. Zero strictly lower part of Rcopy.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, R, & Rcopy );
  FLA_Copy( R, Rcopy );
  MyFLA_Zero_strict_lower_triangular( Rcopy );
  //// FLA_Obj_show( " Rcopy = [ ", Rcopy, "%le", " ] ; " );

  // Create auxiliary objects.
  min_mn_R  = FLA_Obj_min_dim( Rcopy );
  FLA_Obj_create( FLA_Obj_datatype( R ), min_mn_R, 1, 0, 0, & d );
  FLA_Obj_create( FLA_Obj_datatype( R ),        1, 1, 0, 0, & nrm );

  // Compute svd.
  compute_svd( Rcopy, d );

  // Compute residuals.
  FLA_Norm_frob( vs, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_vs );

  FLA_Axpy( FLA_MINUS_ONE, vs, d );

  FLA_Norm_frob( d, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );

  if( dnorm_vs == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_vs;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & d );
  FLA_Obj_free( & nrm );

  // Remove object Rcopy.
  FLA_Obj_free( & Rcopy );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_sv_of_dense_matrix( FLA_Obj A, FLA_Obj vs, FLA_Obj resid ) {
// Check singular values of dense matrix.
// Compute: || singular values of A - vs ||_F, where A is m-by-n.
  FLA_Obj  d, nrm;
  double   dnorm_vs, dnorm_dif;
  int      min_mn_A;

  // Create auxiliary objects.
  min_mn_A  = FLA_Obj_min_dim( A );
  FLA_Obj_create( FLA_Obj_datatype( A ), min_mn_A, 1, 0, 0, & d );
  FLA_Obj_create( FLA_Obj_datatype( A ),        1, 1, 0, 0, & nrm );

  // Compute svd of Acopy.
  compute_svd( A, d );
  //// FLA_Obj_show( " d = [ ", d, "%le", " ];" );

  // Compute residuals.
  FLA_Norm_frob( vs, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_vs );

  FLA_Axpy( FLA_MINUS_ONE, vs, d );

  FLA_Norm_frob( d, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );

  if( dnorm_vs == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_vs;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & d );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_svd_factorization( FLA_Obj A, FLA_Obj U, FLA_Obj S, FLA_Obj VT,
    FLA_Obj resid ) {
// Check SVD factorization.
// Compute: || A - U * S * VT ||_F, where A is m-by-m.
  FLA_Obj  US, USVT, nrm;
  double   dnorm_dif, dnorm_A;

  // Create auxiliary objects.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, S, & US );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, S, & USVT );
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );

  // Compute: || A  - U * S * VT ||_F.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, U, S, FLA_ZERO, US );
  FLA_Copy( A, USVT );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, US, VT, FLA_MINUS_ONE, USVT );
  FLA_Norm_frob( USVT, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );

  // Compute: || A ||_F.
  FLA_Norm_frob( A, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_A;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & US );
  FLA_Obj_free( & USVT );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_utv( FLA_Obj A, FLA_Obj U, FLA_Obj T, FLA_Obj V,
    FLA_Obj resid ) {
// Check UTV factorization.
// Compute: || A - U * T * V' ||_F, where A is m-by-m.
  FLA_Obj  UT, UTVt, nrm;
  double   dnorm_dif, dnorm_A;

  // Create auxiliary objects.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, T, & UT );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, T, & UTVt );
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );

  // Compute: || A - U * T * V' ||_F.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_ONE, U, T, FLA_ZERO, UT );
  FLA_Copy( A, UTVt );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
            FLA_ONE, UT, V, FLA_MINUS_ONE, UTVt );
  FLA_Norm_frob( UTVt, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );

  // Compute: || A ||_F.
  FLA_Norm_frob( A, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_A;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & UT );
  FLA_Obj_free( & UTVt );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error check_solution_of_system( FLA_Obj A, FLA_Obj X, FLA_Obj B,
    FLA_Obj resid ) {
// Compute: || A * X - B ||_F.
  FLA_Obj  Bcopy, nrm;
  double   dnorm_b, dnorm_dif;
  FLA_Obj  X_top, none;

  // Create a copy of B into Bcopy.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, & Bcopy );
  FLA_Copy( B, Bcopy );
  //// FLA_Obj_show( " Bcopy = [ ", Bcopy, "%le", " ] ; " );

  // Create auxiliary objects.
  FLA_Obj_create( FLA_Obj_datatype( B ), 1, 1, 0, 0, & nrm );

  // Compute residual.

  FLA_Norm_frob( B, nrm );
  //// FLA_Obj_extract_real_scalar( nrm, & dnorm_b );
  dnorm_b = *( ( double * ) FLA_Obj_buffer_at_view( nrm ) );

  // Bcopy = A * top_of_X - Bcopy.
  FLA_Part_2x1( X,    & X_top,
                      & none,          FLA_Obj_width( A ), FLA_TOP );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_MINUS_ONE, A, X_top, FLA_ONE, Bcopy );

  //// FLA_Obj_show( " diff = [ ", Bcopy, "%le", " ];" );

  FLA_Norm_frob( Bcopy, nrm );
  //// FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );
  dnorm_dif = *( ( double * ) FLA_Obj_buffer_at_view( nrm ) );

  if( dnorm_b == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_b;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & nrm );

  // Remove object Bcopy.
  FLA_Obj_free( & Bcopy );

  return FLA_SUCCESS;
}

