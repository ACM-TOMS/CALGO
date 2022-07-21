// ============================================================================
// matrix_generate:
//   Version:  0.05
//   Date:     2018-05-07
// ============================================================================
#include "MyFLA_Utils.h"
#include "compute_svd.h"
#include "matrix_generate.h"

#include <mkl.h>

void dgemg_( int*, int*, int*, int*, double*, double*, int*, double*, double* );

// ============================================================================
void matrix_generate( int mat_type, int * check_result, int * seeds,
    double rthresh, FLA_Obj A, FLA_Obj sv ) {

  int     ldim_A, m_A, n_A, mn_A, lwork, i, j;
  double  * buff_A, * buff_sv, * buff_w, num, rnum;
  FLA_Obj d, w;

  // Matrix types 15, 17, and 19 should not be used because of an error in
  // lapack dlatms.f for rectangular matrices.
  if( ( mat_type == 15 )||
      ( mat_type == 17 )||
      ( mat_type == 19 ) ) {
    mat_type = 14;
  }

  //// printf( "Received seeds: %d %d %d %d\n",
  ////         *(seeds), *(seeds+1), *(seeds+2), *(seeds+3) );
  //// printf( " mat_type: %d\n",  mat_type );
  //// printf( " Rthresh: %le\n", rthresh );

  // Some initializations.
  buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
  ldim_A  = FLA_Obj_col_stride( A );
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width( A );
  mn_A    = min( m_A, n_A );
  buff_sv = ( double * ) FLA_Obj_buffer_at_view( sv );

  // Create matrix.
  if( mat_type == 0 ) {
    //
    // Matrix with integer values.
    // ---------------------------
    //
    num = 1;
    for ( j = 0; j < n_A; j++ ) {
      for ( i = (j % m_A); i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
      for ( i = 0; i < (j % m_A); i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
    if( ( m_A >= 1 )&&( n_A >= 1 ) ) {
      buff_A[ 0 + 0 * ldim_A ] = 1.2;
    }
    // Scale down matrix.
    if( num == 0.0 ) {
      rnum = 1.0;
    } else {
      rnum = 1.0 / num;
    }
    for ( j = 0; j < n_A; j++ ) {
      for ( i = 0; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] *= rnum;
      }
    }
#if 0
#endif
    // Compute svd of matrix.
    compute_svd( A, sv );

  } else if( ( 0 < mat_type )&&( mat_type <= 19 ) ) {
    //
    // Matrices generated with subroutine DGEMG.
    // SUBROUTINE DGEMG( MATTYP, M, N, seeds, RTHRESH, A, LDA, WORK )
    // --------------------------------------------------------------
    //

    // Create workspace need by dgemg.
    lwork = max( max( m_A*n_A + 3*m_A + max( m_A, n_A ),
                      max( 3*min( m_A, n_A ), 8 ) ),
                      3*max( m_A, n_A ) ) + 10000 ;
    //// printf( "lwork: %d \n", lwork );
    FLA_Obj_create( FLA_Obj_datatype( A ), lwork, 1, 0, 0, & w );
    buff_w = ( double * ) FLA_Obj_buffer_at_view( w );

    dgemg_( & mat_type, & m_A, & n_A, seeds, & rthresh, buff_A, & ldim_A,
            buff_sv, buff_w );

    // Remove workspace need by dgemg.
    FLA_Obj_free( & w );

    //// FLA_Obj part, none;
    ////  
    //// // Show the first and the last RETURNED singular values.
    //// FLA_Part_2x1( sv, & part,
    ////                   & none,   min( 5, mn_A ), FLA_TOP );
    //// FLA_Obj_show( " top_returned_sv    = [ ", part, "%24.16le", " ];" );
    //// FLA_Part_2x1( sv, & none,
    ////                   & part,   min( 5, mn_A ), FLA_BOTTOM );
    //// FLA_Obj_show( " bottom_returned_sv = [ ", part, "%24.16le", " ];" );

    // Compute singular values of matrix because LAPACK routines inside 
    // "dgemg" routine do not return accurate singular values of A.
    MyFLA_Obj_set_to_one( sv );
    compute_svd( A, sv );

    //// // Show the first and the last COMPUTED singular values.
    //// FLA_Part_2x1( sv, & part,
    ////                   & none,   min( 5, mn_A ), FLA_TOP );
    //// FLA_Obj_show( " top_computed_sv    = [ ", part, "%24.16le", " ];" );
    //// FLA_Part_2x1( sv, & none,
    ////                   & part,   min( 5, mn_A ), FLA_BOTTOM );
    //// FLA_Obj_show( " bottom_computed_sv = [ ", part, "%24.16le", " ];" );

  } else if( mat_type == 20 ) {
    //
    // Matrix with low rank:  rank = 0.1 * mn_A.
    // -----------------------------------------
    //
    // Create auxiliary objects.
    FLA_Obj_create( FLA_Obj_datatype( A ), min( n_A, m_A ),     1, 0, 0, & d );
    FLA_Obj_create( FLA_Obj_datatype( A ), 3 * max( n_A, m_A ), 1, 0, 0, & w );

    int    seeds[ 4 ] = { 20, 21, 22, 23 };
    int    i_mode     = 0;
    double dmax       = -1.0;
    int    kl         = m_A - 1;
    int    ku         = n_A - 1;
    double * buff_d   = ( double * ) FLA_Obj_buffer_at_view( d );
    double * buff_w   = ( double * ) FLA_Obj_buffer_at_view( w );
    int    irank, info;

    // Create vector with singular values.
    irank = max( 1, mn_A / 10 );
    for( i = 0; i < irank; i++ ) {
      buff_d[ i ] = 1.0;
    }
    for( i = irank; i < mn_A; i++ ) {
      buff_d[ i ] = rthresh * 0.1;
    }
    //// FLA_Obj_show( " dini = [ ", d, "%le", " ]; " );

    // Generate matrix.
    dlatms_( & m_A, & n_A, "U", seeds, "Non-symmetric",
             buff_d, & i_mode, & rthresh, & dmax, & kl, & ku, "no packing",
             buff_A, & ldim_A, buff_w, & info );

    // Remove auxiliary objects.
    FLA_Obj_free( & d );
    FLA_Obj_free( & w );

    // Compute svd of matrix.
    compute_svd( A, sv );
    //// FLA_Obj_show( " Aini = [ ", A, "%le", " ]; " );
    //// FLA_Obj_show( " sv = [ ", sv, "%le", " ]; " );

  } else if( mat_type == 100 ) {
    //
    // Random matrix.
    // --------------
    //
    // This is a special case: No singular values are computed to accelerate
    // the matrix generation. Therefore, no residual checking can be done
    // for this matrix type.
    // This matrix type must be used just for timing, and not for checking.
    //
    FLA_Random_matrix( A );
    MyFLA_Obj_set_to_zero( sv );

    // Check no residual checking is asked.
    if( * check_result == 1 ) {
      printf( "%%\n" );
      printf( "%% WARNING: This matrix type does not allow residual " );
      printf( "checking.\n" );
      printf( "%%          It must be used only for fast generation and " );
      printf( "timings.\n" );
      printf( "%%\n" );
      * check_result = 0;
    }

  } else {
    //
    // Zero matrix.
    //
    MyFLA_Obj_set_to_zero( A );
    MyFLA_Obj_set_to_zero( sv );

  }

  //////FLA_Obj_show( " Aini = [ ", A, "%le", " ]; " );
  //////FLA_Obj_show( " sv = [ ", sv, "%le", " ]; " );
}

