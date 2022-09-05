// ============================================================================
// compute_svd:
//   Version:  0.02
//   Date:     2018-04-20
// ============================================================================
#include <math.h>
#include "FLAME.h"
#include "compute_svd.h"


// ============================================================================
FLA_Error compute_svd( FLA_Obj A, FLA_Obj sv ) {
// Compute singular values of A.
  FLA_Obj  B, Iwork, Work;
  double   * buff_B, * buff_Work, * buff_U, * buff_VT, * buff_sv, dwork;
  int      * buff_Iwork, info, m_B, n_B, ldim_B, ldim_U, ldim_VT, lwork;
  char     job;

  //// printf( "compute_svd\n" );

  // Make a temporal copy of A into B.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & B );
  FLA_Copy( A, B );
  //// FLA_Obj_show( " A = [ ", A, "%le", " ] ; " );
  //// FLA_Obj_show( " B = [ ", B, "%le", " ] ; " );

  // Some initializations.
  buff_B   = ( double * ) FLA_Obj_buffer_at_view( B );
  ldim_B   = FLA_Obj_col_stride( B );
  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  buff_sv  = FLA_Obj_buffer_at_view( sv );

  job      = 'N';
  buff_U   = NULL;
  ldim_U   = m_B;
  buff_VT  = NULL;
  ldim_VT  = n_B;

  // Create integer workspace.
  FLA_Obj_create( FLA_INT, 8 * min( m_B, n_B ), 1, 0, 0, & Iwork );
  buff_Iwork = ( int * ) FLA_Obj_buffer_at_view( Iwork );

  // Compute optimal length of real workspace.
  lwork = -1;
  //// dgesvd_( & jobu, & jobvt, & m_B, & n_B, buff_B, & ldim_B,
  ////     buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
  ////     & dwork, & lwork, & info );
  dgesdd_( & job, & m_B, & n_B, buff_B, & ldim_B,
      buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
      & dwork, & lwork, buff_Iwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_: %d \n", info );
  }
  lwork = ( int ) dwork;
  //// printf( "  Optimal lwork: %d\n", lwork );

  // Create real workspace.
  FLA_Obj_create( FLA_Obj_datatype( A ), lwork, 1, 0, 0, & Work );
  buff_Work = ( double * ) FLA_Obj_buffer_at_view( Work );

  // Call to SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
  //                            WORK, LWORK, INFO )
  //// dgesvd_( & jobu, & jobvt, & m_B, & n_B, buff_B, & ldim_B,
  ////     buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
  ////     buff_Work, & lwork, & info );
  dgesdd_( & job, & m_B, & n_B, buff_B, & ldim_B,
      buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
      buff_Work, & lwork, buff_Iwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_: %d \n", info );
  }

  // Remove real workspace.
  FLA_Obj_free( & Work );

  // Remove integer workspace.
  FLA_Obj_free( & Iwork );

  // Remove object B.
  FLA_Obj_free( & B );

  return FLA_SUCCESS;
}

