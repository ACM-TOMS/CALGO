#include <stdio.h>
#include "FLAME.h"
#include "MyFLA_Utils.h"
#include "MyFLA_Tools_QR_UT_var102.h"


// ============================================================================
// Declaration of local prototypes.

static FLA_Error FLA_Apply_left_Qt_of_dense_QR_UT_var102_with_loop( FLA_Obj U, 
    FLA_Obj S, FLA_Obj C, int nb_alg );

static FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_var102_with_loop( FLA_Obj U, 
    FLA_Obj S, FLA_Obj C, int nb_alg );


static FLA_Error FLA_Compute_dense_QR_UT_unb_var102( FLA_Obj A, FLA_Obj t, 
    FLA_Obj S );

static FLA_Error FLA_Compute_td_QR_UT_unb_var102( FLA_Obj U, FLA_Obj D, 
    FLA_Obj t, FLA_Obj S );

static FLA_Error FLA_Update_block_with_left_Qt_of_dense_QR(
    FLA_Obj W, FLA_Obj S, FLA_Obj U1, FLA_Obj U2, FLA_Obj C1, FLA_Obj C2 );

static FLA_Error FLA_Update_block_with_right_Q_of_dense_QR(
    FLA_Obj W, FLA_Obj S, FLA_Obj U1, FLA_Obj U2, FLA_Obj C1, FLA_Obj C2 );

static FLA_Error MyFLA_Apply_H2_UT_l_opd(
    int m_u2_A2,
    int n_a1t,
    double * tau,
    double * u2, int inc_u2,
    double * a1t, int inc_a1t,
    double * A2, int ldim_A2,
    double * workspace );












// ============================================================================
FLA_Error FLA_Compute_dense_QR_UT_var102( FLA_Obj A, FLA_Obj t, FLA_Obj S,
    int nb_alg ) {
//
// Compute QR factorization of dense matrix A.
//
  FLA_Obj  ATL,   ATR,      A00, A01, A02,
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;
  FLA_Obj  tT,              t0,
           tB,              t1,
                            t2;
  FLA_Obj  SL,    SR,       S0,  S1,  W12;
  FLA_Obj  ABR_l, S1_tl, W12_t, none, none2, none3;
  int      b, m_A, n_A;

  //// printf( "FLA_Compute_dense_QR_UT_var102_with_loop\n" );
  //// printf( "  A dims: %d x %d\n", ( int ) A.m, ( int ) A.n );
  //// printf( "  t dims: %d x %d\n", ( int ) t.m, ( int ) t.n );
  //// printf( "  S dims: %d x %d\n", ( int ) S.m, ( int ) S.n );
  //// printf( "    nb_alg: %d \n", nb_alg );
  //// FLA_Obj_show( " Ai = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " ti = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Si = [ ", S, "%le", " ];" );

  // Some initializations.
  m_A = FLA_Obj_length( A );
  n_A = FLA_Obj_width ( A );

  // Initial partitionings.
  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x1( t,    &tT,
                      &tB,            0, FLA_TOP );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );

  while ( ( FLA_Obj_width( ATL ) < n_A )&&( FLA_Obj_length( ATL ) < m_A ) ) {
    b = min( nb_alg, min( FLA_Obj_length( ABR ), FLA_Obj_width( ABR ) ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    FLA_Repart_2x1_to_3x1( tT,                  &t0,
                        /* ** */              /* ****** */
                                                &t1,
                           tB,                  &t2,        b, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &W12,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    // [ U1, t1, S1 ] = QR( [ A11; A21 ], t1, S1 );
    FLA_Part_1x2( ABR, &ABR_l, &none, b, FLA_LEFT );
    FLA_Compute_dense_QR_UT_unb_var102( ABR_l, t1, S1 );

    // S1_tl = FLA_Top_Left_part( S1, b, b );
    FLA_Part_2x2( S1,  & S1_tl, & none,
                       & none2, & none3,   b, b, FLA_TL );

    // W12_t = FLA_Top_part( W12, b );
    FLA_Part_2x1( W12,  & W12_t,
                        & none,   b, FLA_TOP );

    if ( FLA_Obj_width( A12 ) > 0 ) {
      FLA_Update_block_with_left_Qt_of_dense_QR(
          W12_t, S1_tl, A11, A21, A12, A22 );
    }

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &tT,                   t0,
                                                     t1,
                            /* ** */              /* ****** */
                              &tB,                   t2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ W12,
                              FLA_LEFT );
  }

  //// printf( "  End of FLA_Compute_dense_QR_UT_var102_with_loop\n" );
  //// FLA_Obj_show( " Af = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " tf = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Sf = [ ", S, "%le", " ];" );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error FLA_Apply_left_Qt_of_dense_QR_UT_var102( FLA_Obj U, FLA_Obj S, 
    FLA_Obj C, int nb_alg ) {
//
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from left side:
//   for i = 1, n_U/nb_alg
//     C := ( I - U_i * S_i * U_i' ) * C.
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
//
  FLA_Obj  Utl, Sl, none1, none2, none3;
  int      m_Utl, n_Utl;

  //// printf( "  Starting FLA_Apply_left_Qt_of_dense_QR_UT_var102\n" );

  // Extract the top left part of U and the left part of S: 
  // In some cases (when a copy of U is created and used to increase 
  // parallelism), dimensions of U and S can be larger than the right ones, 
  // so the top left part of U and the left part of S must be extracted.
  m_Utl = FLA_Obj_length( C );
  n_Utl = min( m_Utl, FLA_Obj_width( U ) );
  FLA_Part_2x2( U,  & Utl,   & none1,
                    & none2, & none3,  m_Utl, n_Utl, FLA_TL );
  FLA_Part_1x2( S,  & Sl,  & none1,   n_Utl, FLA_LEFT );

  // Apply orthogonal transformations in Utl,Sl to C.
  FLA_Apply_left_Qt_of_dense_QR_UT_var102_with_loop( Utl, Sl, C, nb_alg );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_left_Qt_of_dense_QR_UT_var102_with_loop( FLA_Obj U, 
    FLA_Obj S, FLA_Obj C, int nb_alg ) {
//
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from left side:
//   for i = 1, n_U/nb_alg
//     C := ( I - U_i * S_i * U_i' ) * C.
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
//
  FLA_Obj  UTL, UTR,  U00, U01, U02,
           UBL, UBR,  U10, U11, U12,
                      U20, U21, U22;
  FLA_Obj  SL,  SR,   S0,  S1,  S2;
  FLA_Obj  CT,        C0,
           CB,        C1,
                      C2;
  FLA_Obj  W, UBR_l, S1_tl, W12_t, none, none2, none3;
  int      b, m_U, n_U;

  //// printf( "  Starting FLA_Apply_left_Qt_of_dense_QR_UT_var102_with_loop\n" );

  // Some initializations.
  m_U = FLA_Obj_length( U );
  n_U = FLA_Obj_width ( U );

  // Create object W.
  FLA_Obj_create( FLA_Obj_datatype( U ), nb_alg, FLA_Obj_width( C ), 0, 0,
      & W );
  MyFLA_Obj_set_to_one( W );

  //// printf( "    U:      %d x %d \n", U.m, U.n );
  //// printf( "    t:      %d x %d \n", t.m, t.n );
  //// printf( "    S:      %d x %d \n", S.m, S.n );
  //// printf( "    nb_alg: %d \n", nb_alg );

  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );
  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( ( FLA_Obj_width( UTL ) < n_U )&&( FLA_Obj_length( UTL ) < m_U ) ) {

    b = min( nb_alg, min( FLA_Obj_length( UBR ), FLA_Obj_width( UBR ) ) );

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( CT,                  &C0,
                        /* ** */              /* ** */
                                                &C1,
                           CB,                  &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Some partitionings.
    FLA_Part_1x2( UBR,  & UBR_l, & none,     b, FLA_LEFT );
    FLA_Part_2x2( S1,   & S1_tl, & none,
                        & none2, & none3,    b, b, FLA_TL );
    FLA_Part_2x1( W,    & W12_t,
                        & none,  b, FLA_TOP );

    FLA_Update_block_with_left_Qt_of_dense_QR(
        W12_t, S1_tl, U11, U21, C1, C2 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &CT,                   C0,
                                                     C1,
                            /* ** */              /* ** */
                              &CB,                   C2,     FLA_TOP );
  }

  // Remove object W.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_var102( FLA_Obj U, FLA_Obj S, 
    FLA_Obj C, int nb_alg ) {
//
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from the right side:
//   for i = 1, n_U/nb_alg
//     C := C * ( I - U_i * S_i' * U_i' ).
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
//
  FLA_Obj  Utl, Sl, none1, none2, none3;
  int      m_Utl, n_Utl;

  //// printf( "  Starting FLA_Apply_right_Q_of_dense_QR_UT_var102\n" );

  // Extract the top left part of U and the left part of S: 
  // In some cases (when a copy of U is created and used to increase 
  // parallelism), dimensions of U and S can be larger than the right ones, 
  // so the top left part of U and the left part of S must be extracted.
  m_Utl = FLA_Obj_width( C );
  n_Utl = min( m_Utl, FLA_Obj_length( U ) );
  FLA_Part_2x2( U,  & Utl,   & none1,
                    & none2, & none3,  m_Utl, n_Utl, FLA_TL );
  FLA_Part_1x2( S,  & Sl,  & none1,   n_Utl, FLA_LEFT );

  // Apply orthogonal transformations in U,S to C.
  FLA_Apply_right_Q_of_dense_QR_UT_var102_with_loop( Utl, Sl, C, nb_alg );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_var102_with_loop( FLA_Obj U, 
    FLA_Obj S, FLA_Obj C, int nb_alg ) {
//
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from the right side:
//   for i = 1, n_U/nb_alg
//     C := C * ( I - U_i * S_i' * U_i' ).
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   S    in      nb-by-n  matrix with factors S from factorization.
//   C    in/out  m-by-n   matrix to be updated.
//
  FLA_Obj  UTL, UTR,  U00, U01, U02,
           UBL, UBR,  U10, U11, U12,
                      U20, U21, U22;
  FLA_Obj  SL,  SR,   S0,  S1,  S2;
  FLA_Obj  CL,  CR,   C0,  C1,  C2;
  FLA_Obj  W, UBR_l, S1_tl, W12_t, Ut, none, none2, none3;
  int      b, m_Ut, n_Ut, m_C1;

  //// printf( "  Starting FLA_Apply_right_Q_of_dense_QR_UT_var102_with_loop\n" );

  // Create object W.
  FLA_Obj_create( FLA_Obj_datatype( U ), FLA_Obj_length( C ), nb_alg, 0, 0,
      & W );
  MyFLA_Obj_set_to_one( W );

  //// printf( "    ui:     %d x %d \n", U.m, U.n );
  //// printf( "    si:     %d x %d \n", S.m, S.n );
  //// printf( "    ci:     %d x %d \n", C.m, C.n );
  //// printf( "    nb_alg: %d \n", nb_alg );
  //// printf( "    wi:     %d x %d \n", W.m, W.n );

  //// FLA_Obj_show( " ui = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " si = [ ", S, "%le", " ];" );
  //// FLA_Obj_show( " ci = [ ", C, "%le", " ];" );

  // Partition U.
  FLA_Part_2x1( U,    &Ut,
                      &none,      FLA_Obj_width( C ), FLA_TOP );
  m_Ut = FLA_Obj_length( Ut );
  n_Ut = FLA_Obj_width ( Ut );

  //// printf( "    ut:     %d x %d \n", Ut.m, Ut.n );

  FLA_Part_2x2( Ut,   &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );
  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  while ( ( FLA_Obj_width( UTL ) < n_Ut )&&( FLA_Obj_length( UTL ) < m_Ut ) ) {
    b = min( nb_alg, min( FLA_Obj_length( UBR ), FLA_Obj_width( UBR ) ) );

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    //// printf( "\nIteration:  %d \n", FLA_Obj_width( UTL ) );

    // Some partitionings.
    m_C1 = FLA_Obj_length( C1 );
    FLA_Part_1x2( UBR,  & UBR_l, & none,     b, FLA_LEFT );
    FLA_Part_2x2( S1,   & S1_tl, & none,
                        & none2, & none3,    b, b, FLA_TL );
    FLA_Part_2x2( W,    & W12_t, & none,
                        & none2, & none3,    m_C1, b, FLA_TL );

    FLA_Update_block_with_right_Q_of_dense_QR(
        W12_t, S1_tl, U11, U21, C1, C2 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );
  }

  // Remove object W.
  FLA_Obj_free( & W );

  //// FLA_Obj_show( " uf = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " sf = [ ", S, "%le", " ];" );
  //// FLA_Obj_show( " cf = [ ", C, "%le", " ];" );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error FLA_Compute_td_QR_UT_var102( FLA_Obj U, FLA_Obj D, FLA_Obj t, 
    FLA_Obj S, int nb_alg ) {
//
// Compute:
//   [ U; D ] = QR( [ U; D ] );
// where U is upper triangular.
//
  FLA_Obj  UTL,   UTR,      U00, U01, U02,
           UBL,   UBR,      U10, U11, U12,
                            U20, U21, U22;
  FLA_Obj  DL,    DR,       D0,  D1,  D2;
  FLA_Obj  tT,              t0,
           tB,              t1,
                            t2;
  FLA_Obj  SL,    SR,       S0,  S1,  S2;
  FLA_Obj  S1TL, W12, none1, none2, none3;
  int      b;

  //// printf( "FLA_Compute_td_QR_UT_var102\n" );
  //// FLA_Obj_show( " Uini = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " Dini = [ ", D, "%le", " ];" );

  // Initialize object.
  MyFLA_Obj_set_to_one( S );

  FLA_Part_2x2( U,    &UTL, &UTR,
                      &UBL, &UBR,     0, 0, FLA_TL );
  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );
  FLA_Part_2x1( t,    &tT,
                      &tB,            0, FLA_TOP );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );

  while ( ( FLA_Obj_width ( UTL ) < FLA_Obj_width ( U ) )&&
          ( FLA_Obj_length( UTL ) < FLA_Obj_length( U ) ) ) {

    b = min( nb_alg, min( FLA_Obj_width( UBR ), FLA_Obj_length( UBR ) ) );

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( tT,                &t0,
                        /* ** */            /* ** */
                                              &t1,
                           tB,                &t2,        b, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    // Factorize current block [ U11; D1 ] and save values into t1 and S1.
    FLA_Compute_td_QR_UT_unb_var102( U11, D1, t1, S1 );

    // S1TL = FLA_Top_Left_part( S1, b, b );
    FLA_Part_2x2( S1,  &S1TL,  &none1,
                       &none2, &none3,   b, b, FLA_TL );

    // W12 = FLA_Top_part( S2, b );
    FLA_Part_2x1( S2,  &W12,
                       &none1,  b, FLA_TOP );

    // Update rest of matrix: [ U12; D2 ] with previous transformations.

    // W12 = triu( S1TL )' * ( U12 + U21' * D2 );

    FLA_Copy( U12, W12 );

    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, D1, D2, FLA_ONE, W12 );

    FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
              FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, S1TL, W12 );

    // U12 = U12 -      W12;
    // D2  = D2  - D1 * W12;

    FLA_Axpy( FLA_MINUS_ONE, W12, U12 );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_MINUS_ONE, D1, W12, FLA_ONE, D2 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &tT,                t0,
                                                  t1,
                            /* ** */           /* ** */
                              &tB,                t2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
                              FLA_LEFT );
  }

  //// FLA_Obj_show( " Ufin = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " Dfin = [ ", D, "%le", " ];" );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error FLA_Apply_left_Qt_of_td_QR_UT_var102( FLA_Obj D, FLA_Obj S, 
    FLA_Obj F, FLA_Obj G, int nb_alg ) {
//
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F; G ] from left side.
//   for i = 1, n_U/nb_alg
//     [ F; G ] := ( I - [ I; D_i ] * S_i * [ I; D_i ]' ) * [ F; G ].
//
  FLA_Obj  DL,    DR,       D0,  D1,  D2;
  FLA_Obj  SL,    SR,       S0,  S1,  S2;
  FLA_Obj  FT,              F0,
           FB,              F1,
                            F2;
  FLA_Obj  W, W12_t, S1_tl, none1, none2, none3;
  int      b;

  //// printf( "FLA_Apply_left_Qt_of_td_QR_UT_var102\n" );

  // Create temporal objects.
  FLA_Obj_create( FLA_Obj_datatype( D ), nb_alg, FLA_Obj_width( F ), 0, 0,
      & W );

  // Initial partitionings.
  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );
  FLA_Part_2x1( F,    &FT,
                      &FB,            0, FLA_TOP );

  // Main loop.
  while ( FLA_Obj_width( DL ) < FLA_Obj_width( D ) ){

    b = min( FLA_Obj_width( DR ), nb_alg );
    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( FT,                &F0,
                        /* ** */            /* ** */
                                              &F1,
                           FB,                &F2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // S1_tl = FLA_Top_Left_part( S1, b, b );
    FLA_Part_2x2( S1,  &S1_tl,  &none1,
                       &none2,  &none3,   b, b, FLA_TL );

    // W12_t = FLA_Top_part( W, b );
    FLA_Part_2x1( W,   &W12_t,
                       &none1,  b, FLA_TOP );

    // Update rest of matrix: [ U12; D2 ] with previous transformations.

    // W12_t = triu( S1_tl )' * ( F1 + D1' * G );

    FLA_Copy( F1, W12_t );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, D1, G, FLA_ONE, W12_t );
    FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
              FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, S1_tl, W12_t );

    // F1 = F1 -      W12_t;
    // G  = G  - D1 * W12_t;

    FLA_Axpy( FLA_MINUS_ONE, W12_t, F1 );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_MINUS_ONE, D1, W12_t, FLA_ONE, G );

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &FT,                F0,
                                                  F1,
                            /* ** */           /* ** */
                              &FB,                F2,     FLA_TOP );
  }

  // Remove temporal objects.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

// ============================================================================
FLA_Error FLA_Apply_right_Q_of_td_QR_UT_var102( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, int nb_alg ) {
//
// Apply several block Householder transformations stored
// in ( [ I; D ], S ) to matrix [ F, G ] from the right side.
//   for i = 1, n_U/nb_alg
//     [ F G ] := [ F G ] * ( I - [ I; D_i ] * S_i * [ I; D_i ]' ).
//
  FLA_Obj  DL,    DR,       D0,  D1,  D2;
  FLA_Obj  SL,    SR,       S0,  S1,  S2;
  FLA_Obj  FL,    FR,       F0,  F1,  F2;
  FLA_Obj  W, W12_t, S1_tl, none1, none2, none3;
  int      b;

  //// printf( "FLA_Apply_right_Q_of_td_QR_UT_var102.\n" );

  //// printf( "D dims: %d x %d\n", D.m, D.n );
  //// printf( "S dims: %d x %d\n", S.m, S.n );
  //// printf( "F dims: %d x %d\n", F.m, F.n );
  //// printf( "G dims: %d x %d\n", G.m, G.n );

  //// FLA_Obj_show( " Di = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( " Si = [ ", S, "%le", " ];" );
  //// FLA_Obj_show( " Fi = [ ", F, "%le", " ];" );
  //// FLA_Obj_show( " Gi = [ ", G, "%le", " ];" );


  // Create object for temporal data.
  FLA_Obj_create( FLA_Obj_datatype( D ), FLA_Obj_length( F ), nb_alg, 0, 0,
      & W );

  //// printf( "W dims: %d x %d\n", W.m, W.n );

  // Initial partitionings.
  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );
  FLA_Part_1x2( F,    &FL,  &FR,      0, FLA_LEFT );

  // Main loop.
  while ( FLA_Obj_width( DL ) < FLA_Obj_width( D ) ){
    b = min( FLA_Obj_width( DR ), nb_alg );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( FL,  /**/ FR,        &F0, /**/ &F1, &F2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    // S1_tl = FLA_Top_Left_part( S1, b, b );
    FLA_Part_2x2( S1,  &S1_tl,  &none1,
                       &none2,  &none3,   b, b, FLA_TL );

    // W12_t = FLA_Top_part( W, b );
    FLA_Part_1x2( W,   &W12_t, &none1,  b, FLA_LEFT );

    //// printf( "b:  %d\n", b );
    //// printf( "F1 dims: %d x %d\n", F1.m, F1.n );
    //// printf( "S1_tl dims: %d x %d\n", S1_tl.m, S1_tl.n );
    //// printf( "W12_t dims: %d x %d\n", W12_t.m, W12_t.n );

    // Update rest of matrix: [ U12; D2 ] with previous transformations.

    // W12_t = ( F1 + G * D1 ) * triu( S1_tl );

    FLA_Copy( F1, W12_t );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, G, D1, FLA_ONE, W12_t );
    FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
              FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, S1_tl, W12_t );

    // F1 = F1 - W12_t;
    // G  = G  - W12_t * D1';

    FLA_Axpy( FLA_MINUS_ONE, W12_t, F1 );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
              FLA_MINUS_ONE, W12_t, D1, FLA_ONE, G );

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &FL,  /**/ &FR,        F0, F1, /**/ F2,
                              FLA_LEFT );
  }

  // Remove object for temporal data.
  FLA_Obj_free( & W );

  //// FLA_Obj_show( " Df = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( " Sf = [ ", S, "%le", " ];" );
  //// FLA_Obj_show( " Ff = [ ", F, "%le", " ];" );
  //// FLA_Obj_show( " Gf = [ ", G, "%le", " ];" );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Compute_dense_QR_UT_unb_var102( FLA_Obj A, FLA_Obj t, 
    FLA_Obj S ) {
//
// Compute QR of dense matrix A with unblocked algorithm. 
// Compute and return both t and S.
//
  FLA_Obj  workspace;
  int      j, m_A, n_A, mn_A, dtype_A, ldim_A,
           ldim_S, m_A20, n_A20, m_a21, m_A22, n_A22, n_dB, idx_max_col,
           i_one = 1;
  double   * buff_A, * buff_t, * buff_S,
           * buff_workspace, d_one = 1.0;

  int      idamax_();

  //// printf( "  FLA_Compute_dense_QR_UT_unb_var102\n" );
  //// FLA_Obj_show( " Ain = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " tin = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Sin = [ ", S, "%le", " ];" );

  // Some initializations.
  dtype_A = FLA_Obj_datatype( A );
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width ( A );
  mn_A    = min( m_A, n_A );
  ldim_A  = FLA_Obj_col_stride( A );
  buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
  buff_t  = ( double * ) FLA_Obj_buffer_at_view( t );
  buff_S  = ( double * ) FLA_Obj_buffer_at_view( S );
  ldim_S  = FLA_Obj_col_stride( S );

  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & workspace );
  buff_workspace = ( double * ) FLA_Obj_buffer_at_view( workspace );

  // Main Loop.
  for( j = 0; j < mn_A; j++ ) {
    n_dB  = n_A - j;
    m_a21 = m_A - j - 1;
    m_A22 = m_A - j - 1;
    n_A22 = n_A - j - 1;
    m_A20 = m_A - j - 1;
    n_A20 = j;

    // Compute tau1 and u21 from alpha11 and a21 such that tau1 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    FLA_Househ2_UT_l_opd( m_a21,
                          & buff_A[ j + j * ldim_A ],
                          & buff_A[ ( j+1 ) + j * ldim_A ], i_one,
                          & buff_t[ j ] );

    // / a12t \ =  H / a12t \
    // \ A22  /      \ A22  /
    //
    // where H is formed from tau1 and u21.
    MyFLA_Apply_H2_UT_l_opd(
        m_A22, 
        n_A22, 
        & buff_t[ j ], 
        & buff_A[ ( j+1 ) + j * ldim_A ], 1, 
        & buff_A[ j + ( j+1 ) * ldim_A ], ldim_A, 
        & buff_A[ ( j+1 ) + ( j+1 ) * ldim_A ], ldim_A,
        buff_workspace );

    // Update current column of S.
    // rho11 = tau1;
    buff_S[ j + j * ldim_S ] = buff_t[ j ];
    // t01 = a10t' + A20' * u21;
    dcopy_( & j, & buff_A[ j + 0 * ldim_A ], & ldim_A,
                 & buff_S[ 0 + j * ldim_S ], & i_one );
    dgemv_( "Transpose", & m_A20, & n_A20,
        & d_one,
        & buff_A[ ( j+1 ) + 0 * ldim_A ], & ldim_A,
        & buff_A[ ( j+1 ) + j * ldim_A ], & i_one,
        & d_one,
        & buff_S[ 0 + j * ldim_S ], & i_one );
  }

  FLA_Obj_free( & workspace );

  //// FLA_Obj_show( " Afi = [ ", A, "%le", " ];" );
  //// FLA_Obj_show( " tfi = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Sfi = [ ", S, "%le", " ];" );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Compute_td_QR_UT_unb_var102( FLA_Obj U, FLA_Obj D, 
    FLA_Obj t, FLA_Obj S ) {
//
// Compute QR factorization of [ U; D ], where:
//   U is an upper triangular matrix, and
//   D is a dense matrix.
// Strictly lower triangular part of U will not be accessed.
//
  FLA_Obj  workspace;
  int      dtype_U, m_U, n_U, minmn_U, m_D, ldim_U, ldim_D, ldim_S, i_one = 1,
           m_D_plus_1, j, n_wt2;
  double   * buff_U, * buff_D, * buff_t, * buff_S, * buff_workspace, 
           d_one = 1.0, d_zero = 0.0, d_minus_one = -1.0;

  //// printf( "    FLA_Compute_td_QR_UT_unb_var102\n" );
  //// FLA_Obj_show( " Uini = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " Dini = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( " tini = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Sini = [ ", S, "%le", " ];" );

  // Some initializations.
  dtype_U = FLA_Obj_datatype( U );
  m_U     = FLA_Obj_length( U );
  n_U     = FLA_Obj_width( U );
  buff_U  = ( double * ) FLA_Obj_buffer_at_view( U );
  ldim_U  = FLA_Obj_col_stride( U );
  minmn_U = min( m_U, n_U );
  m_D     = FLA_Obj_length( D );
  buff_D  = ( double * ) FLA_Obj_buffer_at_view( D );
  ldim_D  = FLA_Obj_col_stride( D );
  buff_t  = ( double * ) FLA_Obj_buffer_at_view( t );
  buff_S  = ( double * ) FLA_Obj_buffer_at_view( S );
  ldim_S  = FLA_Obj_col_stride( S );

  // Create temporal object for workspace.
  FLA_Obj_create( dtype_U, n_U, 1, 0, 0, & workspace );
  buff_workspace = ( double * ) FLA_Obj_buffer_at_view( workspace );

  m_D_plus_1 = m_D + 1;

  // Main loop.
  for( j = 0; j < minmn_U; j++ ) {
    n_wt2 = n_U - j - 1;

    // Compute tau1 and u21 from alpha11 and a21 such that tau1 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    FLA_Househ2_UT_l_opd( m_D,
                          & buff_U[ j + j * ldim_U ],
                          & buff_D[ 0 + j * ldim_D ], i_one,
                          & buff_t[ j ] );

    // Update the rest of U and D by applying the Householder transform:
    MyFLA_Apply_H2_UT_l_opd(
        m_D,
        n_wt2, 
        & buff_t[ j ], 
        & buff_D[ 0 + j * ldim_D ], 1, 
        & buff_U[ j + ( j+1 ) * ldim_U ], ldim_U, 
        & buff_D[ 0 + ( j+1 ) * ldim_D ], ldim_D,
        buff_workspace ); 

    // Update current columm of S.
    // rho11 = tau1;
    buff_S[ j + j * ldim_S ] = buff_t[ j ];
    // t01 = A20' * u21;
    dgemv_( "Transpose", & m_D, & j,
        & d_one,
        & buff_D[ 0 + 0 * ldim_D ], & ldim_D,
        & buff_D[ 0 + j * ldim_D ], & i_one,
        & d_zero,
        & buff_S[ 0 + j * ldim_S ], & i_one );
  }

  // Remove temporal object for workspace.
  FLA_Obj_free( & workspace );

  //// FLA_Obj_show( " Ufin = [ ", U, "%le", " ];" );
  //// FLA_Obj_show( " Dfin = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( " tfin = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Sfin = [ ", S, "%le", " ];" );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Update_block_with_left_Qt_of_dense_QR(
    FLA_Obj W, FLA_Obj S, FLA_Obj U1, FLA_Obj U2, FLA_Obj C1, FLA_Obj C2 ) {
//
// Update rest of matrix: [ C1; C2 ] with transformations stored into S,U1,U2.
//
  if ( FLA_Obj_width( C1 ) > 0 ) {
    // W  = triu( S )' * ( U1' * C1 + U2' * C2 );
    FLA_Copy( C1, W );
    FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_ONE, U1, W );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, U2, C2, FLA_ONE, W );
    FLA_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, 
              FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, S, W );

    // C2 = C2 - U2 * W;
    // C1 = C1 - U1 * W;
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_MINUS_ONE, U2, W, FLA_ONE, C2 );
    FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_MINUS_ONE, U1, W );
    FLA_Axpy( FLA_ONE, W, C1 );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Update_block_with_right_Q_of_dense_QR(
    FLA_Obj W, FLA_Obj S, FLA_Obj U1, FLA_Obj U2, FLA_Obj C1, FLA_Obj C2 ) {
//
// Update rest of matrix: [ C1; C2 ] with transformations stored into S,U1,U2.
//
  // W = ( C1 * U1 + C2 * U2 ) * triu( S );
  FLA_Copy( C1, W );
  FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_ONE, U1, W );
  if( FLA_Obj_min_dim( C2 ) > 0 ) {
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, C2, U2, FLA_ONE, W );
  }
  FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR, 
            FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
            FLA_ONE, S, W );

  // C2 = C2 - W * U2';
  // C1 = C1 - W * U1';
  if( FLA_Obj_min_dim( C2 ) > 0 ) {
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
              FLA_MINUS_ONE, W, U2, FLA_ONE, C2 );
  }
  FLA_Trmm( FLA_RIGHT, FLA_LOWER_TRIANGULAR, FLA_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_MINUS_ONE, U1, W );
  FLA_Axpy( FLA_ONE, W, C1 );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Apply_H2_UT_l_opd(
    int m_u2_A2,
    int n_a1t,
    double * tau,
    double * u2, int inc_u2,
    double * a1t, int inc_a1t,
    double * A2, int ldim_A2,
    double * workspace ) {
//
//
  double  one_p       = 1.0;
  double  minus_one_p = -1.0;
  double  rtau;
  int     inc_w1t;

  // FLA_Obj w1t;
  double * w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  w1t = workspace;
  inc_w1t = 1;

  // w1t = a1t;
  dcopy_( & n_a1t,
          a1t, & inc_a1t,
          w1t, & inc_w1t );

  // w1t = w1t + u2' * A2;
  // w1t = w1t + A2^T * conj(u2);
  dgemv_( "Transpose",
          & m_u2_A2,
          & n_a1t,
          & one_p,
          A2, & ldim_A2,
          u2, & inc_u2,
          & one_p,
          w1t, & inc_w1t );

  // w1t = w1t / tau;
  if( * tau == 0.0 ) {
    fprintf( stderr, "ERROR in MyFLA_Apply_H2_UT_l_opd: Tau is zero.\n" );
  } else {
    rtau = 1.0 / ( * tau );
  }
  dscal_( & n_a1t,
          & rtau,
          w1t, & inc_w1t );

  // a1t = - w1t + a1t;
  daxpy_( & n_a1t,
          & minus_one_p,
          w1t, & inc_w1t,
          a1t, & inc_a1t );

  // A2 = - u2 * w1t + A2;
  dger_( & m_u2_A2,
         & n_a1t,
         & minus_one_p,
         u2, & inc_u2,
         w1t, & inc_w1t,
         A2, & ldim_A2 );

  return FLA_SUCCESS;
}

