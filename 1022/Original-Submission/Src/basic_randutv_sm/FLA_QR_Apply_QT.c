#include <stdlib.h>
#include <omp.h>
#include "FLAME.h"
#include "MyFLA_Utils.h"
#include "FLA_QR_Apply_QT.h"


// ============================================================================
// Declaration of local prototypes.

static FLA_Error FLA_QR_apply_QT_of_column_block( FLA_Obj U, FLA_Obj Tau,
    FLA_Obj C, int nb_alg );

static FLA_Error FLA_QR_apply_left_Qt_of_dense_QR_WY_to_row_block( FLA_Obj U, 
    FLA_Obj tau, FLA_Obj C, int nb_alg );

static FLA_Error FLA_QR_apply_left_Qt_of_td_QR_WY_to_row_block( FLA_Obj D,
    FLA_Obj t, FLA_Obj F, FLA_Obj G, int nb_alg );

static FLA_Error FLA_Apply_left_Qt_of_dense_QR_WY_var50a( FLA_Obj U, FLA_Obj t, 
    FLA_Obj C, int nb_alg );

static FLA_Error FLA_Apply_left_Qt_of_td_QR_WY_var50a( FLA_Obj D, FLA_Obj t, 
    FLA_Obj F, FLA_Obj G, int nb_alg );

static FLA_Error MyFLA_Compute_dense_S_unb_var101( FLA_Obj A, FLA_Obj t, 
    FLA_Obj S );

static FLA_Error MyFLA_Compute_tdqr_S_unb_var101( FLA_Obj D, FLA_Obj t, 
    FLA_Obj S );



// ============================================================================
FLA_Error FLA_QR_Apply_QT( FLA_Obj U, FLA_Obj Tau, FLA_Obj C, int nb_alg ) {
// Apply QT defined in (U,Tau) to C:   C := QT * C.
  FLA_Obj  UTL, UTR,   U00, U01, U02,
           UBL, UBR,   U10, U11, U12,
                       U20, U21, U22;
  FLA_Obj  TTL, TTR,   T00, T01, T02,
           TBL, TBR,   T10, T11, T12,
                       T20, T21, T22;
  FLA_Obj  CT,         C0,
           CB,         C1,
                       C2;
  FLA_Obj  UBR_l, TBR_l, none;
  int      b;

  //// printf( "\n\n\n" );

  FLA_Part_2x2( U,     &UTL, &UTR,
                       &UBL, &UBR,   0, 0, FLA_TL );
  FLA_Part_2x2( Tau,   &TTL, &TTR,
                       &TBL, &TBR,   0, 0, FLA_TL );
  FLA_Part_2x1( C,     &CT,
                       &CB,          0, FLA_TOP );

  while ( ( FLA_Obj_width ( UTL ) < FLA_Obj_width ( U ) ) &&
          ( FLA_Obj_length( UTL ) < FLA_Obj_length( U ) ) ) {
    b = min( nb_alg, min( FLA_Obj_length( UBR ), FLA_Obj_width( UBR ) ) );
    if ( 0 == b ) break;

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );
    FLA_Repart_2x2_to_3x3( TTL, /**/ TTR,       &T00, /**/ &T01, &T02,
                        /* ************* */   /* ******************** */
                                                &T10, /**/ &T11, &T12,
                           TBL, /**/ TBR,       &T20, /**/ &T21, &T22,
                           b, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                        /* ** */            /* ** */
                                              &C1,
                           CB,                &C2,        b, FLA_BOTTOM );
    /*-----------------------------------------------------------------------*/

    //// printf( " FLA_QR_Apply_QT. Iter: %d \n", UTL.m );

    //// printf( "Iter:  %d \n", FLA_Obj_width( UTL ) );
    // Some partitionings.
    FLA_Part_1x2( UBR,  & UBR_l, & none,    b, FLA_LEFT );
    FLA_Part_1x2( TBR,  & TBR_l, & none,    1, FLA_LEFT );

    // Apply QT from UBR.
    //// printf( " FLA_QR_Apply_QT. Chivato 1 \n" );
    FLA_QR_apply_QT_of_column_block( UBR_l, TBR_l, CB, nb_alg );
    //// printf( " FLA_QR_Apply_QT. Chivato 2 \n" );
    //// printf( "End of Iter:  %d \n", FLA_Obj_width( UTL ) );
    //// FLA_Obj_show( " QTi = [ ", C, "%le", " ];" );

    /*-----------------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &TTL, /**/ &TTR,       T00, T01, /**/ T02,
                                                     T10, T11, /**/ T12,
                            /* ************** */  /* ****************** */
                              &TBL, /**/ &TBR,       T20, T21, /**/ T22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                  C1,
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }
  //// printf( "End of main routine: \n" );

  //// FLA_Obj_show( " QT = [ ", C, "%le", " ];" );

  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_QR_apply_QT_of_column_block( FLA_Obj U, FLA_Obj Tau,
    FLA_Obj C, int nb_alg ) {
// Apply QT defined in (U,Tau) to C:   C := QT * C.
  FLA_Obj  UT,     U0,
           UB,     U1,
                   U2;
  FLA_Obj  TT,     T0,
           TB,     T1,
                   T2;
  FLA_Obj  CT,     C0,
           CB,     C1,
                   C2;
  int      b, first_nb;
  FLA_Obj  Ct;

  //// printf( "  Beginning FLA_QR_apply_QT_of_column_block\n" );
  //// printf( "  >>> U dims:   %d x %d\n", U.m, U.n );
  //// printf( "  >>> Tau dims: %d x %d\n", Tau.m, Tau.n );
  //// printf( "  >>> C dims:   %d x %d\n", C.m, C.n );
  //// printf( "  >>> nb_alg:   %d \n", nb_alg );

  // Quick return.
  if( ( FLA_Obj_length( U ) == 0 )||( FLA_Obj_width( U ) == 0 )||
      ( FLA_Obj_length( Tau ) == 0 )||( FLA_Obj_width( Tau ) == 0 )||
      ( FLA_Obj_length( C ) == 0 )||( FLA_Obj_width( C ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  first_nb = min( nb_alg, FLA_Obj_length( U ) );

  // Initial partitioning.
  FLA_Part_2x1( U,    &UT,
                      &UB,            first_nb, FLA_TOP );
  FLA_Part_2x1( Tau,  &TT,
                      &TB,            first_nb, FLA_TOP );
  FLA_Part_2x1( C,    &CT,
                      &CB,            first_nb, FLA_TOP );

  // Update CT with transformations from factorization of dense UT.
  FLA_QR_apply_left_Qt_of_dense_QR_WY_to_row_block( UT, TT, CT, nb_alg );
  //// FLA_Obj_show( " QT_after_qr = [ ", C, "%le", " ];" );

  Ct = CT;

  while ( FLA_Obj_length( UT ) < FLA_Obj_length( U ) ) {
    b = min( FLA_Obj_length( UB ), nb_alg );
    //// printf( "Row block. UT.m:    %d\n", FLA_Obj_length( UT ) );
    //// printf( "Row block. nb_alg:  %d\n", nb_alg );
    //// printf( "Row block. b:       %d\n", b );
    FLA_Repart_2x1_to_3x1( UT,                &U0,
                                              &U1,
                        /* ** */            /* ** */
                           UB,                &U2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( TT,                &T0,
                                              &T1,
                        /* ** */            /* ** */
                           TB,                &T2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                                              &C1,
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Update [ Ct; C1 ] with transformations from QR factorization of
    // [ I; U1 ], where U1 is dense.

    //// printf( "Row block. U1 dims: %d x %d\n", U1.m, U1.n );
    //// printf( "Row block. T1 dims: %d x %d\n", T1.m, T1.n );
    //// printf( "Row block. Ct dims: %d x %d\n", Ct.m, Ct.n );
    //// printf( "Row block. C1 dims: %d x %d\n", C1.m, C1.n );

    FLA_QR_apply_left_Qt_of_td_QR_WY_to_row_block( U1, T1, Ct, C1, nb_alg );
    //// FLA_Obj_show( " QT_after_tdqr = [ ", C, "%le", " ];" );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &UT,                U0,
                            /* ** */           /* ** */
                                                  U1,
                              &UB,                U2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &TT,                T0,
                            /* ** */           /* ** */
                                                  T1,
                              &TB,                T2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                            /* ** */           /* ** */
                                                  C1,
                              &CB,                C2,     FLA_TOP );
  }
  //// printf( "  End of FLA_QR_apply_QT_of_column_block\n" );
  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_QR_apply_left_Qt_of_dense_QR_WY_to_row_block( FLA_Obj U, 
    FLA_Obj tau, FLA_Obj C, int nb_alg ) {
// Apply QT defined in (U,tau) to C:   C := QT * C.

  // Quick return.
  if( ( FLA_Obj_length( U ) == 0 )||( FLA_Obj_width( U ) == 0 )||
      ( FLA_Obj_length( tau ) == 0 )||( FLA_Obj_width( tau ) == 0 )||
      ( FLA_Obj_length( C ) == 0 )||( FLA_Obj_width( C ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  //// printf( "  Beginning FLA_QR_apply_left_Qt_of_dense_QR_WY_to_row_block\n" );
  FLA_Apply_left_Qt_of_dense_QR_WY_var50a( U, tau, C, nb_alg );
  //// printf( "  End of FLA_QR_apply_left_Qt_of_dense_QR_WY_to_row_block\n" );
  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_QR_apply_left_Qt_of_td_QR_WY_to_row_block( FLA_Obj D,
    FLA_Obj tau, FLA_Obj F, FLA_Obj G, int nb_alg ) {
// Apply QT defined in (D,tau) to [ F; G ]:   [ F; G ] := QT * [ F; G ].

  // Quick return.
  if( ( FLA_Obj_length( D ) == 0 )||( FLA_Obj_width( D ) == 0 )||
      ( FLA_Obj_length( tau ) == 0 )||( FLA_Obj_width( tau ) == 0 )||
      ( FLA_Obj_length( F ) == 0 )||( FLA_Obj_width( F ) == 0 )||
      ( FLA_Obj_length( G ) == 0 )||( FLA_Obj_width( G ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  //// printf( "  Beginning FLA_Apply_left_Qt_of_td_QR_WY_to_row_block \n" );
  FLA_Apply_left_Qt_of_td_QR_WY_var50a( D, tau, F, G, nb_alg );
  //// printf( "  End of FLA_Apply_left_Qt_of_td_QR_WY_to_row_block \n" );
  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_left_Qt_of_dense_QR_WY_var50a( FLA_Obj U, FLA_Obj t, 
    FLA_Obj C, int nb_alg ) {
// Apply several block Householder transformations stored in (U_i,t_i)
// to matrix C from left side:
//   for i = 1, n_U
//     C := ( I - U_i * t_i * U_i' ) * C.
// where:
//   U    in      m-by-k   matrix with Householder vectors.
//   t    in      nb-by-k  matrix with factors tau from factorization.
//   C    in/out  m-by-n   matrix to be updated.
  FLA_Obj  UTL,   UTR,      U00, U01, U02,
           UBL,   UBR,      U10, U11, U12,
                            U20, U21, U22;
  FLA_Obj  tT,              t0,
           tB,              t1,
                            t2;
  FLA_Obj  CT,              C0,
           CB,              C1,
                            C2;
  FLA_Obj  S, W, UBR_l, S1_tl, W12_t, none, none2, none3;
  int      b, m_U, n_U;

  //// printf( "  Starting FLA_Apply_left_Qt_of_dense_QR_WY_var50a\n" );

  // Some initializations.
  m_U = FLA_Obj_length( U );
  n_U = FLA_Obj_width ( U );

  // Create object S.
  FLA_Obj_create( FLA_Obj_datatype( U ), nb_alg, nb_alg, 0, 0, & S );
  MyFLA_Obj_set_to_one( S );

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
  FLA_Part_2x1( t,    &tT,
                      &tB,            0, FLA_TOP );
  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( ( FLA_Obj_width( UTL ) < n_U )&&( FLA_Obj_length( UTL ) < m_U ) ) {

    b = min( nb_alg, min( FLA_Obj_length( UBR ), FLA_Obj_width( UBR ) ) );

    FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                        /* ************* */   /* ******************** */
                                                &U10, /**/ &U11, &U12,
                           UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                           b, b, FLA_BR );
    FLA_Repart_2x1_to_3x1( tT,                  &t0,
                        /* ** */              /* ****** */
                                                &t1,
                           tB,                  &t2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( CT,                  &C0,
                        /* ** */              /* ** */
                                                &C1,
                           CB,                  &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Some partitionings.
    FLA_Part_1x2( UBR,  & UBR_l, & none,     b, FLA_LEFT );
    FLA_Part_2x2( S,    & S1_tl, & none,
                        & none2, & none3,    b, b, FLA_TL );
    FLA_Part_2x1( W,    & W12_t,
                        & none,  b, FLA_TOP );

    // Compute factor S.
    MyFLA_Compute_dense_S_unb_var101( UBR_l, t1, S1_tl );

    // U11 = trilu( U11 );
    // U21 = U21;
    // W12_t = triu( S1_tl )' * ( U11' * C1 + U21' * C2 );
    FLA_Copy( C1, W12_t );
    FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
              FLA_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_ONE, U11, W12_t );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, U21, C2, FLA_ONE, W12_t );
    FLA_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
              FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, S1_tl, W12_t );

    // C2 = C2 - U21 * W12_t;
    // C1 = C1 - U11 * W12_t;
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_MINUS_ONE, U21, W12_t, FLA_ONE, C2 );
    FLA_Trmm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
              FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_MINUS_ONE, U11, W12_t );
    FLA_Axpy( FLA_ONE, W12_t, C1 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                     U10, U11, /**/ U12,
                            /* ************** */  /* ****************** */
                              &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &tT,                   t0,
                                                     t1,
                            /* ** */              /* ****** */
                              &tB,                   t2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &CT,                   C0,
                                                     C1,
                            /* ** */              /* ** */
                              &CB,                   C2,     FLA_TOP );
  }

  // Remove object S.
  FLA_Obj_free( & S );

  // Remove object W.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_left_Qt_of_td_QR_WY_var50a( FLA_Obj D, FLA_Obj t, 
    FLA_Obj F, FLA_Obj G, int nb_alg ) {
// Apply several block Householder transformations stored
// in ( [ U; D ], t, S ) to matrix [ F; G ] from left side.
//   for i = 1, n_U/nb_alg
//     [ F; G ] := ( I - [ U_i; D_i ] * S_i * [ U_i; D_i ]' ) * [ F; G ].
  FLA_Obj  DL,    DR,       D0,  D1,  D2;
  FLA_Obj  tT,              t0,
           tB,              t1,
                            t2;
  FLA_Obj  SL,    SR,       S0,  S1,  S2;
  FLA_Obj  FT,              F0,
           FB,              F1,
                            F2;
  int      b, bt;
  FLA_Obj  W, W12, S, S1TL, none1, none2, none3;

  //// printf( "FLA_Apply_left_Qt_of_td_QR_WY_var50a\n" );

  // Create temporal objects.
  FLA_Obj_create( FLA_Obj_datatype( D ), nb_alg, FLA_Obj_width( F ), 0, 0,
      & W );
  FLA_Obj_create( FLA_Obj_datatype( D ), nb_alg, nb_alg, 0, 0, & S );

  FLA_Part_1x2( D,    &DL,  &DR,      0, FLA_LEFT );
  FLA_Part_2x1( t,    &tT,
                      &tB,            0, FLA_TOP );
  FLA_Part_1x2( S,    &SL,  &SR,      0, FLA_LEFT );
  FLA_Part_2x1( F,    &FT,
                      &FB,            0, FLA_TOP );

  //// printf( "Dims D: %d x %d \n", D.m, D.n );
  //// printf( "Dims t: %d x %d \n", t.m, t.n );
  //// printf( "Dims S: %d x %d \n", S.m, S.n );
  //// printf( "Dims F: %d x %d \n", F.m, F.n );

  //// FLA_Obj_show( "  Di = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( "  ti = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( "  Fi = [ ", F, "%le", " ];" );
  //// FLA_Obj_show( "  Gi = [ ", G, "%le", " ];" );

  while ( FLA_Obj_width( DL ) < FLA_Obj_width( D ) ) {

    b = min( nb_alg, FLA_Obj_width( DR ) );
    bt = min( nb_alg, FLA_Obj_width( tB ) );

    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( tT,                &t0,
                        /* ** */            /* ** */
                                              &t1,
                           tB,                &t2,        bt, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( FT,                &F0,
                        /* ** */            /* ** */
                                              &F1,
                           FB,                &F2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    //// printf( "  Iter of FLA_Apply_left_Qt_of_td_QR_WY_var50a: %d \n", DL.m );
    //// printf( "      D1. Dims: %d x %d\n", D1.m, D1.n );
    //// printf( "      t1. Dims: %d x %d\n", t1.m, t1.n );
    //// printf( "      S1. Dims: %d x %d\n", S1.m, S1.n );
    //// printf( "      W.  Dims: %d x %d\n", W.m, W.n );

    MyFLA_Obj_set_to_one( S1 );
    MyFLA_Compute_tdqr_S_unb_var101( D1, t1, S1 );
    //// FLA_Obj_show( "         Si = [ ", S1, "%le", " ];" );

    // S1TL = FLA_Top_Left_part( S1, b, b );
    FLA_Part_2x2( S1,  &S1TL,  &none1,
                       &none2, &none3,   b, b, FLA_TL );

    // W12 = FLA_Top_part( W, b );
    FLA_Part_2x1( W,   &W12,
                       &none1,  b, FLA_TOP );

    // Update rest of matrix: [ F1; G ] with previous transformations.

    // W12 = triu( S1TL )' * ( F1 + D1' * G );
    FLA_Copy( F1, W12 );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, D1, G, FLA_ONE, W12 );

    FLA_Trmm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
              FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, S1TL, W12 );

    // F1 = F1 -      W12;
    // G  = G  - D1 * W12;
    FLA_Axpy( FLA_MINUS_ONE, W12, F1 );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_MINUS_ONE, D1, W12, FLA_ONE, G );
    //// printf( "    End of iter of FLA_Apply_left_Qt_of_td_QR_WY_var50a\n" );

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &tT,                   t0,
                                                     t1,
                            /* ** */              /* ** */
                              &tB,                   t2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &FT,                   F0,
                                                     F1,
                            /* ** */              /* ** */
                              &FB,                   F2,     FLA_TOP );
  }

  // Remove temporal objects.
  FLA_Obj_free( & W );
  FLA_Obj_free( & S );

  //// printf( "End of FLA_Apply_left_Qt_of_td_QR_WY_var50a\n" );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Compute_dense_S_unb_var101( FLA_Obj A, FLA_Obj t, 
    FLA_Obj S ) {
// Compute matrix S from previously factorized matrix A, and column vector t.
#if 0
  FLA_Obj  ATL,   ATR,     A00,  a01,     A02,
           ABL,   ABR,     a10t, alpha11, a12t,
                           A20,  a21,     A22;
  FLA_Obj  STL,   STR,     S00,  s01,     S02,
           SBL,   SBR,     s10t, sigma11, s12t,
                           S20,  s21,     S22;
  FLA_Obj  tT,             t0,
           tB,             tau1,
                           t2;
#endif
  int      m_A, n_A, ldim_A, ldim_S;
  double   * buff_A, * buff_t, * buff_S;

  //// printf( "    MyFLA_Compute_dense_S_unb_var101\n" );
  //// printf( "    Dims A:  %d x %d\n", A.m, A.n );
  //// printf( "    Dims of A: %d x %d \n", A.m, A.n );
  //// printf( "    Dims of t: %d x %d \n", t.m, t.n );
  //// printf( "    Dims of S: %d x %d \n", S.m, S.n );

  // Some initializations.
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width( A );
  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  ldim_A = FLA_Obj_col_stride( A );
  buff_t = ( double * ) FLA_Obj_buffer_at_view( t );
  buff_S = ( double * ) FLA_Obj_buffer_at_view( S );
  ldim_S = FLA_Obj_col_stride( S );

  // Quick return.
  if( ( m_A == 0 )||( n_A == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Perform operation.
  dlarft_( "Forward", "Columnwise", & m_A, & n_A, buff_A, & ldim_A, buff_t,
           buff_S, & ldim_S );

#if 0
  // The main loop starts here.
  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );
  FLA_Part_2x2( S,    &STL, &STR,
                      &SBL, &SBR,     0, 0, FLA_TL );
  FLA_Part_2x1( t,    &tT,
                      &tB,            0, FLA_TOP );
  while ( ( FLA_Obj_length( ATL ) < m_A )&&( FLA_Obj_width( ATL )  < n_A ) ) {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );
    FLA_Repart_2x2_to_3x3( STL, /**/ STR,       &S00,  /**/ &s01,     &S02,
                        /* ************* */   /* ************************** */
                                                &s10t, /**/ &sigma11, &s12t,
                           SBL, /**/ SBR,       &S20,  /**/ &s21,     &S22,
                           1, 1, FLA_BR );
    FLA_Repart_2x1_to_3x1( tT,                &t0,
                        /* ** */            /* **** */
                                              &tau1,
                           tB,                &t2,        1, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Update S: s01  = - triu( S00 ) * ( a10t' + A20' * a21 ) * sigma11;
    FLA_Copy( tau1, sigma11 );

    FLA_Copyt( FLA_TRANSPOSE, a10t, s01 );
    FLA_Gemv( FLA_TRANSPOSE, FLA_ONE, A20, a21, FLA_ONE, s01 );
    FLA_Trmv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              S00, s01 );
    FLA_Scal( sigma11, s01 );
    FLA_Scal( FLA_MINUS_ONE, s01 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x3_to_2x2( &STL, /**/ &STR,       S00,  s01,     /**/ S02,
                                                     s10t, sigma11, /**/ s12t,
                            /* ************** */  /* ************************ */
                              &SBL, /**/ &SBR,       S20,  s21,     /**/ S22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &tT,                t0,
                                                  tau1,
                            /* ** */           /* **** */
                              &tB,                t2,     FLA_TOP );
  }
#endif

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Compute_tdqr_S_unb_var101( FLA_Obj D, FLA_Obj t, 
    FLA_Obj S ) {
// Compute factor S of QR factorization of [ U; D ], where U is an upper
// triangular matrix, and D is a dense matrix.
  int    m_D, n_D, ldim_D, ldim_S, i_one = 1, j;
  double  * buff_D, * buff_t, * buff_S, s_one = 1.0,
         s_zero = 0.0, s_minus_one = -1.0;

  //// printf( "      Beginning FLA_Compute_tdqr_S_unb_var101: \n" );
  //// FLA_Obj_show( " Di = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( " ti = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Si = [ ", S, "%le", " ];" );

  // Some initializations.
  m_D    = FLA_Obj_length( D );
  n_D    = FLA_Obj_width( D );
  buff_D = ( double * ) FLA_Obj_buffer_at_view( D );
  ldim_D = FLA_Obj_col_stride( D );

  buff_t = ( double * ) FLA_Obj_buffer_at_view( t );

  buff_S = ( double * ) FLA_Obj_buffer_at_view( S );
  ldim_S = FLA_Obj_col_stride( S );

  for( j = 0; j < n_D; j++ ) {
    // Update current column of S:
    //   sigma11 := tau1;
    //   s01     := - tau1 * triu( S00 ) * D0' * d1;

    // FLA_Copy( tau1, sigma11 );
    buff_S[ j * ldim_S + j ] = buff_t[ j ];

    // FLA_Gemv( FLA_TRANSPOSE, FLA_ONE, D0, d1, FLA_ZERO, s01 );
    dgemv_( "Transpose", & m_D, & j,
            & s_one, & buff_D[ 0 * ldim_D + 0 ], & ldim_D,
                     & buff_D[ j * ldim_D + 0 ], & i_one,
            & s_zero, & buff_S[ j * ldim_S + 0 ], & i_one );

    // FLA_Trmv( FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
    //           S00, s01 );
    dtrmv_( "Upper", "No transpose", "Non-unit", & j,
            &( buff_S[ 0 * ldim_S + 0 ] ), & ldim_S,
            &( buff_S[ j * ldim_S + 0 ] ), & i_one );

    // FLA_Scal( sigma11, s01 );
    dscal_( & j, &( buff_S[ j * ldim_S + j ] ),
                 &( buff_S[ j * ldim_S + 0 ] ), & i_one );

    // FLA_Scal( FLA_MINUS_ONE, s01 );
    dscal_( & j, & s_minus_one, &( buff_S[ j * ldim_S + 0 ] ), & i_one );
  }

  //// FLA_Obj_show( " Df = [ ", D, "%le", " ];" );
  //// FLA_Obj_show( " tf = [ ", t, "%le", " ];" );
  //// FLA_Obj_show( " Sf = [ ", S, "%le", " ];" );
  //// printf( "      End of FLA_Compute_tdqr_S_unb_var101: \n" );

  return FLA_SUCCESS;
}

