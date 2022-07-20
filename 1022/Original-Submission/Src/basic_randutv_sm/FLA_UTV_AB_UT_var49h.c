#include <stdlib.h>
#include <omp.h>
#include "FLAME.h"
#include "MyFLA_Utils.h"
#include "MyFLA_Normal_random.h"
#include "FLA_UTV_AB_UT_var49h.h"
#include "FLASH_Queue_macro_defs_extra.h"

#define NANO_BLOCK_SIZE 32

extern fla_gemm_t* fla_gemm_cntl_blas;

void create_block_matrices( 
  dim_t mb_A,     dim_t nb_A,     
  int build_u,
  dim_t mb_U,     dim_t nb_U,     
  int build_v,
  dim_t mb_V,     dim_t nb_V,     
  dim_t mb_G,     dim_t nb_G,     FLA_Obj G,     FLA_Obj * Gblk,
  dim_t mb_Y,     dim_t nb_Y,     FLA_Obj Y,     FLA_Obj * Yblk,
  dim_t mb_SA,    dim_t nb_SA,    FLA_Obj SA,    FLA_Obj * SAblk,
  dim_t mb_DA,    dim_t nb_DA,    FLA_Obj DA,    FLA_Obj * DAblk,
  dim_t mb_SY,    dim_t nb_SY,    FLA_Obj SY,    FLA_Obj * SYblk,
  dim_t mb_DY,    dim_t nb_DY,    FLA_Obj DY,    FLA_Obj * DYblk,
  dim_t mb_svdU,  dim_t nb_svdU,  FLA_Obj svdU,  FLA_Obj * svdUblk,
  dim_t mb_svdVT, dim_t nb_svdVT, FLA_Obj svdVT, FLA_Obj * svdVTblk ) {

  FLASH_Obj_create_hier_copy_of_flat_ext( G, 1, &mb_G, &nb_G, Gblk );
  FLASH_Obj_create_hier_copy_of_flat_ext( Y, 1, &mb_Y, &nb_Y, Yblk );

  FLASH_Obj_create_hier_copy_of_flat_ext( SA, 1, &mb_SA, &nb_SA, SAblk );
  FLASH_Obj_create_hier_copy_of_flat_ext( DA, 1, &mb_DA, &nb_DA, DAblk );

  FLASH_Obj_create_hier_copy_of_flat_ext( SY, 1, &mb_SY, &nb_SY, SYblk );
  FLASH_Obj_create_hier_copy_of_flat_ext( DY, 1, &mb_DY, &nb_DY, DYblk );

  FLASH_Obj_create_hier_copy_of_flat_ext( svdU, 1, &mb_svdU, &nb_svdU, svdUblk );
  FLASH_Obj_create_hier_copy_of_flat_ext( svdVT, 1, &mb_svdVT, &nb_svdVT, svdVTblk );

}

void free_block_matrices(
  int build_u,
  int build_v,
  FLA_Obj * Gblk, FLA_Obj * Yblk,
  FLA_Obj * SAblk, FLA_Obj * DAblk,
  FLA_Obj * SYblk, FLA_Obj * DYblk,
  FLA_Obj * svdUblk, FLA_Obj * svdVTblk ) {

  FLASH_Obj_free( Gblk );
  FLASH_Obj_free( Yblk );

  FLASH_Obj_free( SAblk );
  FLASH_Obj_free( DAblk );

  FLASH_Obj_free( SYblk );
  FLASH_Obj_free( DYblk );

  FLASH_Obj_free( svdUblk );
  FLASH_Obj_free( svdVTblk );

}


int FLAX_Get_nanoblock_size_var49h( void ) {
  // Return nanoblock size.
  return( NANO_BLOCK_SIZE );
}

// ============================================================================
// Declaration of local prototypes.

static FLA_Error FLA_Compute_QR_of_column_panel( FLA_Obj B, FLA_Obj W,
    FLA_Obj S );

static FLA_Error FLA_Apply_left_Qt_of_QR_of_column_panel( FLA_Obj B,
    FLA_Obj W, FLA_Obj S, FLA_Obj C );

static FLA_Error FLA_Apply_left_Qt_of_dense_QR_UT_to_row_panel( FLA_Obj B,
    FLA_Obj S, FLA_Obj C );

static FLA_Error FLA_Apply_left_Qt_of_td_QR_UT_to_row_panels( FLA_Obj D,
    FLA_Obj S, FLA_Obj F, FLA_Obj G );

static FLA_Error FLA_Apply_right_Q_of_QR_of_column_panel(
    FLA_Obj B, FLA_Obj W, FLA_Obj S, FLA_Obj C );

static FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_to_column_panel(
    FLA_Obj B, FLA_Obj S, FLA_Obj C );

static FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_to_column_panels( FLA_Obj D,
    FLA_Obj S, FLA_Obj F, FLA_Obj G );

static FLA_Error FLA_Keep_upper_triangular_part_of_column_panel( FLA_Obj B );
static FLA_Error FLA_Set_column_panel_to_zero( FLA_Obj B );

static FLA_Error FLA_Generate_normal_column_panel( FLA_Obj A );

static FLA_Error MyFLA_Gemm_tn_mp_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C );
static FLA_Error MyFLA_Gemm_tn_pb_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C );
static FLA_Error MyFLA_Gemm_tn_pb_oo( FLA_Obj A, FLA_Obj B, FLA_Obj C );

static FLA_Error MyFLA_Gemm_nn_mp_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C );
static FLA_Error MyFLA_Gemm_nn_pb_oo( FLA_Obj A, FLA_Obj B, FLA_Obj C );
static FLA_Error MyFLA_Gemm_nn_pb_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C );

static FLA_Error MyFLA_Compute_svd( FLA_Obj U, FLA_Obj A, FLA_Obj VT );

static FLA_Error MyFLA_Gemm_abta( FLA_Obj A, FLA_Obj B );
static FLA_Error MyFLA_Gemm_aabt( FLA_Obj A, FLA_Obj B );
static FLA_Error MyFLA_Gemm_aab( FLA_Obj A, FLA_Obj B );


// ============================================================================
FLA_Error FLA_UTV_AB_UT_var49h( int m_A, int n_A, FLA_Obj Ablk,
    int build_u, FLA_Obj Ublk, int build_v, FLA_Obj Vblk,
    dim_t nb_alg, int n_iter, int nt ) {
//
  FLA_Obj  ATL, ATR,    A00, A01, A02,
           ABL, ABR,    A10, A11, A12,
                        A20, A21, A22;
  FLA_Obj  AR;

  FLA_Obj  UL, UR,      U0,  U1,  U2;
  FLA_Obj  VL, VR,      V0,  V1,  V2;

  FLA_Obj  STL, STR,    S00, S01, S02,
           SBL, SBR,    S10, S11, S12,
                        S20, S21, S22;
  FLA_Obj  DL,  DR,     D0,  D1,  D2;

  FLA_Obj  ABR_l, ABR_r, SAblk, DAblk,
           SBR_l, none;
  dim_t    b, nano_nb, m_SA, n_SA, m_DA, n_DA, prev_num_threads, j;
  FLA_Obj  GT,          G0,
           GB,          G1,
                        G2;
  FLA_Obj  YT,          Y0,
           YB,          Y1,
                        Y2;

  dim_t    m_G, n_G, m_Y, n_Y;
  FLA_Obj  G, Y;
  FLA_Obj  Gblk, Yblk;

  FLA_Obj  LT,          L0,
           LB,          L1,
                        L2;
  FLA_Obj  ET,          E0,
           EB,          E1,
                        E2;
  dim_t    m_SY, n_SY, m_DY, n_DY;

  FLA_Obj  SYblk, DYblk;
  FLA_Obj  svdUblk, svdVTblk ;

  dim_t    m_svdU, n_svdU, m_svdVT, n_svdVT ;

  FLASH_Queue_enable();
  FLASH_Queue_set_verbose_output( FLASH_QUEUE_VERBOSE_NONE );
  FLASH_Queue_set_num_threads( nt );

  FLASH_Queue_begin( );

  // Set seed so that the sequence of random numbers always starts with the
  // same number.
  MyFLA_Reset_normal_random_generator( 123 );

  // Create matrix of S factors: SA.
  nano_nb = FLAX_Get_nanoblock_size_var49h();
  m_SA = nano_nb * ( ( m_A + nb_alg - 1 ) / nb_alg );
  n_SA = n_A;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_SA, n_SA, 1, &nano_nb, &nb_alg, & SAblk );
  FLASH_Set( FLA_ONE, SAblk );

  // Create matrix of diagonal blocks of A: DA.
  m_DA = nb_alg;
  n_DA = n_A;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_DA, n_DA, 1, &nb_alg, &nb_alg, & DAblk );
  FLASH_Set( FLA_ONE, DAblk );

  // Create matrix G.
  m_G = m_A;
  n_G = nb_alg;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_G, n_G, 1, &nb_alg, &nb_alg, & Gblk );
  FLASH_Set( FLA_ONE, Gblk );

  // Create matrix Y.
  m_Y = n_A;
  n_Y = nb_alg;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_Y, n_Y, 1, &nb_alg, &nb_alg, & Yblk );
  FLASH_Set( FLA_ONE, Yblk );

  // Create matrix of S factors for Y: SY.
  nano_nb = FLAX_Get_nanoblock_size_var49h();
  m_SY = nano_nb * ( ( n_A + nb_alg - 1 ) / nb_alg );
  n_SY = nb_alg;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_SY, n_SY, 1, &nano_nb, &nb_alg, & SYblk );
  FLASH_Set( FLA_ONE, SYblk );

  // Create matrix of diagonal blocks of Y: DY.
  m_DY = n_A;
  n_DY = nb_alg;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_DY, n_DY, 1, &nb_alg, &nb_alg, & DYblk );
  FLASH_Set( FLA_ONE, DYblk );

  // Create matrix for storing matrix U of the svd: svdU.
  m_svdU = nb_alg;
  n_svdU = nb_alg;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_svdU, n_svdU, 1, &nb_alg, &nb_alg, & svdUblk );
  FLASH_Set( FLA_ONE, svdUblk );

  // Create matrix for storing matrix VT of the svd: svdVT.
  m_svdVT = nb_alg;
  n_svdVT = nb_alg;

  FLASH_Obj_create_ext( FLASH_Obj_datatype( Ablk ), m_svdVT, n_svdVT, 1, &nb_alg, &nb_alg, & svdVTblk );
  FLASH_Set( FLA_ONE, svdVTblk );

  // The usual Flame algorithm starts here.

  FLA_Part_2x2( Ablk,   &ATL, &ATR,
                        &ABL, &ABR,   0, 0, FLA_TL );

  if( build_u == 1 ) {
    FLA_Part_1x2( Ublk,   &UL,  &UR,    0, FLA_LEFT );
  }

  if( build_v == 1 ) {
    FLA_Part_1x2( Vblk,   &VL,  &VR,    0, FLA_LEFT );
  }

  FLA_Part_2x2( SAblk,  &STL, &STR,
                        &SBL, &SBR,   0, 0, FLA_TL );
  FLA_Part_1x2( DAblk,  &DL,  &DR,    0, FLA_LEFT );

  FLA_Part_2x1( Gblk,   &GT,
                        &GB,          0, FLA_TOP );
  FLA_Part_2x1( Yblk,   &YT,
                        &YB,          0, FLA_TOP );

  FLA_Part_2x1( SYblk,  &LT,
                        &LB,          0, FLA_TOP );

  FLA_Part_2x1( DYblk,  &ET,
                        &EB,          0, FLA_TOP );


  // Main loop.
  while ( ( FLA_Obj_width ( ATL ) < FLA_Obj_width ( Ablk ) ) &&
          ( FLA_Obj_length( ATL ) < FLA_Obj_length( Ablk ) ) ) {
    b = min( 1, min( FLA_Obj_length( ABR ), FLA_Obj_width( ABR ) ) );

    if ( 0 == b ) break;

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );
    if( build_u == 1 ) {
      FLA_Repart_1x2_to_1x3( UL,  /**/ UR,        &U0, /**/ &U1, &U2,
                             b, FLA_RIGHT );

    }
    if( build_v == 1 ) {
      FLA_Repart_1x2_to_1x3( VL,  /**/ VR,        &V0, /**/ &V1, &V2,
                             b, FLA_RIGHT );
    }
    FLA_Repart_2x2_to_3x3( STL, /**/ STR,       &S00, /**/ &S01, &S02,
                        /* ************* */   /* ******************** */
                                                &S10, /**/ &S11, &S12,
                           SBL, /**/ SBR,       &S20, /**/ &S21, &S22,
                           b, b, FLA_BR );
    FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                           b, FLA_RIGHT );

    FLA_Repart_2x1_to_3x1( GT,                &G0,
                        /* ** */            /* ** */
                                              &G1,
                           GB,                &G2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( YT,                &Y0,
                        /* ** */            /* ** */
                                              &Y1,
                           YB,                &Y2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( LT,                &L0,
                        /* ** */            /* ** */
                                              &L1,
                           LB,                &L2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( ET,                &E0,
                        /* ** */            /* ** */
                                              &E1,
                           EB,                &E2,        b, FLA_BOTTOM );
    /*-----------------------------------------------------------------------*/


    // Perform the following processing only if there are more iterations.
    if( FLA_Obj_width( A12 ) > 0 ) {

      //
      // Generate normal random matrix G.
      // ================================
      //

      FLA_Generate_normal_column_panel( GB );

      //
      // Compute the sampling matrix Y.
      // ==============================
      //

      MyFLA_Gemm_tn_mp_oz( ABR, GB, YB );

      // %%% Perform "power iteration" if requested.
      for( j = 0; j < n_iter; j++ ) {
        MyFLA_Gemm_nn_mp_oz( ABR, YB, GB );
        MyFLA_Gemm_tn_mp_oz( ABR, GB, YB );
      }

      //
      // Update A from the right side. Update V if asked.
      // ================================================
      //

      // Factorize column block of Y.
      FLA_Compute_QR_of_column_panel( YB, E1, LB );

      // Update matrix A from the right side
      // with transformations from Y, HB, and LB.
      FLA_Merge_2x1( ATR,
                     ABR, & AR );

      FLA_Apply_right_Q_of_QR_of_column_panel( YB, E1, LB, AR );

      // Update matrix V from the right side
      // with transformations from Y, HB, and LB.
      if( build_v == 1 ) {
        FLA_Apply_right_Q_of_QR_of_column_panel( YB, E1, LB, VR );
      }
    }

    // Perform the following processing only if there are rows in A21.
    if( FLA_Obj_length( A21 ) > 0 ) {

      //
      // Update A from the left side. Update U if asked.
      // ===============================================
      //

      // Some partitionings.
      FLA_Part_1x2( ABR,  & ABR_l, & ABR_r,    1, FLA_LEFT );
      FLA_Part_1x2( SBR,  & SBR_l, & none,     1, FLA_LEFT );

      // Factorize column panel of A:  ABR_l = [ A11; A21 ].
      FLA_Compute_QR_of_column_panel( ABR_l, D1, SBR_l );
    }

    // This operation has been moved upward to improve the parallelism.
    // Compute svd of diagonal block.
    if( FLA_Obj_length( A21 ) > 0 ) {
      FLA_Keep_upper_triangular_part_of_column_panel( A11 );
    }

    MyFLA_Compute_svd( svdUblk, A11, svdVTblk );

    if( FLA_Obj_length( A21 ) > 0 ) {
      // Update rest of matrix ABR_r with transformations from ABR_l and
      // SBR_l.
      FLA_Apply_left_Qt_of_QR_of_column_panel( ABR_l, D1, SBR_l, ABR_r );

      // Update matrix U from the right
      // with transformations from ABR_l, and SBR_l.
      if( build_u == 1 ) {
        FLA_Apply_right_Q_of_QR_of_column_panel( ABR_l, D1, SBR_l, UR );
      }

      FLA_Set_column_panel_to_zero( A21 );
    }

    //
    // Compute SVD of A11, and update matrices A, U, and V.
    // ====================================================
    //

    // Apply U of miniSVD to A.

    MyFLA_Gemm_abta( A12, svdUblk );

    // Apply VT of miniSVD to A.
    MyFLA_Gemm_aabt( A01, svdVTblk );

    // Apply U of miniSVD to U.
    if( build_u == 1 ) {
      MyFLA_Gemm_aab( U1, svdUblk );
    }

    // Apply VT of miniSVD to V.
    if( build_v == 1 ) {
      MyFLA_Gemm_aabt( V1, svdVTblk );
    }

    /*-----------------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    if( build_u == 1 ) {
      FLA_Cont_with_1x3_to_1x2( &UL,  /**/ &UR,        U0, U1, /**/ U2,
                                FLA_LEFT );
    }
    if( build_v == 1 ) {
      FLA_Cont_with_1x3_to_1x2( &VL,  /**/ &VR,        V0, V1, /**/ V2,
                                FLA_LEFT );
    }
    FLA_Cont_with_3x3_to_2x2( &STL, /**/ &STR,       S00, S01, /**/ S02,
                                                     S10, S11, /**/ S12,
                            /* ************** */  /* ****************** */
                              &SBL, /**/ &SBR,       S20, S21, /**/ S22,
                              FLA_TL );
    FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &GT,                G0,
                                                  G1,
                            /* ** */           /* ** */
                              &GB,                G2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &YT,                Y0,
                                                  Y1,
                            /* ** */           /* ** */
                              &YB,                Y2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &LT,                L0,
                                                  L1,
                            /* ** */           /* ** */
                              &LB,                L2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &ET,                E0,
                                                  E1,
                            /* ** */           /* ** */
                              &EB,                E2,     FLA_TOP );
  }

  FLASH_Queue_end( );
  
  // Free data structures.
  free_block_matrices( build_u, 
                       build_v, 
                       & Gblk,  & Yblk,
                       & SAblk, & DAblk,
                       & SYblk, & DYblk,
                       & svdUblk, & svdVTblk );

  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_Compute_QR_of_column_panel( FLA_Obj B, FLA_Obj W,
    FLA_Obj S ) {
// Factorize column block B, starting from top to bottom, block by block.
  FLA_Obj  BT,      B0,
           BB,      B1,
                    B2;
  FLA_Obj  ST,      S0,
           SB,      S1,
                    S2;
  int      b;
  FLA_Obj  Bt;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Initial partitioning.
  FLA_Part_2x1( B,    &BT,
                      &BB,            1, FLA_TOP );
  FLA_Part_2x1( S,    &ST,
                      &SB,            1, FLA_TOP );

  Bt = BT;

  // Enqueue task for QR factorization of dense BT.
  // ----------------------------------------------

  ENQUEUE_FLASH_COMP_DENSE_QR_UT( 
		  *FLASH_OBJ_PTR_AT( BT ),
		  *FLASH_OBJ_PTR_AT( ST )
		  );

  // Enqueue task for getting a copy of Householder vectors of previous QR.
  // ----------------------------------------------------------------------

  ENQUEUE_FLASH_MYCOPY(
		  *FLASH_OBJ_PTR_AT( BT ),
		  *FLASH_OBJ_PTR_AT( W )
		  );

  // Loop for factorization of blocks.
  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ) {

    b = min( FLA_Obj_length( BB ), 1 );
    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( ST,                &S0,
                        /* ** */            /* ** */
                                              &S1,
                           SB,                &S2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Enqueue task for QR factorization of [ Bt; B1 ],
    // where Bt is upper triangular, and B1 is dense.
    // -------------------------------------------------

    ENQUEUE_FLASH_COMP_TD_QR_UT(
		  *FLASH_OBJ_PTR_AT( Bt ),
		  *FLASH_OBJ_PTR_AT( B1 ),
		  *FLASH_OBJ_PTR_AT( S1 ) );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &ST,                S0,
                                                  S1,
                            /* ** */           /* ** */
                              &SB,                S2,     FLA_TOP );
  }


  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_right_Q_of_QR_of_column_panel( FLA_Obj B,
    FLA_Obj W, FLA_Obj S, FLA_Obj C ) {
// Update C applying transformations from the right,
// with transformations saved in B, T, and S.
  FLA_Obj  BT,     B0,
           BB,     B1,
                   B2;
  FLA_Obj  ST,     S0,
           SB,     S1,
                   S2;
  FLA_Obj  CL,     CR,       C0,  C1,  C2;
  dim_t    b;
  FLA_Obj  Ct;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 )||
      ( FLA_Obj_length( C ) == 0 )||( FLA_Obj_width( C ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Initial partitioning.
  FLA_Part_2x1( B,    &BT,
                      &BB,            1, FLA_TOP );
  FLA_Part_2x1( S,    &ST,
                      &SB,            1, FLA_TOP );
  FLA_Part_1x2( C,    &CL,  &CR,      1, FLA_LEFT );

  // Update CL with transformations from factorization of dense BB.
  FLA_Apply_right_Q_of_dense_QR_UT_to_column_panel( W, ST, CL );

  // Save object copies of the top row blocks.
  Ct = CL;

  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ){
    b = min( FLA_Obj_length( BB ), 1 );

    FLA_Repart_2x1_to_3x1( BT,                &B0,
                                              &B1,
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( ST,                &S0,
                                              &S1,
                        /* ** */            /* ** */
                           SB,                &S2,        b, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    // Update [ Ct; C1 ] with transformations from QR factorization of
    // [ I; B1 ], where Bt is upper triangular.
    FLA_Apply_right_Q_of_dense_QR_UT_to_column_panels( B1, S1, Ct, C1 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                            /* ** */           /* ** */
                                                  B1,
                              &BB,                B2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &ST,                S0,
                            /* ** */           /* ** */
                                                  S1,
                              &SB,                S2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_to_column_panel(
    FLA_Obj B, FLA_Obj S, FLA_Obj C ) {
// Update row block C with transformations from QR factorization of dense B.
  FLA_Obj  CT,              C0,
           CB,              C1,
                            C2;
  dim_t    b;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 )||
      ( FLA_Obj_length( C ) == 0 )||( FLA_Obj_width( C ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( CT ) < FLA_Obj_length( C ) ){

    b = min( FLA_Obj_length( CB ), 1 );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                        /* ** */            /* ** */
                                              &C1,
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    ENQUEUE_FLASH_APPLY_RIGHT_Q_OF_DENSE_QR_UT(
		    *FLASH_OBJ_PTR_AT(B), 
		    *FLASH_OBJ_PTR_AT(S),
		    *FLASH_OBJ_PTR_AT(C1) );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                  C1,
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_to_column_panels( FLA_Obj D,
    FLA_Obj S, FLA_Obj F, FLA_Obj G ) {
// Update row blocks [ F G ] from the right side
// with transformations from the TD QR factorization.
  FLA_Obj  FT,              F0,
           FB,              F1,
                            F2;
  FLA_Obj  GT,              G0,
           GB,              G1,
                            G2;
  dim_t    b;

  // Quick return.
  if( ( FLA_Obj_length( D ) == 0 )||( FLA_Obj_width( D ) == 0 )||
      ( FLA_Obj_length( F ) == 0 )||( FLA_Obj_width( F ) == 0 )||
      ( FLA_Obj_length( G ) == 0 )||( FLA_Obj_width( G ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  FLA_Part_2x1( F,    &FT,
                      &FB,            0, FLA_TOP );
  FLA_Part_2x1( G,    &GT,
                      &GB,            0, FLA_TOP );

  while ( FLA_Obj_length( FT ) < FLA_Obj_length( F ) ){
    b = min( FLA_Obj_length( FB ), 1 );

    FLA_Repart_2x1_to_3x1( FT,                &F0,
                        /* ** */            /* ** */
                                              &F1,
                           FB,                &F2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( GT,                &G0,
                        /* ** */            /* ** */
                                              &G1,
                           GB,                &G2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    ENQUEUE_FLASH_APPLY_RIGHT_Q_OF_TD_QR_UT(
		    *FLASH_OBJ_PTR_AT(D), 
		    *FLASH_OBJ_PTR_AT(S), 
		    *FLASH_OBJ_PTR_AT(F1), 
		    *FLASH_OBJ_PTR_AT(G1) );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &FT,                F0,
                                                  F1,
                            /* ** */           /* ** */
                              &FB,                F2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &GT,                G0,
                                                  G1,
                            /* ** */           /* ** */
                              &GB,                G2,     FLA_TOP );

  }
  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_Keep_upper_triangular_part_of_column_panel( FLA_Obj B ) {
// Keep upper triangular part of panel B. The rest is set to zero.
  FLA_Obj  BT,      B0,
           BB,      B1,
                    B2;
  dim_t    b;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Initial partitioning.
  FLA_Part_2x1( B,    &BT,
                      &BB,            1, FLA_TOP );

  // Enqueue task for zeroing strictly lower part of BT.
  // ---------------------------------------------------

  ENQUEUE_FLASH_KEEP_UPPER_TRIANG( *FLASH_OBJ_PTR_AT( BT ) );

  // Loop for factorization of blocks.
  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ) {
    b = min( FLA_Obj_length( BB ), 1 );

    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Enqueue task for zeroing B1.
    // ----------------------------

    ENQUEUE_FLASH_SET_TO_ZERO( *FLASH_OBJ_PTR_AT( B1 ) );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Apply_left_Qt_of_QR_of_column_panel( FLA_Obj B,
    FLA_Obj W, FLA_Obj S, FLA_Obj C ) {
// Update C, with transformations saved in B, and taus saved in T.
  FLA_Obj  BT,     B0,
           BB,     B1,
                   B2;
  FLA_Obj  ST,     S0,
           SB,     S1,
                   S2;
  FLA_Obj  CT,     C0,
           CB,     C1,
                   C2;
  dim_t    b;
  FLA_Obj  Ct;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 )||
      ( FLA_Obj_length( C ) == 0 )||( FLA_Obj_width( C ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Initial partitioning.
  FLA_Part_2x1( B,    &BT,
                      &BB,            1, FLA_TOP );
  FLA_Part_2x1( S,    &ST,
                      &SB,            1, FLA_TOP );
  FLA_Part_2x1( C,    &CT,
                      &CB,            1, FLA_TOP );

  // Update CB with transformations from factorization of dense BB.
  FLA_Apply_left_Qt_of_dense_QR_UT_to_row_panel( W, ST, CT );

  // Save object copies of the top row blocks.
  Ct = CT;

  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ){

    b = min( FLA_Obj_length( BB ), 1 );
    FLA_Repart_2x1_to_3x1( BT,                &B0,
                                              &B1,
                        /* ** */            /* ** */
                           BB,                &B2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( ST,                &S0,
                                              &S1,
                        /* ** */            /* ** */
                           SB,                &S2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                                              &C1,
                        /* ** */            /* ** */
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Update [ Ct; C1 ] with transformations from QR factorization of
    // [ I; B1 ], where Bt is upper triangular.
    FLA_Apply_left_Qt_of_td_QR_UT_to_row_panels( B1, S1, Ct, C1 );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                            /* ** */           /* ** */
                                                  B1,
                              &BB,                B2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &ST,                S0,
                            /* ** */           /* ** */
                                                  S1,
                              &SB,                S2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                            /* ** */           /* ** */
                                                  C1,
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_Apply_left_Qt_of_dense_QR_UT_to_row_panel( FLA_Obj B,
    FLA_Obj S, FLA_Obj C ) {
// Update row block C with transformations from QR factorization of dense B.
  FLA_Obj  CL,    CR,       C0,  C1,  C2;
  dim_t    b;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 )||
      ( FLA_Obj_length( C ) == 0 )||( FLA_Obj_width( C ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

  while ( FLA_Obj_width( CL ) < FLA_Obj_width( C ) ){

    b = min( FLA_Obj_width( CR ), 1 );
    FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    ENQUEUE_FLASH_APPLY_LEFT_QT_OF_DENSE_QR(
                                        *FLASH_OBJ_PTR_AT(B),
                                        *FLASH_OBJ_PTR_AT(S),
                                        *FLASH_OBJ_PTR_AT(C1) );

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error FLA_Apply_left_Qt_of_td_QR_UT_to_row_panels( FLA_Obj D,
    FLA_Obj S, FLA_Obj F, FLA_Obj G ) {
// Update row blocks [ F; G ] with transformations from QR factorization
// of matrix ([ I; V ], S).
  FLA_Obj  FL, FR,    F0,  F1,  F2;
  FLA_Obj  GL, GR,    G0,  G1,  G2;
  dim_t    b;

  // Quick return.
  if( ( FLA_Obj_length( D ) == 0 )||( FLA_Obj_width( D ) == 0 )||
      ( FLA_Obj_length( F ) == 0 )||( FLA_Obj_width( F ) == 0 )||
      ( FLA_Obj_length( G ) == 0 )||( FLA_Obj_width( G ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  FLA_Part_1x2( F,    &FL,  &FR,      0, FLA_LEFT );
  FLA_Part_1x2( G,    &GL,  &GR,      0, FLA_LEFT );

  while ( FLA_Obj_width( FL ) < FLA_Obj_width( F ) ){

    b = min( FLA_Obj_width( FR ), 1 );
    FLA_Repart_1x2_to_1x3( FL,  /**/ FR,        &F0, /**/ &F1, &F2,
                           b, FLA_RIGHT );
    FLA_Repart_1x2_to_1x3( GL,  /**/ GR,        &G0, /**/ &G1, &G2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    ENQUEUE_FLASH_APPLY_LEFT_QT_OF_TD_QR( 
		    *FLASH_OBJ_PTR_AT(D), 
		    *FLASH_OBJ_PTR_AT(S), 
		    *FLASH_OBJ_PTR_AT(F1), 
		    *FLASH_OBJ_PTR_AT(G1) );

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &FL,  /**/ &FR,        F0, F1, /**/ F2,
                              FLA_LEFT );
    FLA_Cont_with_1x3_to_1x2( &GL,  /**/ &GR,        G0, G1, /**/ G2,
                              FLA_LEFT );
  }
  return FLA_SUCCESS;
}




// ============================================================================
static FLA_Error FLA_Set_column_panel_to_zero( FLA_Obj B ) {
// Set column panel B to zero.
  FLA_Obj  BT,      B0,
           BB,      B1,
                    B2;
  dim_t    b;

  // Quick return.
  if( ( FLA_Obj_length( B ) == 0 )||( FLA_Obj_width( B ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Initial partitioning.
  FLA_Part_2x1( B,    &BT,
                      &BB,            0, FLA_TOP );

  // Loop for factorization of blocks.
  while ( FLA_Obj_length( BT ) < FLA_Obj_length( B ) ) {
    b = min( FLA_Obj_length( BB ), 1 );

    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Enqueue task for zeroing B1.
    // ----------------------------

    ENQUEUE_FLASH_SET_TO_ZERO(
                  *FLASH_OBJ_PTR_AT(B1)
		);

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error FLA_Generate_normal_column_panel( FLA_Obj A ) {
// Set column panel to a normal random matrix.
  FLA_Obj  AT,              A0,
           AB,              A1,
                            A2;
  dim_t    b;

  //FLASH_Obj_show( "  Aprerand = [ ", A, "%le", " ];" );

  // Quick return.
  if( ( FLA_Obj_length( A ) == 0 )||( FLA_Obj_width( A ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Check that input argument A is a column panel.
  if( FLA_Obj_width( A ) != 1 ) {
    fprintf( stderr, "+++ ERROR in FLA_Generate_normal_column_panel:  " );
    fprintf( stderr, "Input argument is not a column panel\n" );

    return -1;
  }

  FLA_Part_2x1( A,    &AT,
                &AB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){
    b = min( FLA_Obj_length( AB ), 1 );
    ////////////////////printf( " b = %d\n", b );
    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* ** */
                                              &A1,
                           AB,                &A2,        b /*era b*/, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Enqueue task for generating normal random matrix.
    // -------------------------------------------------

    ENQUEUE_FLASH_NORMAL_RANDOM_MATRIX( *FLASH_OBJ_PTR_AT(A1) );

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  A1,
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Compute_svd( FLA_Obj U, FLA_Obj A, FLA_Obj VT ) {

  // Enqueue task for computing the SVD of A.
  // ----------------------------------------

  ENQUEUE_SVD_OF_BLOCK( 
		  *FLASH_OBJ_PTR_AT( U ),
		  *FLASH_OBJ_PTR_AT( A ),
		  *FLASH_OBJ_PTR_AT( VT ) );

  return FLA_SUCCESS;
}



// ============================================================================
static FLA_Error MyFLA_Gemm_tn_mp_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C ) {
// Gemm: Transpose-No transpose, Matrix-Panel, One-Zero.
//
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;
  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;
  int     b;

  FLA_Part_2x1( A,    &AT,
                      &AB,            0, FLA_TOP );
  FLA_Part_2x1( B,    &BT,
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){
    b = min( FLA_Obj_length( AB ), 1 );

    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* ** */
                                              &A1,
                           AB,                &A2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    if( FLA_Obj_length( AT ) == 0 ) {
      // First iteration: Assign.
      MyFLA_Gemm_tn_pb_oz( A1, B1, C );
    } else {
      // Rest of iterations: Accumulate.
      MyFLA_Gemm_tn_pb_oo( A1, B1, C );
    }

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  A1,
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Gemm_tn_pb_oo( FLA_Obj A, FLA_Obj B, FLA_Obj C ) {
// Gemm: Transpose-No transpose, Panel-Block, One-One.
//
  FLA_Obj AL,    AR,       A0,  A1,  A2;
  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;
  int     b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );
  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){
    b = min( FLA_Obj_width( AR ), 1 );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                        /* ** */            /* ** */
                                              &C1,
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // C = C + A1' * B.

    // Enqueue task for C = C + A1 * B.
    // --------------------------------

    ENQUEUE_FLASH_Gemm( 
			FLA_TRANSPOSE,
			FLA_NO_TRANSPOSE,
			FLA_ONE,
		  	*FLASH_OBJ_PTR_AT( A1 ),
		  	*FLASH_OBJ_PTR_AT( B ),
		  	FLA_ONE,
		  	*FLASH_OBJ_PTR_AT( C1 ), fla_gemm_cntl_blas
		  );



    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                  C1,
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}



// ============================================================================
static FLA_Error MyFLA_Gemm_tn_pb_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C ) {
// Gemm: Transpose-No transpose, Panel-Block, One-Zero.
//
  FLA_Obj AL,    AR,       A0,  A1,  A2;
  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;
  int     b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );
  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){
    b = min( FLA_Obj_width( AR ), 1 );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                        /* ** */            /* ** */
                                              &C1,
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // C = A1' * B.

    // Enqueue task for C = A1 * B.
    // ----------------------------

    ENQUEUE_FLASH_Gemm( 
			FLA_TRANSPOSE,
			FLA_NO_TRANSPOSE,
			FLA_ONE,
		  	*FLASH_OBJ_PTR_AT( A1 ),
		  	*FLASH_OBJ_PTR_AT( B ),
		  	FLA_ZERO,
		  	*FLASH_OBJ_PTR_AT( C1 ), fla_gemm_cntl_blas
		  );


    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                  C1,
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Gemm_nn_mp_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C ) {
// Gemm: No transpose-No transpose, Matrix-Panel, One-Zero.
//
  FLA_Obj AL,    AR,       A0,  A1,  A2;
  FLA_Obj BT,              B0,
          BB,              B1,
                           B2;
  int     b;

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );
  FLA_Part_2x1( B,    &BT,
                      &BB,            0, FLA_TOP );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){
    b = min( FLA_Obj_width( AR ), 1 );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );
    FLA_Repart_2x1_to_3x1( BT,                &B0,
                        /* ** */            /* ** */
                                              &B1,
                           BB,                &B2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    if( FLA_Obj_width( AL ) == 0 ) {
      // First iteration: Assign.
      MyFLA_Gemm_nn_pb_oz( A1, B1, C );
    } else {
      // Rest of iterations: Accumulate.
      MyFLA_Gemm_nn_pb_oo( A1, B1, C );
    }

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
    FLA_Cont_with_3x1_to_2x1( &BT,                B0,
                                                  B1,
                            /* ** */           /* ** */
                              &BB,                B2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Gemm_nn_pb_oo( FLA_Obj A, FLA_Obj B, FLA_Obj C ) {
// Gemm: No transpose-No transpose, Panel-Block, One-One.
//
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;
  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;
  int     b;

  FLA_Part_2x1( A,    &AT,
                      &AB,            0, FLA_TOP );
  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){
    b = min( FLA_Obj_length( AB ), 1 );

    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* ** */
                                              &A1,
                           AB,                &A2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                        /* ** */            /* ** */
                                              &C1,
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // C = C + A1 * B.

    // Enqueue task for C = C + A1 * B.
    // --------------------------------

    ENQUEUE_FLASH_Gemm( 
			FLA_NO_TRANSPOSE,
			FLA_NO_TRANSPOSE,
			FLA_ONE,
		  	*FLASH_OBJ_PTR_AT( A1 ),
		  	*FLASH_OBJ_PTR_AT( B ),
		  	FLA_ONE,
		  	*FLASH_OBJ_PTR_AT( C1 ), fla_gemm_cntl_blas
		  );


    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  A1,
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                  C1,
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Gemm_nn_pb_oz( FLA_Obj A, FLA_Obj B, FLA_Obj C ) {
// Gemm: No transpose-No transpose, Panel-Block, One-Zero.
//
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;
  FLA_Obj CT,              C0,
          CB,              C1,
                           C2;
  int     b;

  FLA_Part_2x1( A,    &AT,
                      &AB,            0, FLA_TOP );
  FLA_Part_2x1( C,    &CT,
                      &CB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){
    b = min( FLA_Obj_length( AB ), 1 );

    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* ** */
                                              &A1,
                           AB,                &A2,        b, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( CT,                &C0,
                        /* ** */            /* ** */
                                              &C1,
                           CB,                &C2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // C = A1 * B.

    // Enqueue task for C = A1 * B.
    // ----------------------------

    ENQUEUE_FLASH_Gemm( 
			FLA_NO_TRANSPOSE,
			FLA_NO_TRANSPOSE,
			FLA_ONE,
		  	*FLASH_OBJ_PTR_AT( A1 ),
		  	*FLASH_OBJ_PTR_AT( B ),
		  	FLA_ZERO,
		  	*FLASH_OBJ_PTR_AT( C1 ), fla_gemm_cntl_blas
		  );


    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  A1,
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &CT,                C0,
                                                  C1,
                            /* ** */           /* ** */
                              &CB,                C2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

// ============================================================================
//
static FLA_Error MyFLA_Gemm_abta( FLA_Obj A, FLA_Obj B ) {
// Compute:  A := B' * A, where A is a row panel.
//
  FLA_Obj AL,    AR,       A0,  A1,  A2;
  dim_t   b;

  // Quick return.
  if( ( FLA_Obj_length( A ) == 0 )||( FLA_Obj_width( A ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Check that input argument A is a row panel.
  if( FLA_Obj_length( A ) != 1 ) {
    fprintf( stderr, "+++ ERROR in MyFLA_Gemm_abta:  " );
    fprintf( stderr, "Input argument is not a row panel\n" );

    return -1;
  }

  FLA_Part_1x2( A,    &AL,  &AR,      0, FLA_LEFT );

  while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) ){
    b = min( FLA_Obj_width( AR ), 1 );

    FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
                           b, FLA_RIGHT );
    /*------------------------------------------------------------*/

    // A1 = B' * A1.

    // Enqueue task for A1 = B' * A1.
    // ------------------------------

    ENQUEUE_FLASH_GEMM_ABTA(
		  *FLASH_OBJ_PTR_AT( B ),
		  *FLASH_OBJ_PTR_AT( A1 )
		  );

    /*------------------------------------------------------------*/
    FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
                              FLA_LEFT );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Gemm_aabt( FLA_Obj A, FLA_Obj B ) {
// Compute:  A := A * B', where A is a column panel.
//
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;
  dim_t   b;

  // Quick return.
  if( ( FLA_Obj_length( A ) == 0 )||( FLA_Obj_width( A ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Check that input argument A is a column panel.
  if( FLA_Obj_width( A ) != 1 ) {
    fprintf( stderr, "+++ ERROR in MyFLA_Gemm_aabt:  " );
    fprintf( stderr, "Input argument is not a column panel\n" );

    return -1;
  }

  FLA_Part_2x1( A,    &AT,
                      &AB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){
    b = min( FLA_Obj_length( AB ), 1 );

    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* ** */
                                              &A1,
                           AB,                &A2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Enqueue task for A = A * B'.
    // ----------------------------

    ENQUEUE_FLASH_GEMM_AABT(
		  *FLASH_OBJ_PTR_AT( A1 ),
		  *FLASH_OBJ_PTR_AT( B )
		  );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  A1,
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}


// ============================================================================
static FLA_Error MyFLA_Gemm_aab( FLA_Obj A, FLA_Obj B ) {
// Compute:  A := A * B, where A is a column panel.
//
  FLA_Obj AT,              A0,
          AB,              A1,
                           A2;
  dim_t   b;

  // Quick return.
  if( ( FLA_Obj_length( A ) == 0 )||( FLA_Obj_width( A ) == 0 ) ) {
    return FLA_SUCCESS;
  }

  // Check that input argument A is a column panel.
  if( FLA_Obj_width( A ) != 1 ) {
    fprintf( stderr, "+++ ERROR in MyFLA_Gemm_aab:  " );
    fprintf( stderr, "Input argument is not a column panel\n" );

    return -1;
  }

  FLA_Part_2x1( A,    &AT,
                      &AB,            0, FLA_TOP );

  while ( FLA_Obj_length( AT ) < FLA_Obj_length( A ) ){
    b = min( FLA_Obj_length( AB ), 1 );

    FLA_Repart_2x1_to_3x1( AT,                &A0,
                        /* ** */            /* ** */
                                              &A1,
                           AB,                &A2,        b, FLA_BOTTOM );
    /*------------------------------------------------------------*/

    // Enqueue task for A = A * B.
    // ---------------------------

    ENQUEUE_FLASH_GEMM_AAB(
		  *FLASH_OBJ_PTR_AT( A1 ),
		  *FLASH_OBJ_PTR_AT( B )
		  );

    /*------------------------------------------------------------*/
    FLA_Cont_with_3x1_to_2x1( &AT,                A0,
                                                  A1,
                            /* ** */           /* ** */
                              &AB,                A2,     FLA_TOP );
  }

  return FLA_SUCCESS;
}

