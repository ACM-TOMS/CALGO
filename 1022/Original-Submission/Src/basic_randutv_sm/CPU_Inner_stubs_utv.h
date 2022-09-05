#include "FLAME.h"

FLA_Error CPU_Compute_dense_QR_WY_inner_utv( FLA_Obj A, FLA_Obj S, dim_t nb_alg );

FLA_Error CPU_Apply_left_Qt_of_dense_QR_WY_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, dim_t nb_alg );

FLA_Error CPU_Apply_right_Q_of_dense_QR_WY_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, dim_t nb_alg );

FLA_Error CPU_Compute_td_QR_WY_inner_utv( FLA_Obj U, FLA_Obj D, FLA_Obj S,
    dim_t nb_alg );

FLA_Error CPU_Apply_left_Qt_of_td_QR_WY_inner_utv( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, dim_t nb_alg );

FLA_Error CPU_Apply_right_Q_of_td_QR_WY_inner_utv( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, dim_t nb_alg );


FLA_Error CPU_Compute_dense_QR_UT_inner_utv( FLA_Obj A, FLA_Obj S, dim_t nb_alg );

FLA_Error CPU_Apply_left_Qt_of_dense_QR_UT_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, dim_t nb_alg );

FLA_Error CPU_Apply_right_Q_of_dense_QR_UT_inner_utv( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, dim_t nb_alg );

FLA_Error CPU_Compute_td_QR_UT_inner_utv( FLA_Obj U, FLA_Obj D, FLA_Obj S,
    dim_t nb_alg );

FLA_Error CPU_Apply_left_Qt_of_td_QR_UT_inner_utv( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, dim_t nb_alg );

FLA_Error CPU_Apply_right_Q_of_td_QR_UT_inner_utv( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, dim_t nb_alg );


FLA_Error CPU_Keep_upper_triang_inner_utv( FLA_Obj A );

FLA_Error CPU_Set_to_zero_inner_utv( FLA_Obj A );

FLA_Error CPU_Gemm_tn_inner_utv( FLA_Obj alpha, FLA_Obj A, FLA_Obj B,
    FLA_Obj beta, FLA_Obj C );

FLA_Error CPU_Gemm_nn_inner_utv( FLA_Obj alpha, FLA_Obj A, FLA_Obj B,
    FLA_Obj beta, FLA_Obj C );

FLA_Error CPU_Compute_svd_inner_utv( FLA_Obj U, FLA_Obj A, FLA_Obj VT );

FLA_Error CPU_Gemm_abta_inner_utv( FLA_Obj A, FLA_Obj B );

FLA_Error CPU_Gemm_aabt_inner_utv( FLA_Obj A, FLA_Obj B );

FLA_Error CPU_Gemm_aab_inner_utv( FLA_Obj A, FLA_Obj B );

FLA_Error CPU_Mycopy_inner_utv( FLA_Obj A, FLA_Obj B );

