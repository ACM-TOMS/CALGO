#include "FLAME.h"

FLA_Error CPU_Compute_dense_QR_WY_inner_qr( FLA_Obj A, FLA_Obj t, FLA_Obj S,
    int nb_alg );

FLA_Error CPU_Compute_td_QR_WY_inner_qr( FLA_Obj U, FLA_Obj D, FLA_Obj t,
    FLA_Obj S, int nb_alg );

FLA_Error CPU_Apply_left_Qt_of_dense_QR_WY_inner_qr( FLA_Obj U, FLA_Obj S, 
    FLA_Obj C, int nb_alg );

FLA_Error CPU_Apply_left_Qt_of_td_QR_WY_inner_qr( FLA_Obj D, FLA_Obj S, 
    FLA_Obj F, FLA_Obj G, int nb_alg );

FLA_Error CPU_Mycopy_inner_qr( FLA_Obj A, FLA_Obj B );

