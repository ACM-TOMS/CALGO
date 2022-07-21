#include "FLAME.h"


FLA_Error FLA_Compute_dense_QR_UT_var102( FLA_Obj A, FLA_Obj t, FLA_Obj S,
    int nb_alg );

FLA_Error FLA_Apply_left_Qt_of_dense_QR_UT_var102( FLA_Obj U, FLA_Obj S, 
    FLA_Obj C, int nb_alg );

FLA_Error FLA_Apply_right_Q_of_dense_QR_UT_var102( FLA_Obj U, FLA_Obj S,
    FLA_Obj C, int nb_alg );


FLA_Error FLA_Compute_td_QR_UT_var102( FLA_Obj U, FLA_Obj D, FLA_Obj t, 
    FLA_Obj S, int nb_alg );

FLA_Error FLA_Apply_left_Qt_of_td_QR_UT_var102( FLA_Obj D, FLA_Obj S, 
    FLA_Obj F, FLA_Obj G, int nb_alg );

FLA_Error FLA_Apply_right_Q_of_td_QR_UT_var102( FLA_Obj D, FLA_Obj S,
    FLA_Obj F, FLA_Obj G, int nb_alg );

