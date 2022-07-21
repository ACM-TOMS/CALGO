#include "FLAME.h"


FLA_Error build_q_matrix( FLA_Obj U, FLA_Obj s, FLA_Obj C );

FLA_Error check_ort( FLA_Obj Q, FLA_Obj resid );

FLA_Error check_qr_with_q( FLA_Obj A, FLA_Obj Q, FLA_Obj R, FLA_Obj resid );

FLA_Error check_qr_with_qt( FLA_Obj A, FLA_Obj QT, FLA_Obj R, FLA_Obj resid );

FLA_Error check_qrcp_with_q( FLA_Obj A, FLA_Obj Q, FLA_Obj R, FLA_Obj p,
    FLA_Obj resid );

FLA_Error check_sv_of_upper_triang_matrix( FLA_Obj R, FLA_Obj vs,
    FLA_Obj resid );

FLA_Error check_sv_of_dense_matrix( FLA_Obj A, FLA_Obj vs, FLA_Obj resid );

FLA_Error check_svd_factorization( FLA_Obj A, FLA_Obj U, FLA_Obj S, FLA_Obj VT,
    FLA_Obj resid );

FLA_Error check_utv( FLA_Obj A, FLA_Obj U, FLA_Obj T, FLA_Obj V,
    FLA_Obj resid );

FLA_Error check_solution_of_system( FLA_Obj A, FLA_Obj X, FLA_Obj B,
    FLA_Obj resid );

