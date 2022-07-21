/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLASH_QUEUE_MACRO_DEFS_EXTRA_H
#define FLASH_QUEUE_MACRO_DEFS_EXTRA_H


//#ifdef FLA_ENABLE_SUPERMATRIX
#include "CPU_Inner_stubs_utv.h"
#include "MyFLA_Parallel_normal_random.h"
#include "MyFLA_Utils.h"


#define FLASH_OBJ_PTR_ID( A )  ( A ).base->id

/*
Reminder to create a macro to enqueue when SuperMatrix is configured, and
also to create a macro for when it is not below to return an error code.
*/

// LAPACK-level

#define ENQUEUE_FLASH_LU_piv_macro( A, p, cntl ) \
        FLASH_Queue_push( (void *) FLA_LU_piv_macro_task, \
                          (void *) cntl, \
                          "LU   ", \
                          FALSE, \
                          0, 0, 0, 2, \
                          A, p )

#define ENQUEUE_FLASH_COMP_DENSE_QR_UT( B, T ) \
        FLASH_Queue_push( (void*)CPU_Compute_dense_QR_UT_inner_utv, \
        NULL, \
        "Compute_dense_QR_UT", \
        FALSE, \
        0, 0, 0, 2, \
	B, T )

#define ENQUEUE_FLASH_COMP_TD_QR_UT( U, D, S ) \
        FLASH_Queue_push( (void*)CPU_Compute_td_QR_UT_inner_utv, \
        NULL, \
        "Compute_td_QR", \
        FALSE, \
        0, 0, 0, 3, \
	U, D, S )

#define ENQUEUE_FLASH_KEEP_UPPER_TRIANG( A ) \
        FLASH_Queue_push( (void*)CPU_Keep_upper_triang_inner_utv, \
        NULL, \
        "Keep_upper_triang", \
        FALSE, \
        0, 0, 0, 1, \
	A )

#define ENQUEUE_FLASH_SET_TO_ZERO( A ) \
        FLASH_Queue_push( (void*)CPU_Set_to_zero_inner_utv, \
        NULL, \
        "Set_to_zero", \
        FALSE, \
        0, 0, 0, 1, \
	A )

#define ENQUEUE_FLASH_APPLY_LEFT_QT_OF_DENSE_QR( U, S, C ) \
        FLASH_Queue_push( (void*)CPU_Apply_left_Qt_of_dense_QR_UT_inner_utv, \
        NULL, \
        "Apply_left_qt_of_dense_qr", \
        FALSE, \
        0, 0, 2, 1, \
	U, S, C )

#define ENQUEUE_FLASH_APPLY_LEFT_QT_OF_TD_QR( D, S, F, G ) \
        FLASH_Queue_push( (void*)CPU_Apply_left_Qt_of_td_QR_UT_inner_utv, \
        NULL, \
        "Apply_left_qt_of_td_qr", \
        FALSE, \
        0, 0, 2, 2, \
	D, S, F, G )

#define ENQUEUE_FLASH_APPLY_RIGHT_Q_OF_DENSE_QR_UT( U, S, C ) \
        FLASH_Queue_push( (void*)CPU_Apply_right_Q_of_dense_QR_UT_inner_utv, \
        NULL, \
        "Apply_right_q_of_dense_qr_ut", \
        FALSE, \
        0, 0, 2, 1, \
	U, S, C )

#define ENQUEUE_FLASH_APPLY_RIGHT_Q_OF_TD_QR_UT( D, S, F, G ) \
        FLASH_Queue_push( (void*)CPU_Apply_right_Q_of_td_QR_UT_inner_utv, \
        NULL, \
        "Apply_right_q_of_td_qr_ut", \
        FALSE, \
        0, 0, 2, 2, \
	D, S, F, G )

#define ENQUEUE_FLASH_NORMAL_RANDOM_MATRIX( A ) \
        FLASH_Queue_push( (void*)MyFLA_Normal_random_matrix, \
        NULL, \
        "Normal_random_matrix", \
        FALSE, \
        0, 0, 0, 1, \
	A )

#define ENQUEUE_SVD_OF_BLOCK( U, A, VT ) \
        FLASH_Queue_push( (void*)CPU_Compute_svd_inner_utv, \
        NULL, \
        "SVD_of_block", \
        FALSE, \
        0, 0, 0, 3, \
	U, A, VT )

#define ENQUEUE_FLASH_GEMM_ABTA( A, B ) \
        FLASH_Queue_push( (void*)CPU_Gemm_abta_inner_utv, \
        NULL, \
        "GEMM_abta", \
        FALSE, \
        0, 0, 1, 1, \
	A, B )

#define ENQUEUE_FLASH_GEMM_NN_OZ( C, A, B ) \
        FLASH_Queue_push( (void*)FLA_GEMM, \
        NULL, \
        "GEMM_NN_OZ", \
        FALSE, \
        0, 0, 1, 2, \
	C, A, B )

#define ENQUEUE_FLASH_GEMM_AABT( A, B ) \
        FLASH_Queue_push( (void*)CPU_Gemm_aabt_inner_utv, \
        NULL, \
        "GEMM_aabt", \
        FALSE, \
        0, 0, 1, 1, \
	A, B )

#define ENQUEUE_FLASH_GEMM_AAB( A, B ) \
        FLASH_Queue_push( (void*)CPU_Gemm_aab_inner_utv, \
        NULL, \
        " GEMM_aab", \
        FALSE, \
        0, 0, 1, 1, \
	A, B )

#define ENQUEUE_FLASH_MYCOPY( A, B ) \
        FLASH_Queue_push( (void *) CPU_Mycopy_inner_utv, \
                          NULL, \
                          "MYCOPY   ", \
                          FALSE, \
                          0, 0, 1, 1, \
                          A, B  )


//#endif // FLA_ENABLE_SUPERMATRIX


#endif // FLASH_QUEUE_MACRO_DEFS_EXTRA_H
