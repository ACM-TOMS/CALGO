#include "FLAME.h"


void time_utv( int variant, int num_threads, dim_t nb, int num_execs,
    int print_data, int check_result,
    FLA_Obj A, FLA_Obj Acopy, FLA_Obj sv,
    int build_ort_matrices, int q,
    double * dtime, double * gflops, double * res );

