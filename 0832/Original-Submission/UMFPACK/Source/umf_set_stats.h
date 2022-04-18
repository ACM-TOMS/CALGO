/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.3 (Jan. 16, 2004), Copyright (c) 2004 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL void UMF_set_stats
(
    double Info [ ],
    SymbolicType *Symbolic,
    double max_usage,
    double num_mem_size,
    double flops,
    double lnz,
    double unz,
    double maxfrsize,
    double ulen,
    double npiv,
    double maxnrows,
    double maxncols,
    Int scale,
    Int prefer_diagonal,
    Int what
) ;
