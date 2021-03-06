/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.3 (Jan. 16, 2004), Copyright (c) 2004 by Timothy A.      */
/* Davis.  All Rights Reserved.  See ../README for License.                   */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

GLOBAL Int UMF_analyze
(
    Int n_row,
    Int n_col,
    Int Ai [ ],
    Int Ap [ ],
    Int Up [ ],
    Int fixQ,
    Int Front_ncols [ ],
    Int W [ ],
    Int Link [ ],
    Int Front_nrows [ ],
    Int Front_npivcol [ ],
    Int Front_parent [ ],
    Int *nfr_out,
    Int *p_ncompactions
) ;
