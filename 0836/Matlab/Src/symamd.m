function [p, stats] = symamd (S, knobs)
%SYMAMD Symmetric approximate minimum degree permutation.
%    P = SYMAMD (S) for a symmetric positive definite matrix S, returns the
%    permutation vector p such that S(p,p) tends to have a sparser Cholesky
%    factor than S.  Sometimes SYMAMD works well for symmetric indefinite
%    matrices too.  SYMAMD tends to be faster than SYMMMD and tends to return
%    a better ordering.  The matrix S is assumed to be symmetric; only the
%    strictly lower triangular part is referenced.   S must be square.
%
%    See also COLMMD, COLPERM, SPPARMS, COLAMD, SYMMMD, SYMRCM.
%
%    Usage:  P = symamd (S)
%            P = symamd (S, knobs)
%            [P, stats] = symamd (S)
%            [P, stats] = symamd (S, knobs)
%
%    knobs is an optional input argument.  If S is n-by-n, then rows and
%    columns with more than knobs(1)*n entries are removed prior to ordering,
%    and ordered last in the output permutation P. If the knobs parameter is not
%    present, then the default of 0.5 is used instead.  knobs (2) controls the
%    printing of statistics and error messages.
%
%    stats is an optional 20-element output vector that provides data about the
%    ordering and the validity of the input matrix S.  Ordering statistics are
%    in stats (1:3).  stats (1) = stats (2) is the number of dense or empty
%    rows and columns ignored by SYMAMD and stats (3) is the number of
%    garbage collections performed on the internal data structure used by
%    SYMAMD (roughly of size 8.4*nnz(tril(S,-1)) + 9*n integers).
%
%    MATLAB built-in functions are intended to generate valid sparse matrices,
%    with no duplicate entries, with ascending row indices of the nonzeros
%    in each column, with a non-negative number of entries in each column (!)
%    and so on.  If a matrix is invalid, then SYMAMD may or may not be able
%    to continue.  If there are duplicate entries (a row index appears two or
%    more times in the same column) or if the row indices in a column are out
%    of order, then SYMAMD can correct these errors by ignoring the duplicate
%    entries and sorting each column of its internal copy of the matrix S (the
%    input matrix S is not repaired, however).  If a matrix is invalid in other
%    ways then SYMAMD cannot continue, an error message is printed, and no
%    output arguments (P or stats) are returned.  SYMAMD is thus a simple way
%    to check a sparse matrix to see if it's valid.
%
%    stats (4:7) provide information if SYMAMD was able to continue.  The
%    matrix is OK if stats (4) is zero, or 1 if invalid.  stats (5) is the
%    rightmost column index that is unsorted or contains duplicate entries,
%    or zero if no such column exists.  stats (6) is the last seen duplicate
%    or out-of-order row index in the column index given by stats (5), or zero
%    if no such row index exists.  stats (7) is the number of duplicate or
%    out-of-order row indices.
%
%    stats (8:20) is always zero in the current version of SYMAMD (reserved
%    for future use).
%
%    The ordering is followed by a column elimination tree post-ordering.
%
%    Authors:
%
%	The authors of the code itself are Stefan I. Larimore and Timothy A.
%	Davis (davis@cise.ufl.edu), University of Florida.  The algorithm was
%	developed in collaboration with John Gilbert, Xerox PARC, and Esmond
%	Ng, Oak Ridge National Laboratory.
%
%    Date:
%
%	September 8, 2003.  Version 2.3.
%
%    Acknowledgements:
%
%	This work was supported by the National Science Foundation, under
%	grants DMS-9504974 and DMS-9803599.

%    Notice:
%	Copyright (c) 1998-2003 by the University of Florida.
%	All Rights Reserved.
%
%	See http://www.cise.ufl.edu/research/sparse/colamd (the colamd.c
%	file) for the License.
%
%    Availability:
%
%	The colamd/symamd library is available at
%
%	    http://www.cise.ufl.edu/research/sparse/colamd/
%

%-------------------------------------------------------------------------------
% perform the symamd ordering:
%-------------------------------------------------------------------------------

if (nargout <= 1 & nargin == 1)
    p = symamdmex (S) ;
elseif (nargout <= 1 & nargin == 2)
    p = symamdmex (S, knobs) ;
elseif (nargout == 2 & nargin == 1)
    [p, stats] = symamdmex (S) ;
elseif (nargout == 2 & nargin == 2)
    [p, stats] = symamdmex (S, knobs) ;
else
    error ('symamd:  incorrect number of input and/or output arguments') ;
end

%-------------------------------------------------------------------------------
% symmetric elimination tree post-ordering:
%-------------------------------------------------------------------------------

[ignore, q] = sparsfun ('symetree', S (p,p)) ;
p = p (q) ;

