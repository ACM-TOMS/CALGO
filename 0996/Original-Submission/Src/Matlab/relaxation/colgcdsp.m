%COLGCDSP   Compute the greatest common divisor of each column.
%
%    Usage: gcdvec = COLGCDSP(mat);
%
%    Input:
%       mat: matrix in sparse format.
%
%    Output:
%       gcd: row vector where the ii-th element is the gcd of ii-column.
%
%    Example:
%       >> mat = sparse([2 6 7 0; 4 9 9 0; 8 12 0 0]);
%       >> gcdvec = colgcdsp(mat)
%
%       gcdvec =
%
%              2     3     1     1