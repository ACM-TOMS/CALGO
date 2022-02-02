function sz = size(a,idx)
%SIZE Size of tensor_as_matrix.
%
%   D = SIZE(X) returns the two-element row vector D = [M, N]
%   containing the number of rows and columns in the matrix.
% 
%   M = SIZE(X,DIM) returns the length of the dimension specified by
%   the scalar DIM.  For example, SIZE(X,1) returns the number of
%   rows.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR_AS_MATRIX, TENSOR_AS_MATRIX/TSIZE.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.

if exist('idx', 'var')
    sz = size(a.data, idx);
else
    sz = size(a.data);
end
