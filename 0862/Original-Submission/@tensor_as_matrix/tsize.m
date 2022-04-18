function sz = tsize(a,idx)
%TSIZE Tensor size of tensor_as_matrix.
%
%   D = SIZE(X) returns the size of the tensor being stored as a
%   matrix. 
% 
%   M = SIZE(X,DIM) returns the length of the dimension(s) specified
%   by DIM.  For example, SIZE(X,1) returns the size of the first
%   dimension of the tensor.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR_AS_MATRIX, TENSOR_AS_MATRIX/SIZE.

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
    sz = tsize(idx);
else
    sz = tsize;
end
