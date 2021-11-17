function C = times(A,B)
%TIMES Element-wise multiplication for tensors.
%
%    TIMES(A,B) denotes element-by-element multiplication.  A and B
%    must have the same dimensions unless one is a scalar.
%    A scalar can be multiplied into anything.
% 
%    C = TIMES(A,B) is called for the syntax 'A .* B' when A or B is
%    a tensor.
% 
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR.

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


A = tensor(A);
B = tensor(B);

if (numel(A.data) == 1) || (numel(B.data) == 1) 
    C = tensor(A.data * B.data);
    return;
end


if ~issamesize(A,B)
  error('Tensor order and size must agree.');
end

C = A.data .* B.data;
C = tensor(C, size(A));

