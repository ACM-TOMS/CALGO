function C = mtimes(A,B)
%MTIMES Scalar or inner product for tensors.
%
%   C = MTIMES(A,B) is the inner product of two, same-sized tensors A
%   and B.
%
%   C = MTIMES(A,B) multiplies the tensor A by the scalar B (or
%   vice versa).
%
%   C = MTIMES(A,B) is called for the syntax 'A * B' when A or B is a
%   tensor.
% 
%   Examples
%      A = tensor(rand(3,4,2));
%      B = tensor(rand(3,4,2));
%      C = tensor(rand(4,3,2));
%      D = A * B; <-- Computes the inner product of tensors A & B.
%      D = A * C; <-- ERROR! (Tensors are not of the same size);
%      D = 5 * A;
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/TTT.

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

if (numel(B.data) == 1) 
    C = tensor(A.data * B.data, size(A));    
    return;
elseif (numel(A.data) == 1)
    C = tensor(A.data * B.data, size(B));    
    return;
end

if ~issamesize(A,B)
  error('Tensor dimensions must agree.');
end

C = ttt(A,B,1:ndims(A));
C = C.data; % extract the scalar

