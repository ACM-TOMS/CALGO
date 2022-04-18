function C = mtimes(A,B)
%MTIMES Implement A*B (scalar multiply) for cp_tensor.
%
%   C = mtimes(A,B) computes A * B where A is a CP tensor and B is
%   a scalar (or vice versa). The result C is the same size as A.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also CP_TENSOR.

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


% Note: We can do scalar times a tensor, but anything more complex is
% an error.

if isa(B,'numeric') && isequal(size(B),[1 1])
    C = cp_tensor(B * A.lambda, A.u);
elseif isa(A,'numeric') && isequal(size(A),[1 1])
    C = cp_tensor(A * B.lambda, B.u);
else
    error('Use mtimes(full(A),full(B)).');
end
