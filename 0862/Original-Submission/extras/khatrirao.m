function P = khatrirao(varargin)
%KHATRIRAO Khatri-Rao product
%
%   KHATRIRAO(A,B) computes the Khatri-Rao product of matrices A and
%   B that have the same number of columns.  The result is the
%   column-wise Kronecker product
%   [KRON(A(:,1),B(:,1)) ... KRON(A(:,n),B(:,n))]
%
%   KHATRIRAO(A1,A2,...) computes the Khatri-Rao product of
%   multiple matrices that have the same number of columns.
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

%$Id: khatrirao.m,v 1.9 2006/04/04 23:55:52 tgkolda Exp $

if nargin==1
    A = varargin{1};
else
    A = varargin;
end

M = size(A{1},1);
N = size(A{1},2);

for i = 2 : length(A)
    M = M * size(A{i},1);
    if (N ~= size(A{i},2))
	error('All matrices must have the same number of columns.')
    end
end

P = zeros(M,N);
for n = 1:N
    ab = A{1}(:,n);
    for i = 2 : length(A)
       ab = A{i}(:,n) * ab(:).';  % Compute outer product of nth columns
    end
    P(:,n) = ab(:);          % Fill nth column of P with reshaped result
end

