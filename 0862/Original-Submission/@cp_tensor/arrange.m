function A = arrange(A)
%ARRANGE Arranges the rank-1 terms of a CP tensor.
%
%   ARRANGE(A) normalizes the columns of each matrix, absorbing the
%   excess weight into lambda and then sorts everything so that the
%   lambda values are in decreasing order.
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

% Ensure that matrices are normalized
for r = 1:length(A.lambda)
    for n = 1:ndims(A)
    tmp = norm(A.u{n}(:,r));
    A.u{n}(:,r) = A.u{n}(:,r) / tmp;
    A.lambda(r) = A.lambda(r) * tmp;
    end
end

% Sort
[A.lambda, idx] = sort(A.lambda, 1, 'descend');
for i = 1 : ndims(A)
    A.u{i} = A.u{i}(:,idx);
end

return;

% s.type = '.';
% s.subs = 'U';
% factors = subsref(A,s);
% seq = indx(end:-1:1);
% for i = 1:ndims(A)
%   tmp = factors{i};
%   factors{i} = tmp(:,seq);
% end

