function t = full(t)
%FULL Convert a cp_tensor to a (dense) tensor.
%
%   T = FULL(C) converts a CP tensor to a (dense) tensor.
%
%   Examples
%      A = CP_TENSOR([3; 2], rand{4,2), rand(5,2), rand(3,2));
%      B = full(A) is a dense tensor.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also CP_TENSOR, TENSOR.

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


K = length(t.lambda);
M = order(t);
I = size(t);

for k = 1 : K           

    % Add in rank-1 matrix corresponding to
    % lambda(k)

    tmp = 1;
    for m = 1 : M
	tmp = tmp * t.u{m}(:,k)';
	tmp = reshape(tmp, prod(I(1:m)), 1);
    end

    if length(I) == 1
	tmpI = [I 1];
    else
	tmpI = I;
    end
    tmp = reshape(tmp, tmpI);

    if k == 1
	a = t.lambda(k) * tmp;
    else
	a = a + t.lambda(k) * tmp;
    end            
end

t = tensor(a, I);

return;
