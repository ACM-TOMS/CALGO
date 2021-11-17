function A = full(B)
%FULL Convert a tucker_tensor to a (dense) tensor.
%
%   A = FULL(B) converts Tucker tensor B to (dense) tensor A.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TUCKER_TENSOR, TENSOR.

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

M = order(B);
tmp = ttm(B.lambda, B.u{1}, 1);
for m = 2 : M
   tmp = ttm(tmp, B.u{m}, m);
end
A = tensor(tmp, size(B));
return;
