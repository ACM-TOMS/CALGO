function C = xor(A,B)
%XOR Logical EXCLUSIVE OR.
%
%   XOR(A,B) is the logical symmetric difference of elements A and B.
%   The result is one where either A or B, but not both, is nonzero.
%   The result is zero where A and B are both zero or nonzero.  A and
%   B must have the same dimensions (or one can be a scalar).
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

if ~( issamesize(A,B)  ||  (numel(A.data) == 1)  ||  (numel(B.data) == 1) )
    error('Tensor size mismatch.')
end

C = multiarrayop(@xor,A,B);
