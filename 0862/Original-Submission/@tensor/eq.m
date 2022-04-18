function C = eq(A,B)
%EQ Equal for tensors.
%
%   A == B does element by element comparisons between A and B and
%   returns a matrix of the same size with elements set to one where
%   the relation is true and elements set to zero where it is not.  A
%   and B must have the same dimensions unless one is a scalar. A
%   scalar can be compared with anything.
% 
%   C = EQ(A,B) is called for the syntax 'A == B' when A or B is a
%   tensor.
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

C = multiarrayop(@eq,A,B);
