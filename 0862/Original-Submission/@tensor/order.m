function n = order(t)
%ORDER Return the order of a tensor.
%
%   ORDER(T) returns the mathematical order of tensor T.  In most
%   cases, ORDER is equal to NDIMS, the number of dimensions of a
%   tensor.  The single difference is that the ORDER of a scalar is 0
%   and not 1.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, NDIMS.

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


n = length(t.size);
