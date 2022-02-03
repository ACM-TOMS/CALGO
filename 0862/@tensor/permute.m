function T = permute(T,order)
%PERMUTE Permute tensor dimensions.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The result has the
%   same values of A, but the order of the subscripts needed to access
%   any particular element are rearranged as specified by ORDER.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/SIZE, TENSOR/NDIMS, PERMUTE.

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

N = ndims(T);

if length(order) == 1
  if (order == 1)
    return;
  else
    error('Invalid Order.');
  end
end

if N ~= length(order)
  error('ORDER must have at least N elements.')
elseif ~isempty(setxor(1:N,order))
  error('ORDER cannot contain repeated or invalid permutation indices.');
end

T.data = permute(T.data,order);
T.size = T.size(order);

return;
