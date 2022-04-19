function B = shiftdim(varargin)
%SHIFTDIM Shift dimensions for a tensor.
%
%   B = SHIFTDIM(X,N) shifts the dimensions of tensor X by N.  When N 
%   is positive, SHIFTDIM shifts the dimensions to the left and wraps 
%   the N leading dimensions to the end.  When N is negative, SHIFTDIM
%   shifts the dimensions to the right and pads with singletons.
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


A = varargin{1};
n = varargin{2};

asize = size(A);

if (n >= 0)
  m = mod(n-1,length(size(A))) + 1;
  bsize = [asize(m+1:end) asize(1:m)];
  n = m - sum(asize(1:m) == 1);  % Don't shiftdim for leading singleton dims
else
  bsize = [ones(1,-n) asize];
end

if (n >= 0) && (mod(n,length(asize)) == 0)
  B = A;
else
  B = feval(@shiftdim, A.data, n);
  B = reshape(B,bsize);
  B = tensor(B, bsize);
end

return;
