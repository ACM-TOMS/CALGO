function c = nck(n, k)

%NCK Binomial coefficient or all combinations.
%   NCK(N,K) where N and K are non-negative integers returns N!/K!(N-K)!.
%   This is the number of combinations of N things taken K at a time.
%   When a coefficient is large, the result is only accurate  to 15 digits
%   for double-precision inputs, or 8 digits for single-precision inputs.
%
%   Class support for inputs N,K,V:
%      float: double, single
%
%   See also PERMS.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.21.4.7 $  $Date: 2006/10/02 16:33:01 $

if k > n/2,
    k = n-k;
end

if k <= 1,
    c = n^k;
else
    nums = (n-k+1):n;
    dens = 1:k;
    nums = nums./dens;
    c = round(prod(nums));
end