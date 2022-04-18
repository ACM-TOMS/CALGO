%ENUMCOMBWITHREP   Enumerate combinations with repetition.
%
%   This mex function returns two dimensional array
%   which enumerates k-combinations of {0,1,...,n-1} with repetition.
%   Each column vector represents a k-combination.
%
%   Usage:
%      comb = ENUMCOMBWITHREP(n,k);
%
%   Example:
%      comb = ENUMCOMBWITHREP(4,2)
%
%      comb =
%
%           0     0     0     0     1     1     1     2     2     3
%           0     1     2     3     1     2     3     2     3     3
