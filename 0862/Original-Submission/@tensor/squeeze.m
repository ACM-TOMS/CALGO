function B = squeeze(A)
%SQUEEZE Remove singleton dimensions from a tensor.
%
%   B = SQUEEZE(A) returns a tensor B with the same elements as
%   A but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(A,dim)==1.  
%
%   Examples
%      squeeze( tensor(rand(2,1,3)) ) is a 2-by-3 tensor
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

sizeA = size(A);
B = squeeze(A.data);
B = tensor(B, sizeA(sizeA > 1));
