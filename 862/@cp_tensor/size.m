function m = size(t,idx)
%SIZE Size of cp_tensor.
%  
%   D = SIZE(T) returns the size of the tensor. 
%
%   I = SIZE(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also CP_TENSOR, CP_TENSOR/ORDER, CP_TENSOR/NDIMS.

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

if exist('idx','var')

    m = size(t.u{idx}, 1);

else

    for i = 1 : order(t)
	m(i) = size(t.u{i}, 1);
    end

end
