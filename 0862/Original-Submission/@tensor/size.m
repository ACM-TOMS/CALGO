function m = size(t,idx)
%SIZE Size of tensor.
%  
%   D = SIZE(T) returns the size of the tensor.  Note that tensor
%   objects explicitly store trailing singleton dimensions, so they
%   will be reported.
%
%   I = size(T,DIM) returns the size of the dimension specified by
%   the scalar DIM.
%
%   Examples
%      T = tensor(rand(3,4,2,1),[3 4 2 1]);
%      size(T) = [ 3 4 2 1 ]; <-- Tensor explicitly tracks trailing
%                                 singleton dimesions
%      size(double(T)) = [3 4 2]]; <--- Multidimensional arrays do
%                                       not.
%      size(T,5) <--- ERROR!
%      size(double(T),5) = 1 <--- Okay because the number of
%                                 dimensions is not tracked
%                                 explicitly as with tensors. 
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/ORDER, TENSOR/NDIMS, SIZE.

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

    m = t.size(idx);

else

    m = t.size;

end
