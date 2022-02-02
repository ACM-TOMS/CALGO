function t = cp_tensor(varargin)
%CP_TENSOR Tensor stored in CANDECOMP/PARAFAC form.
%
%   A = CP_TENSOR(lambda,U1,U2,...,UM) creates a CP tensor from its
%   constituent parts. Here lambda is a k-vector and each Um is a
%   matrix with k columns.
%
%   A = CP_TENSOR(lambda, U) is the same as above except that U is a
%   cell array containing matrix Um in cell m.
%
%   A = CP_TENSOR(U) assumes U is a cell array containt matrix Um
%   in cell m and assigns the weight of each factor to be one.
%
%   A = CP_TENSOR(T) creates a CP tensor by copying an existing CP
%   tensor.
%
%   Examples
%      A = CP_TENSOR([3; 2], rand{4,2), rand(5,2), rand(3,2))
%      creates a 4 x 5 x 2  CP tensor.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TUCKER_TENSOR, CP_TENSOR/FULL.

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


% Copy CONSTRUCTOR
if (nargin == 1) && isa(varargin{1}, 'cp_tensor')
    t.lambda = varargin{1}.lambda;
    t.u = varargin{1}.u;
    t = class(t, 'cp_tensor');
    return;
end

if isa(varargin{1},'cell')

    u = varargin{1};
    t.lambda = ones(size(u{1},2),1);
    t.u = u;
    
else

    t.lambda = varargin{1};
    if ~isa(t.lambda,'numeric') || ndims(t.lambda) ~=2 || size(t.lambda,2) ~= 1
	error('LAMBDA must be a column vector.');
    end
    
    if isa(varargin{2},'cell')
	t.u = varargin{2};
    else
	for i = 2 : nargin
	    t.u{i-1} = varargin{i};
	end
    end

end
    
    
% Check that each Um is indeed a matrix
for i = 1 : length(t.u)
    if ndims(t.u{i}) ~= 2
	error(['Matrix U' int2str(i) ' is not a matrix!']);
    end
end

% Size error checking			     
k = length(t.lambda); 
for i = 1 : length(t.u)            
    if  size(t.u{i},2) ~= k
       error(['Matrix U' int2str(i) ' does not have ' int2str(k) ' columns.']);
    end
end

t = class(t, 'cp_tensor');
return;
