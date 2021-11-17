function t = tucker_tensor(varargin)
%TUCKER_TENSOR Tensor stored in Tucker form.
%
%   TUCKER_TENSOR(T) creates a TUCKER tensor by copying an existing
%   TUCKER tensor.
%
%   TUCKER_TENSOR(lambda,U1,U2,...,UM) creates a TUCKER tensor from
%   its constituent parts. Here lambda is a dense tensor of size
%   K1 x K2 x ... x KM and each Um is a matrix with Km columns. 
%
%   TUCKER_TENSOR(lambda, U) is the same as above except that U is a
%   cell array containing matrix Um in cell m.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, CP_TENSOR

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
if (nargin == 1) && isa(varargin{1}, 'tucker_tensor')
    t.lambda = varargin{1}.lambda;
    t.u = varargin{1}.u;
    t = class(t, 'tucker_tensor');
    return;
end

t.lambda = varargin{1};
if ~isa(t.lambda,'tensor')
    error('LAMBDA must be a tensor.');
end

if isa(varargin{2},'cell')
    t.u = varargin{2};
else
    for i = 2 : nargin
	t.u{i-1} = varargin{i};
    end
end

% Check that each Um is indeed a matrix
for i = 1 : length(t.u)
    if ndims(t.u{i}) ~= 2
	error(['Matrix U' int2str(i) ' is not a matrix!']);
    end
end

% Size error checking			     
k = size(t.lambda); 

if length(k) ~= length(t.u)
    error(['LAMBDA has order ', int2str(length(k)), ...
	   ' but there are ', int2str(length(t.u)), ' matrices.']);
end

for i = 1 : length(t.u)            
    if  size(t.u{i},2) ~= k(i)
	error(['Matrix U' int2str(i) ' does not have ' int2str(k(i)) ...
	       ' columns.']);
    end
end

t = class(t, 'tucker_tensor');
return;
