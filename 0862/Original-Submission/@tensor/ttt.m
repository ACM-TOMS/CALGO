function c = ttt(varargin)
%TTT Tensor mulitplication (tensor times tensor).
% 
%   TTT(A,B) computes the outer product of tensors A and B.
%
%   TTT(A,B,ADIMS,BDIMS) computes the contracted product of tensors 
%   A and B in the dimensions specified by the row vectors ADIMS and 
%   BDIMS.  The sizes of the dimensions specified by ADIMS and BDIMS 
%   must match; that is, size(A,ADIMS) must equal size(B,BDIMS). 
%
%   TTT(A,B,DIMS) computes the inner product of tensors A and B in the
%   dimensions specified by the vector DIMS.  The sizes of the
%   dimensions specified by DIMS must match; that is, size(A,DIMS) must
%   equal size(B,DIMS). 
%
%   Examples
%      A = tensor(rand(4,2,3));
%      B = tensor(rand(3,4,2));
%      C = ttt(A,B) computes the outer produce of A and B, and the
%      result is of size 4 x 2 x 3 x 3 x 4 x 2.
%      D = ttt(A,A,1:3) computes the inner product of A with itself.
%      E = ttt(A,B,[1 2 3],[2 3 1]) computes the inner (or contracted)
%      product between A and F = permute(B,[2 3 1]);
%      G = ttt(A,B,[1 3],[2 1]) computes the product of A and B
%      along the specified dimenions, and the result is a 2 x 2
%      tensor. 
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/TTM, TENSOR/TTV.

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

%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

% Check the number of arguments
if (nargin < 2)
    error('TTT requires at least two arguments.');
end

% Check the first argument
if ~isa(varargin{1}, 'tensor')
    error('First argument must be a tensor.');
else
    a = varargin{1};
end

% Check the second argument
if ~isa(varargin{2}, 'tensor')
    error('Second argument must be a tensor.');
else
    b = varargin{2};
end

% Optional 3rd argument
if nargin >= 3
    adims = varargin{3};
else
    adims = [];
end

% Optional 4th argument
if nargin >= 4
    bdims = varargin{4};
else
    bdims = adims;
end

if ~isequal(size(a,adims),size(b,bdims))
    error('Dimensions specified by ADIMS and BDIMS do not match.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE THE PRODUCT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Avoid transpose by reshaping A and computing C = A * B
amatrix = tensor_as_matrix(a,setdiff(1:ndims(a),adims));
bmatrix = tensor_as_matrix(b,bdims);
cmatrix = amatrix * bmatrix;

c = tensor(cmatrix);
