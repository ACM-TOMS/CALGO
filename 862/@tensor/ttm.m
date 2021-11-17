function c = ttm(a,v,n)
%TTM Tensor times matrix.
%
%   B = TTM(A,V,N) computes the n-mode product of tensor A with a
%   matrix V; i.e., A x_N V.  The integer N specifies the dimension
%   (or mode) of A along which V should be multiplied.  If size(V) =
%   [J,I], then A must have size(A,N) = I.  The result will be the
%   same order and size as A except that size(B,N) = J.
%
%   B = TTM(A,U) computes the n-mode product of tensor A with a
%   sequence of matrices in the cell array U.  The n-mode products are
%   computed sequentially along all dimensions (or modes) of A. The
%   cell array U contains ndims(A) matrices.
%
%   B = TTM(A,U,DIMS) computes the sequence tensor-matrix products
%   along the dimensions specified by DIMS.
%
%   Examples
%        A = tensor(rand(5,3,4,2));
%        V = rand(4,5); X = rand(4,3); Y = rand(3,4); Z = rand(3,2);
%        U = {V,X,Y,Z};
%     The following equivalent statements compute B = A x_1 V, the
%     result of which is a 4 x 3 x 4 x 2 tensor.
%        B = ttm(A, V, 1) 
%        B = ttm(A, U, 1) 
%     The following equivalent statements compute 
%     B = A x_1 V x_2 X x_3 Y x_4 Z, the result of which is a 
%     4 x 4 x 3 x 3 tensor.
%        B = ttm(A, {V,X,Y,Z}, [1 2 3 4])
%        B = ttm(A, U, [1 2 3 4])
%        B = ttm(A, U) 
%     The following equivalent statements compute B = A x_3 Y x_4 Z,
%     the result of which is a 5 x 3 x 3 x 3 tensor.
%        B = ttm(A, {Y,Z}, [3 4]) 
%        B = ttm(A, U, [3 4])
%     The following equivalent statements compute 
%     B = A x_1 V x_2 X x_4 Z, the result of which is a 
%     4 x 4 x 4 x 3 tensor.
%        B = ttm(A, {V,X,Z}, [1 2 4]) 
%        B = ttm(A, U, [1 2 4]) 
%        B = ttm(A, {V,X,Z}, -3)
%        B = ttm(A, U, -3) 
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/TTT, TENSOR/TTV.

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
    error('TTM requires at least two arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHECK FOR CELL ARRAY %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(v,'cell')
    if (nargin < 3)
        n = 1:ndims(a);
    end
    c = ttms(a,v,n);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%
%%% ERROR CHECKING %%%
%%%%%%%%%%%%%%%%%%%%%%

% Check the first argument
if ~isa(a, 'tensor')
    error('A must be a tensor.');
end

% Check the second argument
if ndims(v) ~= 2
    error('V must be a matrix.');
end

% Check the third argument
if (nargin < 3) || ~isequal(size(n),[1 1]) || (n < 0) || (n > ndims(a))
    error('Dimension N must be between 1 and NDIMS(A).');
end

% Check that sizes match!
if size(a,n) ~= size(v,2)
    error('Nth dimension of A is not equal to number of columns in V.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTE THE PRODUCT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert v to a tensor
v = tensor(v);

% Compute product
c = ttt(a, v, n, 2);

% Permute rows to proper order in c
dims = [setdiff(1:ndims(a),n) n];
[sdims sidx] = sort(dims);
c = permute(c,sidx);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = ttms(a,v,dims)

% Check validity of DIMS and then check for "minus" case
if (max(abs(dims)) > ndims(a))
    error('An entry in DIMS exceeds order of A.');
elseif (max(dims) < 0)
    dims = setdiff(1:ndims(a), -dims);
end
    
% Check validity of parameters passed to TTMS
n = length(dims);
if (n > ndims(a)) || (n > length(v))
    error('DIMS is too long.');
elseif (n < length(v)) && (length(v) < ndims(a))
    error('If length(DIMS) < length(V), then length(V) must equal ndims(A).');
elseif (length(v) > ndims(a))
    error('Length of V greater than order of A.');
end

% Check sizes of V and DIMS to determine version to use.
if (n == length(v))
    vidx = 1:n;   % index V by order 
else 
    vidx = dims;  % index V by dimension 
end

% Calculate individual products
c = ttm(a, v{vidx(1)}, dims(1));
for i = 2 : n
    c = ttm(c, v{vidx(i)}, dims(i));
end
