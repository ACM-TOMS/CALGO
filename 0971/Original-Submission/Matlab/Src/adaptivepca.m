function [U,S,V] = adaptivepca(A,epsilon,b,p)
%ADAPTIVEPCA  Adaptive Randomized Range Finder (Algorithm randQB_pb of 
%Martinsson and Voronin (2015))
%
%
%   [V,D] = ADAPTIVEPCA(A) constructs an Orthonormal matrix Q such that
%           |(I-QQ')A| < 1.0e-7 using a block size of 10 and 1 power
%           iteration.  Then, it computes an approximate
%           factorization A=USV', where S is a nonnegative diagonal matrix,
%           and U and V are orthonormal.
%
%   [V,D] = ADAPTIVEPCA(A, epsilon) constructs an Orthonormal matrix Q such that
%           |(I-QQ')A| < epsilon using a block size of 10 and 1 power
%           iteration.  Then, it computes an approximate
%           factorization A=USV', where S is a nonnegative diagonal matrix,
%           and U and V are orthonormal.
%
%   [V,D] = ADAPTIVEPCA(A, epsilon,b) constructs an Orthonormal matrix Q such that
%           |(I-QQ')A| < epsilon using a block size of b and 1 power
%           iteration.  Then, it computes an approximate
%           factorization A=USV', where S is a nonnegative diagonal matrix,
%           and U and V are orthonormal.
%
%   [V,D] = ADAPTIVEPCA(A, epsilon,b,p) constructs an Orthonormal matrix Q such that
%           |(I-QQ')A| < epsilon using a block size of b and p power
%           iterations.  Then, it computes an approximate
%           factorization A=USV', where S is a nonnegative diagonal matrix,
%           and U and V are orthonormal.
%
%   inputs (the first is required):
%   A -- matrix being approximated
%   epsilon -- the precision of the approximation to the range of A:
%               |(I-QQ')A| < epsilon
%   b -- block size
%   p -- number of power iterations

%
%   outputs (all three are required):
%   U -- m x k matrix in the approximation USV' to A ,
%        where A is m x n; the columns of U are orthonormal
%   S -- k x k matrix in the rank-k approximation USV' to A, where A is 
%         m x n; the entries of S are all nonnegative, and all nonzero
%        entries appear in nonincreasing order on the diagonal
%   V -- n x k matrix in the rank-k approximation USV' to A,
%        where A is m x n; the columns of V are orthonormal
%
%   Please note that the present version of this code does not center the
%   matrix.
%
%   References:
%   [1] P. Martinsson and S. Voronin, A randomized blocked algorithm
%       for efficiently computing rank-revealing factorizations of 
%       matrices, pp. 1-12, 2015.

if(nargin < 1)
    error('MATLAB:adaptivepca:TooFewIn',...
        'There must be at least 1 input.')
end
if(nargin < 2)
    epsilon = 1.0e-7;
end
if(nargin < 3)
    b = 10;
end
if(nargin < 4)
    p = 1;
end

%
% Check the number of outputs.
%
if(nargout ~= 3)
    error('MATLAB:adaptivepca:WrongNumOut',...
        'There must be exactly 3 outputs.')
end

%
% Check the first input argument.
%
if(~isfloat(A))
    error('MATLAB:adaptivepca:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
    error('MATLAB:adaptivepca:In1Empty',...
        'Input 1 must not be empty.')
end

%
% Retrieve the dimensions of A.
%
A_org = A;
[m,n] = size(A);

%
% Check the remaining input arguments.
%
if(~isscalar(epsilon))
    error('MATLAB:adaptivepca:In2Not1x1',...
        'Input 2 must be a scalar.')
end

if(~isscalar(b))
    error('MATLAB:adaptivepca:In3Not1x1',...
        'Input 3 must be a scalar.')
end

if(~(b > 0))
    error('MATLAB:adaptivepca:In3NonPos',...
        'Input 3 must be > 0.')
end

if(~(p > 0))
    error('MATLAB:adaptivepca:In4NonPos',...
        'Input 4 must be > 0.')
end

if(~isscalar(p))
    error('MATLAB:adaptivepca:In4Not1x1',...
        'Input 4 must be a scalar.')
end

if(~(b <= min(m,n)))
    error('MATLAB:adaptivepca:In3TooBig',...
        'Input 3 must be <= min(size(A)).')
end

i = 0;
while normest(A,epsilon/3) > epsilon
    i = i + 1;
    if m < n
        % Draw standard Gaussian matrix W
        W = 2*rand(b,m)-1 + ~isreal(A)*1i*(2*rand(b,m)-1);
        % Compute Y = A*W
        Y = A'*W';
    else
        % Draw standard Gaussian matrix W
        W = 2*randn(n,b)-1 + ~isreal(A)*1i*(2*rand(n,b)-1);
        % Compute Y = A*W
        Y = A*W;
    end
    
    % Form a matrix Q whose columns constitute a well-conditioned basis
    % for the columns of the earlier Y
    [Q{i},~] = lu(Y);
    
    % Power scheme to enhance accuracy
    if m < n
        for j = 1:p
            [Q{i},~] = lu(A*Q{i});
            [Q{i},~] = qr(A'*Q{i},0);
        end
    else
        for j = 1:p
            [Q{i},~] = lu(A'*Q{i});
            [Q{i},~] = qr(A*Q{i},0);
        end
    end
    
    % Double Gram Schmidt
    if m < n
        if i > 1
            temp = cellfun(@(x,y)x*y',Q(1:i-1),Q(1:i-1),'un',0);
            temp = sum(cat(3, temp{:}), 3);
            [Q{i},~] = qr(Q{i}-temp*Q{i},0);
            
        end
    else
        if i > 1
            temp = cellfun(@(x,y)x*y',Q(1:i-1),Q(1:i-1),'un',0);
            temp = sum(cat(3, temp{:}), 3);
            [Q{i},~] = qr(Q{i}-temp*Q{i},0);
        end
    end
    
    % remove the orthogonal projection
    if m < n
        A = A - (Q{i}*(A*Q{i})')';
    else
        A = A - Q{i}*(Q{i}'*A);
    end
end

Q = [Q{1:end}];

if m < n
    % SVD the A applied to Q to obtain approximations
    % to the singular values and left singular vectors of the A;
    [U,S,R] = svd(A_org*Q,'econ');
    
    % Adjust the right singular vectors to approximate
    % the right singular vectors of A
    V = Q*R;
else
    % SVD Q' apploed to the A to obtain approximations
    % to the singular values and right singular vectors of the A;
    R = Q';
    [R,S,V] = svd(R*A_org,'econ');
    
    % Adjust the left singular vectors to approximate
    % the left singular vectors of A.
    U = Q*R;
end

% Retain only the leftmost k columns of U,
% the leftmost k columns of V,
% and the uppermost leftmost k * k block of S.
[k,~] = size(S);
U = U(:,1:k);
V = V(:,1:k);
S = S(1:k,1:k);
end

