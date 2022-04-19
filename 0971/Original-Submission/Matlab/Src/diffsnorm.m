function snorm = diffsnorm(A,U,S,V,raw,its)
%DIFFSNORM  2-norm accuracy of an approx. to a matrix.
%
%
%   snorm = DIFFSNORM(A,U,S,V)  computes an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of C(A)-USV', using 20 iterations of the power
%           method started with a random vector.  Where C(A) refers
%           to the matrix A from the input, after centering its columns.
%
%   snorm = DIFFSNORM(A,U,S,V,raw)  finds an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of A - USV' or C(A) - USV, using 20 iterations of the
%           power method started with a random vector; where C(A) refers
%           to the matrix A from the input, after centering its columns.
%
%   snorm = DIFFSNORM(A,U,S,V,raw, its)   finds an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of A - USV' or c(A) - USV, using its iterations of the
%           power method started with a random vector; where C(A) refers
%           to the matrix A from the input, after centering its columns.
%           its must be a positive integer.
%
%
%   Increasing its improves the accuracy of the estimate snorm of the
%   spectral norm of A-USV'.
%
%
%   Note: DIFFSNORM invokes RANDN. To obtain repeatable results, reset
%         the seed for the pseudorandom number generator.
%
%
%   inputs (the first four are required):
%   A -- first matrix in c(A)-USV' or A-USV' whose spectral norm is being
%        estimated; C(A) refers to A after centering its columns
%   U -- second matrix in c(A)-USV' or A-USV' whose spectral norm is being
%        estimated; C(A) refers to A after centering its columns
%   S -- third matrix in c(A)-USV' or A-USV' whose spectral norm is being
%        estimated; C(A) refers to A after centering its columns
%   V -- fourth matrix in c(A)-USV' or A-USV' whose spectral norm is being
%        estimated; C(A) refers to A after centering its columns
%   raw -- centers A when raw is false but does not when raw is true;
%          raw must be a Boolean and defaults to false
%   its -- number of iterations of the power method to conduct;
%          its must be a positive integer, and defaults to 20
%
%   output:
%   snorm -- an estimate of the spectral norm of A-USV' or c(A) - USV' (the
%            estimate fails to be accurate with exponentially low probability
%            as its increases; see references 1 and 2 below)
%
%
%   Example:
%     A = rand(1000,2)*rand(2,1000);
%     A = A/normest(A);
%     [U,S,V] = pcafast(A,2,true);
%     diffsnorm(A,U,S,V,true)
%
%     This example produces a rank-2 approximation USV' to A such that
%     the columns of U are orthonormal, as are the columns of V, and the
%     entries of S are all nonnegative and are zero off the diagonal.
%     diffsnorm(A,U,S,V) outputs an estimate of the spectral norm
%     of A-USV', which should be close to the machine precision.
%
%
%   References:
%   [1] Jacek Kuczynski and Henryk Wozniakowski, Estimating the largest
%       eigenvalues by the power and Lanczos methods with a random
%       start, SIAM Journal on Matrix Analysis and Applications, 13 (4):
%       1094-1122, 1992.
%   [2] Edo Liberty, Franco Woolfe, Per-Gunnar Martinsson, Vladimir
%       Rokhlin, and Mark Tygert, Randomized algorithms for the low-rank
%       approximation of matrices, Proceedings of the National Academy
%       of Sciences (USA), 104 (51): 20167-20172, 2007. (See the
%       appendix.)
%   [3] Franco Woolfe, Edo Liberty, Vladimir Rokhlin, and Mark Tygert,
%       A fast randomized algorithm for the approximation of matrices,
%       Applied and Computational Harmonic Analysis, 25 (3): 335-366,
%       2008. (See Section 3.4.)
%
%

%   Copyright 2016 Huamin Li, George C. Linderman, Kelly Stanton

%
% Check the number of inputs.
%
if(nargin < 4)
    error('MATLAB:diffsnorm:TooFewIn',...
        'There must be at least 4 inputs.')
end

%
% Set the input its to its default value 20, if necessary.
%
if(nargin < 5 )
    raw = false;
end

%
% Set the input its to its default value 20, if necessary.
%
if(nargin < 6 )
    its = 20;
end


%
% Check the input arguments.
%
if(~isfloat(A))
    error('MATLAB:diffsnorm:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
    error('MATLAB:diffsnorm:In1Empty',...
        'Input 1 must not be empty.')
end

if(~isfloat(U))
    error('MATLAB:diffsnorm:In2NotFloat',...
        'Input 2 must be a floating-point matrix.')
end

if(~isfloat(S))
    error('MATLAB:diffsnorm:In3NotFloat',...
        'Input 3 must be a floating-point matrix.')
end

if(isempty(S))
    error('MATLAB:diffsnorm:In3Empty',...
        'Input 3 must not be empty.')
end

if(~isfloat(V))
    error('MATLAB:diffsnorm:In4NotFloat',...
        'Input 4 must be a floating-point matrix.')
end

if(~isscalar(its))
    error('MATLAB:diffsnorm:In5Not1x1',...
        'Input 5 must be a scalar.')
end

if(~(its > 0))
    error('MATLAB:diffsnorm:In5NonPos',...
        'Input 5 must be > 0.')
end

%
% Retrieve the dimensions of A, U, S, and V.
%
[m n] = size(A);
[m2 k] = size(U);
[k2 l] = size(S);
[n2 l2] = size(V);

%
% Make sure that the dimensions of A, U, S, and V are commensurate.
%
if(m ~= m2)
    error('MATLAB:diffsnorm:In1In2BadDim',...
        'The 1st dims. of Inputs 1 and 2 must be equal.')
end

if(k ~= k2)
    error('MATLAB:diffsnorm:In2In3BadDim',...
        'The 2nd dim. of Input 2 must equal the 1st dim. of Input 3.')
end

if(l ~= l2)
    error('MATLAB:diffsnorm:In3In4BadDim',...
        'The 2nd dims. of Inputs 3 and 4 must be equal.')
end

if(n ~= n2)
    error('MATLAB:diffsnorm:In1In4BadDim',...
        'The 2nd dim. of Input 1 must equal the 1st dim. of Input 4.')
end

if raw
    % Do not normalize
    applyA = @(X) A*X;
    applyAT = @(X) A'*X;
else
    % Normalize by the average of the entries in every column.
    c = mean(A,1);
    applyA = @(X) bsxfun(@minus,A*X,c*X);
    applyAT = @(X) A'*X-c'*sum(X,1);
end

if(m >= n)
    
    %
    % Generate a random vector x.
    %
    x = randn(n,1) + ...
        (~isreal(A) || ~isreal(U) || ~isreal(S) || ~isreal(V)) * 1i*randn(n,1);
    x = x/norm(x);
    
    %
    % Run its iterations of the power method, starting with the random x.
    %
    for it = 1:its
        %
        %   Set y = (A-USV')x.
        %
        y = V'*x;
        y = S*y;
        y = U*y;
        y = applyA(x) -y;
        %
        %   Set x = (A'-c'*ones(1,m)-VS'U')y.
        %
        x = U'*y;
        x = S'*x;
        x = V*x;
        x = applyAT(y) -x;
        %
        %   Normalize x, memorizing its Euclidean norm.
        %
        snorm = norm(x);
        if(snorm == 0)
            break;
        end
        x = x/snorm;
    end
    
end


if(m < n)
    
    %
    % Generate a random vector y.
    %
    
    y = randn(1,m) + ...
        (~isreal(A) || ~isreal(U) || ~isreal(S) || ~isreal(V)) * 1i*randn(1,m);
    y = y/norm(y);
    
    %
    % Run its iterations of the power method, starting with the random y.
    %
    for it = 1:its
        %
        %   Set x = y(A-ones(m,1)*c-USV').
        %
        x = y*U;
        x = x*S;
        x = x*V';
        x = applyAT(y')' - x;
        
        %
        %   Set y = x(A'-c'*ones(1,m)-VS'U').
        %
        y = x*V;
        y = y*S';
        y = y*U';
        y = applyA(x')' - y;
        
        %
        %   Normalize y, memorizing its Euclidean norm.
        %
        snorm = norm(y);
        if(snorm == 0)
            break;
        end
        y = y/snorm;
    end
    
end

snorm = sqrt(snorm);
