function snorm = diffsnormschur(A,V,T,its)
%DIFFSNORMSCHUR  2-norm accuracy of a Schur decomp. of a matrix
%
%
%   snorm = DIFFSNORMSCHUR(A,V,T)  computes an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of A-VTV', using 20 iterations of the power
%           method started with a random vector.
%
%   snorm = DIFFSNORMSCHUR(A,V,T,its)  finds an estimate snorm of the
%           spectral norm (the operator norm induced by the Euclidean
%           vector norm) of A-VTV', using its iterations of the power
%           method started with a random vector; its must be a positive
%           integer.
%
%
%   Increasing its improves the accuracy of the estimate snorm of the
%   spectral norm of A-VTV'.
%
%
%   Note: DIFFSNORMSCHUR invokes RANDN. To obtain repeatable results, reset
%         the seed for the pseudorandom number generator.
%
%
%   inputs (the first three are required):
%   A -- first matrix in A-VTV' whose spectral norm is being estimated
%   V -- second matrix in A-VTV' whose spectral norm is being estimated
%   T -- third matrix in A-VTV' whose spectral norm is being estimated
%   its -- number of iterations of the power method to conduct;
%          its must be a positive integer, and defaults to 20
%
%   output:
%   snorm -- an estimate of the spectral norm of A-VTV' (the estimate
%            fails to be accurate with exponentially low probability
%            as its increases; see references 1 and 2 below)
%
%
%   Example:
%     A = rand(1000,2);
%     A = A*A';
%     A = A/normest(A);
%     [V,D] = eigen(A,2,4,4,true);
%     diffsnormschur(A,V,D)
%
%     This example produces a rank-2 approximation VDV' to A such that
%     the columns of V are orthonormal and D is a diagonal matrix whose
%     diagonal entries are nonnegative and nonincreasing.
%     diffsnormschur(A,V,D) outputs an estimate of the spectral norm of
%     A-VDV', which should be close to the machine precision.
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
%

%   Copyright 2016 Huamin Li, George C. Linderman, Kelly Stanton

%
% Check the number of inputs.
%
if(nargin < 3)
  error('MATLAB:diffsnormschur:TooFewIn',...
        'There must be at least 3 inputs.')
end

if(nargin > 4)
  error('MATLAB:diffsnormschur:TooManyIn',...
        'There must be at most 4 inputs.')
end

%
% Check the number of outputs.
%
if(nargout > 1)
  error('MATLAB:diffsnormschur:TooManyOut',...
        'There must be at most 1 output.')
end

%
% Set the input its to its default value 20, if necessary.
%
if(nargin == 3)
  its = 20;
end

%
% Check the input arguments.
%
if(~isfloat(A))
  error('MATLAB:diffsnormschur:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
  error('MATLAB:diffsnormschur:In1Empty',...
        'Input 1 must not be empty.')
end

if(~isfloat(V))
  error('MATLAB:diffsnormschur:In2NotFloat',...
        'Input 2 must be a floating-point matrix.')
end

if(~isfloat(T))
  error('MATLAB:diffsnormschur:In3NotFloat',...
        'Input 3 must be a floating-point matrix.')
end

if(isempty(T))
  error('MATLAB:diffsnormschur:In3Empty',...
        'Input 3 must not be empty.')
end

if(~isscalar(its))
  error('MATLAB:diffsnormschur:In4Not1x1',...
        'Input 4 must be a scalar.')
end

if(~(its > 0))
  error('MATLAB:diffsnormschur:In4NonPos',...
        'Input 4 must be > 0.')
end

%
% Retrieve the dimensions of A, V, and T.
%
[m n] = size(A);
[m2 k] = size(V);
[k2 k3] = size(T);

%
% Make sure that the dimensions of A, V, and T are commensurate.
%
if(m ~= n)
  error('MATLAB:diffsnormschur:In1BadDim',...
        'Input 1 must be square.')
end

if(m ~= m2)
  error('MATLAB:diffsnormschur:In1In2BadDim',...
        'The 1st dims. of Inputs 1 and 2 must be equal.')
end

if(k ~= k2)
  error('MATLAB:diffsnormschur:In2In3BadDim',...
        'The 2nd dim. of Input 2 must equal the 1st dim. of Input 3.')
end

if(k2 ~= k3)
  error('MATLAB:diffsnormschur:In3BadDim',...
        'Input 3 must be square.')
end


%
% Generate a random vector x.
%
x = randn(n,1) + (~isreal(A) || ~isreal(V) || ~isreal(T))*1i*randn(n,1);

x = x/norm(x);

%
% Run its iterations of the power method, starting with the random x.
%
for it = 1:its
%
% Set y = (A-VTV')x.
%
  y = (x'*V)';
  y = T*y;
  y = V*y;
  y = A*x-y;
%
% Set x = (A'-VT'V')y.
%
  x = V'*y;
  x = T'*x;
  x = V*x;
  x = A'*y-x;
%
% Normalize x, memorizing its Euclidean norm.
%
  snorm = norm(x);
  if(snorm == 0)
    break;
  end
  x = x/snorm;
end

snorm = sqrt(snorm);
