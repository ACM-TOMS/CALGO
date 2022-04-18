function [varargout] = eigen(A,k,its,l,nnd_flag)
%EIGEN  Low-rank eigendecomposition of a SELF-ADJOINT matrix
%
%
%   [V,D] = EIGEN(A)  constructs a nearly optimal rank-6 approximation
%           VDV' to A, using 4 normalized power iterations, with block
%           size 6+2=8, started with an n x 8 random matrix, when A is
%           n x n; the reference below explains "nearly optimal." When
%           A is the only input to EIGEN, the dimension n of A must be
%           >= 6.
%
%   [V,D] = EIGEN(A,k)  constructs a nearly optimal rank-k approx.
%           VDV' to A, using 4 normalized power iterations, with block
%           size k+2, started with an n x (k+2) random matrix, when A
%           is n x n; the reference below explains "nearly optimal."
%           k must be a positive integer <= the dimension n of A.
%
%   [V,D] = EIGEN(A,k,its)  constructs a nearly optimal rank-k approx.
%           VDV' to A, using its normalized power iterations, with block
%           size k+2, started with an n x (k+2) random matrix, when A
%           is n x n; the reference below explains "nearly optimal."
%           k must be a positive integer <= the dimension n of A, and
%           its must be a nonnegative integer.
%
%   [V,D] = EIGEN(A,k,its,l)  constructs a nearly optimal rank-k
%           approximation VDV' to A, using its normalized power
%           iterations, with block size l, started with an n x l random
%           matrix, when A is n x n; the reference below explains
%           "nearly optimal." k must be a positive integer <= the
%           dimension n of A, its must be a nonnegative integer, and l
%           must be a positive integer >= k.
%
%   [V,D] = EIGEN(A,k,its,l,nnd_flag)  constructs a nearly optimal rank-k
%           approximation VDV' to A, using its normalized power
%           iterations, with block size l, started with an n x l random
%           matrix, when A is n x n; the reference below explains
%           "nearly optimal." k must be a positive integer <= the
%           dimension n of A, its must be a nonnegative integer, and l
%           must be a positive integer >= k.  The nnd flag is a boolean
%           variable signifying the input matrix A is nonnegative
%           definite.
%           
%
%   Replacing "[V,D]" with "lambda" in any of the above produces the
%   vector lambda = diag(D) rather than both V and D.
%
%
%   The low-rank approximation VDV' comes in the form of an
%   eigendecomposition -- the columns of V are orthonormal and D is a
%   real diagonal matrix whose diagonal entries are nonnegative and
%   nonincreasing. V is n x k and D is k x k, when A is n x n.
%
%   Increasing its or l improves the accuracy of the approximation VDV';
%   the reference below describes how the accuracy depends on its and l.
%   Please note that even its=1 guarantees superb accuracy, whether or
%   not there is any gap in the singular values of the matrix A being
%   approximated, at least when measuring accuracy as the spectral norm
%   ||A-VDV'|| of A-VDV' (relative to the spectral norm ||A|| of A).
%
%
%   Note: THE MATRIX A MUST BE SELF-ADJOINT
%
%   Note: EIGEN invokes RAND. To obtain repeatable results, reset the
%         seed for the pseudorandom number generator.
%
%   Note: The user may ascertain the accuracy of the approximation VDV'
%         to A by invoking DIFFSNORMS(A,V,D).
%
%
%   inputs (the first is required):
%   A -- matrix being approximated
%   k -- rank of the approximation being constructed;
%        k must be a positive integer <= the dimension of A, and
%        defaults to 6
%   its -- number of normalized power iterations to conduct;
%          its must be a nonnegative integer, and defaults to 4
%   l -- block size of the normalized power iterations;
%        l must be a positive integer >= k, and defaults to k+2
%   nnd_flag -- set to true if the matrix is nonnegative definite,
%               but defaults to false
%
%   outputs (produces either the first two, V and D, or just lambda):
%   V -- n x k matrix in the rank-k approximation VDV' to A, where A is
%        n x n
%   D -- k x k diagonal matrix in the rank-k approximation VDV' to A,
%        with nonnegative and nonincreasing diagonal entries
%   lambda -- k x 1 vector equal to diag(D)
%
%
%   Example:
%     A = rand(1000,2);
%     A = A*A';
%     A = A/normest(A);
%     [V,D] = eigen(A,2);
%     diffsnormschur(A,V,D)
%
%     This example produces a rank-2 approximation VDV' to A such that
%     the columns of V are orthonormal and D is a diagonal matrix whose
%     diagonal entries are real and their absolute nonincreasing.
%     diffsnormschur(A,V,D) outputs an estimate of the spectral norm of
%     A-VDV', which should be close to the machine precision.
%
%
%   Reference:
%   Nathan Halko, Per-Gunnar Martinsson, and Joel Tropp,
%   Finding structure with randomness: probabilistic algorithms
%   for constructing approximate matrix decompositions,
%   arXiv:0909.4061 [math.NA; math.PR], 2009
%   (available at http://arxiv.org).
%

%   Copyright 2016 Huamin Li, George C. Linderman, Kelly Stanton

%
% Check the number of inputs.
%
if(nargin < 1)
    error('MATLAB:eigen:TooFewIn',...
        'There must be at least 1 input.')
end

%
% Check the number of outputs.
%
if(nargout > 2)
    error('MATLAB:eigen:TooManyOut',...
        'There must be at most 2 outputs.')
end

%
% Set the inputs k, its, and l to default values, if necessary.
%
if (nargin < 2)
    k=6;
end
if (nargin <3)
    its = 4;
end

if (nargin <4)
    l = k+2;
end

if (nargin <5)
    nnd_flag = false;
end


%
% Check the first input argument.
%
if(~isfloat(A))
    error('MATLAB:eigen:In1NotFloat',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
    error('MATLAB:eigen:In1Empty',...
        'Input 1 must not be empty.')
end

%
% Retrieve and check the dimensions of A.
%
[m, n] = size(A);

if(m ~= n)
    error('MATLAB:eigen:In1NotSymm',...
        'Input 1 must be square.')
end

%
% Check the remaining input arguments.
%
if(~isscalar(k))
    error('MATLAB:eigen:In2Not1x1',...
        'Input 2 must be a scalar.')
end

if(~isscalar(its))
    error('MATLAB:eigen:In3Not1x1',...
        'Input 3 must be a scalar.')
end

if(~isscalar(4))
    error('MATLAB:eigen:In4Not1x1',...
        'Input 4 must be a scalar.')
end

if(~(k > 0))
    error('MATLAB:eigen:In2NonPos',...
        'Input 2 must be > 0.')
end

if((nargin > 1) && (k > n))
    error('MATLAB:eigen:In2TooBig',...
        'Input 2 must be <= the dimension of Input 1.')
end

if((nargin == 1) && (k > n))
    error('MATLAB:eigen:InTooTiny',...
        'The dimension of the input must be >= 6.')
end

if(~(its >= 0))
    error('MATLAB:eigen:In3Bad',...
        'Input 3 must be >= 0.')
end

if(l < k)
    error('MATLAB:eigen:In4ltIn2',...
        'Input 4 must be >= Input 2.')
end

%
% Check whether A is self-adjoint to nearly the machine precision
% and warn the user if not (but do NOT report an error).
%
x = rand(n,1);
y = A*x;
z = A'*x;

if((norm(y-z) > .1d-11*norm(y)) || (norm(y-z) > .1d-11*norm(z)))
    warning('The input matrix is not exactly self-adjoint.')
    if(nnd_flag)
        warning('The matrix must also be nonnegative definite.')
        warning('This routine does not check for definiteness.')
    end
end



%
% Eigendecompose A directly if l >= n/1.25.
%
if(l >= n/1.25)
    [V,D] = eig(full(A));
    
    %
    % Rearrange the decomposition so that the absolute values
    % of the eigenvalues are nonincreasing.
    %
    [E,P] = sort(abs(diag(D)),'descend');
    D = diag(D);
    D = diag(real(D(P)));
    V = V(:,P);
    
    
    %
    % Retain only the leftmost k columns of V and
    % the uppermost leftmost k x k block of D.
    %
    V = V(:,1:k);
    
    D = D(1:k,1:k);
    
    if(nnd_flag)
        D = abs(D);
    end
    
    %
    % Fill the output array.
    %
    if((nargout == 0) || (nargout == 1))
        varargout{1} = diag(D);
    end
    
    if(nargout == 2)
        varargout{1} = V;
        varargout{2} = D;
    end
    
    return
    
end


%
% Apply A to a random matrix, obtaining Q.
%
R = (2*rand(n,l)-ones(n,l)) + ~isreal(A)*1i*(2*rand(n,l)-ones(n,l));

Q = A*R;

%
% Form a matrix Q whose columns constitute a well-conditioned basis
% for the columns of the earlier Q.
%
if(its == 0)
    if(nnd_flag)
        anorm = 0;
        for j = 1:l
            anorm = max(anorm,norm(Q(:,j))/norm(R(:,j)));
        end
    end
    [Q,~,~] = qr(Q,0);
else
    [Q,~] = lu(Q);
end

%
% Conduct normalized power iterations.
%
for it = 1:its
    
    if(nnd_flag)
        cnorm = zeros(l,1);
        for j = 1:l
            cnorm(j) = norm(Q(:,j));
        end
    end
    Q = A*Q;
    
    if(it < its)
        [Q,R] = lu(Q);
    else
        if(nnd_flag)
            anorm = 0;
            for j = 1:l
                anorm = max(anorm,norm(Q(:,j))/cnorm(j));
            end
        end
        [Q,R,E] = qr(Q,0);
        
    end
    
end


%
% Use the Nystrom method to obtain approximations to the eigenvalues
% and eigenvectors of A (shifting A on the subspace spanned by the
% columns of Q in order to make the shifted A be positive definite).
% An alternative is to use the (symmetric) square root in place of the
% Cholesky factor of the shift.
%

if(nnd_flag)
    anorm = .1d-6*anorm*sqrt(n);
    E = A*Q+anorm*Q;
    R = Q'*E;
    R = (R+R')/2;
    R = chol(R);
    [~,D,V] = svd(R'\E','econ');
    D = D*D-anorm*eye(l,l);
else
    R = Q'*A*Q;
    R = (R+R')/2;
    [V,D] = eig(R);
    V = Q*V;
    
    [~,P] = sort(abs(diag(D)),'descend');
    D = diag(D);
    D = diag(real(D(P)));
    V = V(:,P);
end

%
% Retain only the leftmost k columns of V and
% the uppermost leftmost k x k block of D.
%
V = V(:,1:k);
D = D(1:k,1:k);

if(nnd_flag)
    D = abs(D);
end

%
% Fill the output array.
%
if((nargout == 0) || (nargout == 1))
    varargout{1} = diag(D);
end

if(nargout == 2)
    varargout{1} = V;
    varargout{2} = D;
end
