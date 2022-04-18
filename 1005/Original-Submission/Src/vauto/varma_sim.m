% VARMA_SIM  Simulate an ARMA or a VARMA time series model
% 
%   ARMA MODELS:
%     X = VARMA_SIM(A,B,sigma,n) with scalar sigma generates a zero-mean time
%     series of length n using the ARMA(p,q) model:
% 
%              x(t) = A1·x(t-1) + ... + Ap·x(t-p) + y(t)                    (1a)
%     where:
%              y(t) = eps(t) + B1·eps(t-1) + ... + Bq·eps(t-q)              (1b)
% 
%     and eps(t) is N(0,sigma). A should be the row vector [A1,...,Ap] and B
%     should be [B1,...,Bq]. X returns the series as a row vector. To start the
%     simulation, x(1),...,x(g),eps(1),...,eps(g) are drawn from N(0,SCE) where
%     SCE = [SS CC; CC' EE] is obtained by solving the modified Yule-Walker
%     equations (see vyw_factorize, vyw_solve and S_build), so there is no need
%     to throw away the first values to avoid spin-up effects. The model must be
%     stationary.
% 
%     X = VARMA_SIM(A,B,sigma,n,mu) returns a time series with mean mu using
%     y(t) as above, and x(t) given by:
% 
%             x(t) - mu = A1·(x(t-1) - mu) + ... + Ap·(x(t-p) - mu) + y(t)   (2)
% 
%     In this case x and eps are drawn from N([mu; 0], SCE).
% 
%     X = VARMA_SIM(A,B,sigma,n,mu,M) generates M sequences simultaneously and
%     returns the i-th one in the i-th row of X. Use mu=[] to obtain zero-mean.
% 
%     X = VARMA_SIM(A,B,sigma,n,mu,M,X0) sets x(1)...x(g) to X0(end-p+1:end) -
%     mu and draws eps(1)...eps(g) from the conditional distribution of eps(1),
%     ..., eps(g) given x(1),...,x(g)). For this option the model need not be
%     stationary.
% 
%     [X,EPS] = VARMA_SIM(...) returns also the shock series, eps(t), in the
%     t-th row of EPS.
% 
%   VARMA MODELS:
%     X = VARMA_SIM(A,B,Cov,n), where Cov is an r×r matrix uses a VARMA(p,q)
%     model given by (1), but now x(t) is r-dimensional, eps(t) is r-variate
%     normal with mean 0 and covariance Cov, and Cov and the Ai's and Bi's are
%     r×r matrices. A and B should contain A = [A1 A2...Ap] (an r × r·p matrix)
%     and B = [B1 B2 ... Bq] (an r × r·q matrix). The r×n matrix X returns x(t)
%     in its t-th column. As in the scalar case the simulated series is spin-up
%     free, the starting values x(1),...,x(g) and eps(1),...,eps(g) being drawn
%     from the correct distribution. The model must be stationary.
%  
%     X = VARMA_SIM(...,mu) uses (2) for x(t) instead of (1a).
% 
%     X = VARMA_SIM(...,mu,M) returns M sequences in an r×n×M multidimensional X
%     (when r > 1) with the i-th sequence in X(1:r, 1:n, i). Use mu = [] to
%     obtain zero-mean.
% 
%     X = VARMA_SIM(...,mu,M,X0) initializes x(1),...,x(h) with the last h
%     columns of X0 instead of drawing from N(mu,SS). The model need not be
%     stationary.
% 
%     [X,EPS] = VARMA_SIM(...) returns also the shock series; the i-th eps(t) is
%     returned in EPS(:, t, i).
% 
%   For both ARMA and VARMA, use VARMA_SIM(A,[],...) for a pure autoregressive 
%   model, and VARMA_SIM([],B,...) for a pure moving average model.
% 
%   The method used is described in [1] and [2].
%
%   [1] K Jonasson and SE Ferrando 2006. Efficient likelihood evaluation for
%       VARMA processes with missing values. Report VHI-01-2006, Engineering
%       Research Institute, University of Iceland.
%
%   [2] K Jonasson 2006. Matlab programs for complete and incomplete data exact
%       VARMA likelihood and its gradient. Report VHI-02-2006, Engineering
%       Research Institute, University of Iceland.
%
%   Kristján Jónasson, Dept. of Computer Science, University of Iceland, 2006.
%   jonasson@hi.is.

function [x, eps] = varma_sim(A, B, Sig, n, mu, M, x0)
  r = size(Sig, 1);
  if isempty(A), A = zeros(r,0); end
  if isempty(B), B = zeros(r,0); end  
  p = size(A, 2)/r;
  q = size(B, 2)/r;
  Aflp = reshape(flipdim(reshape(A,r,r,p), 3), r, r*p);  %  [Ap...A2 A1]
  Bflp = reshape(flipdim(reshape(B,r,r,q), 3), r, r*q);  %  [Bq...B2 B1]
  h = max(p,q);
  if n<h, error('Too short series'); end
  if nargin < 5 || isempty(mu), mu = zeros(r,1); else mu = mu(:); end
  if nargin < 6 || isempty(M), M=1; end
  if nargin < 7, x0 = []; end
  I = 1:r*h;
  x = zeros(n*r,M);
  eps = zeros(n*r,M);
  mup = repmat(mu,h,1);
  [C, G] = find_CGW(A, B, Sig);
  PLU = vyw_factorize(A);
  if ~isempty(PLU) && ~isempty(PLU{1}) && PLU{1}(1) == 0,
    error('Non-stationary model: X0 must be specified.'),
  end
  S = vyw_solve(A, PLU, G);
  SS = S_build(S, A, G, h);
  CC = CC_build(A, C, h);
  [R,pp] = chol(SS);
  if pp~=0, error('Non-stationary model: X0 must be specified.'), end
  EE = mat2cell(zeros(r*h,r*h), repmat(r,1,h), repmat(r,1,h));
  for i=1:h, EE{i,i} = Sig; end
  EE = cell2mat(EE);
  if isempty(x0) %  Start the sequence from scratch
    x(I,:) = randnm(M, SS)';
  else %  Initialize the series with x0 (every generated series starts the same)
    ascertain(length(x0)>=h);
    if r == 1, x0 = x0(:)'; end
    x(I,:) = repmat(reshape(x0(:,end-h+1:end),[],1) - mup, 1, M);
  end
  eps(I,:) = randnm(M, EE - CC'*(SS\CC))' + CC'*(SS\x(I,:));
  eps(r*h+1:end,:) = reshape(randnm((n-h)*M,Sig)', [], M); %  shocks
  % eps = reshape(randnm(n*M,Sig)', [], M); %  for testing
  Ip = r*(h-p)+1:r*h;
  Iq = r*(h-q)+1:r*h;
  I = r*h+1:r*h+r;
  for t=h+1:n %  Generate the series
    x(I,:) = eps(I,:) + Aflp*x(Ip,:) + Bflp*eps(Iq,:);
    I = I+r;
    Ip = Ip+r;
    Iq = Iq+r;
  end
  if r==1 && M==1  %  one ARMA sequence:
    x = reshape(x,1,n) + mu;
    eps = reshape(eps,1,n);
  elseif r==1      %  several ARMA sequences:
    x = reshape(x,n,M) + mu;
    eps = reshape(eps,n,M);
  else             %   one or more VARMA sequences in r×n×M array:
    x = reshape(x,r,n,M) + repmat(mu,[1,n,M]);
    eps = reshape(eps,r,n,M); 
  end
end

function x = randnm(n,Sig,mu)
  %  RANDNM  Multivariate normal random vectors
  %    X = RANDNM(N,SIG,MU) generates N random vectors from an r-variate normal
  %    distribution with mean MU and covariance matrix SIG. MU should be an
  %    r-vector. X will be N×r. X = RANDNM(N,S) uses mean 0.
  %
  %    NOTE: The formula for del was found empirically to be approximately the
  %    minimum necessary to ensure that Cholesky factorization of Sig+del*I
  %    works. The while loop is to further ensure that.
  r=size(Sig,1);
  if nargin<3, mu=zeros(1,r); end
  mu = mu(:)';
  [R,p] = chol(Sig);
  if p~=0, %  Add delta to diagonal of Sig to handle postive semidefinite Sig
    del = max(1,norm(diag(Sig),inf))*(3+r/50)*eps; %  eps is machine epsilon
    kdel = 0;
    while p~=0
      Sig = Sig + del*eye(r);
      kdel = kdel + 1;
      [R,p] = chol(Sig);
    end
    ascertain(kdel<10);
  end
  x = randn(n,r)*R + repmat(mu,n,1);
end
