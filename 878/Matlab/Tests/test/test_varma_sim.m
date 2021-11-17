% TEST_VARMA_SIM  Check that varma_sim and CC_build work as documented
%
%   TEST_VARMA_SIM prints theoretical covariance, data covariance of simulated
%     series, and the max relative difference between the two, for 8 testcases.
%   TEST_VARMA_SIM QUIET prints only the max difference and only for 4 cases.
%   TEST_VARMA_SIM(N) tests only case N
%
%   The routine CC_build is also checked.

function test_varma_sim(varargin)
  if nargin > 0
    if isequal(lower(varargin{1}),'quiet'), tcase = 1:2:7; quiet = true;
    else tcase = varargin{1}; quiet = false; end
  else
    quiet = false;
    tcase = 1:1:8;
  end
  n = 20000;
  M = 20000;
  dmax = 0;
  fprintf('TESTING VARMA_SIM... ');
  if quiet, step=2; else, step=1; end
  for i=tcase
    [A, B, Sig, p, q, r, name] = testcase(i);
    h = max(p,q);
    %
    % TEST DIMENSIONS OF RESULTS
    X1 = varma_sim(A,B,Sig,10);
    X2 = varma_sim(A,B,Sig,10,[],1);
    X3 = varma_sim(A,B,Sig,10,[],5); 
    [X4,eps4] = varma_sim(A,B,Sig,10);
    [X5,eps5] = varma_sim(A,B,Sig,10,[],5); 
    assertsize(X1,r,10);
    assertsize(X2,r,10);
    assertsize(X4,r,10); assertsize(eps4,r,10); 
    if r==1,
      assertsize(X3,10,5);
      assertsize(X5,10,5); assertsize(eps5,10,5); 
    else
      assertsize(X3,r,10,5);
      assertsize(X5,r,10,5); assertsize(eps5,r,10,5); 
    end
    %
    % TEST MEAN
    mu = (1:r)';
    X = varma_sim(A, B, Sig, 200, mu, 200);
    if r==1, X=reshape(X,1,200,200); end
    mud = mean(mean(X,3),2);
    d = max(abs(mud-mu))/max(mu);
    if ~quiet
      fprintf('\nTESTCASE %s', name);
      if is_stationary(A, Sig) fprintf(' (model is stationary)'); 
      else fprintf(' *** WARNING--NONSTATIONARY MODEL--'); end
      fprintf(':\n');
      fprintf('  Difference between theoretical and data mean %.1f%%\n', d*100);
    end
    dmax = max(dmax,d);
    %
    % START FROM SCRATCH; COMPARE DATA COVARIANCE WITH THEORETICAL
    X = varma_sim(A, B, Sig, n, [])';     % simulate a long single series
    SSd = cov([X(1:end-1,:) X(2:end,:)]); % data covariance, lag 0 and 1
    [C,G,W,S] = find_CGWS(A, B, Sig);
    SS = S_build(S, A, G, 2); % theoretical covariance
    d = max(abs(SSd(:)-SS(:)))/max(abs(SS(:)));
    dmax = max(dmax,d);
    if ~quiet
      disp('  Theoretical covariance at lags 0 and 1:')
      for j=1:r, fprintf(1,'  %6.3f',SS(j,:)); disp ' '; end
      disp('  Data covariance at lags 0 and 1:')
      for j=1:r, fprintf(1,'  %6.3f',SSd(j,:)); disp ' '; end
      fprintf('  Max relative difference: %.1f%%\n', d*100);
    end
    X = varma_sim(A, B, Sig, h+3, [], M); % simulate multiple short series
    X = reshape(X,r*(h+3),M);
    SSd = cov(X');                % data covariance
    SS = S_build(S, A, G, h+3);  % theoretical covariance
    d = max(abs(SSd(:)-SS(:)))/max(abs(SS(:)));
    dmax = max(dmax,d);
    if ~quiet
      disp('  Data covariance for multiple series of length h+3:')
      for j=1:r, fprintf(1,'  %6.3f', SSd(j,:)); disp(' '); end
      fprintf('  Max relative difference: %.1f%%\n', d*100);
    end
    % CHECK CC_BUILD ON A SERIES OF LENGTH 3*h; COMPARE CC WITH DATA COV(X,EPS)
    CC = CC_build(A, C, h);
    [X,eps] = varma_sim(A, B, Sig, 2*h, [], M);
    if r==1, X=shiftdim(X,-1); end
    if r==1, eps=shiftdim(eps,-1); end
    X1 =reshape(X  (:,1:h,:),r*h,[]); X2 =reshape(X  (:,end-h+1:end,:),r*h,[]);
    ep1=reshape(eps(:,1:h,:),r*h,[]); ep2=reshape(eps(:,end-h+1:end,:),r*h,[]);
    CC1=cov([X1' ep1']);              CC2=cov([X2', ep2']);
    CC1=CC1(1:r*h,r*h+1:end);         CC2=CC2(1:r*h,r*h+1:end);
    d1 = max(abs(CC1(:)-CC(:)))/max(abs(CC(:)));
    d2 = max(abs(CC2(:)-CC(:)))/max(abs(CC(:)));
    dmax = max([dmax,d1,d2]);
    if ~quiet
      disp('  Theoretical cov(x,eps):')
      for j=1:r*h, fprintf(1,'  %6.3f',CC(j,:)); disp ' '; end
      disp('  CC1:')
      for j=1:r*h, fprintf(1,'  %6.3f',CC1(j,:)); disp ' '; end
      disp('  CC2:')
      for j=1:r*h, fprintf(1,'  %6.3f',CC2(j,:)); disp ' '; end
      fprintf('  Max relative difference: %.1f%%\n', max(d1,d2)*100);
    end
    %
    % TEST STARTING FROM GIVEN X'S
    X0 = 2+repmat((1:r)', 1, h);
    X = varma_sim(A, B, Sig, h+2, mu, M, X0);          % simulated forecast
    if r==1, X=shiftdim(X,-1); end
    X1 = reshape(X(:,h+1,:),r,M);
    X2 = reshape(X(:,h+2,:),r,M);
    [Ex1,Ex2,Vx1] = expX(A,B,C,X0,mu,SS(1:r*h,1:r*h),Sig); % theoretical E & V
    d1 = max(abs(Ex1 - mean(X1,2)))/max(abs(Ex1)); % difference one step ahead
    d2 = max(abs(Ex2 - mean(X2,2)))/max(abs(Ex2)); % difference two steps ahead
    d3 = max(max(abs(Vx1-cov(X1'))))/max(abs(Vx1(:)));  % x1 covar difference
    if ~quiet
      fprintf('  Max one-step-ahead forecast difference: %.1f%%\n', d1*100);
      fprintf('  Max two-step-ahead forecast difference: %.1f%%\n', d2*100);
      fprintf('  Max one-step-ahead forecast var difference: %.1f%%\n', d3*100);
    end
    dmax = max([dmax,d1,d2,d3]);
  end
  fprintf('Max difference: %.1f%%\n', dmax*100);
end

function [Ex1,Ex2,Vx1] = expX(A,B,C,X,mu,SS,Sig)
  % Calculate theoretical one and two step expected forecast & one step variance
  r=size(C{1},1);
  t = size(X,2)+1;
  CC = CC_build(A, C, t-1);
  A=makecell(A);
  B=makecell(B);
  X=X-repmat(mu,1,t-1);
  p=length(A);
  q=length(B);
  Ex1 = zeros(r,1); Ex2 = zeros(r,1);
  xp = reshape(X,[],1);
  for i=1:p
    Ex1 = Ex1 + A{i}*X(:,t-i);
  end
  Eeps = reshape(CC'*(SS\xp), r, t-1);
  BB = zeros(r,0);
  EE = [];
  for i=1:q
    Ex1 = Ex1 + B{i}*Eeps(:,t-i);
    BB = [B{i} BB];
    EE = blkdiag(EE, Sig);
  end
  Vx1 = CC'*(SS\CC);
  Vx1 = Sig + BB*(EE - Vx1(end-r*q+1:end,end-r*q+1:end))*BB';
  for i=1:p
    if i==1, Ex2 = Ex2 + A{i}*Ex1;        end
    if i> 1, Ex2 = Ex2 + A{i}*X(:,t-i+1); end
  end
  for i=2:q
    Ex2 = Ex2 + B{i}*Eeps(:,t-i+1);
  end
  Ex1 = Ex1 + mu;
  Ex2 = Ex2 + mu;
end

function assertsize(A,n,m,k)
  if nargin == 3 || k==1
    ascertain(isequal(size(A), [n,m]));
  else
    ascertain(isequal(size(A), [n,m,k]));
  end
end
