function [z,c,qdat] = deparam(w,beta,z0,options)
%DEPARAM Schwarz-Christoffel exterior parameter problem.
%       [Z,C,QDAT] = DEPARAM(W,BETA) solves the Schwarz-Christoffel
%       mapping parameter problem with a disk as fundamental domain and
%       the exterior of the polygon specified by W as the target.  W
%       must be a vector of the vertices of the polygon, specified in
%       clockwise order, and BETA should be a vector of the turning
%       angles of the polygon; see SCANGLES for details.  If successful,
%       DEPARAM will return Z, a vector of the pre-images of W; C, the
%       multiplicative constant of the conformal map; and QDAT, a matrix
%       of quadrature data used by some of the other S-C routines.
%
%       [Z,C,QDAT] = DEPARAM(W,BETA,Z0) uses Z0 as an initial guess for
%       Z.
%
%       [Z,C,QDAT] = DEPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
%       parameters.  See SCPARMOPT.
%	
%	See also SCPARMOPT, DRAWPOLY, DEDISP, DEPLOT, DEMAP, DEINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w); 				% no. of vertices
w = w(:);
beta = beta(:);

% Set up defaults for missing args
if nargin < 4
  options = [];
  if nargin < 3
    z0 = [];
  end
end

err = sccheck('de',w,beta);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

[trace,tol] = scparmopt(options);
nqpts = max(ceil(-log10(tol)),2);
qdat = scqdata(beta,nqpts); 		% quadrature data
  
%%if length(beta)~=n
%%  error('Mismatched angles and vertices')
%%elseif any(beta > 1) | any(beta <= -1)
%%  error('Each entry of beta must be in (-1,1]')
%%elseif abs(sum(beta)-2) > tol
%%  disp('Warning: angles do not sum to +2')
%%  if abs(sum(beta)+2) < tol
%%    disp('Vertices were probably specified in the wrong order.')
%%    disp('Use flipud and scangle to reverse ordering.')
%%    return
%%  end
%%elseif (beta(n)==0 | beta(n)==1) & (n > 2)
%%  error('Sides adjacent to w(n) must not be collinear')
%%elseif n < 2
%%  error('Polygon must have at least two vertices')
%%end

if n==2					% it's a slit
  z = [-1;1];

else
  % Set up normalized lengths for nonlinear equations
  len = abs(diff(w([n,1:n])));
  nmlen = abs(len(3:n-1)/len(2));
  
  % Set up initial guess
  if isempty(z0)
    y0 = zeros(n-1,1);
  else
    th = angle(z0(:));
    th(th<=0) = th(th<=0) + 2*pi;
    dt = diff([0;th(1:n-1);2*pi]);
    y0 = log(dt(1:n-1)./dt(2:n));
  end
  
  % Solve nonlinear system of equations:

  % package data
  nrow = max([n,nqpts,3]);
  ncol = 3+2*(n+1);
  fdat = zeros(nrow,ncol);
  fdat(1:3,1) = [n;nqpts;ncol];
  fdat(1:n,2) = beta;
  if n > 3
    fdat(1:n-3,3) = nmlen(:);
  end
  fdat(1:nqpts,4:ncol) = qdat;
  % set options
  opt = zeros(16,1);
  opt(1) = trace;
  opt(6) = 100*(n-3);
  opt(8) = tol;
  opt(9) = tol/10;
  opt(12) = nqpts;
  % do it
  [y,termcode] = nesolve('depfun',y0,opt,fdat);
  if termcode~=1
    disp('Warning: Nonlinear equations solver did not terminate normally')
  end
  
  % Convert y values to z
  cs = cumsum(cumprod([1;exp(-y)]));
  theta = 2*pi*cs/cs(n);
  z = ones(n,1);
  z(1:n-1) = [exp(i*theta(1:n-1))];
end

% Determine scaling constant
mid = exp(i*mean(angle(z(n-1:n))));
c = (w(n) - w(n-1)) / (dequad(z(n-1),mid,n-1,z,beta,qdat)-...
    dequad(z(n),mid,n,z,beta,qdat));


