function [x,c,qdat] = hpparam(w,beta,x0,options);
%HPPARAM Schwarz-Christoffel half-plane parameter problem.
%	[X,C,QDAT] = HPPARAM(W,BETA) solves the Schwarz-Christoffel
%	parameter problem with the upper half-plane as fundamental
%	domain and interior of the specified polygon as the target.  W
%	must be a vector of the vertices of the polygon, specified in
%	counterclockwise order.  BETA is a vector of turning angles; see
%	SCANGLES.  If successful, HPPARAM will return X, a vector of the
%	pre-images of W; C, the multiplicative constant of the conformal
%	map; and QDAT, a matrix of quadrature data used by some of
%	the other S-C routines.
%
%       [X,C,QDAT] = HPPARAM(W,BETA,X0) uses X0 as an initial guess for
%       X.
%
%       [X,C,QDAT] = HPPARAM(W,BETA,X0,OPTIONS) uses a vector of control
%       parameters.  See SCPARMOPT.
%	
%	See also SCPARMOPT, DRAWPOLY, HPDISP, HPPLOT, HPMAP, HPINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/26/95.

n = length(w); 				% no. of vertices
w = w(:);
beta = beta(:);

% Set up defaults for missing args
if nargin < 4
  options = [];
  if nargin < 3
    x0 = [];
  end
end

err = sccheck('hp',w,beta);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

[trace,tol] = scparmopt(options);
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta(1:n-1),nqpts); 	% quadrature data

atinf = (beta <= -1);

% Find prevertices (solve param problem)
if n==3
  x = [-1;0;Inf];

else

  % Set up normalized lengths for nonlinear equations:

  % indices of left and right integration endpoints
  left = 1:n-2;				
  right = 2:n-1;				
  % delete indices corresponding to vertices at Inf
  left(find(atinf)) = [];
  right(find(atinf) - 1) = [];
  cmplx = ((right-left) == 2);
  % normalize lengths by w(2)-w(1)
  nmlen = (w(right)-w(left))/(w(2)-w(1));
  % abs value for finite ones; Re/Im for infinite ones
  nmlen = [abs(nmlen(~cmplx));real(nmlen(cmplx));imag(nmlen(cmplx))];
  % first entry is useless (=1)
  nmlen(1) = [];
  
  % Set up initial guess
  if isempty(x0)
    y0 = zeros(n-3,1);
  else
    x0 = x0(:);
    x0 = (x0-x0(2))/(x0(2)-x0(1));
    y0 = log(diff(x0(2:n-1)));
  end

  % Solve nonlinear system of equations:

  % package data
  nrow = max([n-1,nqpts,4]);
  ncol = 6+2*n;
  fdat = zeros(nrow,ncol);
  fdat(1:4,1) = [n;length(left);nqpts;ncol];
  fdat(1:n-1,2) = beta(1:n-1);
  fdat(1:n-3,3) = nmlen(:);
  fdat(1:fdat(2,1),4:6) = [left(:),right(:),cmplx(:)];
  fdat(1:nqpts,7:ncol) = qdat;
  % set options
  opt = zeros(16,1);
  opt(1) = 2*trace;
  opt(6) = 100*(n-3);
  opt(8) = tol;
  opt(9) = tol/10;
  opt(12) = nqpts;
  % do it
  [y,termcode] = nesolve('hppfun',y0,opt,fdat);
  if termcode~=1
    disp('Warning: Nonlinear equations solver did not terminate normally')
  end

  % Convert y values to x
  x = [-1;cumsum([0;exp(y)]);Inf];
end

% Determine multiplicative constant
mid = mean(x(1:2));
g = hpquad(x(2),mid,2,x(1:n-1),beta(1:n-1),qdat) -...
    hpquad(x(1),mid,1,x(1:n-1),beta(1:n-1),qdat);
c = (w(1)-w(2))/g;

