function [z,c,qdat] = dparam(w,beta,z0,options);
%DPARAM Schwarz-Christoffel disk parameter problem.
%	[Z,C,QDAT] = DPARAM(W,BETA) solves the Schwarz-Christoffel
%	mapping parameter problem with the disk as fundamental domain
%	and the polygon specified by W as the target.  W must be a
%	vector of the vertices of the polygon, specified in
%	counterclockwise order, and BETA should be a vector of the
%	turning angles of the polygon; see SCANGLE for details.  If
%	successful, DPARAM will return Z, a vector of the pre-images of
%	W; C, the multiplicative constant of the conformal map; and
%	QDAT, a matrix of quadrature data used by some of the other
%	S-C routines.
%
%       [Z,C,QDAT] = DPARAM(W,BETA,Z0) uses Z0 as an initial guess for
%       Z.
%
%       [Z,C,QDAT] = DPARAM(W,BETA,Z0,OPTIONS) uses a vector of control
%       parameters.  See SCPARMOPT.
%	
%	See also SCPARMOPT, DRAWPOLY, DDISP, DPLOT, DMAP, DINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/26/95.

n = length(w);				% no. of vertices
w = w(:);
beta = beta(:);

% Set up defaults for missing args
if nargin < 4
  options = [];
  if nargin < 3
    z0 = [];
  end
end

err = sccheck('d',w,beta);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

[trace,tol] = scparmopt(options);
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta,nqpts); 		% quadrature data

atinf = (beta <= -1);

if n==3
  % Trivial solution
  z = [-i;(1-i)/sqrt(2);1];

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
  if isempty(z0)
    y0 = zeros(n-3,1);
  else
    z0 = z0(:)./abs(z0(:));
    % Moebius to make th(n-2:n)=[1,1.5,2]*pi;
    Am = moebius(z0(n-2:n),[-1;-i;1]);
    z0 = (Am(1)*z0+Am(2))./(Am(3)*z0+Am(4));
    th = angle(z0);
    th(th<=0) = th(th<=0) + 2*pi;
    dt = diff([0;th(1:n-2)]);
    y0 = log(dt(1:n-3)./dt(2:n-2));
  end
  
  % Solve nonlinear system of equations:
  
  % package data
  nrow = max([n,nqpts,4]);
  ncol = 6+2*(n+1);
  fdat = zeros(nrow,ncol);
  fdat(1:4,1) = [n;length(left);nqpts;ncol];
  fdat(1:n,2) = beta;
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
  [y,termcode] = nesolve('dpfun',y0,opt,fdat);
  if termcode~=1
    disp('Warning: Nonlinear equations solver did not terminate normally')
  end

  % Convert y values to z
  cs = cumsum(cumprod([1;exp(-y)]));
  theta = pi*cs(1:n-3)/cs(n-2);
  z = ones(n,1);
  z([1:n-3]) = exp(i*theta);
  z(n-2:n-1) = [-1;-i];
end

% Determine scaling constant
mid = (z(1)+z(2))/2;
c = (w(1) - w(2))/...
    (dquad(z(2),mid,2,z,beta,qdat) - dquad(z(1),mid,1,z,beta,qdat));


