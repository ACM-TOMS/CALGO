function [z,c,L,qdat] = rparam(w,beta,cnr,z0,options);
%RPARAM Schwarz-Christoffel rectangle parameter problem.
%       [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS) solves the
%       Schwarz-Christoffel parameter problem with a rectangle as
%       fundamental domain and interior of the specified polygon as the
%       target.  W must be a vector of the vertices of the polygon,
%       specified in counterclockwise order.  BETA is a vector of
%       turning angles; see SCANGLES.  CORNERS is a 4-component vector
%       specifying the indices of the vertices which are the images of
%       the corners of the rectangle.  *BE SURE* the first two entries
%       describe the LONG sides of the rectangle, and go in
%       counterclockwise order.  If CORNERS is omitted, the user is
%       requested to select these vertices using the mouse.
%
%       If successful, RPARAM will return Z, a vector of the
%       prevertices; C, the multiplicative constant of the conformal
%       map; L, a parameter related to aspect ratio; and QDAT, a matrix
%       of quadrature data used by some of the other SC routines.
%
%       [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,Z0) uses Z0 as an initial
%       guess for Z.  In this case, Z0 represents the image of
%       prevertices on the strip 0 <= Im z <= 1.  You can use R2STRIP to
%       transform prevertices from the rectangle to the strip.
%
%       [Z,C,L,QDAT] = RPARAM(W,BETA,CORNERS,Z0,OPTIONS) uses a vector of
%       control parameters.  See SCPARMOPT.  
%	
%	See also SCPARMOPT, DRAWPOLY, RDISP, RPLOT, RMAP, RINVMAP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w); 				% no. of vertices
w = w(:);
beta = beta(:);

% Set up defaults for missing args
if nargin < 5
  options = [];
  if nargin < 4
    z0 = [];
    if nargin < 3
      cnr = [];
    end
  end
end

if isempty(cnr)
  disp('Use mouse to select images of rectangle corners.')
  disp('Go in counterclockwise order and select a long edge first.')
  figure(gcf)
  cnr = scselect(w,beta,4);
end

% Renumber vertices so that cnr(1)=1
renum = [cnr(1):n,1:cnr(1)-1];
w = w(renum);
beta = beta(renum);
cnr = rem(cnr-cnr(1)+1+n-1,n)+1;

[trace,tol] = scparmopt(options);
nqpts = max(ceil(-log10(tol)),4);
qdat = scqdata(beta,nqpts); 		% quadrature data

% Check input data.
err = sccheck('r',w,beta,cnr);
if err==1
  fprintf('Use SCFIX to make polygon obey requirements\n')
  error(' ')
end

atinf = (beta <= -1);

if isempty(z0)
  % Try to find a reasonable initial guess.
  dw = abs(diff(w([1:n,1])));		% side lengths
  dw(isinf(dw)) = mean(dw(~isinf(dw)))*ones(size(dw(isinf(dw))));
  % Estimate length and width, and thus conformal modulus
  len = mean([sum(dw(cnr(1):cnr(2)-1)), sum(dw(cnr(3):cnr(4)-1))]);
  wid = mean([sum(dw(cnr(2):cnr(3)-1)), sum(dw([cnr(4):n,1:cnr(1)-1]))]);
  modest = len/wid;
  % Evenly space prevertices to match this conformal modulus
  z0(cnr(1):cnr(2)) = linspace(0,modest,diff(cnr(1:2))+1);
  dx = z0(cnr(1)+1)-z0(cnr(1));
  z0(cnr(1)-1:-1:1) = z0(cnr(1))-dx*(1:cnr(1)-1);
  z0(cnr(2)+1:cnr(3)-1) = z0(cnr(2)) + dx*(1:diff(cnr(2:3))-1);
  z0(cnr(4):-1:cnr(3)) = linspace(0,modest,diff(cnr(3:4))+1);
  dx = z0(cnr(4)-1)-z0(cnr(4));
  z0(cnr(4)+1:n) = z0(cnr(4))-dx*(1:n-cnr(4));

else
  if length(z0)~=n
    error('Initial guess has wrong number of prevertices')
  end
  z0 = z0(renum);
  if any(imag(z0(1:cnr(3)-1))) | any(~imag(z0(cnr(3):n)))
    error('Initial guess has prevertices on wrong side of strip')
  end
end

% Convert z0 to unconstrained vars
y0 = zeros(n-3,1);
dz = diff(z0);
dz(cnr(3):n-1) = -dz(cnr(3):n-1);
y0(1:cnr(2)-2) = log(dz(1:cnr(2)-2));
y0(cnr(2)-1) = mean(log(dz([cnr(2)-1,cnr(3)])));
y0(cnr(2):n-3) = log(dz([cnr(2):cnr(3)-2,cnr(3)+1:n-1]));

% Find prevertices (solve param problem)

% Set up normalized lengths for nonlinear equations:
% indices of left and right integration endpoints
left = 1:n-2;				
right = 2:n-1;				
% delete indices corresponding to vertices at Inf
left(find(atinf)) = [];
right(find(atinf) - 1) = [];
if atinf(n-1)
  right = [right,n];
end
cmplx = ((right-left) == 2);
% normalize lengths by w(2)-w(1)
nmlen = (w(right)-w(left))/(w(2)-w(1));
% abs value for finite ones
nmlen(~cmplx) = abs(nmlen(~cmplx));
% first entry is useless (=1)
nmlen(1) = [];

beta = [0;beta(1:cnr(3)-1);0;beta(cnr(3):n)];

% Solve nonlinear system of equations:
% package data
nrow = max([n+2,nqpts,8]);
ncol = 6+2*n+2;
fdat = zeros(nrow,ncol);
fdat(1:4,1) = [n;length(left);nqpts;ncol];
fdat(5:8,1) = cnr(:);
fdat(1:n+2,2) = beta;
fdat(1:fdat(2,1)-1,3) = nmlen(:);
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
[y,termcode] = nesolve('rpfun',y0,opt,fdat);
if termcode~=1
  disp('Warning: Nonlinear equations solver did not terminate normally')
end

% Convert y values to z on strip
zs = rptrnsfm(y,cnr);
nb = cnr(3)-1;

% Determine multiplicative constant
mid = mean(zs(1:2));
g = stquad(zs(2),mid,2,zs,beta,qdat) -...
    stquad(zs(1),mid,1,zs,beta,qdat);
c = (w(1)-w(2))/g;

% Find corners of rectangle
L = zs(cnr(2))-zs(cnr(1));
if L < 5.9
  m = exp(-2*pi*L);
  K = ellipke(m);
  Kp = ellipke(1-m);
else
  K = pi/2;
  Kp = pi*L + log(4);
end
rect = [K;K+i*Kp;-K+i*Kp;-K];
bounds = [-K,K,0,Kp];

% Find prevertices on the rectangle:
%   initial values evenly spaced on the rectangle
z = zeros(size(zs));
z(cnr) = rect;
z(1:cnr(2)) = linspace(rect(1),rect(2),diff(cnr(1:2))+1).';
tmp = linspace(rect(2),i*imag(rect(3)),diff(cnr(2:3))+1).';
z(cnr(2):cnr(3)-1) = tmp(1:diff(cnr(2:3)));
z(cnr(3):cnr(4)) = linspace(rect(3),rect(4),diff(cnr(3:4))+1).';
tmp = linspace(rect(4),0,n-cnr(4)+2).';
z(cnr(4):n) = tmp(1:n-cnr(4)+1);
zn = z(:);

% Which are on left/right sides?
lr = zeros(n,1);
lr([1:cnr(2),cnr(3):cnr(4)]) = ones(cnr(2)+cnr(4)-cnr(3)+1,1);

%   Newton iteration
maxiter = 50;
done = zeros(size(zn));
done(cnr) = ones(4,1);
k = 0;
while ~all(done) & k < maxiter
  [F,dF] = r2strip(zn(~done),z(cnr),L);
  F = zs(~done) - F;
  % Must keep points from leaving the rectangle:
  %   pure real/imaginary, and not too big
  step = 2*F./dF;
  step(lr(~done)) = i*imag(step(lr(~done)));
  step(~lr(~done)) = real(step(~lr(~done)));
  bad = ones(size(step));
  while any(bad)
    step(bad) = step(bad)/2;
    znew = zn(~done) + step;
    bad = real(znew) < bounds(1) | real(znew) > bounds(2) |...
	imag(znew) < bounds(3) | imag(znew) > bounds(4);
  end
  % Update
  zn(~done) = znew;
  done(~done) =  (abs(F) < tol);
  k = k + 1;
end
if any(abs(F)> tol)
  disp('Warning in rparam: Iteration for rectangle prevertices DNC')
end
z(:) = zn; 

% Undo renumbering
z(renum) = z;



