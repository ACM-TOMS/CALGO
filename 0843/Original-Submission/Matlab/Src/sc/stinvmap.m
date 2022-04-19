function zp = stinvmap(wp,w,beta,z,c,qdat,z0,options)
%STINVMAP Schwarz-Christoffel strip inverse map.
%   STINVMAP(WP,W,BETA,Z,C,TOL) computes the inverse of the
%   Schwarz-Christoffel strip map (i.e., from the polygon to the strip)
%   at the points given in vector WP. The other arguments are as in
%   STPARAM. TOL is a scalar tolerance, or a quadrature-data matrix QDAT
%   as returned by SCQDATA, or may be omitted.
%
%   The default algorithm is to solve an ODE in order to obtain a fair
%   approximation for ZP, and then improve ZP with Newton
%   iterations. The ODE solution at WP requires a vector Z0 whose
%   forward image W0 is such that for each j, the line segment
%   connecting WP(j) and W0(j) lies inside the polygon. By default Z0 is
%   chosen by a fairly robust automatic process. Using a parameter (see
%   below), you can choose to use either an ODE solution or Newton
%   iterations exclusively.
%
%   STINVMAP(WP,W,BETA,Z,C,TOL,Z0) has two interpretations. If the ODE
%   solution is being used, Z0 overrides the automatic selection of
%   initial points. (This can be handy in convex polygons, where the
%   choice of Z0 is trivial.) Otherwise, Z0 is taken as an initial guess
%   to ZP. In either case, if length(Z0)==1, the value Z0 is used for
%   all elements of WP; otherwise, length(Z0) should equal length(WP).
%
%   STINVMAP(WP,W,BETA,Z,C,TOL,Z0,OPTIONS) uses a vector of parameters
%   that control the algorithm. See SCINVOPT.
%
%   See also SCINVOPT, STPARAM, STMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: stinvmap.m,v 2.4 2002/04/03 15:27:05 driscoll Exp $

n = length(w);
w = w(:);
beta = beta(:);
z = z(:);
zp = zeros(size(wp));
wp = wp(:);
lenwp = length(wp);

if nargin < 8
  options = [];
  if nargin < 7
    z0 = [];
    if nargin < 6
      qdat = [];
    end
  end
end

[ode,newton,tol,maxiter] = scinvopt(options);

if isempty(qdat)
  qdat = tol;
end

if length(qdat)==1
  qdat = scqdata(beta,max(ceil(-log10(qdat)),2));
end

done = zeros(size(wp));
% First, trap all points indistinguishable from vertices, or they will cause
% trouble.
for j=1:n
  idx = find(abs(wp-w(j)) < 3*eps);
  zp(idx) = z(j)*ones(length(idx),1);
  done(idx) = ones(length(idx),1);
end
lenwp = lenwp - sum(done);
if lenwp==0, return, end

% ODE
if ode
  if isempty(z0)
    % Pick a value z0 (not a singularity) and compute the map there.
    [z0,w0] = scimapz0('st',wp(~done),w,beta,z,c,qdat);
  else
    w0 = stmap(z0,w,beta,z,c,qdat);
    if length(z0)==1 & lenwp > 1
      z0 = z0(:,ones(lenwp,1)).';
      w0 = w0(:,ones(lenwp,1)).';
    end
    w0 = w0(~done);
    z0 = z0(~done);
  end

  % Use relaxed ODE tol if improving with Newton.
  odetol = max(tol,1e-2*(newton));

  % Rescale dependent coordinate
  scale = (wp(~done) - w0(:));

  % Solve ODE
  z0 = [real(z0);imag(z0)];
  [t,y] = ode23('stimapfun',[0,0.5,1],z0,odeset('abstol',odetol),...
      scale,z,beta,c);
  [m,leny] = size(y);
  zp(~done) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
end

% Newton iterations
if newton
  if ~ode
    zn = z0(:);
    if length(z0)==1 & lenwp > 1
      zn = zn(:,ones(lenwp,1));
    end
    zn(done) = zp(done);
  else
    zn = zp(:);
  end
    
  wp = wp(:);
  k = 0;
  while ~all(done) & k < maxiter
    F = wp(~done) - stmap(zn(~done),w,beta,z,c,qdat);
    dF = c*stderiv(zn(~done),z,beta);
    zn(~done) = zn(~done) + F(:)./dF(:);
    done(~done) =  (abs(F) < tol);
    k = k + 1;
  end
  if any(abs(F)> tol)
    str = sprintf('Check solution; maximum residual = %.3g\n',max(abs(F)));
    warning(str)
  end
  zp(:) = zn; 
end;

