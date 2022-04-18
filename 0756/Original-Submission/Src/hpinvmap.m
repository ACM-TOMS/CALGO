function zp = hpinvmap(wp,w,beta,x,c,qdat,z0,options)
%HPINVMAP Schwarz-Christoffel half-plane inverse map.
%	HPINVMAP(WP,W,BETA,X,C,QDAT) computes the inverse of the
%	Schwarz-Christoffel half-plane map (i.e., from the polygon
%	to the upper half-plane ) at the points given in vector WP. The
%	other arguments are as in HPPARAM.  QDAT may be omitted.
%
%	The default algorithm is to solve an ODE in order to obtain a fair
%	approximation for ZP, and then improve ZP with Newton iterations.
%	The ODE solution at WP requires a vector Z0 whose forward image W0
%	is such that for each j, the line segment connecting WP(j) and W0(j)
%	lies inside the polygon.  By default Z0 is chosen by a fairly robust
%	automatic process.  Using a parameter (see below), you can choose to
%	use either an ODE solution or Newton iterations exclusively.
%
%	HPINVMAP(WP,W,BETA,X,C,QDAT,Z0) has two interpretations.  If the ODE
%	solution is being used, Z0 overrides the automatic selection of
%	initial points.  (This can be handy in convex polygons, where the
%	choice of Z0 is trivial.)  Otherwise, Z0 is taken as an initial
%	guess to ZP.  In either case, if length(Z0)==1, the value Z0 is used
%	for all elements of WP; otherwise, length(Z0) should equal
%	length(WP).
%
%       HPINVMAP(WP,W,BETA,X,C,QDAT,Z0,OPTIONS) uses a vector of parameters
%       that control the algorithm.  See SCIMAPOPT.
%
%	See also SCIMAPOPT, HPPARAM, HPMAP.
%
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
w = w(:);
beta = beta(:);
x = x(:);
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

[ode,newton,tol,maxiter] = scimapopt(options);

nfin = n - isinf(x(n));
if isempty(qdat)
  qdat = scqdata(beta(1:nfin),max(ceil(-log10(tol)),2));
end

% ODE
if ode
  if isempty(z0)
    % Pick a value z0 (not a singularity) and compute the map there.
    [z0,w0] = scimapz0('hp',wp,w,beta,x,c,qdat);
  else
    w0 = hpmap(z0,w,beta,x,c,qdat);
    if length(z0)==1 & lenwp > 1
      z0 = z0(:,ones(lenwp,1)).';
      w0 = w0(:,ones(lenwp,1)).';
    end
  end

  % Use relaxed ODE tol if improving with Newton.
  odetol = max(tol,1e-3*(newton));

  % Set up data for the ode function.
  global HPIMDATA	
  HPIMDATA = (wp - w0(:))/c; 		% adjusts "time" interval
  HPIMDATA(1:nfin,2:3) = [x(1:nfin), beta(1:nfin)];
  HPIMDATA(1,4) = nfin;

  z0 = [real(z0);imag(z0)];
  [t,y] = ode45('hpimapf1',0,1,z0,odetol);
  [m,leny] = size(y);
  zp(:) = y(m,1:lenwp)+sqrt(-1)*y(m,lenwp+1:leny);
end

% Newton iterations
if newton
  if ~ode
    zn = z0(:);
    if length(z0)==1 & lenwp > 1
      zn = zn(:,ones(lenwp,1));
    end
  else
    zn = zp(:);
  end
    
  wp = wp(:);
  done = zeros(size(zn));
  k = 0;
  while ~all(done) & k < maxiter
    F = wp(~done) - hpmap(zn(~done),w,beta,x,c,qdat);
    m = length(F);
    dF = c*exp(sum(beta(1:nfin,ones(m,1)).*...
         log(zn(~done,ones(nfin,1)).'-x(1:nfin,ones(m,1)))));
    zn(~done) = zn(~done) + F(:)./dF(:);
    done(~done) =  (abs(F) < tol);
    k = k + 1;
  end
  if any(abs(F)> tol)
    disp('Warning in hpinvmap: Solution may be inaccurate')
    fprintf('Maximum residual = %.3g\n',max(abs(F)))
  end
  zp(:) = zn; 
end;

