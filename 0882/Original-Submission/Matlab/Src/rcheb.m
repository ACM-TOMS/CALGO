function [xi,lambda,err] = rcheb(alpha,w)

%RCHEB Rational Gauss-Chebyshev quadrature and interpolation
%   [XI,LAMBDA] = RCHEB(ALPHA,W) computes the nodes XI and weights
%   LAMBDA of the rational Gauss-Chebyshev quadrature formula with poles
%   in ALPHA, corresponding to the Chebyshev weight function W. The
%   poles must be outside the interval [-1,1] and may be complex.
%
%   W is optional (default W == 1) and indicates the weight function:
%    W == 1 ==> w(x) = (1-x^2)^(-1/2)
%    W == 2 ==> w(x) = ((1-x)/(1+x))^(1/2)
%    W == 3 ==> w(x) = (1-x^2)^(1/2)
%
%   XI = RCHEB(ALPHA,W) only returns the nodes, which are almost optimal
%   for rational interpolation with prescribed poles.
%
%   [XI,LAMBDA,ERR] = RCHEB(ALPHA,W) also returns an accuracy estimate
%   for each of the nodes in XI. This is only useful for verification
%   purposes in case some or all poles are extremely close to the
%   interval [-1,1].

% -------------------------------------------------------------------------
% Authors: Joris Van Deun & Karl Deckers & Adhemar Bultheel
%
% Reference: Algorithm XXX: Near best fixed pole rational interpolation
%            with applications in spectral methods
%
% Software revision date: September 10, 2007
% -------------------------------------------------------------------------

% constants
pinf = 5 / eps;	% poles >= pinf considered equal to inf
ptol = 5 * eps; % poles closer than ptol apart considered equal
xtol = 50 * eps; % tolerance on accuracy of zeros

% check input
if nargin < 1,
    error('No input.');
elseif nargin > 2,
    error('Too many input arguments.');
end

if isempty(alpha),
    error('List of poles is empty.');
elseif any(isnan([real(alpha),imag(alpha)])),
    error('NaN is not a valid pole.');
end
alpha = alpha(:).';
n = length(alpha);
alpha = real(alpha) + sqrt(-1) * abs(imag(alpha)); % to upper half plane
if any((imag(alpha) <= ptol) & (abs(real(alpha)) <= 1 + ptol)),
    error('Poles too close to the interval [-1,1].');
end

if nargin < 2,
    w = [0 0];
elseif w == 1 | w == 2 | w == 3,
    w = [floor(w / 2), w - 1];
else
    error('Wrong weight.');
end

% initialise output variables
xi = zeros(1,n);
lambda = xi;

% process poles
alpha(abs(alpha) >= pinf) = inf;

an = alpha(end);
if isinf(an),
    bn = 0;
else
    bn = real(1 / (an + (2 * (real(an)>=0) - 1) * sqrt(an^2 - 1)));  
    if abs(bn) <= (1 / (2 * pinf)),
	bn = 0;
	alpha(end) = inf;
    else
	alpha(end) = (bn + 1 / bn) / 2;
    end
end
			 
% the case n = 1 is explicitly known
if n == 1, 
    xi = (bn + w(1,2) - 2 * w(1,1)) / (1 + w(1,1));
    if nargout > 1,
	lambda = pi / (1 + w(1,2) - w(1,1));
    end
    if nargout > 2,
	err = eps; % no digits lost in computing bn
    end
    return
end

% determine the number of distinct poles and sort
a = psort([alpha, alpha(1:n-1)], ptol);
m = size(a,2);
	
% points for pchip interpolation
ttk = pi * ([n:-1:1] - (1 - w(1,1)) / 2) / (n + w(1,2)/2);            

if m == 1, % all poles equal 
    if bn == 0, % polynomial case (all poles equal to infinity)
	xi = cos(ttk);
	if nargout > 1,
	    lambda = 2 * pi * (1 - w(1,1) * xi.^w(1,2)) ./ (2 * n + w(1,2));
	end
	if nargout > 2,
	    err = eps * ones(1,n); % no digits lost in this case
	end
	return
    end
    beta = [bn; n];
else % not all poles equal
    beta = 1 ./ (a(1,:) + (2*(real(a(1,:))>=0) - 1) .* sqrt(a(1,:).^2 - 1));
    if isinf(a(1,end)),
	beta(end) = 0;
    end
    beta = [beta; a(2,:)/2];
end

% compute abs and arg of beta for later use
r = abs(beta(1,:)); phi = angle(beta(1,:));
b = (1-r) .* exp(sqrt(-1)*phi);
beta = [beta; r; phi; b]; 

% construct and evaluate piecewise cubic hermite interpolating polynomial
%  to obtain initial values for the nodes
f = ceval(ttk,n,beta,w,[],[0,0]);
ttk = [pi, ttk, 0]; % values
f = [(n + w(1,2) / 2) * pi, f ,0]; % abscissae
tk0 = pchip(f, ttk, (pi * ([n:-1:1] - (1 - w(1,1)) / 2)));

% compute the nodes using newton-raphson iteration
bnd = repmat([0; pi],1,n);
[tk,k1,bnd] = newton(tk0,bnd,n,beta,w,[n:-1:1],xtol);

% find additional initial values if divergence occurs and iterate
if ~isempty(k1),
    tk0 = newinit(n,beta,w,k1,ttk,f,ptol); % new initial values
    if ~isempty(tk0),
	[tk(n-k1+1),k2,bnd] = newton(tk0(1,:),bnd,n,beta,w,k1,xtol);
	% use third set of initial values and/or bisection if still
	% diverges
	if ~isempty(k2),
	    if size(tk0,1) == 2, % third set available
		i = 1; j = 1; s = zeros(size(k2));
		% find indices of nodes which have not converged yet
		while i <= length(k2),
		    if k1(j) == k2(i),
			s(i) = j;
			i = i + 1; 
		    end
		    j = j + 1;
		end
		[tk(n-k2+1),k3,bnd] = newton(tk0(2,s),bnd,n,beta,w,k2,xtol);
		if ~isempty(k3), % still diverges, use bisection
		    bnd = bound(bnd,n,k3,tk,ttk,f);
		    tk(n-k3+1) = bisect(bnd,n,beta,w,k3,xtol);
		end
	    else % no third set available, use bisection
		bnd = bound(bnd,n,k2,tk,ttk,f);
		tk(n-k2+1) = bisect(bnd,n,beta,w,k2,xtol);
	    end
	end
    else % use bisection if no new initial values could be found
	bnd = bound(bnd,n,k1,tk,ttk,f);
	tk(n-k1+1) = bisect(bnd,n,beta,w,k1,xtol);
    end
end

% finally compute the zeros xi from tk
xi = cos(tk);

% return the quadrature weights if requested
if nargout > 1,
    lambda = weight(xi,tk,n,beta,w);
end

% return accuracy estimate if requested
if nargout > 2,
    err = ceval(tk,n,beta,w,[n:-1:1],[1,1]);
end

% -------------------------------------------------------------------------

function a = psort(a,tol)
% sort values in a and consider equal if less than tol apart
%  if a has two rows, the values in the second row are multiplicities
%  (this only occurs when psort is called from newinit)

[p,n] = size(a);

if n == 1,
    at = [a(1);1];
    a2 = []; i = 1;
else
    [a1,i] = sort(a(1,:));
    warning('off'); % division by zero is allowed in the next line
    a2 = find(abs(diff(a1)./a1(1:n-1))>tol);
    warning('on');
    m = length(a2);
    if m == 0,
        at = [a1(end);n];
    else
        a3 = diff(a2); % multiplicities
        at = [a1(a2) a1(end);a2(1) a3 n-a2(end)];
    end
end

if p == 2,
    b = cumsum(a(2,i) - 1);
    b = [b(a2), b(end)];
    b = [b(1), diff(b)];
    at(2,:) = at(2,:) + b;
end

a = at;


% -------------------------------------------------------------------------

function f = ceval(theta,n,beta,c,k,e)
% evaluate F(theta) if e = [0,0]
%  or dF(theta)/dtheta if e = [0,1]
%  or F(theta) - (k-d)*pi if e = [1,0]
%  or (F(theta) - (k-d)*pi)/(dF(theta)/dtheta) if e = [1,1]
%  where d = (1-c(1))/2

m = size(beta,2); % number of distinct poles
mb = beta(2,:)'; phi = beta(4,:)'; b = beta(5,:).';
n2 = length(theta); % number of evaluation points
t = ((m > n2) | (10 * m > (n2 + 780))) + 2 * ((m <= 108) & ... 
    (n2 <= 300) & (30 * m > (n2 + 120)) & (10 * n2 > (m + 36)));

if ((e(1) == 1) | (e(2) == 0)),
    f = -(n - 1 - c(2)/2) * theta; % depends on the weight function
    if (e(1) == 1)
            f = f - (pi * (k - (1 - c(1))/2));
    end
    if (m == 1) % number of different poles == 1
        f = f + (2 * n - 1) * atan2(sin(theta),cos(theta)-beta(1));
    elseif (n2 == 1) % number of zeros to compute == 1
        f = f + ftheta(theta,phi,b,mb);
    else
        switch t % most efficient computation based on t
            case 0
                for j = 1:n2
                    f(j) = f(j) + ftheta(theta(j),phi,b,mb);
                end
            case 1
                for j = m:-1:1
                    f = f + ftheta(theta,phi(j),b(j),mb(j));
                end
            otherwise
                thetam = repmat(theta,m,1);
                bn = repmat(b,1,n2); phin = repmat(phi,1,n2);
                mbn = repmat(mb,1,n2);
                f = f + ftheta(thetam,phin,bn,mbn);
        end
    end
end

if (e(2) == 1), % derivative also needed
    df = -(n - 1 - c(2)/2) * ones(1,n2);
    if (m == 1) % number of different poles == 1
        df = df + (2 * n - 1) * (1 - beta(1) * cos(theta)) ./ ... 
	    (1 - 2 * beta(1) * cos(theta) + beta(1)^2);
    elseif (n2 == 1) % number of zeros to compute == 1
        df = df + fdtheta(theta,beta(1,:).',mb);
    else
        switch t
            case 0
                for j = 1:n2
                    df(j) = df(j) + fdtheta(theta(j),beta(1,:).',mb);
                end
            case 1
                for j = m:-1:1
                    df = df + fdtheta(theta,beta(1,j),mb(j));
                end
            otherwise
                if (e(1)==0)
                    thetam = repmat(theta,m,1);
                    mbn = repmat(mb,1,n2);
                end
		betan = repmat(beta(1,:).',1,n2);
                df = df + fdtheta(thetam,betan,mbn);
        end
    end
    if (e(1) == 1)
        f = f ./ df;
    else
        f = df;
    end
end

% -------------------------------------------------------------------------

function ft = ftheta(theta,phi,b,m)
% compute arg(z-beta) + arg(z - conj(beta)) where z = exp(i*theta)
%  and using the right branch cut

s = imag(b); c = real(b);
sin1 = sin((theta+phi)/2); sin2 = sin((theta-phi)/2);
y1 = 2*sin2.*cos((theta+phi)/2) + s;
y2 = 2*sin1.*cos((theta-phi)/2) - s;
x = -2*sin1.*sin2 + c;
ft = m .* (atan2(y1,x) + atan2(y2,x) + 2 * pi * ((x < 0) & (y2 < 0)));
if (size(b,1) > 1),
    ft = sum(ft);
end

% -------------------------------------------------------------------------

function fdt = fdtheta(theta,beta,m)
% derivative of ft

z = exp(sqrt(-1) * theta);
fdt = m .* real(z./(z-beta) + z./(z-conj(beta)));
if (size(beta,1) > 1),
    fdt = sum(fdt);
end

% -------------------------------------------------------------------------

function [xi,k,bnd] = newton(xi,bnd,n,beta,w,k,xtol)
% perform Newton-Raphson iterations starting from initial values in xi
%  and confidence interval in bnd

maxiter = 10;

% initialisation parameters
j = 0; % number of iterations
c = 0; % keep track of nodes which have converged
n2 = length(xi); % number of nodes remaining
xc = [];
kc = [];
cc = [];
bndc = [];

% iteration
while any(~c) & j < maxiter & n2 > 0, 
    j = j + 1;
    corr = ceval(xi,n,beta,w,k,[1,1]);
    xo = xi;
    xi = xi - corr;
    
    c = (abs(corr) <= xtol);
    
    % shrink confidence interval
    i = find(corr >= 0);
    bnd(2,i) = xo(i); % upper bound
    i = find(corr <= 0);
    bnd(1,i) = xo(i); % lower bound
    
    % check for diverging nodes
    conv = (xi > bnd(1,:) & xi < bnd(2,:)) | (xi == xo);
    if (n2 > 1)
        e = conv & ~c;
	% converged or diverging
        i = find(~e);
        xc = [xc, xi(i)]; kc = [kc, k(i)]; cc = [cc, c(i) & conv(i)];
	bndc = [bndc, bnd(:,i)];
	% continue with remaining (converging) nodes
        i = find(e);
        xi = xi(i); k = k(i); c = c(i);
	bnd = bnd(:,i);
        n2 = length(xi);
    elseif conv == 0               
        j = 2 * maxiter; % divergence
        c = 0;
    end
end

xi = [xi,xc]; c = [c,cc]; k = [k,kc];
bnd = [bnd, bndc];

% put back in order
% [k,j] = sort(k,2,'descend'); % use this in recent versions of Matlab
[k,j] = sort(k); k = k(end:-1:1); j = j(end:-1:1); % and not this
xi = xi(j); c = c(j);
bnd = bnd(:,j);

% nodes which have not converged
i = find(~c);
k = k(i);
bnd = bnd(:,i);

% -------------------------------------------------------------------------

function tk0 = newinit(n,beta,w,k,xi,f,tol)
% obtain new initial values

% find poles close to the boundary 
b = beta(:,3 * beta(3,:) >= 1);
if isempty(b),
    tk0 = [];
    return;
end

% estimate maxima of d(F_n)/d(theta)
t(1,:) = real(b(1,:)) .* ((1 + b(3,:).^2).^2 + ...
    4*imag(b(1,:)).^2) ./ (2 * b(3,:) .* (b(3,:) .* ...
    (1 + b(3,:).^2) + abs(imag(b(1,:))) .* sqrt( ...
    (1 + b(3,:).^2).^2 - 4*real(b(1,:)).^2)));
t(1,t(1,:) > 1) = 1; t(1,t(1,:) < -1) = -1;
t(1,:) = acos(t(1,:)); % estimates 
t = psort([t; b(2,:)],tol); % remove duplicates but keep multiplicities
t = t(:,end:-1:1); % store in reverse order (like the nodes)

dft = ceval(t(1,:),n,beta,w,[],[0,1]); % d(F_n)/d(theta)
df = diff(f); % f contains F_n evaluated in the interpolation points xi
ds = df./diff(xi);
% for each local maximum determine whether it satisfies the assumption
%  of 'peak shaped' and only retain it if it does
u = []; v = []; taux = []; dftaux = [];
for j = 1:size(t,2),
    % [xi(j2),xi(j1)] is smallest interval containing t(1,j)
    % j1 = find(xi >= t(1,j),1,'last'); % use this in recent versions
    % j2 = find(xi <= t(1,j),1); % of Matlab instead of the following
    j1 = find(xi >= t(1,j)); j1 = j1(end);
    j2 = find(xi <= t(1,j)); j2 = j2(1);
    if j1 ~= j2,
	s = dft(j)/ds(j1);
    elseif j2 == 1,
	s = dft(j)/ds(1);
    elseif j2 == n + 2,
	s = dft(j)/ds(n+1);
    else
	s = dft(j)*max(1/ds(j2),1/ds(j2-1));
    end
    if s > 2, % keep only peak shaped local maxima in t
	u = [u, j1]; v = [v, j2];
	taux = [taux, t(:,j)]; dftaux = [dftaux, dft(j)];
    end
end        
mt = length(u); % final number of local maxima we will use
if mt == 0,
    tk0 = [];
    return
else
    t = taux; dft = dftaux;
end

c = (k - (1 - w(1))/2) * pi;
ft = ceval(t(1,:),n,beta,w,[],[0,0]);
if mt == 1,
    tk0 = t(1,:) + (c - ft) / dft; % perform one newton step already
else
    i = find(c <= ft(end));
    tk0(i) = t(1,end) + (c(i) - ft(end)) ./ dft(end);
    i = find(c >= ft(1));
    tk0(i) = t(1,1) + (c(i) - ft(1)) ./ dft(1);
    i = find((c > ft(end)) & (c < ft(1)));
    if isempty(i),
	return
    end
    for j = 1:length(i),
	% p = find(ft >= c(i(j)),1,'last'); % use this in recent Matlab
	p = find(ft >= c(i(j))); p = p(end);
	theta = (t(1,p) + t(1,p+1)) / 2 + pi/4 * (t(2,p)/dft(p) - ...
	    t(2,p+1)/dft(p+1));
	ftheta = ceval(theta,n,beta,w,[],[0,0]);
	if c(i(j)) >= ftheta,
	    itk(:,j) = p + [0; 1];
	else
	    itk(:,j) = p + [1; 0];
	end
    end
    tk0(i) = t(1,itk(1,:)) + (c(i) - ft(itk(1,:))) ./ dft(itk(1,:));
    tk0 = [tk0; tk0];
    tk0(2,i) = t(1,itk(2,:)) + (c(i) - ft(itk(2,:))) ./ dft(itk(2,:));
end

% -------------------------------------------------------------------------

function bnd = bound(bnd,n,k,tk,ttk,f)
% improve upper and lower bounds for nodes that have not converged
%  based on those that have

m = length(k);

tk = [pi tk(1,:) 0];
kd = -diff(k);
i = find(kd > 1);

ku = zeros(1,m);
kl = ku;

ku(1) = k(1) + 1;
ku(1+i) = k(1+i) + 1;
for j = 2:m,
    if ku(j) == 0,
	ku(j) = ku(j-1);
    end
end

kl(m) = k(m) - 1;
kl(i) = k(i) - 1;
for j = m-1:-1:1,
    if kl(j) == 0,
	kl(j) = kl(j+1);
    end
end

bnd(1,:) = max(bnd(1,:),tk(n-kl+2));
bnd(2,:) = min(bnd(2,:),tk(n-ku+2));

% -------------------------------------------------------------------------

function xi = bisect(bnd,n,beta,w,k,xtol)
% use bisection to compute the nodes in the confidence intervals in bnd

maxiter = 52;
l = bnd(1,:); u = bnd(2,:);

% initialisation
fl = ceval(l,n,beta,w,k,[1,0]);
fu = ceval(u,n,beta,w,k,[1,0]);

d = u - l;

i = find(d <= xtol); % nodes which are already accurate enough
uc = u(i); lc = l(i); kc = k(i);

i = find(d > xtol); % need bisection to improve accuracy
yi = l(i); n2 = length(i);

% iteration
j = 0;
while j < maxiter & n2 > 0,  
    j = j + 1;
    l = u(i); fl = fu(i);
    k = k(i);
    u = l - d(i) / 2; % bisection
    fu = ceval(u,n,beta,w,k,[1,0]);
    i = find((fl < 0 & fu > 0) | (fl > 0 & fu < 0));
    yi(i) = l(i);
    
    % checking the criterium of accuracy
    d = u - yi;

    i = find(abs(d) <= xtol); % accurate enough
    uc = [uc,u(i)]; lc = [lc,yi(i)]; kc = [kc,k(i)];

    i = find(abs(d) > xtol); % accuracy not reached yet
    yi = yi(i); n2 = length(i);
end

% output
uc = (uc + lc) / 2;
u(i) = (u(i) + yi) / 2;
u = [uc,u(i)]; k = [kc,k(i)];
% [k,j] = sort(k,2,'descend'); % use this in recent versions of Matlab
[k,j] = sort(k); k = k(end:-1:1); j = j(end:-1:1);
xi = u(j);

% -------------------------------------------------------------------------

function lambda = weight(xi,tk,n,beta,w)
% compute the quadrature weights

m = size(beta,2);
mb = beta(2,:); r = beta(3,:); phi = beta(4,:);

lambda = (1 + w(1,2)) * ones(1,n);
if (m == 1)
    lambda = lambda + (2*n - 1)*(1 + r)./((1 - r) + ...
        4*r*sin((tk-phi)/2).^2/(1-r));
else
    t = (10 * m > (n + 780)) + 2 * ((m <= 108) & (n <= 300) & ...
        (30 * m > (n + 120))); % parameter which determines efficiency
    switch t
        case 0
            for j = 1:n
                lambda(j) = lambda(j) + sum(mb.*(1 + r).* ...
		(1./((1 - r) + 4*r.*sin((tk(j)-phi)/2).^2./(1-r)) + ...
		1 ./ ((1 - r) + 4*r.*sin((tk(j)+phi)/2).^2./(1-r))));
            end
        case 1
            for k = 1:1:m
                lambda = lambda + mb(k)*(1 + r(k)) * ...
		(1./((1-r(k))+4*r(k)*sin((tk-phi(k))/2).^2/(1-r(k))) ...
	       	+ 1./((1-r(k))+4*r(k)*sin((tk+phi(k))/2).^2/(1-r(k))));
            end
        otherwise
            tkm = repmat(tk,m,1);
            rn = repmat(r',1,n); mbn = repmat(mb',1,n);
	    phin = repmat(phi',1,n);
            lambda = lambda + sum(mbn.*(1 + rn) .* (1 ./ ...
	        ((1 - rn)+4*rn.*sin((tkm-phin)/2).^2./(1-rn)) + 1./ ...
		((1 - rn)+4*rn.*sin((tkm+phin)/2).^2./(1-rn))));
    end
end
lambda = 2 * pi * (1 - w(1,1) * xi.^w(1,2)) ./ lambda;


% END ---------------------------------------------------------------------

