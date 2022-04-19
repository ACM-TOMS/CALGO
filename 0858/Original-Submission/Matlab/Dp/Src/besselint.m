function [f,err] = besselint(a,nu,m,varargin)
%BESSELINT Compute the integral of a product of Bessel functions
%   F = BESSELINT(A,NU,M) tries to compute the integral of the product
%
%       X^M * J_NU(1)(A(1)*X) * ... * J_NU(k)(A(k)*X)
%
%   for X from 0 to Inf, where A is a vector of length k whose
%   entries are strictly positive real numbers, NU is a vector
%   of length k whose entries are real numbers (or alternatively,
%   NU can be a scalar when all Bessel functions are of the same
%   order), and M is a real number. J_NU(i) is the Bessel function
%   of the first kind and order NU(i).
%
%   F = BESSELINT(A,NU,M,RELTOL) uses a relative error tolerance of
%   RELTOL instead of the default, which is 50*eps.
%
%   F = BESSELINT(A,NU,M,RELTOL,ABSTOL) also uses an absolute error
%   tolerance of ABSTOL instead of the default, which is 0. Pass an
%   empty matrix for RELTOL if this is not required. If both are
%   present, the algorithm takes whichever tolerance requires less
%   computational effort.
%
%   [F,ERR] = BESSELINT(...) also returns an error estimate in ERR.
%   This is a 1-by-2 vector whose first component contains an
%   estimate of the relative error and whose second component
%   contains an estimate of the absolute error.
%
%   If the program detects that the integral is a discontinuous
%   function of its parameters (as explained in section 2 of the
%   reference article) a warning message is given.
%
%   Occasionally the program detects that the numerical quadrature
%   rule used in the finite part integration does not reach the
%   requested precision. Then a warning message is given.
%   We encountered this mainly in cases where the exact value is zero
%   and the routine aims for maximal precision.

% ----------------------------------------------------------------------
%   Authors: Joris Van Deun & Ronald Cools,
%            Dept. of Computer Science, K.U.Leuven, Belgium.
%
%   Reference: Algorithm 858: Computing infinite range integrals of
%              an arbitrary product of Bessel functions
%              ACM Transactions on Mathematical Software, Vol. 32, Nr. 4, 2006.
%
%   Software revision date: November 6, 2006
% ----------------------------------------------------------------------

% check validity of input arguments
if nargin < 3
    error('BESSELINT requires at least three input arguments.')
elseif nargin > 5
    error('Too many input arguments.')
end

if any(a <= 0) | any(imag(a) ~= 0),
    error('Coefficients must be strictly positive.');
elseif imag(m) ~= 0,
    error('Power of x must be real.');
elseif any(imag(nu) ~= 0),
    error('Order must be real');
end

% normalize coefficients such that min(a)=1
a = a(:)';
k = length(a);
am = min(a);
a = a/am;

nu = nu(:)';
if length(nu) == 1,
    nu = nu*ones(1,k);
elseif length(nu) ~= k,
    error('Wrong length for NU');
end

% transform negative integer orders
n = find(nu < 0 & fix(nu) == nu);
nu(n) = -nu(n);
si = (-1)^sum(nu(n));

% check if function is integrable near zero
if sum(nu) + m <= -1,
    error('sum(NU)+M must be greater than -1');
end

% process input arguments to determine precision
tol = [0 0];
i = 1; v = varargin;
while ~isempty(v) & i,
    if isempty(v{end}),
        v = v(1:end-1); % remove empty last arguments
    else
        i = 0;
    end
end
i = length(v);
if i ~= 0,
    if ~isempty(varargin{1}), tol(1) = varargin{1}; end
    if i > 1, tol(2) = varargin{2}; end
end
tol = abs(tol);
if isempty(find(tol)), tol(1) = 50*eps; end
i = find(tol);
if length(i) == 2;
    [x1,c1] = gop(a,nu,m,tol(1),0); % c1 is value of cost function
    [x2,c2] = gop(a,nu,m,tol(2),1); % idem for c2
    if c1 <= c2, x = x1; i = 1; else x = x2; i = 2; end
else
    [x,dummy] = gop(a,nu,m,tol(i),i-1); % dummy not used
end

% compute infinite range approximation
[f1,x,param,s,sflag] = ira(a,nu,m,x,k,tol(i),i-1);
if sflag, % exceptional case, display warning message
    disp(sprintf(s));
    warning('Possible discontinuous case; result and error estimate may be inaccurate');
end

% perform finite range integration
[f2,err2,eflag] = fri(a,nu,m,x(1),tol);

% renormalize to get final result
f = (f1+f2)*si/am^(m+1);

% error estimate
if nargout > 1,
    abserr1 = (2/pi)^(k/2) * (param(1)*x(1)+param(2)) *...
        x(1)^(m-k/2-2*x(2)-2) / sqrt(prod(a));
    err(2) = (abs(abserr1)+err2) / am^(m+1); % renormalized!
    err(1) = err(2)/abs(f);
elseif eflag & ~sflag,
    warning('Requested precision not reached');
end

% ----------------------------------------------------------------------

function [x,chi] = gop(a,nu,m,tol,errflag)
%GOP Get optimal parameters x0 (breakpoint) and n (order)

k = length(a);
% on average, computation of incomplete gamma function takes
%  tGJ times longer than that of bessel function
tGJ = 15;
N = (15+19); % average number of points in quadrature formula
c = 2^(k+1)*tGJ*pi/(N*sum(a));

% compute LambertW(c) (need not be very accurate)
if c < 2, w1 = 0; else w1 = log(c) - log(log(c)); end
w2 = (c*exp(-w1)+w1^2) / (w1+1);
while abs(w2-w1) > abs(w2)*1e-5, % Newton-iteration
    w1 = w2;
    w2 = (c*exp(-w1)+w1^2)/(w1+1);
end
cl = w2;
if errflag, % absolute error required
    ce = k-2*m;
else
    ce = 0;
end

% determine degree n of asymptotic expansion by finding zero of
%  optimfun using Dekker-Brent algorithm
% nmin assures that assumptions in error analysis are satisfied
nmin = max([1,ceil((m-k/2)/2),abs(nu)/2]);
fmin = optimfun(nmin,a,nu,m,tol,errflag,c,cl,ce,k);
if fmin < 0,
    n = nmin; % smallest n suffices (optimfun is decreasing with n)
else
    nmax = nmin + 10;
    fmax = optimfun(nmax,a,nu,m,tol,errflag,c,cl,ce,k);
    while (fmax > 0) & (nmax < 35)
        nmax = nmax + 5;
        fmax = optimfun(nmax,a,nu,m,tol,errflag,c,cl,ce,k);
    end
    if ~(fmax < 0), % required precision unreasonably small
        n = nmax;
    else % Dekker-Brent to find zero in interval [nmin,nmax]
        n1 = nmin; nm1 = n1 + 1; n0 = n1 + 0.5;
        while abs(nmax-n1) > 0.5,
            if n1 ~= nm1,
                d = fmax*(nmax-nmin)/(fmax-fmin);
                if (sign(nmax-n1) ~= sign(d)) | (abs(d)>abs(n1-nmax)),
                    d = (nmax-n1)/2;
                end
            else
                d = (nmax-n1)/2;
            end
            nmin = nmax; fmin = fmax; nm1 = n0; n0 = n1;
            nmax = nmax - d;
            fmax = optimfun(nmax,a,nu,m,tol,errflag,c,cl,ce,k);
            if sign(fmax) ~= sign(fmin), n1 = nmin; end
        end
        n = (nmax+n1)/2;
    end
end
x = [c/cl*(n+(3+ce-(1+ce)/(1+cl))/4),n];
x = [x(1),ceil(x(2))]; % n should be integer
chi = 2^k*x(2)*tGJ + x(1)/pi*N*sum(a); % value of cost function

% ----------------------------------------------------------------------

function f = optimfun(n,a,nu,m,tol,errflag,c,cl,ce,k)
%OPTIMFUN Implements the precision constraint
x = c/cl*(n+(3+ce-(1+ce)/(1+cl))/4);

t = auxfun([x,n],a,nu,m,k); % auxiliary function
if t ~= 0,
    t = log(t);
    if errflag, % absolute error required
        f = (m-k/2-2*n-2)*log(x) + t + ...
           k/2*log(2/pi) - sum(log(a))/2 - log(tol);
    else
        f = t - (2*n+2)*log(x) + log(sum(a)) - log(tol);
    end
else
    f = -realmax;
end

% ----------------------------------------------------------------------

function f = auxfun(x,a,nu,m,k)
%AUXFUN Auxiliary function

x0 = x(1); n = x(2);
nu_shift = rem(nu,2); % take cosine mod 2*pi
ch = abs(cos(pi*nu_shift));
ch(find(abs(rem(nu_shift,1)-0.5) < eps)) = 0;

% coefficients in asymptotic expansion
cn = ch .* exp(gammaln(2*n+2.5-nu) - gammaln(2*n+3) +...
    gammaln(2*n+2.5+nu) - log(pi) - (n+1)*log(4));
dn = cn .* (2*n+2.5-nu) .* (2*n+2.5+nu) / (4*n+6);

A = sum(cn./a.^(2*n+2));
B = sum(dn./a.^(2*n+3));

f = A*x0/(2*n+1+k/2-m) + B/(2*n+2+k/2-m);

% ----------------------------------------------------------------------

function [f,x,param,s,sflag] = ira(a,nu,m,x,k,tol,errflag)
%IRA Infinite range approximation

iu = sqrt(-1);    % imaginary unit
x0 = x(1); n = x(2);

% coefficients of asymptotic expansion for J_nu
mu = 4*nu.*nu;
c(1,:) = ones(1,k);
d(1,:) = (mu-1)/8;
i = 1;
while ((max(abs(c(i,:))) < 1e3/eps) & (max(abs(d(i,:))) <...
        1e3/eps) & (i <= n)) | (i <= m/2 - k/4 + 1), % test for very large coefficients
    c(i+1,:) = -(mu - (4*i-3)^2) .* (mu - (4*i-1)^2) ./...
        ((2*i-1)*i*128).*c(i,:);
    d(i+1,:) = -(mu - (4*i-1)^2) .* (mu - (4*i+1)^2) ./...
        (i*(2*i+1)*128).*d(i,:);
    i = i + 1;
end

% parameters used in error estimate
A = sum(abs(c(end,:))./a.^(2*i)) / (2*i-1+k/2-m);
B = sum(abs(d(end,:))./a.^(2*i+1)) / (2*i+k/2-m);
param = [A,B];

% if value of n obtained from gop was too large, determine new
%  corresponding x0 such that precision constraint is satisfied
if i <= n,
    n = i - 1; % new value of n
    if errflag, % absolute error required
        ch = (2/pi)^(k/2)/sqrt(prod(a));
        A = ch*A; B = ch*B; p = m - k/2 - 2*n - 2;
    else
        ch = sum(a);
        A = A/ch; B = B/ch; p = -2*n - 2;
    end
    x1 = -1; x2 = x0; % initial value is original x0
    while abs(x1-x2) > abs(x2)*1e-2, % Newton-iteration
        x1 = x2;
        x2 = x1 - (log(A*x1+B) + p*log(x1) - log(tol))/...
            (A/(A*x1+B) + p/x1);
    end
    x0 = x2; x = [x0,n];
end

P = zeros(2*n+2,k);
P(1:2:2*n+1,:) = c;
P(2:2:2*n+2,:) = iu*d;
P = P(end:-1:1,:);
P = (1-iu)/(2*sqrt(pi))*repmat(exp(-pi/2*nu*iu),2*n+2,1).*P;

% sum over 2^(k-1) terms
If = 0; % integral of infinite part
s = ''; sflag = 0; % used to store warning in discontinuous case
Fe = 0; Fea = 0; % keep track of 'infinity' in discontinuous case
for j = 0:2^(k-1)-1,
    if k==1,
        b = 0;
    else % convert j to binary representation
        [f,e] = log2(j);
        b = rem(floor(j*pow2(2 - k:0)),2);
        b = b(end:-1:1);
    end
    % product of k laurent-"polynomials"
    F = a(1).^(-[2*n + 1:-1:0]).*P(:,1).'/sqrt(a(1));
    for i = 2:k,
        if b(i-1), Q = conj(P(:,i).'); else Q = P(:,i).'; end
        Ph = a(i).^(-[2*n + 1:-1:0]).*Q/sqrt(a(i));
        F = conv(F,Ph); % convolution gives coeffs of product
    end
    F = F(end:-1:1);
    % exponent
    e = a(1) + sum((-1).^b.*a(2:end));
    eabs = sum(a);
    I1 = 0; I2 = 1; n0 = k/2 - m - 1; i = 1;
    if (abs(e) <= eabs*eps) % special case
        nu_rem = mod(nu(1) + 0.5 + sum((-1).^b.*(nu(2:end) +...
            0.5)) , 2);
        % avoid cancellation by setting coefficients to zero
        nua = sum(abs(nu) + k/2);
        if abs(nu_rem - 1) <= max(5,nua)*eps,
            F(1:2:end) = 0;
        elseif abs(nu_rem) <= max(5,nua)*eps,
            F(2:2:end) = 0;
        end
        if n0 <= 0, % possible discontinuous problem, return warning
            s1 = ' a(1)'; sflag = 1;
            for j = 2:k, % make string stating which combination is zero
                if b(j-1), op = '-'; else op = '+'; end
                s1 = strcat(s1,op,'a(',num2str(j),')');
            end
            s = strcat(s,'Assuming ',s1,'=0;');
            Fe = Fe + F(1:fix(-n0)+1); % gather coeffs of infinite integrals...
            Fea = Fea + abs(F(1:fix(-n0)+1));
            i = i + 1 + fix(-n0); % ...then skip infinities
            n_temp = mod(n0,1);
            if n_temp <= max(5,(-n0))*eps, % n_temp probably zero
                n0 = 1;
            else
                n0 = n_temp;
            end
        end
    end
    % compute sum of integrals, until constant value or out of coefficients
    while (abs(I2-I1) >= abs(I1)*eps) & (i <= length(F)),
        if abs(F(i)) > eps,
            if abs(e) <= eabs*eps,
                Ih = -x0^(-n0)/(-n0);
            else % compute integrals using incomplete Gamma function
                Ih = (-iu*e)^n0*igamma(-n0,-iu*e*x0);
            end
            I2 = I1;
            I1 = I2 + F(i)*Ih;
        end
        n0 = n0 + 1; i = i + 1;
    end
    If = If + I1;
end
% check for non-zero coeffs in front of infinite integrals
tmp = find(abs(Fe) > max(5,Fea)*eps );
if ~isempty(tmp),
    f = real(Fe(tmp(1)))*inf; % first non-zero coeff is dominant for inf
else
    f = 2*real(If);
end

% ----------------------------------------------------------------------

function [f,err,eflag] = fri(a,nu,m,x0,tol)
%FRI Finite range integration

% get rough approximation to number of zeros in integrand
n = max(ceil(sum(x0*a/pi+1/4-nu/2)),0)+1;
h = x0/n; % length of subintervals
z = [0:h:x0];

f = 0; err1 = 0;
% treat possible algebraic singularity in 0 (first interval)
p = sum(nu) + m; % degree of singularity
if abs(fix(p)-p) > abs(p)*50*eps, % fractional degree -> extrapolate
    xs = z(2); z = z(2:end);
    rho = 1/4; f1 = 0; i = 0;
    [f2,err,eflag] = nqf(@fun,[xs/2,xs],tol,a,nu,m);
    I1(1) = f2; xs = xs/2;
    while ~eflag & abs(f2-f1) > tol(1)*abs(f2) & abs(f2-f1) > tol(2),
        i = i + 1;
        [I2(1),err,eflag] = nqf(@fun,[xs*rho,xs],tol,a,nu,m);
        I2(1) = I1(1) + I2(1);
        for j = 1:i, % extrapolate in sense of Richardson
            I2(j+1) = (I2(j) - rho^(p+j)*I1(j)) / (1 - rho^(p+j));
        end;
        f1 = I1(i); f2 = I2(i+1); I1 = I2; xs = xs*rho;
    end
    f = f2;
    err1 = abs(f-f1);
    if length(z) == 1,
        err = err1;
        return
    end
end

z = [z(1:2:end-1) x0]; % take 2 intervals together to save work

% split integrand at equidistant points to integrate
V = [z(1:end-1)' z(2:end)'];
[f2,err,eflag] = nqf(@fun,V,tol,a,nu,m);
f = f + f2;
err = err + err1;

% ----------------------------------------------------------------------

function f = fun(x,a,nu,m)
%FUN Integrand

f = x.^m;
for i = 1:length(a),
    f = f.*besselj(nu(i),a(i)*x);
end

% ----------------------------------------------------------------------

function [I,err,eflag] = nqf(fun,X,tol,varargin)
%NQF Numerical quadrature formula

% some hard coded Gauss-Legendre rules
xw5=[0,                         0.28444444444444444444;
     0.53846931010568309104,    0.23931433524968323402;
     0.90617984593866399280,    0.11846344252809454375;];

xw7=[0,                         0.20897959183673469388;
     0.40584515137739716691,    0.19091502525255947247;
     0.74153118559939443986,    0.13985269574463833395;
     0.94910791234275852453,    0.064742483084434846630;];

xw11=[0,                        0.13646254338895031536;
      0.26954315595234497233,   0.13140227225512333109;
      0.51909612920681181593,   0.11659688229599523996;
      0.73015200557404932409,   0.093145105463867125710;
      0.88706259976809529908,   0.062790184732452312314;
      0.97822865814605699280,   0.027834283558086833246;];

xw15=[0,                        0.10128912096278063644;
      0.20119409399743452230,   0.099215742663555788230;
      0.39415134707756336990,   0.093080500007781105516;
      0.57097217260853884754,   0.083134602908496966776;
      0.72441773136017004742,   0.069785338963077157225;
      0.84820658341042721620,   0.053579610233585967505;
      0.93727339240070590431,   0.035183023744054062354;
      0.98799251802048542849,   0.015376620998058634177;];

xw19=[0,                        0.080527224924391847988;
      0.16035864564022537587,   0.079484421696977173823;
      0.31656409996362983199,   0.076383021032929833389;
      0.46457074137596094572,   0.071303351086803305889;
      0.60054530466168102347,   0.064376981269668113838;
      0.72096617733522937862,   0.055783322773666997358;
      0.82271465653714282498,   0.045745010811224999731;
      0.90315590361481790164,   0.034522271368820613291;
      0.96020815213483003085,   0.022407113382849800168;
      0.99240684384358440319,   0.0097308941148632385170;];

xw23=[0,                        0.066827286093053087676;
      0.13325682429846611093,   0.066231019702348308687;
      0.26413568097034493053,   0.064452861094041074987;
      0.39030103803029083142,   0.061524542153364765234;
      0.50950147784600754969,   0.057498320111205682472;
      0.61960987576364615639,   0.052446045732270705037;
      0.71866136313195019446,   0.046457883030017573738;
      0.80488840161883989215,   0.039640705888359477462;
      0.87675235827044166738,   0.032116210704262926063;
      0.93297108682601610235,   0.024018835865542334286;
      0.97254247121811523196,   0.015494002928489722153;
      0.99476933499755212352,   0.0067059297435708860456;];

xw={xw5,xw7,xw11,xw15,xw19,xw23};

if (tol(1) <= 50*eps) & (tol(2) <= 50*eps),
    k = 4; % avoid taking not enough points
else
    k = 1;
end
N = size(X,1); eflag = 0;
I1 = 0; % dummy value
for j = k:length(xw),
    If = 0;
    w = [xw{j}(end:-1:2,2); xw{j}(:,2)];
    for i = 1:N, % integrate over different subintervals
        L = (X(i,2)-X(i,1))/2; M = (X(i,2)+X(i,1))/2;
        x = L*[-xw{j}(end:-1:2,1);xw{j}(:,1)]' + M;
        If = If + 2*L*feval(fun,x,varargin{:})*w;
    end
    I2 = I1;
    I1 = If;
    err = abs(I2-I1);
    if (j > k) & ((err <= tol(1)*abs(I1)) | (err <= tol(2))),
        I = I1; % precision reached
        return
    end
end
% requested precision not reached
I = I1; eflag = 1;

% --- END BESSELINT ----------------------------------------------------
