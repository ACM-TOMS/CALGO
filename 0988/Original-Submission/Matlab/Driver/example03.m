% Example script to test improper integrals and other features of AMGKQ.
% Copyright (C) 2014 Robert W. Johnson
%
% This file is part of AMGKQ.  See AMGKQ.M for details.
% Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
% This is free software; see GPLv3 or greater for copying conditions.
% There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.

% preparations
close all
clear all
tabcell = cell(1,9);
ntab = 0;
disp(' ')

fx = @(z) [1./(1+z.^2).^2; exp(i*z)./(1+z.^2)]
a = -1
b = -1
c = [1, 2*i]
exact = pi*[1/2; exp(-1)]
[RES, ERR, NSUB, FL] = amgkq(fx,a,b,c), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) 1./sqrt(abs(x))
a = 0
b = 10
exact = 2*sqrt(10)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) 1./sqrt(abs(x))
a = -10
b = 10
exact = 4*sqrt(10)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) 1./(sqrt (x).*(1+x))
a = 0
b = inf
exact = pi
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) log(x)./(1-x.^2)
a = 0
b = 1
exact = -pi^2/8
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x).*x./(1-exp(-2*x))
a = 0
b = inf
exact = pi^2/8
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

%fx = @(x) bsxfun(@times,exp(-x),[x; x.^2; x.^3; x.^4; x.^5])
fx = @(x) [exp(-x).*x; exp(-x).*x.^2; exp(-x).*x.^3; exp(-x).*x.^4; exp(-x).*x.^5]
a = 0
b = inf
exact = gamma(2:6)'
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x.^2)
a = -inf
b = inf
exact = sqrt(pi)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x.^2).*cos(x)
a = 0
b = inf
exact = sqrt(pi)/2/exp(1/4)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x.^2)./(1+x.^2)
a = 0
b = 1
exact = pi/4*exp(1)*(1-erf(1)^2)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2)
a = [-inf, -inf]
b = [inf, inf]
exact = sqrt(2*pi)*pi
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2)
a = [-10, -10]
b = [10, 10]
exact = sqrt(2*pi)*erf(b(1)/sqrt(2))*2*atan(b(2))
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) [exp(-x(1,:).^2/2); 1./(1+x(2,:).^2)]
a = [-10, -10]
b = [10, 10]
exact = [sqrt(2*pi)*erf(b(1)/sqrt(2))*(b(2)-a(2));2*atan(b(2))*(b(1)-a(1))]
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2).*x(3,:).^10.*(1-x(3,:)).^10
a = [-10, -10, 0]
b = [10, 10, 1]
exact = sqrt(2*pi)*erf(b(1)/sqrt(2))*2*atan(b(2))*beta(11,11)
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) [exp(-x(1,:).^2/2); 1./(1+x(2,:).^2); x(3,:).^10.*(1-x(3,:)).^10]
a = [-10, -10, 0]
b = [10, 10, 1]
exact = [sqrt(2*pi)*erf(b(1)/sqrt(2))*(b(2)-a(2));2*atan(b(2))*(b(1)-a(1));beta(11,11)*(b(1)-a(1))*(b(2)-a(2))]
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

p = 1/2; q = 1/2;
fx = @(x) x.^(-1/2).*(1-x).^(-1/2)
a = 0
b = 1
exact = beta(p,q)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

p = 1/3; q = 1/3;
fx = @(x) x.^(-2/3).*(1-x).^(-2/3)
a = 0
b = 1
exact = beta(p,q)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

p = 1/4; q = 1/4;
fx = @(x) x.^(-3/4).*(1-x).^(-3/4)
a = 0
b = 1
exact = beta(p,q)
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) (sin(x)./x).^2
a = 0
b = inf
exact = pi/2
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) (sin(x)./x).^3
a = 0
b = inf
exact = 3*pi/8
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) (sin(x)./x).^4
a = 0
b = inf
exact = pi/3
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) [(sum(x.^2,1) < 1); (sum(x.^2,1) > 1)]
a = [-1, -1]
b = [1, 1]
exact = [pi; 4-pi]
[RES, ERR, NSUB, FL] = amgkq(fx,a',b'), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) sin(x)./x
a = 0
b = inf
exact = pi/2
longarg{6} = 0; longarg{3} = 1e3
[RES, ERR, NSUB, FL] = amgkq(fx,a,b,longarg{:}), ACC = RES - exact, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,ACC};

fx = @(x) sin(3*x).*cosh(x).*sinh(x)
a = 10
b = 15
exact = 2.588424538641647e+10
[RES, ERR, NSUB, FL] = amgkq(fx,a,b), RELACC = RES/exact - 1, disp(' ')
ntab = ntab + 1;
tabcell(ntab,:) = {fx,a,b,exact,RES,ERR,NSUB,FL,RELACC};

