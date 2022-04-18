% Example script to test functionality of AMGKQ.
% Copyright (C) 2014 Robert W. Johnson
%
% This file is part of AMGKQ.  See AMGKQ.M for details.
% Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
% This is free software; see GPLv3 or greater for copying conditions.
% There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.

% AMGKQ: [RES, ERR, NSUB, FL] = amgkq(F, A, B, C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, CFLAG, VERB, P1, P2, ...)

disp(' ')
[RES, ERR, NSUB, FL] = amgkq('sin',-pi,pi,[],[],[],7)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(inline('sin'),-pi,pi,[],[],[],8)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@sin,-pi,pi,[],[],[],9)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) sin(x),-pi,pi,[],[],[],10)

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x.^2),-inf,inf), ACC = RES - sqrt(pi)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x.^2).*cos(x),0,inf), ACC = RES - sqrt(pi)/2/exp(1/4)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x.^2)./(1+x.^2),0,1), ACC = RES - pi/4*exp(1)*(1-erf(1)^2)

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x.^2),-inf,inf,[],[],[],[],2), ACC = RES - sqrt(pi)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./(sqrt(x.*(1-x))),0,1,[],[],[],[],2), ACC = RES - pi

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x).*x./(1-exp(-2*x)),0,inf), ACC = RES - pi^2/8
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x).*x./(1-exp(-2*x)),0,inf,[],0,1e3), ACC = RES - pi^2/8
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) log(x)./(1-x.^2),0,1), RES + pi^2/8
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) log(x)./(1-x.^2),0,1,[],0,1e3), RES + pi^2/8

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./sqrt(abs(x)),0,10), ACC = RES - 2*sqrt(10)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./sqrt(abs(x)),-10,10,[],0,1e3), ACC = RES - 4*sqrt(10)

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./(sqrt (x).*(1+x)),0,inf), ACC = RES - pi
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./(sqrt (-x).*(1-x)),-inf,0), ACC = RES - pi

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) bsxfun(@times,exp(-x),[x; x.^2; x.^3; x.^4; x.^5]),0,inf), ACC = RES - gamma(2:6)'
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) sin(3*x).*cosh(x).*sinh(x),10,15,[],[0 1e-14]), RELACC = RES / 2.588424538641647e+10 - 1

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./sqrt(abs(x(1,:))).*sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:)),[0;10],[10;15],[],[0 1e-14]), RELACC = RES ./ (2*sqrt(10)*2.588424538641647e+10) - 1
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) 1./sqrt(abs(x(1,:))).*sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:)),[-10;10],[10;15],[],[0 1e-14],2e3), RELACC = RES ./ (4*sqrt(10)*2.588424538641647e+10) - 1
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) [1./sqrt(abs(x(1,:))); sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:))],[0;10],[10;15],[],[0 1e-14]), RELACC = RES ./ [2*sqrt(10)*5;10*2.588424538641647e+10] - 1
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) [1./sqrt(abs(x(1,:))); sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:))],[-10;10],[10;15],[],[],2e3), RELACC = RES ./ [4*sqrt(10)*5;20*2.588424538641647e+10] - 1

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(z) [1./(1+z.^2).^2; exp(i*z)./(1+z.^2)],-1,-1,[1,2*i]), ACC = RES - pi*[1/2; exp(-1)]

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) prod(exp(-x.^2),1),[-inf;inf],[inf;-inf],[],0), RES + pi
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2),[-inf;-inf],[inf;inf]), ACC = RES - sqrt(2*pi)*pi

p = 1/2; q = 1/2;
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) x.^(p-1).*(1-x).^(q-1),0,1), ACC = RES - beta(p,q)
p = 1/3; q = 1/3;
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) x.^(p-1).*(1-x).^(q-1),0,1), ACC = RES - beta(p,q)
p = 1/4; q = 1/4;
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) x.^(p-1).*(1-x).^(q-1),0,1), ACC = RES - beta(p,q)

disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) sin(x)./x,0,inf,[],[],1e3,[],[],0), ACC = RES - pi/2
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) (sin(x)./x).^2,0,inf,[],[],1e3), ACC = RES - pi/2
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) (sin(x)./x).^3,0,inf,[],[],1e3), ACC = RES - 3*pi/8
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) (sin(x)./x).^4,0,inf,[],[],1e3), ACC = RES - pi/3

AB = [-1 1; -1 1];
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) [(sum(x.^2,1) < 1); (sum(x.^2,1) > 1)],AB(:,1),AB(:,2),[],[],2e3), ACC = RES - [pi; 4-pi]
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) [(sum(x.^2,1) < 1); (sum(x.^2,1) > 1)],AB(:,1),AB(:,2),[],[],1e3,40), ACC = RES - [pi; 4-pi]

AB = [-1 1; -1 1]*10;
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2),AB(:,1),AB(:,2)), ACC = RES - sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*2*atan(AB(2,2))
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) [exp(-x(1,:).^2/2); 1./(1+x(2,:).^2)],AB(:,1),AB(:,2)), ACC = RES - [sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*diff(AB(2,:));2*atan(AB(2,2))*diff(AB(1,:))]

AB = [[-1 1; -1 1]*10; 0 1];
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2).*x(3,:).^10.*(1-x(3,:)).^10,AB(:,1),AB(:,2),[],0), ACC = RES - sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*2*atan(AB(2,2))*beta(11,11)
disp(' ')
[RES, ERR, NSUB, FL] = amgkq(@(x) [exp(-x(1,:).^2/2); 1./(1+x(2,:).^2); x(3,:).^10.*(1-x(3,:)).^10],AB(:,1),AB(:,2)), ACC = RES - [sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*diff(AB(2,:));2*atan(AB(2,2))*diff(AB(1,:));beta(11,11)*prod(diff(AB(1:2,:),[],2))]

