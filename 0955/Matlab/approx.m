%
% this routine produces various polynomial approximations
%

function approx()

close all;
clear all;

%
% set extent of central region in which approximations are defined
%

rhi = 3.25;
rlo = 0.4;

slo = sqrt(2*(1-rlo + rlo.*log(rlo))) .* sign(rlo-1);
shi = sqrt(2*(1-rhi + rhi.*log(rhi))) .* sign(rhi-1);

file = ' ';
file = char(file,'//  use polynomial approximations in central region',' ');
file = char(file,sprintf('    if ( (s>%10.7gf) && (s<%8.7gf) ) {',slo,shi));
disp(file)

%
% plot f(r)
%

r = linspace(0.0001, rhi);
y = sqrt(2*(1-r + r.*log(r))) .* sign(r-1);
figure(1)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
plot(r,y)
xlabel('r'); ylabel('f(r)')

print('-deps2','approx1.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% first, approximate f^{-1}
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th  = pi*(0.5 + 0:199)'/200;
r   = 0.5*(rhi+rlo) - 0.5*(rhi-rlo)*cos(th);
y   = sqrt(2*(1-r + r.*log(r))) .* sign(r-1);
s   = (r-1)./y;

%
% least-squares polynomial fit with diagonal weighting
%

degree = 14;

D = diag(1./s);
A = zeros(length(y),degree+1);
for deg = 0:degree
  A(:,degree+1-deg) = y.^deg;
end

A = A(:,2:end-2);
p1 = ((D*A)\(D*(s-1-y/6)))'; % subtract out leading order behaviour
p1 = [p1 1/6 1];             % add back in leading order behaviour
s2 = polyval(p1,y);
p1 = [p1 0];

%
% plot error
%

figure(2)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);

subplot(1,2,1)
plot(y,y.*(s-s2))
xlabel('s'); ylabel('error'); title('    f^{-1} approximation')

%
% output polynomial approximation
%

file = ' ';
file = char(file,'//  polynomial approximation to f^{-1}(s) - 1',' ');
file = char(file,sprintf('      rm =%16.9gf;',p1(1)));

%for deg = 1:degree
%  file = char(file,sprintf('      rm =%16.9gf + rm*s;',p1(deg+1)));
%end

for deg = 1:degree-2
  file = char(file,sprintf('      rm =%16.9gf + rm*s;',p1(deg+1)));
end
file = char(file,'      rm =             s + s*(rm*s);');

disp(file)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% second, approximate c_0(r)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

th  = pi*(0.5 + 0:199)'/200;
r   = 0.5*(rhi+rlo) - 0.5*(rhi-rlo)*cos(th);
y   = sqrt(2*(1-r + r.*log(r))) .* sign(r-1);
c   = log(y.*sqrt(r)./(r-1)) ./ log(r);

%
% least-squares polynomial fit with diagonal weighting
%

degree = 12;

D = diag(1./r.^0.5);
A = zeros(length(r),degree);
for deg = 1:degree
  A(:,degree+1-deg) = (r-1).^deg;
end

p2 = ((D*A)\(D*(c-1/3)))'; % subtract out leading order behaviour
p2 = [p2 1/3];             % add back in leading order behaviour

c2 = polyval(p2,r-1);

%
% plot error
%

subplot(1,2,2)
plot(r,c-c2); 
xlabel('r'); ylabel('error'); title('    c_0 approximation')

print('-deps2','approx2.eps')

%
% output polynomial approximation
%

file = ' ';
file = char(file,'//  polynomial approximation to correction c0(r)',' ');
% file = char(file,sprintf('      rm = r - 1f;'));
file = char(file,sprintf('      t  = %16.9gf;',p2(1)));

for deg = 1:degree
  file = char(file,sprintf('      t  = %16.9gf + t*rm;',p2(deg+1)));
end

disp(file)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% third, additional O(1/lam) correction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xs = [10 10.5 11 12 13 16 20 50 80 100];

degree = 10;

D = eye(1000*length(xs));
A = zeros(1000*length(xs),degree+1);
c = zeros(1000*length(xs),1);

for k = 1:length(xs)
  x = xs(k);

  lam_min = x/rhi;
  lam_max = x/rlo;

  lam = linspace(lam_min,lam_max,1000);
  w   = zeros(size(lam));
  x3  = zeros(size(lam));
  
%
% for accuracy in tails, need to select correct variant
%

  upper = find(x<lam);
  w(upper) = norminv(gammainc(lam(upper),x,'upper'));
  lower = find(x>=lam);
  w(lower) = -norminv(gammainc(lam(lower),x,'lower'));

%
% use Newton iteration to find rm
%

  s  = w./sqrt(lam);
  rm = finv(s);
  c0 = log(s.*sqrt(1+rm)./rm) ./ log(1+rm);
  x2 = lam + lam.*rm;

%
% full correction of leading order error
%
  t = abs(rm)<0.0001;
  x3(t)  = x2(t) + 1/3 - rm(t)/36;
  x3(~t) = x2(~t) + c0(~t);

  c(1000*(k-1)+1:1000*k) = (x - x3)';
  D(1000*(k-1)+1:1000*k,1000*(k-1)+1:1000*k) = ...
       diag((1:1000).^(-0.75).*(1000:-1:1).^(-0.25));

  for deg = 0:degree
    A(1000*(k-1)+1:1000*k,degree+1-deg) = rm.^deg ./ lam;
  end
end

%
% compute polynomial approximation
%

p3 = ((D*A)\(D*c))';

file = ' ';
file = char(file,'//  O(1/lam) correction',' ');
% file = char(file,sprintf('      rm = r - 1f;'));
file = char(file,sprintf('      x  = %16.9gf;',p3(1)));

for deg = 1:degree
  file = char(file,sprintf('      x  = %16.9gf + x*rm;',p3(deg+1)));
end

disp(file)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% final validation -- part 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);

xs  = [10 100 exp(linspace(log(10),log(1000),2000))];
err = zeros(size(xs));

for k = 1:length(xs)
  x = xs(k);
  
% double precision range
  lam_min = max([gammaincinv(1e-16, x,'lower') x/rhi 4]);
  lam_max = min([gammaincinv(1e-300,x,'upper') x/rlo]);

% single precision range
  lam_min = max([gammaincinv(1e-7, x,'lower') x/rhi 4]);
  lam_max = min([gammaincinv(1e-38,x,'upper') x/rlo]);

  lam = linspace(lam_min,lam_max,1000);
  w   = zeros(size(lam));

  upper = find(x<lam);
  w(upper) = norminv(gammainc(lam(upper),x,'upper'));
  lower = find(x>=lam);
  w(lower) = -norminv(gammainc(lam(lower),x,'lower'));

  rm = polyval(p1,w./sqrt(lam));

  x2 = lam.*(1+rm) + polyval(p2,rm) + polyval(p3,rm)./lam;

  if k<=2
    subplot(1,2,k)
    plot(lam,x2-x)
    title(['$x = ' num2str(x) '$'],'interpreter','latex')
    xlabel('$\lambda$','interpreter','latex'); ylabel('error')
  else
    err(k) = max(abs(x2-x));
  end
end

print('-deps2','approx3.eps')

figure(4)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
semilogx(xs(3:end),err(3:end));
xlabel('$x$','interpreter','latex')
ylabel('maximum error')
%axis([1e1 4e3 0 1e-4])

print('-deps2','approx4.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% final validation -- part 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(5)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.0]; set(gcf,'pos',pos);

xs  = [10 100 exp(linspace(log(10),log(1000),1000))];
err = zeros(size(xs));

for k = 1:length(xs)
  x  = xs(k);

% double precision range

% brute force search for lam_min because gammaincinv fails

  lam_min = x;
  fac     = 2;
  for kount = 1:20
     while(gammainc(lam_min/fac,x,'lower')>1e-300)
       lam_min = lam_min/fac;
     end
     fac = sqrt(fac);
  end
  lam_min = max(4,lam_min);

%  lam_min = max(4, gammaincinv(eps(0.5),x,'lower'));
  lam_max =        gammaincinv(1e-300,x,'upper');

% single precision range
%  lam_min = max(4, gammaincinv(6.e-8,x,'lower'));
%  lam_max =        gammaincinv(1e-38,x,'upper');

  lam = linspace(lam_min,lam_max,1000);
  w   = zeros(size(lam));

  upper = find(x<lam);
  w(upper) = norminv(gammainc(lam(upper),x,'upper'));
  lower = find(x>=lam);
  w(lower) = -norminv(gammainc(lam(lower),x,'lower'));

  s  = w./sqrt(lam);
  rm = finv(s);
  c0 = log(s.*sqrt(1+rm)./rm) ./ log(1+rm);
  x1 = lam+lam.*rm;

  t = abs(rm)<0.0001;
  x1(t)  = x1(t) + 1/3 - rm(t)/36;
  x1(~t) = x1(~t) + c0(~t);

%  x2 = x1 - (8/405)*(1+0.025)./(x1+0.025*lam);

  x2 = x1 - 0.0218./(x1+0.065*lam);

  if k<=2
    subplot(2,2,k)
    plot(lam,x1-x)
    title(['x = ' num2str(x)])
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 1')
    subplot(2,2,k+2)
    plot(lam,x2-x)
    xlabel('$\lambda$','interpreter','latex'); ylabel('error 2')
  else
    err1(k) = max(abs(x1-x));
    err2(k) = max(abs(x2-x));
  end
end

print('-deps2','approx5.eps')

figure(6)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.6]; set(gcf,'pos',pos);
loglog(xs(3:end),err1(3:end),'-k',xs(3:end),err2(3:end),'-.k');
xlabel('$x$','interpreter','latex')
ylabel('maximum error')
legend('error 1','error 2')
%axis([100 1e4 0 0.002])
%axis([10 1e4 0 0.0025])

print('-deps2','approx6.eps')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Newton iteration to calculate f^{-1}(s) - 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rm = finv(s)

t     = (s > -sqrt(2));
r(~t) = 0;
r(t)  = max(0.1,1+s(t));
t     = t & (s~=0);

while (sum(t)>0)
  r_old = r;
  logr  = log(r(t));
  f = sqrt(2*(1 - r(t) + r(t).*logr)) .* sign(r(t)-1);
  r(t) = r_old(t) - (f-s(t)).*f./logr;
  r(t) = max(r(t), 0.1*r_old(t));
  t = abs(r-r_old)>1e-10;
end

rm = r-1;

end
