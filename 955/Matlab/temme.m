%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% implements and tests the Temme asymptotic 
% approximation to the incomplete Gamma function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function temme()

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute coefficients f_n for Temme approximation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 12;
b = zeros(1,2*N+2);
f = zeros(1,2*N+2);

f(1) = -1/3;
f(2) =  1/12;
f(3) = -2/135;

for k=4:2*N+2
  j = 3:k-1;
  f(k) = - (k+1)/(k+2) * ( (k-1)*f(k-1)/(3*k) + ...
                sum( f(j-1).*f(k+1-j)./(k+2-j) ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute and plot error in Temme approximation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 1.5]; set(gcf,'pos',pos);

xs = [10 100 10:1000];

for k = 1:length(xs)
  x = round(xs(k));
  lam_min = gammaincinv(eps(0.5),x,'lower');
  lam_max = gammaincinv(1e-300,x,'upper');
  lam = exp(linspace(log(lam_min),log(lam_max),1000));

  y   = zeros(size(lam));
  y2  = zeros(size(y));
  y2p = zeros(size(y));

  upper = find(x<lam);
  y(upper) =  gammainc(lam(upper),x,'upper');
  lower = find(x>=lam);
  y(lower) = -gammainc(lam(lower),x,'lower');

  for i = 1:length(y)
    [y2(i), y2p(i)] = temme_eval(lam(i),x,f);
  end

  if k<3
    subplot(2,2,k)
%    plot(lam,(y-y2),'-k')
    semilogx(lam,(y-y2),'-k')
    title(['x = ' num2str(x)])
    xlabel('$\lambda$','interpreter','latex'); ylabel('absolute error')
    if k==1
      axis([0.1 1000 -1e-15 1e-15]);
      set(gca,'xtick',[0.1,10,1000])
    else
      axis([10 1200 -1e-15 1e-15]);
      set(gca,'xtick',[10,100,1000])
    end
    subplot(2,2,k+2)
%    plot(lam,(y-y2)./y2p,'-k')
    semilogx(lam,(y-y2)./y2p,'-k')
    xlabel('$\lambda$','interpreter','latex'); ylabel('relative error')
    if k==1
      axis([0.1 1000 -1e-13 1e-13]);
      set(gca,'xtick',[0.1,10,1000])
    else
      axis([10 1200 -1e-13 1e-13]);
      set(gca,'xtick',[10,100,1000])
    end
  else
    re(k-2) = max(abs(y-y2)./y2p);
  end
end

print('-deps2','temme1.eps');

figure(2)
%pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.2 0.8]; set(gcf,'pos',pos);
pos=get(gcf,'pos'); pos(3:4)=pos(3:4).*[1.0 0.65]; set(gcf,'pos',pos);
loglog(xs(3:end),re,'-k')
xlabel('x'); ylabel('maximum relative error')
print('-deps2','temme2.eps');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% generate code to evaluate Temme approximation for S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file = ' ';
file = strvcat(file,'  Double precision Temme approximation to function S',' ');
file = strvcat(file,sprintf('        B1 =%24.17g;              S = B1;',f(2*N+2)));
file = strvcat(file,sprintf('        B0 =%24.17g;              S = B0 + S*Eta;',f(2*N+1)));

for k = 2*N-1:-2:1
  file = strvcat(file,sprintf('        B1 =%24.17g + %2d.0*B1*Xi; S = B1 + S*Eta;',...
                                                   f(k+1),k+2.0) );
  file = strvcat(file,sprintf('        B0 =%24.17g + %2d.0*B0*Xi; S = B0 + S*Eta;',...
                                                   f(k  ),k+1.0) );
end

file = strvcat(file,sprintf('      S  = S / (1.0 + B1*Xi);'));

disp(file)

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function to evaluate Temme approximation for S
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q, Qp] = temme_eval(lam,x,f)

N = length(f)/2 - 1;

b(2*N+1:2*N+2) = f(2*N+1:2*N+2);

for k=2*N:-1:1
  b(k) = f(k) + (k+1)*b(k+2)/x;
end

r   = x/lam;
eta = sqrt(2*(1-r+r*log(r))/r) * sign(1-r);

S  = sum(b.*eta.^(0:2*N+1)) / (1+b(2)/x);
if x<lam
  Q  = ncf(-eta*sqrt(x)) + exp(-0.5*x*eta^2)*S/sqrt(2*pi*x);
else
  Q  = -ncf(eta*sqrt(x)) + exp(-0.5*x*eta^2)*S/sqrt(2*pi*x);
end
Qp = exp(-0.5*x*eta^2)/sqrt(2*pi*x);

% sum downwards

if (r<0.5)
  xi = 1/x;
% coefficients from http://web.mit.edu/kenta/www/three/stirling.txt
  s  = xi/12 - xi^3/360 + xi^5/1260 - xi^7/1680 + xi^9/1188 - xi^11*691/360360;
  t  = exp ( (x-lam) - x*log(x/lam) - s ) * sqrt(x / (2*pi*lam*lam));

%  t = exp(-lam-gammaln(x)+(x-1)*log(lam));

  Q = t;
  for i=1:50
    t = t*(x-i)/lam;
    Q = Q + t;
  end

% sum upwards

elseif (r>2)
  xi = 1/x;
  s  = xi/12 - xi^3/360 + xi^5/1260 - xi^7/1680 + xi^9/1188 - xi^11*691/360360;
  t  = exp ( (x-lam) - x*log(x/lam) - s ) / sqrt(2*pi*x);

%  t = exp(-lam-gammaln(x)+(x-1)*log(lam)) * lam/x;

  Q = -t;
  for i=1:50
    t = t*lam/(x+i);
    Q = Q - t;
  end
%  end
end

end

%
% ncf function based on erfc
%

function P = ncf(z)
  P = 0.5 * erfc(-z/sqrt(2));
end
