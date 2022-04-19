clear all

syms x t mu eps positive
syms w y z

% gam_l == log( sqrt(2*pi) * Gamma(x) )
% arg_l == log( integrand )

gam_l = (x-1/2)*log(x) - x + 1/(12*x);   %  + O(1/x^2) error
arg_l = - gam_l - t + (x-1)*log(t) + log(mu)/2;

% switch from t,x,mu to z,y,eps, and expand in eps

arg_l = subs(arg_l,t,x -sqrt(mu)*z);
arg_l = subs(arg_l,x,mu+sqrt(mu)*y);
arg_l = subs(arg_l,mu,eps^(-2));
arg_l = simple(taylor(arg_l,eps));
term0 = simple(subs(arg_l,eps,0))     % leading order term
arg = taylor(exp(arg_l-term0)-1,eps); % error

% for each power of epsilon, decompose integrand error
% into powers of z, and then integrate by parts symbolically

for m=1:3
  arg  = diff(arg,'eps')/m;
  c    = simple(subs(arg,eps,0));

  nmax = 0;
  while (diff(c,'z',nmax+1)~=0)
    nmax = nmax+1;
  end

  C{1} = subs(c,z,0);
  for n = 1:nmax-1
    C{n+1} = subs(diff(c,'z',n),z,0)/factorial(n);
  end
  C{nmax+1} = diff(c,'z',nmax)/factorial(nmax);

  f{m} = 0;
  for n = nmax:-1:1
    f{m} = f{m} - C{n+1}*y^(n-1);
    if (n>1)
      C{n-1} = C{n-1} + (n-1)*C{n+1};
    end
  end
end

f1 = simple(f{1})
f2 = simple(f{2})
f3 = simple(f{3})

% convert into asymptotic expansion for inverse

f1 = subs(f1,y,w);
f2 = subs(f2,y,w);
f3 = subs(f3,y,w);

g0_1 = 1;
g0_2 = w;
g0_3 = 2*w*g0_2 + diff(g0_2,'w');

g1   = - f1
g1_1 = diff(g1,'w');
g1_2 = w*g1_1 + diff(g1_1,'w');

g2   = - g1_1*f1 - g0_1*f2 - g0_2*(f1^2)/2;
g2   = simple(g2)
g2_1 = diff(g2,'w');

g3   = - g0_1*f3 - g0_2*f1*f2 - g0_3*f1^3/6 ...
       - g1_1*f2 - g1_2*f1^2/2 - g2_1*f1;
g3   = simple(g3)

%
% verification that Temme expansion leads to same result
%

syms r rm positive
syms s w

f  = sqrt(2*(1-r + r.*log(r)));
c0 = log(f*sqrt(r)/(r-1)) / log(r);

f  = subs(f,r,1+rm);
ft = taylor(f,rm)

fi = s;
for n = 1:6
  fi = taylor(s-subs(ft-rm,rm,fi),s);
end
fi = taylor(fi,s)

c0 = subs(c0,r,1+rm);
c0 = taylor(c0,rm)

c1 = sym('-8/405');

temme = 1+rm + c0*eps^2 + c1*eps^4;
temme = subs(temme,rm,fi);
temme = subs(temme,s,eps*w);
g1 = subs(diff(temme,eps,2),eps,0)/2
g2 = subs(diff(temme,eps,3),eps,0)/6
g3 = subs(diff(temme,eps,4),eps,0)/24
