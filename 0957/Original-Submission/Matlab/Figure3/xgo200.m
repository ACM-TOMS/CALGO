function x=xgo200(n)
%XGO200 Lower boundary for a natural domain.
%XGO200(N) determines, for given N, the (negative) value 
%   X=X^{**}(N), and an approximation for it, below which
%   Gaussian quadrature requires more than 200 points to
%   evaluate f_N(X) to an accuracy of 12 decimal digits.
 
global ab ab0
N0=195; dN=1;
eps0=.5e-12; tol=.5e-5;
for x=-1:-1:-30
  [N,y,go200]=N_inerfc(n,x,eps0,N0,dN);
  if go200==1
    a=x; b=x+1;
    break
  end
end
%[a b]
%[N,y,go200]=N_inerfc(n,a,eps0,N0,dN); go200
%[N,y,go200]=N_inerfc(n,b,eps0,N0,dN); go200
ktol=ceil(log((b-a)/tol)/log(2));
for k=1:ktol
%  [a b]
%  [N,y,go200]=N_inerfc(n,a,eps0,N0,dN); go200
%  [N,y,go200]=N_inerfc(n,b,eps0,N0,dN); go200
  x=(a+b)/2;
  [N,y,go200]=N_inerfc(n,x,eps0,N0,dN);
  if go200==0, b=x; else a=x; end
end
