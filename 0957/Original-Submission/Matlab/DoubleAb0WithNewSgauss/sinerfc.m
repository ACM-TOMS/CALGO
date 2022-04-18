function [y,N]=sinerfc(dig,n,x)
%SINERFC Symbolic dig-digit precision counterpart of inerfc.m.

global ab
digits(dig); x=vpa(x);
if n<0 | n>150, error('n out of range'); end
if subs(x)<-14.6 | subs(x)>27, error('x out of range'); end
if subs(x)>=0
  if n<=25, N=ceil(43+3.92*subs(x));
  elseif n<=50, N=ceil(43+4.042*subs(x));
  elseif n<=75, N=ceil(39+4.143*subs(x));
  elseif n<=100, N=ceil(49+3.333*subs(x));
  elseif n<=125, N=ceil(58+3.333*subs(x));
  else N=ceil(67+2.958*subs(x));
  end
else
  if n>50, error('n out of range'); end
  if subs(x)<-5, error('x out of range'); end
  if abs(subs(x))<=4, N=floor((47+37*abs(subs(x)))/4+.46*n);
    else N=floor(57+2.36*n);
  end
end
disp(' Printing whos information from sinerfc just before call to sgauss ')
disp(' Have changed N is call to sgauss to int32(N) ')
disp(' to prevent problem with zeros call in sgauss ')
whos
uv=sgauss(dig,int32(N),ab); u=uv(:,1); v=uv(:,2);
if subs(x)>=0
  y0=sum(v.*u.^n.*exp(-2*x*u));
  y=2*exp(-x^2)*y0/(sqrt(vpa(pi))*gamma(vpa(n+1)));
else
  y0=sum(v.*u.^n.*exp(-(x^2+2*x*u)));
  y=2*y0/(sqrt(vpa(pi))*gamma(vpa(n+1)));
end


