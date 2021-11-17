function [N,y,go200]=N_inerfc(n,x,eps0,N0,dN)
%N_INERFC Number of Gaussian quadrature points needed.
%N_INERFC(n,X,EPS0,N0,DN) determines the number N of Gaussian 
%   quadrature points needed to obtain y=f_n(X) to a relative 
%   accuracy EPS0 and outputs go200=1 if that number is greater
%   than 200. It does so by successively incrementing an initial 
%   estimate N0 by DN units until the required accuracy is
%   achieved. 

go200=0;
N=N0-dN;
y=1; y0=0;
while abs((y-y0)/y)>eps0
  N=N+dN; y0=y;
  if N>200, go200=1; return; end
  y=quad_inerfc(N,n,x);
end
