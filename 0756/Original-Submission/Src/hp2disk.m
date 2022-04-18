function [z,c] = hp2disk(w,beta,x,a)
%HP2DISK Convert solution from the half-plane to one from the disk.
%	[Z,C] = HP2DISK(W,BETA,X,C) quickly transforms the solution X,C
%	of the Schwarz-Christoffel half-plane mapping parameter problem
%	to the solution Z,C of the disk problem.
%	
%	See also DISK2HP, HPPARAM, DPARAM.
%	
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
z = zeros(size(x));
if isinf(x(n))
  z(n) = 1;
  z(1:n-1) = (x(1:n-1)-i)./(x(1:n-1)+i);
else
  z = (x-i)./(x+i);
  z = z/z(n);
end
z = z./abs(z);  % enforce exactly abs(z)==1

% Recalculate constant from scratch.
mid = (z(1)+z(2))/2;
qdat = scqdata(beta,10);
c = (w(1) - w(2))/...
    (dquad(z(2),mid,2,z,beta,qdat) - dquad(z(1),mid,1,z,beta,qdat));

