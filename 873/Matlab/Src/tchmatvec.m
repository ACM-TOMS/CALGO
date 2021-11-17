%
%   File:              tchmatvec.m
%
%   Description:       T(A)*v where T(A) is polynomial in A.
%                      T(w) is chosen a the Tchebyshev polynomials.
%
%   Authors:           Chao Yang (10/30/95), Marielba Rojas (10/23/2004)
%
%   Version:           1.2
%
%   System:            MATLAB 6.0 or higher 
%
%   Date:              original: 10/30/1995 modified: 15 March 2007
%
%   Functions called:  matvec
%
%   Call: y = tchmatvec(v,Balpha,a,b,maxdegree);
%

%
%  if the eigenvalues of the matrix lie in (a0,b)
%  then a0 <= a <= b; a is the 'cut' to the right of the wanted eigenvalues
%       b is an estimate of the largest eigenvalue; has to be upper bound though
%

function [y] = tchmatvec(v,Balpha,a,b,maxdegree)

%
% ignore polynomials of degree lower than 3
%

T(:,1) = v;
w = matvec(v,Balpha);

T(:,2) = ((a+b)/(a-b)*v+2/(b-a)*w); 

for k = 3:maxdegree
   w = matvec(T(:,2),Balpha);
   T(:,3) = 2*( (a+b)/(a-b)*T(:,2)+2/(b-a)*w ) - T(:,1);
   T(:,1) = T(:,2);
   T(:,2) = T(:,3);
end;

y = T(:,3);
 
