function [x,a] = disk2hp(w,beta,z,c)
%DISK2HP Convert solution from the disk to one from the half-plane.
%       [X,C] = DISK2HP(W,BETA,Z,C) quickly transforms the solution Z,C
%       of the Schwarz-Christoffel disk mapping parameter problem to the
%       solution X,C of the half-plane problem.
%	
%	See also HP2DISK, DPARAM, HPPARAM.
%       
%	Written by Toby Driscoll.  Last updated 5/24/95.

n = length(w);
x = zeros(size(z));
x(n) = Inf;
x(1:n-1) = -i*(z(1:n-1)+1)./(z(1:n-1)-1); % Mobius transfmn
x = real(x);				  % enforce exactly imag(x)==0

% Recalculate constant from scratch.
mid = mean(x(1:2));
qdat = scqdata(beta(1:n-1),10);
a = (w(1)-w(2))/(hpquad(x(2),mid,2,x(1:n-1),beta(1:n-1),qdat) - ...
    hpquad(x(1),mid,1,x(1:n-1),beta(1:n-1),qdat));


