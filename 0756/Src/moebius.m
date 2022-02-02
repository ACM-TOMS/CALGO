function A = moebius(z,w)
%MOEBIUS Moebius transformation parameters.
%	A = MOEBIUS(Z,W) computes the coefficients of the Moebius
%	transformation taking the 3-vector Z to W, so that
%	
%	         W = (A(1)*Z + A(2))./(A(3)*Z + A(4)).
%	
%	Infinities are not recognized and will not work.
%	
%	Written by Toby Driscoll.  Last updated 5/26/95.

t1 = -diff(z(1:2))*diff(w(2:3));
t2 = -diff(z(2:3))*diff(w(1:2));

A(1) = w(1)*t1 - w(3)*t2;
A(2) = w(3)*z(1)*t2 - w(1)*z(3)*t1;
A(3) = t1 - t2;
A(4) = z(1)*t2 - z(3)*t1;
