function [v,ier] = tsval1(x,y,yp,sigma,iflag,te)
% tsval1:  Values or derivatives of a tension spline
%
% USAGE:  [v,ier] = tsval1(x,y,yp,sigma,iflag,te);
%
%   This function evaluates a Hermite interpolatory ten-
% sion spline H or its first, second, or third derivative
% at a set of points TE.
%
% On input:
%
%       X = Vector of length N containing the abscissae.
%           These must be in strictly increasing order:
%           X(I) < X(I+1) for I = 1 to N-1.  N >= 2.
%
%       Y = Vector of length N containing data values or
%           function values returned by Function SMCRV.
%           Y(I) = H(X(I)) for I = 1 to N.
%
%       YP = Vector of length N containing first deriva-
%            tives.  YP(I) = HP(X(I)) for I = 1 to N, where
%            HP denotes the derivative of H.
%
%       SIGMA = Vector of length N-1 containing tension fac-
%               tors whose absolute values determine the
%               balance between cubic and linear in each
%               interval.  SIGMA(I) is associated with int-
%               erval (I,I+1) for I = 1 to N-1.
%
%       IFLAG = Output option indicator:
%               IFLAG = 0 if values of H are to be computed.
%               IFLAG = 1 if first derivative values are to
%                         be computed.
%               IFLAG = 2 if second derivative values are to
%                         be computed.
%               IFLAG = 3 if third derivative values are to 
%                         be computed.
%
%       TE = Vector of length NE containing the evaluation
%            points.  The sequence should be strictly in-
%            creasing for maximum efficiency.  Extrapolation
%            is performed if a point is not in the interval
%            [X(1),X(N)].  NE > 0.
%
% On output:
%
%       V = Vector of size(TE) containing function, first 
%           derivative, second derivative, or third deriva-
%           tive values at the evaluation points unless IER 
%           < 0.  If IER = -2, V is zeros.  If IER = -1, V 
%           may be only partially defined.
%
%       IER = Error indicator:
%             IER = 0  if no errors were encountered and
%                      no extrapolation occurred.
%             IER > 0  if no errors were encountered but
%                      extrapolation was required at IER
%                      points.
%             IER = -1 if the abscissae are not in strictly
%                      increasing order.  (This error will
%                      not necessarily be detected.)
%             IER = -2 if N < 2, IFLAG < 0, IFLAG > 3, or
%                      NE < 1.
%
% Modules required by TSVAL1:  HPPPVAL, HPPVAL, HPVAL, HVAL,
%                                SNHCSH
%
%***********************************************************

n = length(x);
ne = length(te);

% Test for invalid input.

if (n < 2  ||  iflag < 0  ||  iflag > 3  ||  ne < 1)
   v = zeros(size(te));
   ier = -2;
   return;
end

if (iflag == 0)
   [v,ier] = hval(te,x,y,yp,sigma);
elseif (iflag == 1) 
   [v,ier] = hpval(te,x,y,yp,sigma);
elseif (iflag == 2)
   [v,ier] = hppval(te,x,y,yp,sigma);
else
   [v,ier] = hpppval(te,x,y,yp,sigma);
end

% Convert v from a column to a row if te is a row vector.

if size(te,1) == 1
  v = v';
end

return;
end  % tsval1
