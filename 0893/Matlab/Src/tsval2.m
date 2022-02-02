function [vx,vy,ier] = tsval2 (t,x,y,xp,yp,sigma,iflag,te)
% tsval2:  Values or derivative vectors of planar curve
%
% USAGE:  [vx,vy,ier] = tsval2(t,x,y,xp,yp,sigma,iflag,te);
%
%   This function returns values or derivatives of a pair
% of Hermite interpolatory tension splines H1 and H2 that
% form the components of a parametric planar curve C(t) =
% (H1(t),H2(t)).  Refer to Functions TSPBP and TSPSP.
%
% On input:
%
%       T = Vector of length N containing a strictly in-
%           creasing sequence of abscissae (parameter
%           values).  N >= 2.  Refer to Function ARCL2D.
%
%       X = Vector of length N containing data values or
%           function values returned by Function SMCRV.
%           X(I) = H1(T(I)) for I = 1 to N.
%
%       Y = Vector of length N containing data values or
%           function values returned by Function SMCRV.
%           Y(I) = H2(T(I)) for I = 1 to N.
%
%       XP = Vector of length N containing first deriva-
%            tives.  XP(I) = H1P(T(I)) for I = 1 to N,
%            where H1P denotes the derivative of H1.
%
%       YP = Vector of length N containing first deriva-
%            tives.  YP(I) = H2P(T(I)) for I = 1 to N,
%            where H2P denotes the derivative of H2.
%
%   Note that C(T(I)) = (X(I),Y(I)) and CP(T(I)) = (XP(I),
% YP(I)), I = 1 to N, are data (control) points and deriva-
% tive (velocity) vectors, respectively.
%
%       SIGMA = Vector of length N-1 containing tension fac-
%               tors whose absolute values determine the
%               balance between cubic and linear in each
%               interval.  SIGMA(I) is associated with int-
%               erval (I,I+1) for I = 1 to N-1.
%
%       IFLAG = Output option indicator:
%               IFLAG = 0 if values of H1 and H2 (points on
%                         the curve) are to be computed.
%               IFLAG = 1 if first derivative vectors are to
%                         be computed.  Unit tangent vectors
%                         can be obtained by normalizing
%                         these to unit vectors.
%               IFLAG = 2 if second derivative (accelera-
%                         tion) vectors are to be computed.
%                         Given a velocity vector U and
%                         acceleration vector V, the
%                         corresponding signed curvature
%                         can be computed as (U X V)/|U|^3.
%               IFLAG = 3 if third derivative vectors are to
%                         be computed.
%
%       TE = Vector of length NE containing the evaluation
%            points.  The sequence should be strictly in-
%            creasing for maximum efficiency.  Extrapolation
%            is performed if a point is not in the interval
%            [T(1),T(N)].  NE > 0.
%
% On output:
%
%       VX,VY = Vectors of size(TE) containing values, first
%               derivatives, second derivatives, or third 
%               derivatives of H1 and H2, respectively, at 
%               the evaluation points (unless IER < 0).  If 
%               IER = -2, VX and VY are zeros.  If IER = -1, 
%               VX and VY may be only partially defined.
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
% Modules required by TSVAL2:  HPPPVAL, HPPVAL, HPVAL, HVAL, 
%                                SNHCSH
%
%***********************************************************

n = length(t);
ne = length(te);

% Test for invalid input.

if (n < 2  ||  iflag < 0  ||  iflag > 3  ||  ne < 1)
   vx = zeros(size(te));
   vy = zeros(size(te));
   ier = -2;
   return;
end

if (iflag == 0)
   [vx,ierx] = hval(te,t,x,xp,sigma);
   [vy,iery] = hval(te,t,y,yp,sigma);
elseif (iflag == 1) 
   [vx,ierx] = hpval(te,t,x,xp,sigma);
   [vy,iery] = hpval(te,t,y,yp,sigma);
elseif (iflag == 2)
   [vx,ierx] = hppval(te,t,x,xp,sigma);
   [vy,iery] = hppval(te,t,y,yp,sigma);
else
   [vx,ierx] = hpppval(te,t,x,xp,sigma);
   [vy,iery] = hpppval(te,t,y,yp,sigma);
end

% Convert vx and vy from columns to rows if te is a row vector.

if size(te,1) == 1
  vx = vx';
  vy = vy';
end

if (ierx > 0  ||  iery > 0)
   ier = 1;
end   
if (ierx < 0  ||  iery < 0)
   ier = -1;
end   
return;
end  % tsval2
