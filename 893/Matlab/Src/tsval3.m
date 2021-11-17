function [vx,vy,vz,ier] = tsval3 (t,x,y,z,xp,yp,zp,sigma, ...
                                  iflag,te)
% tsval3:  Values or derivative vectors of space curve
%
% USAGE:  [vx,vy,vz,ier] = tsval3(t,x,y,z,xp,yp,zp,sigma,iflag,te);
%
%   This function returns values or derivatives of three
% Hermite interpolatory tension splines H1, H2, and H3 that
% form the components of a parametric space curve C(t) =
% (H1(t),H2(t),H3(t)).  Refer to Functions TSPBP and
% TSPSP.
%
% On input:
%
%       T = Vector of length N containing a strictly in-
%           creasing sequence of abscissae (parameter
%           values).  N >= 2.  Refer to Function ARCL3D.
%
%       X = Vector of length N containing data values or
%           function values returned by Function SMCRV.
%           X(I) = H1(T(I)) for I = 1 to N.
%
%       Y = Vector of length N containing data values or
%           function values returned by Function SMCRV.
%           Y(I) = H2(T(I)) for I = 1 to N.
%
%       Z = Vector of length N containing data values or
%           function values returned by Function SMCRV.
%           Z(I) = H3(T(I)) for I = 1 to N.
%
%       XP = Vector of length N containing first deriva-
%            tives.  XP(I) = H1P(T(I)) for I = 1 to N,
%            where H1P denotes the derivative of H1.
%
%       YP = Vector of length N containing first deriva-
%            tives.  YP(I) = H2P(T(I)) for I = 1 to N,
%            where H2P denotes the derivative of H2.
%
%       ZP = Vector of length N containing first deriva-
%            tives.  ZP(I) = H3P(T(I)) for I = 1 to N,
%            where H3P denotes the derivative of H3.
%
%   Note that C(T(I)) = (X(I),Y(I),Z(I)) and CP(T(I)) =
% (XP(I),YP(I),ZP(I)), I = 1 to N, are data (control)
% points and derivative (velocity) vectors, respectively.
%
%       SIGMA = Vector of length N-1 containing tension fac-
%               tors whose absolute values determine the
%               balance between cubic and linear in each
%               interval.  SIGMA(I) is associated with int-
%               erval (I,I+1) for I = 1 to N-1.
%
%       IFLAG = Output option indicator:
%               IFLAG = 0 if values of H1, H2, and H3
%                         (points on the curve) are to be
%                         computed.
%               IFLAG = 1 if first derivative vectors are to
%                         be computed.  Unit tangent vectors
%                         can be obtained by normalizing
%                         these to unit vectors.
%               IFLAG = 2 if second derivative (accelera-
%                         tion) vectors are to be computed.
%                         Given a velocity vector U and
%                         acceleration vector V, the
%                         corresponding curvature vector
%                         can be computed as (U X V X U)/
%                         |U|^4 = (|U|^2*V - <U,V>*U)/|U|^4.
%               IFLAG = 3 if third derivative vectors are to
%                         be computed.  For velocity U,
%                         acceleration V, and third deriva-
%                         tive W, the torsion is det(U,V,W)/
%                         |U X V|^2.
%
%       TE = Vector of length NE containing the evaluation
%            points.  The sequence should be strictly in-
%            creasing for maximum efficiency.  Extrapolation
%            is performed if a point is not in the interval
%            [T(1),T(N)].  NE > 0.
%
% On output:
%
%       VX,VY,VZ = Vectors of size(TE) containing values,
%                  first derivatives, second derivatives,
%                  or third derivatives of H1, H2, and H3, 
%                  respectively, at the evaluation points 
%                  (unless IER < 0).  If IER = -2, VX, VY,
%                  and VZ are zeros.  If IER = -1, VX, VY, 
%                  and VZ may be only partially defined.
%
%       IER = Error indicator:
%             IER = 0  if no errors were encountered and
%                      no extrapolation occurred.
%             IER = 1  if no errors were encountered but
%                      extrapolation was required at one
%                      or more points.
%             IER = -1 if the abscissae are not in strictly
%                      increasing order.  (This error will
%                      not necessarily be detected.)
%             IER = -2 if N < 2, IFLAG < 0, IFLAG > 3, or
%                      NE < 1.
%
% Modules required by TSVAL3:  HPPVAL, HPVAL, HVAL, SNHCSH
%
%***********************************************************

n = length(t);
ne = length(te);

% Test for invalid input.

if (n < 2  ||  iflag < 0  ||  iflag > 3  ||  ne < 1)
   vx = zeros(size(te));
   vy = zeros(size(te));
   vz = zeros(size(te));
   ier = -2;
   return;
end

if (iflag == 0)
   [vx,ierx] = hval(te,t,x,xp,sigma);
   [vy,iery] = hval(te,t,y,yp,sigma);
   [vz,ierz] = hval(te,t,z,zp,sigma);
elseif (iflag == 1) 
   [vx,ierx] = hpval(te,t,x,xp,sigma);
   [vy,iery] = hpval(te,t,y,yp,sigma);
   [vz,ierz] = hpval(te,t,z,zp,sigma);
elseif (iflag == 2)
   [vx,ierx] = hppval(te,t,x,xp,sigma);
   [vy,iery] = hppval(te,t,y,yp,sigma);
   [vz,ierz] = hppval(te,t,z,zp,sigma);
else
   [vx,ierx] = hpppval(te,t,x,xp,sigma);
   [vy,iery] = hpppval(te,t,y,yp,sigma);
   [vz,ierz] = hpppval(te,t,z,zp,sigma);
end

% Convert vx, vy, and vz from columns to rows if te is a row vector.

if size(te,1) == 1
  vx = vx';
  vy = vy';
  vz = vz';
end

if (ierx > 0  ||  iery > 0  ||  ierz > 0)
   ier = 1;
end   
if (ierx < 0  ||  iery < 0  ||  ierz < 0)
   ier = -1;
end   
return;
end  % tsval3
