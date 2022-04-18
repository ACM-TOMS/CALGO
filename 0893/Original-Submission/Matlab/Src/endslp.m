function yp = endslp(x1,x2,x3,y1,y2,y3,sigma)
% endslp:  Endpoint first derivative estimate
%
% USAGE:  yp = endslp(x1,x2,x3,y1,y2,y3,sigma);
%
%   Given data values associated with a strictly increasing
% or decreasing sequence of three abscissae X1, X2, and X3,
% this function returns a derivative estimate at X1 based on
% the tension spline H(x) that interpolates the data points
% and has third derivative equal to zero at X1.  Letting S1
% denote the slope defined by the first two points, the est-
% mate is obtained by constraining the derivative of H at X1
% so that it has the same sign as S1 and its magnitude is
% at most 3*abs(S1).  If SIGMA = 0, H(x) is quadratic and
% the derivative estimate is identical to the value computed
% by Function YPC1 at the first point (or the last point
% if the abscissae are decreasing).
%
% On input:
%
%       X1,X2,X3 = Abscissae satisfying either X1 < X2 < X3
%                  or X1 > X2 > X3.
%
%       Y1,Y2,Y3 = Data values associated with the abscis-
%                  sae.  H(X1) = Y1, H(X2) = Y2, and H(X3)
%                  = Y3.
%
%       SIGMA = Tension factor associated with H in inter-
%               val (X1,X2) or (X2,X1).
%
% On output:
%
%       YP = (Constrained) derivative of H at X1, or zero
%            if the abscissae are not strictly monotonic.
%
% Module required by ENDSLP:  SNHCSH
%
%***********************************************************

dx1 = x2 - x1;
dxs = x3 - x1;
if (dx1*(dxs-dx1) <= 0)
   yp = 0;
   return;
end
sg1 = abs(sigma);
if (sg1 < 1.e-9) 

% SIGMA = 0:  H is the quadratic interpolant.

   t = (dx1/dxs)^2;
else
   sigs = sg1*dxs/dx1;
   if (sigs <= 0.5) 
%
% 0 < SIG1 < SIGS <= .5:  compute approximations to
%   COSHM1 = COSH(SIG1)-1 and COSHMS = COSH(SIGS)-1.
%
      [dummy,coshm1] = snhcsh(sg1);
      [dummy,coshms] = snhcsh(sigs);
      t = coshm1/coshms;
   else
%
% SIGS > .5:  compute T = COSHM1/COSHMS.
%
      t = exp(sg1-sigs)*((1.0-exp(-sg1))/ ...
                         (1.0-exp(-sigs)))^2;
   end
end

% The derivative of H at X1 is
%   T = ((Y3-Y1)*COSHM1-(Y2-Y1)*COSHMS)/
%       (DXS*COSHM1-DX1*COSHMS).
%
% ENDSLP = T unless T*S1 < 0 or abs(T) > 3*abs(S1).

t = ((y3-y1)*t-y2+y1)/(dxs*t-dx1);
s1 = (y2-y1)/dx1;
if (s1 >= 0) 
  yp = min([max([0.0, t]), 3.0*s1]);
else
  yp = max([min([0.0, t]), 3.0*s1]);
end
return;

end  % endslp
