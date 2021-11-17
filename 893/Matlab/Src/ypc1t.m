function [u,v,w,ier] = ypc1t(x,y,z)
% ypc1t:  Local derivative estimates as unit tangent vectors
%
% USAGE:  [u,v,ier] = ypc1t(x,y) 
%     or  [u,v,w,ier] = ypc1t(x,y,z)
%
%   Given an ordered sequence of N control points (X,Y) or
% (X,Y,Z) defining a parametric curve, this function returns
% a sequence of unit tangent vectors (U,V) or (U,V,W) at the
% control points.  The tangent direction at an interior
% point is taken to be a weighted sum of the incident chord
% directions, where the weights are inverse chord lengths.
% If the curve is open, each endpoint tangent is taken to be
% the reflection about the end chord of the tangent at the
% other end of the chord.  If N = 2, both endpoint tangents
% coincide with the chord direction.
%
% On input:
%
%       X,Y,Z = Vectors of length N containing control point
%               coordinates, where N >= 2 and N >= 3 if the 
%               curve is closed.  Adjacent control points 
%               must be distinct.  The curve is assumed to
%               be closed if and only if the first and last
%               points coincide.  Z may be omitted (for a
%               planar curve), in which case W is not 
%               returned.
%
% On output:
%
%       U,V,W = Column vectors of length N containing the 
%               components of unit tangent vectors at the 
%               control points unless IER > 0.
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = 1 if N < 2 or (N < 3 and the curve is
%                     closed).
%             IER = I if control points I-1 and I coincide 
%                     for some I in the range 2 to N.
%
% Modules required by YPC1T:  None
%
%***********************************************************

n = length(x);

% In order to allow the caller to effectively omit w from the 
% output parameter list when z is omitted from the input list,
% w is assigned the same value as ier when nargin = 2.

if n < 2 
   u = 0;
   v = 0;
   w = 1;
   ier = 1;
   return;
end

% Convert input row vectors to column vectors.

x = x(:);
y = y(:);
if nargin == 2
   closed = x(1) == x(n)  &&  y(1) == y(n);
else
   z = z(:);
   closed = x(1) == x(n)  &&  y(1) == y(n)  &&  z(1) == z(n);
end

if n < 3  &&  closed
   u = 0;
   v = 0;
   w = 1;
   ier = 1;
   return;
end

i = 2:n;
im1 = i-1;
if nargin == 2
   wt = sqrt((x(i)-x(im1)).^2 + (y(i)-y(im1)).^2);
else
   wt = sqrt((x(i)-x(im1)).^2 + (y(i)-y(im1)).^2 + (z(i)-z(im1)).^2);
end
if (any(wt <= 0))
   ier = find(wt <= 0, 1) + 1;
   u = 0;
   v = 0;
   w = ier;
   return;
end
w = 0;
ier = 0;

% Treat the case N = 2.

if n == 2
   s = 1./wt(1);
   u(1) = s*(x(2)-x(1));
   v(1) = s*(y(2)-y(1));
   u(2) = u(1);
   v(2) = v(1);
   if nargin == 3
      w(1) = s*(z(2)-z(1));
      w(2) = w(1);
   end
   return;
end

% Compute interior tangent directions.  The weights wt
% are reciprocals of squared chord lengths.

wt = 1./wt.^2;
i = 2:n-1;
u(i) = wt(i-1).*(x(i)-x(i-1)) + wt(i).*(x(i+1)-x(i));
v(i) = wt(i-1).*(y(i)-y(i-1)) + wt(i).*(y(i+1)-y(i));
if nargin == 3
   w(i) = wt(i-1).*(z(i)-z(i-1)) + wt(i).*(z(i+1)-z(i));
end
if closed
   u(1) = wt(n-1).*(x(n)-x(n-1)) + wt(1).*(x(2)-x(1));
   u(n) = u(1);
   v(1) = wt(n-1).*(y(n)-y(n-1)) + wt(1).*(y(2)-y(1));
   v(n) = v(1);
   if nargin == 3
      w(1) = wt(n-1).*(z(n)-z(n-1)) + wt(1).*(z(2)-z(1));
      w(n) = w(1);
   end
else

% Open curve:  compute endpoint tangent vectors.

   if nargin == 2
      s = 2*wt(1)*(u(2)*(x(2)-x(1)) + v(2)*(y(2)-y(1)));
      u(1) = s*(x(2)-x(1)) - u(2);
      v(1) = s*(y(2)-y(1)) - v(2);
      s = 2*wt(n-1)*(u(n-1)*(x(n)-x(n-1)) + ...
                     v(n-1)*(y(n)-y(n-1)));
      u(n) = s*(x(n)-x(n-1)) - u(n-1);
      v(n) = s*(y(n)-y(n-1)) - v(n-1);
   else
      s = 2*wt(1)*(u(2)*(x(2)-x(1)) + v(2)*(y(2)-y(1)) + ...
                   w(2)*(z(2)-z(1)));
      u(1) = s*(x(2)-x(1)) - u(2);
      v(1) = s*(y(2)-y(1)) - v(2);
      w(1) = s*(z(2)-z(1)) - w(2);
      s = 2*wt(n-1)*(u(n-1)*(x(n)-x(n-1)) + ...
                     v(n-1)*(y(n)-y(n-1)) + ...
                     w(n-1)*(z(n)-z(n-1)));
      u(n) = s*(x(n)-x(n-1)) - u(n-1);
      v(n) = s*(y(n)-y(n-1)) - v(n-1);
      w(n) = s*(z(n)-z(n-1)) - w(n-1);
   end
end

% Normalize tangents to unit vectors.

if nargin == 2
   wt = 1./sqrt(u.^2 + v.^2);
   u = wt.*u;
   v = wt.*v;
else
   wt = 1./sqrt(u.^2 + v.^2 + w.^2);
   u = wt.*u;
   v = wt.*v;
   w = wt.*w;
end
   
return;

end  % ypc1t
