function [x,y,z,ier] = h2bsp(nd,t,xk,yk,zk,sigma,per,bv1,bvn)
% h2bsp:  Convert a tension spline from Hermite to B-spline form
%
% USAGE:  [x,y,z,ier] = h2bsp(nd,t,xk,yk,zk,sigma,per,bv1,bvn);
%
%   Given a sequence of knots, knot function values, and 
% tension factors defining a parametric tension spline curve,
% this function computes a sequence of control points for 
% which the corresponding B-spline curve interpolates the
% knot function values.  The control points p_j, j = 1:n are 
% obtained by solving the linear equations
%
%    B_{i-3}(t_i)*p_{i-1} + B_{i-2}(t_i)*p_i + 
%                           B_{i-1}(t_i)*p_{i+1} = C(t_i)
%
% for i = 3:N-1 (periodic end conditions) or i = 3:N-2 
% (interpolatory end conditions), along with two end
% end conditions, where B_j is a C^2 tension B-spline basis 
% function with support on the open interval (t_j,t_{j+4})
% and triple zeros at the endpoints, and C(t_i) is the knot
% function value.  
%
%   In the case of periodic end conditions (a closed curve
% with C(t_1) = C(t_N)), p_1 = p_N and the first and last
% equations are
%
%    B_0(t_2)*p_2 + B_1(t_2)*p_3 + 
%                       B_{-1}(t_2)*p_N = C(t_2), and 
% 
%    B_{N-1}(t_N)*p_2 + B_{N-3}(t_N)*p_{N-1} +
%                       B_{N-2}(t_N)*p_N = C(t_N).
%
% This is an order-(N-1) nonsymmetric almost tridiagonal
% linear system (with nonzero elements in the lower left
% and upper right corners) for each component of the N-1
% control points.
%
%   In the case of interpolatory end conditions, we define
%
%    p_1 = C(t_1) + bv1*(t_2-t_1)/3  and 
%    p_N = C(t_N) - bvn*(t_N-t_{N-1})/3 
%
% for user-specified endpoint derivative vectors bv1 and 
% bvn, and the first and last equations are
%
%    B_0(t_2)*p_2 + B_1(t_2)*p_3 = C(t_2) - 
%                                  B_{-1}(t_2)*p_1, and
%
%    B_{N-4}(t_{N-1})*p_{N-2} + B_{N-3}(t_{N-1})*p_{N-1} =
%                     C(t_{N-1}) - B_{N-2}(t_{N-1})*p_{N+1}.
%
% This is an order-(N-2) nonsymmetric tridiagonal linear
% system for each component of the N-2 unknown control 
% points.
%
%   Refer to Function BSP2H for a means of converting a
% B-spline curve to Hermite form.
%
% On input:
%
%       ND = Number of dimensions:
%            ND = 2 if a planar curve is to be converted.
%            ND = 3 if a space curve is to be converted.
%
%       T = Vector of length N containing a non-decreasing
%           sequence of knots.  These might be computed as
%           cumulative polygonal chord lengths between the
%           knot function values:  Functions ARCL2D and 
%           ARCL3D.  N >= 2 and N >= 4 if PER = TRUE.
%
%       XK,YK,ZK = Vectors of length N containing the 
%                  Cartesian coordinates of an ordered 
%                  sequence of knot function values C(t_i)
%                  for t_i = T(i), i = 1 to N.  ZK is an 
%                  unused dummy parameter if ND = 2.  In the
%                  case of a closed curve (PER = TRUE), the 
%                  first and last points should coincide. 
%
%       SIGMA = Vector of length N-1 containing tension
%               factors.  SIGMA(i) is associated with inter-
%               val (T(i),T(i+1)) for i = 1 to N-1.  If
%               SIGMA(i) = 0, C(t) is cubic, and as SIGMA 
%               increases, C(t) approaches linear on the 
%               interval.
%
%       PER = Logical variable with value TRUE if and only
%             C(t) is a parametric periodic function with 
%             period T(N)-T(1) corresponding to a closed 
%             curve.  It is assumed without a test that the
%             first and last knot function values coincide 
%             in this case.  On output, X(1) = X(N), Y(1) =
%             Y(N), and if ND = 3, Z(1) = Z(N).
%
%       BV1,BVN = Vectors of length ND containing endpoint
%                 derivative vectors if PER = FALSE, or
%                 unused dummy parameters otherwise.
%
% On output:
%
%       X,Y,Z = Row vectors of length N containing the 
%               control points p(t_i) for t_i = T(i), 
%               i = 1 to N.  (Z = 0 if ND = 2.)
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = 1 if ND or N is outside its valid range.
%             IER = 2 if T is not non-decreasing.
%
% Modules required by H2BSP:  SNHCSH, TRISOLVE, TRISOLVP
%
%***********************************************************

n = length(t);
x = zeros(1,n);
y = x;
z = x;

% Test for errors, and compute knot interval lengths dt
%   associated with interpolatory end conditions.

if nd < 2  ||  nd > 3  ||  n < 2  ||  (n < 4  &&  per) 
   ier = 1;
   return;
end

t = t(:)';
dt = [0 diff(t) 0 0];
if any(dt < 0)
   ier = 2;
   return;
end
ier = 0;

% Store tension factors sig in 1-1 correspondence with the
% n+2 knot intervals, again assuming nonperiodic end 
% conditions.

sig = [sigma(1) sigma(:)' sigma(n-1) sigma(n-1)];

% Store knot function values as the rows of an n-1 by nd 
% or n-2 by nd array d (right hand side vectors).

xk = xk(:);
yk = yk(:);
if nd == 3
  zk = zk(:);
end

if ~per
   if nd == 2
      d = [xk(2:n-1) yk(2:n-1)];
   else
      d = [xk(2:n-1) yk(2:n-1) zk(2:n-1)];
   end
else
   if nd == 2
      d = [xk(2:n) yk(2:n)];
   else
      d = [xk(2:n) yk(2:n) zk(2:n)];
   end

% Alter the endpoint values of dt and sig for periodic
% end conditions.

   dt(1) = dt(n);
   dt(n+1) = dt(2);
   dt(n+2) = dt(3);
   sig(1) = sig(n);
   sig(n+1) = sig(2);
   sig(n+2) = sig(3);
end

% Store sdt = sig.*dt and compute functions 
%
%   q = coshm(sdt)/(sig.*sinh(sdt)), 
%   c = sinhm(sdt)/(sig.^2.*sinh(sdt)),
%
% on the n+2 intervals.

sdt = sig.*dt;
q = zeros(1,n+2);
c = q;

% Components with sig = 0 or dt = 0:

k = find(sdt < 1.e-9);
q(k) = 0.5*dt(k);
c(k) = dt(k).^2/6;

% Components with 0 < sdt < 0.5:

k = find(sdt >= 1.e-9  &  sdt <= 0.5);
[sinhm,coshm] = snhcsh(sdt(k));
sinh = sinhm + sdt(k);
q(k) = coshm./(sig(k).*sinh);
c(k) = sinhm./(sig(k).*sig(k).*sinh);

% Components with sdt > 0.5:  scale sinhm, coshm, and sinh
% by 2*exp(-sdt) to avoid overflow.

k = find(sdt > 0.5);
ems = exp(-sdt(k));
sinh = 1.0 - ems.*ems;
sinhm = sinh - 2.0*sdt(k).*ems;
coshm = (1.0-ems).^2;
q(k) = coshm./(sig(k).*sinh);
c(k) = sinhm./(sig(k).*sig(k).*sinh);

% Compute s2(i) = q(i)+q(i+1) and b(i) = (c(i)-c(i+1))/s2(i)
% for i = 1:n+1 with b(i) = 0 if s2(i) = 0, and compute
% s3(i) = b(i)-b(i+1)+dt(i+1) for i = 1:n.

i = 1:n+1;
ip1 = i+1;
s2 = q(i)+q(ip1);
b = zeros(1,n+1);
k = find(s2);
b(k) = (c(k)-c(k+1))./s2(k);
i = 1:n;
ip1 = i+1;
s3 = b(i)-b(ip1)+dt(ip1);
if n > 2

% Store the matrix diagonals in a, b, and e.
% The nonzero corner elements (periodic case) are l and u.

   a = zeros(n-2,1);
   e = a;
   l = 0;
   u = 0;
   i = 3:n;
   k = find(s2(i));
   a(k) = c(k+3)./(s2(k+2).*s3(k+1));
   i = 2:n-1;
   k = find(s2(i));
   e(k) = c(k+1)./(s2(k+1).*s3(k+1));
   if s2(2)
      u = c(3)/(s2(2)*s3(1));
   end
   i = 2:n-2;
   b = [1-u-e(1); 1-a(i-1)-e(i)];
end

if per

% Solve the order-(n-1) system for p_i, i = 2:n, stored in
% the n-1 by nd array d, and store the solution in x, y, 
% and z with p_1 = p_n.

   if s2(n)
      l = c(n)/(s2(n)*s3(n));
   end
   b(n-1) = 1-a(n-2)-l;
   d = trisolvp(a,b,e,d,l,u);
   x = [d(n-1,1); d(:,1)]';
   y = [d(n-1,2); d(:,2)]';
   if nd == 3
      z = [d(n-1,3); d(:,3)]';
   end
else
   if n > 2

% Solve the order-(n-2) system for p(i), i = 2:n-1, stored
% in the n-2 by nd array d.

      d(1,1) = d(1,1) - u*(xk(1) + bv1(1)*dt(2)/3);
      d(1,2) = d(1,2) - u*(yk(1) + bv1(2)*dt(2)/3);
      d(n-2,1) = d(n-2,1) - e(n-2)*(xk(n) - bvn(1)*dt(n)/3);
      d(n-2,2) = d(n-2,2) - e(n-2)*(yk(n) - bvn(2)*dt(n)/3);
      if nd == 3
         d(1,3) = d(1,3) - u*(zk(1) + bv1(3)*dt(2)/3);
         d(n-2,3) = d(n-2,3) - e(n-2)*(zk(n) - bvn(3)*dt(n)/3);
      end
      a(n-2) = [];
      e(n-2) = [];
      d = trisolve(a,b,e,d);
   end

% Store the solution d, along with p_1 and p_n, in x, y, and z.

   x = [xk(1) + bv1(1)*dt(2)/3; d(:,1); xk(n) - bvn(1)*dt(n)/3]';
   y = [yk(1) + bv1(2)*dt(2)/3; d(:,2); yk(n) - bvn(2)*dt(n)/3]';
   if nd == 3
      z = [zk(1) + bv1(3)*dt(2)/3; d(:,3); zk(n) - bvn(3)*dt(n)/3]';
   end
end
return;

end  % h2bsp
