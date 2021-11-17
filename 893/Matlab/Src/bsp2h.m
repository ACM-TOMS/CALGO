function [xk,yk,zk,xp,yp,zp,ier] = bsp2h(nd,t,x,y,z, ...
         sigma,per,bv1,bvn)
% bsp2h:  Convert a tension spline from B-spline to Hermite form
%
% USAGE:  [xk,yk,zk,xp,yp,zp,ier] = bsp2h(nd,t,x,y,z, ...
%         sigma,per,bv1,bvn);
%
%   Given a sequence of knots, control points, and tension 
% factors defining a parametric tension spline curve in 
% B-spline form, this function converts the curve to Hermite
% form, returning knot function values and derivative 
% vectors.  The B-spline form is
%
%    C(t) = Sum {B_{j-2}(t)*p_j},
%
% where the sum is over j = 0 to N+1, p_j is a control point
% in two or three dimensions, and B_j is a C^2 tension 
% B-spline basis function with support on the open interval
% (t_j,t_{j+4}), triple zeros at the endpoints; i.e., B_j is
% a generalization of the cubic B-spline.  The outputs are
% associated with knots t_i, i = 1 to N.  
%
%   In the case of a closed curve with p_1 = p_N and periodic
% end conditions we have p_0 = p_{N-1} and p_{N+1} = p_2.
% In the non-periodic case the first and last control points
% are defined by endpoint derivative vectors C'(t_1) = bv1 
% and C'(t_N) = bvn, and we use quadruple knots at the end-
% points, resulting in endpoint interpolation:  
%
%    C(t_1) = p_0 = p_1 - bv1*(t_2-t_1)/3 and 
%    C(t_N) = p_{N+1} = p_N + bvn*(t_N-t_{N-1})/3. 
%
% Note that, with the exception of the endpoints, duplicating 
% control points reduces geometric continuity, while duplica-
% ting knots reduces parametric continuity.  Thus, with no
% duplicated knots, if two adjacent control points coincide, 
% continuity of curvature is lost at the control point, and 
% three coincident control points in sequence results in a 
% discontinuous tangent direction (a corner).  On the other
% hand, if knots are computed as cumulative chord length,
% then a pair of duplicate control points implies duplicate
% knots and a C^1/G^0 corner with a pair of duplicate knot
% function values each with a zero derivative vector.  The
% derivative vectors are perturbed away from zero in this
% case to avoid warning messages when depicting them with
% arrows (quiver3 objects).
%
%   Refer to Function H2BSP for a means of computing a 
% sequence of control points for which the corresponding
% B-spline curve interpolates a specified sequence of
% knot function values.
%
%   The tension splines may be evaluated by Function
% TSVAL2 (or TSVAL3) or Functions HVAL (values), HPVAL
% (first derivatives), HPPVAL (second derivatives), 
% HPPPVAL (third derivatives), and TSINTL (integrals).
%
% On input:
%
%       ND = Number of dimensions:
%            ND = 2 if a planar curve is to be converted.
%            ND = 3 if a space curve is to be converted.
%
%       T = Vector of length N containing a non-decreasing
%           sequence of knots.  These might be computed as
%           cumulative polygonal chord length in the control
%           polygon:  Functions ARCL2D and ARCL3D.  N >= 2
%           and N >= 4 if PER = TRUE.
%
%       X,Y,Z = Vectors of length N containing the Cartesian
%               coordinates of an ordered sequence of control
%               points p(i), i = 1 to N.  Z is an unused 
%               dummy parameter if ND = 2.  In the case of a
%               closed curve (PER = TRUE), the first and 
%               last points should coincide. 
%
%       SIGMA = Vector of length N-1 containing tension
%               factors.  SIGMA(i) is associated with inter-
%               val (T(i),T(i+1)) for i = 1 to N-1.  If
%               SIGMA(i) = 0, C(t) is cubic, and as SIGMA 
%               increases, C(t) approaches linear on the 
%               interval.
%
%       PER = Logical variable with value TRUE if and only
%             C(t) is to be a periodic function with period
%             T(N)-T(1) corresponding to a closed curve.  It
%             is assumed without a test that the first and 
%             last control points coincide in this case.  On
%             output, XK(1) = XK(N), YK(1) = YK(N), XP(1) = 
%             XP(N), YP(1) = YP(N), and, if ND = 3, then
%             ZK(1) = ZK(N) and ZP(1) = ZP(N).
%
%       BV1,BVN = Endpoint derivative vectors with ND 
%                 components if PER = FALSE, or unused
%                 dummy parameters otherwise.
%
% On output:
%
%       XK,YK,ZK = Row vectors of length N containing the 
%                  knot function values C(t_i) for t_i = 
%                  T(i), i = 1 to N.  (ZK = 0 if ND = 2.)
%
%       XP,YP,ZP = Row vectors of length N containing the 
%                  knot derivative vectors C'(t_i) for 
%                  t_i = T(i), i = 1 to N.  (ZP = 0 if
%                  ND = 2.)
%
%       IER = Error indicator:
%             IER = 0 if no errors were encountered.
%             IER = 1 if ND or N is outside its valid range.
%             IER = 2 if T is not non-decreasing.
%
% Module required by BSP2H:  SNHCSH
%
%***********************************************************

n = length(t);
xk = zeros(1,n);
yk = xk;
zk = xk;
xp = xk;
yp = xk;
zp = xk;

% Test for errors, and compute knot interval lengths dt
% associated with interpolatory end conditions.

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
if per

% Alter the endpoint values of dt and sig for periodic
% end conditions.

   dt(1) = dt(n);
   dt(n+1) = dt(2);
   dt(n+2) = dt(3);
   sig(1) = sig(n);
   sig(n+1) = sig(2);
   sig(n+2) = sig(3);
end

% Store control points as the columns of an nd by n+1 
% array p.  The last column is used only if per = true.

if nd == 2
   p = [x(:)' x(2);
        y(:)' y(2)];
else
   p = [x(:)' x(2);
        y(:)' y(2);
        z(:)' z(2)];
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

% Compute s2(i) = q(i)+q(i+1) and d(i) = (c(i)-c(i+1))/s2(i)
% for i = 1:n+1 with d(i) = 0 if s2(i) = 0, and compute
% s3(i) = d(i)-d(i+1)+dt(i+1) for i = 1:n.

i = 1:n+1;
ip1 = i+1;
s2 = q(i)+q(ip1);
d = zeros(1,n+1);
k = find(s2);
d(k) = (c(k)-c(k+1))./s2(k);
i = 1:n;
ip1 = i+1;
s3 = d(i)-d(ip1)+dt(ip1);

% Compute the coefficients d(i) of p(i) and e(i) of p(i+2)
% in the expression for the knot function values at t(i)
% for i = 1:n-1.  The coefficient of p(i+1) is 1-d(i)-e(i).

i = 1:n-1;
ip1 = i+1;
ip2 = i+2;
d = zeros(1,n-1);
e = d;
k = find(s2(ip1));
d(k) = c(k+2)./(s2(k+1).*s3(k));
e(k) = c(k+1)./(s2(k+1).*s3(k+1));

% Compute the knot function values at t(i), i = 2:n.

xk(ip1) = p(1,ip1) - d.*(p(1,ip1)-p(1,i)) + ...
                     e.*(p(1,ip2)-p(1,ip1));
yk(ip1) = p(2,ip1) - d.*(p(2,ip1)-p(2,i)) + ...
                     e.*(p(2,ip2)-p(2,ip1));
if nd == 3
   zk(ip1) = p(3,ip1) - d.*(p(3,ip1)-p(3,i)) + ...
                        e.*(p(3,ip2)-p(3,ip1));
end

% Compute the derivative vectors at t(i), i = 2:n.

d(k) = q(k+2)./(s2(k+1).*s3(k));
e(k) = q(k+1)./(s2(k+1).*s3(k+1));
xp(ip1) = d.*(p(1,ip1)-p(1,i)) + e.*(p(1,ip2)-p(1,ip1));
yp(ip1) = d.*(p(2,ip1)-p(2,i)) + e.*(p(2,ip2)-p(2,ip1));
if nd == 3
   zp(ip1) = d.*(p(3,ip1)-p(3,i)) + e.*(p(3,ip2)-p(3,ip1));
end

% Compute the endpoint knot values.

if per
   xk(1) = xk(n);
   yk(1) = yk(n);
   xp(1) = xp(n);
   yp(1) = yp(n);
   if nd == 3
      zk(1) = zk(n);
      zp(1) = zp(n);
   end
else
   xk(1) = x(1) - bv1(1)*dt(2)/3;
   yk(1) = y(1) - bv1(2)*dt(2)/3;
   xp(1) = bv1(1);
   yp(1) = bv1(2);
   xk(n) = x(n) + bvn(1)*dt(n)/3;
   yk(n) = y(n) + bvn(2)*dt(n)/3;
   xp(n) = bvn(1);
   yp(n) = bvn(2);   
   if nd == 3
      zk(1) = z(1) - bv1(3)*dt(2)/3;
      zp(1) = bv1(3);
      zk(n) = z(n) + bvn(3)*dt(n)/3;
      zp(n) = bvn(3);
   end
end

% Perturb xp so that all components are significant relative
% to the corresponding components of xk.  This is necessary
% to avoid divide-by-zero warnings in calls to quiver3.

k = find(xp == 0);
xp(k) = 4*eps(xk(k));
return;

end  % bsp2h
