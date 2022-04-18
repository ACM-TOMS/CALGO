function tspack
% tspack:  Tension Spline Package
%
% R. Renka
% 02/01/2008
%
%   This program provides a quick easily-verified test of
% all TSPACK modules except the high-level interface
% functions TSPxx and TSVALx.  The primary test uses data
% values taken from the quadratic function f(x) = x^2 for
% which knot-derivative estimation and interpolation with no
% tension is exact.  Since this test function is positive,
% strictly increasing, and convex, zero tension is suffici-
% ent to preserve these properties in the data.  The
% abscissae are taken to be a set of N points uniformly
% distributed in [0,1]
%
%   The maximum absolute error in the derivative estimates
% is printed for each of the following.
%
%      YPC1 (local approximation of knot-derivatives)
%      YPC2 with end conditions from YPC1 estimates
%      YPC2 with specified endpoint first derivatives
%      YPC2 with specified endpoint second derivatives
%      YPC2 with end conditions computed by ENDSLP
%
%   The maximum error in a tension factor is printed for
% each of the following.  These should all be zero since no
% tension is required in any of the cases.
%
%      SIGBI (minimum tension factor required to satisfy
%             constraints defined by array B)
%      SIGS (minimum tension factor required to preserve
%            monotonicity and convexity on each interval)
%      SIG0 with a lower bound of Y(I) on (X(I),X(I+1))
%      SIG1 with a lower bound of YP(I) on (X(I),X(I+1))
%      SIG2 (minimum tension factor required to preserve
%            convexity on each interval)
%
%   The maximum absolute error on a set of 3*(N-1)+1 points
% uniformly distributed in [0,1] is printed for each of the
% following evaluation functions.
%
%      HVAL (function values)
%      HPVAL (first derivative values)
%      HPPVAL (second derivative values)
%      HPPPVAL (third derivative values)
%      TSINTL (integral from 0 to T for each evaluation
%              point T)
%
%   The following six functions are used to construct para-
% metric tension spline fits to N points uniformly
% distributed on the unit circle.  In the case of SMCRV,
% periodic end conditions are used to fit cos(A) and natural
% end conditions are used for sin(A), where A is the parame-
% ter value (angle).  Thus, both B2TRI and B2TRIP are exer-
% cized.  The weights W are based on a standard deviation of
% EPS, for machine precision EPS, so that the smoothing
% curves are close to the interpolatory curves.  The maximum
% distance from the circle to a set of 3*(N-1)+1 evaluation
% points (angles uniformly distributed in [0,2*PI]) is
% printed for each function.
%
%      YPC1T (local approximations to knot-derivative vec-
%             tors with periodic end conditions)
%      YPC1P (local approximations to knot-derivatives with
%             periodic end conditions)
%      YPC2P (global approximations to knot-derivatives with
%             periodic end conditions)
%      SIGBP (minimum tension factor required to satisfy
%             bounds on distance between line segments and
%             planar curve segments)
%      SIGSP (minimum tension factor required to preserve
%             local convexity of the data)
%      SMCRV (global approximations to both function values
%             and first derivatives at the knots)
%
%   A tension spline interpolant of the N points on the unit
% circle is converted to a B-spline representation by H2BSP
% and back to its interpolatory representation by BSP2H, and
% the maximum change in a knot function value is printed.
%
%      BSP2H (converts a B-spline curve to a Hermite inter-
%             polatory spline curve)
%      H2BSP (converts a parametric interpolatory tension
%             spline curve to a B-spline curve)
%
%   The final test consists of computing a set of parameter
% values associated with the N points on the unit circle.
% These are computed by both ARCL2D and ARCL3D with constant
% Z values, and the maximum difference is printed.
%
%      ARCL2D (cumulative arc lengths for a planar curve)
%      ARCL3D (cumulative arc lengths for a space curve)

% Global variable defining maximum tension factor:

global SBIG
SBIG = 85.0;

n = 33;
nm1 = n-1;

% Print a heading.

fprintf(1, ['\n\n\n',repmat(' ',1,21),'TSPACK Test Program', ...     
        '\n\n\n  The name of each tested module is followed by the ', ...
        'maximum\n  absolute error.  This should be a small ', ...
        'multiple of the\n  machine precision except in the ', ...
        'case of YPC1T, YPC1P, YPC2P,\n  and SMCRV which may ', ...
        'have errors as large as 1.e-4.\n\n\n']);

% Compute abscissae X uniformly distributed in [0,1],
%   ordinates Y from Y(I) = f(X(I)), for f(x) = x^2,
%   true derivative values YPTRUE from f'(x) = 2*x,
%   zero tension factors SIGMA, bounds B for SIGBI, and
%   uniform weights W for SMCRV.

x = (0:1/nm1:1)';
y = x.^2;
yptrue = 2.0*x;
sigma = zeros(1,nm1);
b = zeros(5,n);
b(1,:) = 1.0;
b(2,:) = y';
b(3,:) = 2.0;
b(4,:) = yptrue';
b(5,:) = 1.0;
w = (1/eps^2)*ones(n,1);

% Test YPC1 (ISL=-1) and YPC2 with all four types of end
%   conditions.  The knot-derivative estimates should agree
%   with YPTRUE in all four cases.

bv1 = 0;
bvn = 2.0;
for isl = -1:3
   if (isl < 0) 
%
% * YPC1 test:
%
      [yp,ier] = ypc1(x,y);
   else
      if (isl == 2), bv1 = 2.0; end
%
% * YPC2 test:
%
      [yp,ier] = ypc2(x,y,sigma,isl,isl,bv1,bvn);
   end
   if (ier ~= 0) 

% Error in YPC1 or YPC2.

      if (isl < 0) 
         fprintf(1,'\n\n\n          *** Error in YPC1:  IER = %0.0f\n', ier);   
      else
         fprintf(1,['\n\n\n          *** Error in YPC2 with ', ...
                 'ISL1 = ISLN = %0.0f:  IER = %0.0f\n'], isl, ier);   
      end
      return;
   end
   err = max(abs(yp-yptrue));
   if (isl < 0) 
      fprintf(1,[repmat(' ',1,15),'YPC1',repmat(' ',1,18),'%8.2e\n'], err);
   else
      fprintf(1,[repmat(' ',1,15),'YPC2, ISL1=ISL2=%0.0f', ...
              repmat(' ',1,5),'%8.2e\n'], isl, err);
   end
end 

% Test SIGBI, SIGS, SIG0, SIG1, and SIG2 with bounds
%   which are satisfied by the cubic (zero tension).

tol = 0;

% * SIGBI test:

bmax = 3.0;
[sigma,icflg,dsmax,ier] = sigbi(x,y,yptrue,tol,b,bmax,sigma);
if (ier < 0) 
   fprintf(1,[repmat(' ',1,10),'*** Error in SIGBI, IER = %0.0f\n'], ier);
   return;
end
if (any(icflg))
   i = find(icflg,1);
   fprintf(1,[repmat(' ',1,10),'*** Error in SIGBI.  Constraint %0.0f', ...
           'in interval %0.0f is invalid.\n'], icflg(i), i);
   return; 
end   
fprintf(1,[repmat(' ',1,15),'SIGBI',repmat(' ',1,17),'%8.2e\n'], dsmax);

% * SIGS test:

[sigma,dsmax,ier] = sigs(x,y,yptrue,tol,sigma);
if (ier < 0)
   fprintf(1,[repmat(' ',1,10),'*** Error in SIGS.  IER = %0.0f\n'], ier);
   return;
end
fprintf(1,[repmat(' ',1,15),'SIGS',repmat(' ',1,18),'%8.2e\n'], dsmax); 

% * SIG0 test:

ifl = -1;
err = 0;
for i = 1:nm1
   ip1 = i+1;
   bnd = b(2,i);
   [sig,ier] = sig0(x(i),x(ip1),y(i),y(ip1),yptrue(i), ...
                    yptrue(ip1),ifl,bnd,tol);
   if (ier < 0) 
      fprintf(1,[repmat(' ',1,10),'*** Error in SIG0, Interval %0.0f', ...
              ', IER = %0.0f\n'], i, ier);
      return;
   end
   err = max([err, sig]);
end
fprintf(1,[repmat(' ',1,15),'SIG0',repmat(' ',1,18),'%8.2e\n'], err);

% * SIG1 test:

err = 0;
for i = 1:nm1
   ip1 = i+1;
   bnd = b(4,i);
   [sig,ier] = sig1(x(i),x(ip1),y(i),y(ip1),yptrue(i), ...
                    yptrue(ip1),ifl,bnd,tol);
   if (ier ~= 0) 
      fprintf(1,[repmat(' ',1,10),'*** Error in SIG1, Interval %0.0f', ...
              ', IER = %0.0f\n'], i, ier);
      return;
   end
   err = max([err, sig]);
end
fprintf(1,[repmat(' ',1,15),'SIG1',repmat(' ',1,18),'%8.2e\n'], err);

% * SIG2 test:

err = 0;
ifl = 1;
for i = 1:nm1
   ip1 = i+1;
   [sig,ier] = sig2(x(i),x(ip1),y(i),y(ip1),yptrue(i), ...
                    yptrue(ip1),ifl,tol);
   if (ier ~= 0) 
      fprintf(1,[repmat(' ',1,10),'*** Error in SIG2, Interval %0.0f', ...
              ', IER = %0.0f\n'], i, ier);
      return;
   end
   err = max([err, sig]);
end
fprintf(1,[repmat(' ',1,15),'SIG2',repmat(' ',1,18),'%8.2e\n'], err);

% Test the evaluation functions HVAL, HPVAL, HPPVAL, HPPPVAL, and
%   TSINTL with SIGMA = 0.

sigma = zeros(1,nm1);

%   The number of evaluation points NPTS is taken to be
%     three per subinterval uniformly distributed in [0,1].

npts = 3*nm1 + 1;
dt = 1.0/(npts-1);
t = (0:dt:1)';

% * HVAL test:

[h,ier] = hval(t,x,y,yptrue,sigma);
if (ier < 0) 
   fprintf(1,[repmat(' ',1,10),'*** Error in HVAL:  IER = %0.0f\n'], ier);
end
ht = t.*t;
err0 = max(abs(ht-h));

% * HPVAL test:

[hp,ier] = hpval(t,x,y,yptrue,sigma);
if (ier < 0) 
   fprintf(1,[repmat(' ',1,10),'*** Error in HPVAL:  IER = %0.0f\n'], ier);
end
hpt = 2.0*t;
err1 = max(abs(hpt-hp));

% * HPPVAL test:

[hpp,ier] = hppval(t,x,y,yptrue,sigma);
if (ier < 0) 
   fprintf(1,[repmat(' ',1,10),'*** Error in HPPVAL:  IER = %0.0f\n'], ier);
end
hppt = 2.0;
err2 = max(abs(hppt-hpp));

% * HPPPVAL test:

[hppp,ier] = hpppval(t,x,y,yptrue,sigma);
if (ier < 0) 
   fprintf(1,[repmat(' ',1,10),'*** Error in HPPPVAL:  IER = %0.0f\n'], ier);
end
hpppt = 0.0;
err3 = max(abs(hpppt-hppp));

% * TSINTL test:

hi = zeros(npts,1);
for k = 1:npts
   [hi(k),ier] = tsintl(0.0,t(k),x,y,yptrue,sigma);
   if (ier < 0) 
      fprintf(1,[repmat(' ',1,10),'*** Error in TSINTL:  IER = %0.0f\n'], ier);
      return;
   end
end
hit = t.^3./3.0;
erri = max(abs(hit-hi));

fprintf(1,[repmat(' ',1,15),'HVAL',repmat(' ',1,18),'%8.2e\n'], err0);
fprintf(1,[repmat(' ',1,15),'HPVAL',repmat(' ',1,17),'%8.2e\n'], err1);
fprintf(1,[repmat(' ',1,15),'HPPVAL',repmat(' ',1,16),'%8.2e\n'], err2);
fprintf(1,[repmat(' ',1,15),'HPPPVAL',repmat(' ',1,15),'%8.2e\n'], err3);
fprintf(1,[repmat(' ',1,15),'TSINTL',repmat(' ',1,16),'%8.2e\n'], erri);

% Test YPC1T, YPC1P, YPC2P, SIGBP, SIGSP, and SMCRV (with 
%   both periodic and natural end conditions) by approxima-
%   ting a circle with parametric tension splines.
%
%   Compute N points on the unit circle along with upper and
%     lower bounds B for SIGBP.  Parameter values are angles
%     A uniformly distributed in [0,2*PI].  In each interval
%     the upper bound (distance toward the center) is EPS,
%     and the lower bound is the orthogonal distance from
%     the midpoint of the line segment to the circle.

da = 2*pi/nm1;
a = (0:da:2*pi)';
x = cos(a);
x(n) = x(1);
y = sin(a);
y(n) = y(1);
t = (a(1:nm1) + a(2:n))/2;
bu = ones(nm1,1)*eps;
bl = -sqrt( (cos(t)-(x(1:nm1)+x(2:n))/2).^2 + ...
            (sin(t)-(y(1:nm1)+y(2:n))/2).^2 );

% * YPC1T, YPC1P, YPC2P, SIGBP, SIGSP test:

for ic = 0:5
   if (ic == 0  ||  ic == 3)
      [xp,yp,ier] = ypc1t(x,y);
      if (ier ~= 0) 
         fprintf(1,[repmat(' ',1,10),'*** Error in YPC1T.  IER = %0.0f\n'], ...
                 ier); 
         return;
      end
   elseif (ic == 1  ||  ic == 4)
      [xp,ier] = ypc1p(a,x);
      if (ier ~= 0) 
         fprintf(1,[repmat(' ',1,10),'*** Error in YPC1P.  IER = %0.0f\n'], ...
                 ier); 
         return;
      end
      [yp,ier] = ypc1p(a,y);
      if (ier ~= 0) 
         fprintf(1,[repmat(' ',1,10),'*** Error in YPC1P.  IER = %0.0f\n'], ...
                 ier); 
         return;
      end
   else
      [xp,ier] = ypc2p(a,x,sigma);
      if (ier ~= 0)
         fprintf(1,[repmat(' ',1,10),'*** Error in YPC2P.  IER = %0.0f\n'], ...
                 ier); 
         return;
      end
      [yp,ier] = ypc2p(a,y,sigma);
      if (ier ~= 0)
         fprintf(1,[repmat(' ',1,10),'*** Error in YPC2P.  IER = %0.0f\n'], ...
                 ier); 
         return;
      end
   end

%   Compute tension factors SIGMA.
%     Note that, in the case of YPC2P, XP and YP are
%     not recomputed following the change in SIGMA,
%     and the splines are therefore not C-2.

   if (ic < 3)

%   Compute bounds-constrained tension factors SIGMA.

      bmax = 1.0;
      [sigma,dsmax,ier] = sigbp(x,y,xp,yp,tol,bl,bu,bmax,sigma);
      if (ier < 0)
         fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in SIGBP.  ', ...
                 'IER = %0.0f\n'], ier); 
         return;
      end
   else

%   Compute shape-preserving tension factors SIGMA.

      z = 0;   % Unused dummy variables
      zp = 0;
      [sigma,dsmax,ier] = sigsp(2,a,x,y,z,xp,yp,zp,tol,sigma);
      if (ier < 0)
         fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in SIGSP.  ', ...
                 'IER = %0.0f\n'], ier); 
         return;
      end
   end

%   Compute evaluation points.

   npts = 3*nm1 + 1;
   dt = 2*pi/(npts-1);
   t = (0:dt:2*pi)';
   [xt,ier] = hval(t,a,x,xp,sigma);
   if (ier < 0) 
      fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in HVAL.  ', ...
              'IER = %0.0f\n'], ier); 
      return;
   end
   [yt,ier] = hval(t,a,y,yp,sigma);
   if (ier < 0)
      fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in HVAL.  ', ...
              'IER = %0.0f\n'], ier); 
      return;
   end
   err = max(sqrt((cos(t)-xt).^2 + (sin(t)-yt).^2));

   if (ic == 0) 
      fprintf(1,[repmat(' ',1,15),'YPC1T/SIGBP',repmat(' ',1,11), ...
              '%8.2e\n'], err);
   elseif (ic == 1) 
      fprintf(1,[repmat(' ',1,15),'YPC1P/SIGBP',repmat(' ',1,11), ...
              '%8.2e\n'], err);
   elseif (ic == 2)
      fprintf(1,[repmat(' ',1,15),'YPC2P/SIGBP',repmat(' ',1,11), ...
              '%8.2e\n'], err);
   elseif (ic == 3) 
      fprintf(1,[repmat(' ',1,15),'YPC1T/SIGSP',repmat(' ',1,11), ...
              '%8.2e\n'], err);
   elseif (ic == 4) 
      fprintf(1,[repmat(' ',1,15),'YPC1P/SIGSP',repmat(' ',1,11), ...
              '%8.2e\n'], err);
   else
      fprintf(1,[repmat(' ',1,15),'YPC2P/SIGSP',repmat(' ',1,11), ...
              '%8.2e\n'], err);
   end
end

% * SMCRV test:

sm = n;
smtol = sqrt(2.0/n);
period = true;
[xs,xp,ier] = smcrv(a,x,sigma,period,w,sm,smtol);
if (ier ~= 0)
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in SMCRV.  ', ...
              'IER = %0.0f\n'], ier); 
   return;
end
period = false;
[ys,yp,ier] = smcrv(a,y,sigma,period,w,sm,smtol);
if (ier ~= 0)
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in SMCRV.  ', ...
              'IER = %0.0f\n'], ier); 
   return;
end
t = (0:dt:2*pi)';
[xt,ier] = hval(t,a,xs,xp,sigma);
if (ier < 0) 
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in HVAL.  ', ...
           'IER = %0.0f\n'], ier); 
   return;
end
[yt,ier] = hval(t,a,ys,yp,sigma);
if (ier < 0)
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in HVAL.  ', ...
           'IER = %0.0f\n'], ier); 
   return;
end
err = max(sqrt((cos(t)-xt).^2 + (sin(t)-yt).^2));
fprintf(1,[repmat(' ',1,15),'SMCRV',repmat(' ',1,17),'%8.2e\n'], err);

% Copy (x,y) into (xk,yk), compute random tension factors 
% sigma, use H2BSP to compute control points (x,y) such that
% the B-spline curve (with B-spline basis functions defined 
% by a and sigma) interpolates (xk,yk), use BSP2H to compute 
% knot function values (xs,ys), and compare (xs,ys) to 
% (xk,yk).

nd = 2;
xk = x;
yk = y;
zk = 0;    % Dummy variables z, bv1, bvn
bv1 = 0;
bvn = 0;
sigma = rand*SBIG*ones(1,nm1);
period = true;

% * BSP2H/H2BSP test:

[x,y,z,ier] = h2bsp(nd,a,xk,yk,zk,sigma,period,bv1,bvn);
if (ier < 0) 
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in H2BSP.  ', ...
           'IER = %0.0f\n'], ier); 
   return;
end
[xs,ys,zs,xp,yp,zp,ier] = bsp2h(nd,a,x,y,z,sigma,period,bv1,bvn);
if (ier < 0) 
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in BSP2H.  ', ...
           'IER = %0.0f\n'], ier); 
   return;
end
x = x';     % Convert row vectors to column vectors.
y = y';
xs = xs';
ys = ys';
err = max([norm(xs-xk,inf),norm(ys-yk,inf)]);
fprintf(1,[repmat(' ',1,15),'BSP2H/H2BSP',repmat(' ',1,11),'%8.2e\n'], err);

% Test ARCL2D and ARCL3D by computing the differences
%   between cumulative arc lengths for the data sets
%   (X,Y) and (X,Y,W), where W(I) = 1 for all I.

w = ones(n,1);

% * ARCL2D/ARCL3D test:

[a,ier] = arcl2d(x,y);
if (ier ~= 0)
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in ARCL2D.  ', ...
           'IER = %0.0f\n'], ier); 
   return;
end
[wk,ier] = arcl3d(x,y,w);
if (ier ~= 0)
   fprintf(1,['\n\n\n',repmat(' ',1,10),'*** Error in ARCL3D.  ', ...
           'IER = %0.0f\n'], ier); 
   return;
end
err = max(abs(a-wk));
fprintf(1,[repmat(' ',1,15),'ARCL2D/ARCL3D',repmat(' ',1,9),'%8.2e\n'], err);
return;

end  % tspack
