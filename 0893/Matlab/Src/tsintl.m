function [hi,ier] = tsintl(a,b,x,y,yp,sigma)
% tsintl:  Integral of a Hermite interpolatory tension spline
%
% USAGE:  [hi,ier] = tsintl(a,b,x,y,yp,sigma);
%
%   This function computes the integral from A to B of a
% Hermite interpolatory tension spline H.
%
% On input:
%
%       A,B = Lower and upper limits of integration, respec-
%             tively.  Note that -TSINTL(B,A,...) = 
%             TSINTL(A,B,...).
%
%       X = Vector of length N containing the abscissae.
%           These must be in strictly increasing order:
%           X(I) < X(I+1) for I = 1 to N-1.  N >= 2.
%
%       Y = Vector of length N containing data values.
%           H(X(I)) = Y(I) for I = 1 to N.
%
%       YP = Vector of length N containing first deriva-
%            tives.  HP(X(I)) = YP(I) for I = 1 to N, where
%            HP denotes the derivative of H.
%
%       SIGMA = Vector of length N-1 containing tension fac-
%               tors whose absolute values determine the
%               balance between cubic and linear in each
%               interval.  SIGMA(I) is associated with int-
%               erval (I,I+1) for I = 1 to N-1.
%
% On output:
%
%       HI = Integral of H from A to B, or zero if IER < 0.
%
%       IER = Optional error indicator:
%             IER = 0  if no errors were encountered and
%                      X(1) <= T <= X(N) for T = A and
%                      T = B, or A = B.
%             IER = 1  if no errors were encountered but
%                      extrapolation was necessary:  A or B
%                      is not in the interval (X(1),X(N)).
%             IER = -1 if the abscissae are not in strictly
%                      increasing order.  Only those in or
%                      adjacent to an interval of integra-
%                      tion are tested.
%
% Modules required by TSINTL:  INTRVL, SNHCSH
%
%***********************************************************

global SBIG

n = length(x);

% Accumulate the integral from XL to XU in SUM.

xl = min([a,b]);
xu = max([a,b]);
sum = 0;

% Find left-end indices of intervals containing XL and XU.
%   If XL < X(1) or XU > X(N), extrapolation is performed
%   using the leftmost or rightmost interval.

il = intrvl(xl,x);
iu = intrvl(xu,x);
if (nargout > 1)
   ier = (xl < x(1)  ||  xu > x(n));
end
if (xl == xu)
   hi = 0;
   return;
end
ilp1 = il + 1;

% Compute the integral from XL to X(IL+1).

dx = x(ilp1) - x(il);
if (dx <= 0)
   hi = 0;
   ier = -1;
   return;
end
u = x(ilp1) - xl;
b1 = u/dx;
y2 = y(ilp1);
s = (y2-y(il))/dx;
s2 = yp(ilp1);
d1 = s - yp(il);
d2 = s2 - s;
sig = abs(sigma(il));
if (sig < 1.e-9) 

%   SIG = 0.

   sum = sum + u*(y2 - u*(6.0*s2 - b1*(4.0*d2 + ...
               (3.0*b1-4.0)*(d1-d2)))/12.0);
elseif (sig <= 0.5) 

%   0 < SIG <= .5.

   sb1 = sig*b1;
   [sm,cm,cmm] = snhcsh(sig);
   [sm1,cm1,cmm1] = snhcsh(sb1);
   e = sig*sm - cmm - cmm;
   sum = sum + u*(y2 - s2*u/2.0) + ((cm*cmm1-sm*sm1)* ...
               (d1+d2) + sig*(cm*sm1-(sm+sig)*cmm1)*d2)/ ...
               ((sig/dx)^2*e);
else

%   SIG > .5.

   sb1 = sig*b1;
   sb2 = sig - sb1;
   if (-sb1 > SBIG  ||  -sb2 > SBIG) 
      sum = sum + u*(y2 - s*u/2.0);
   else
      e1 = exp(-sb1);
      e2 = exp(-sb2);
      ems = e1*e2;
      tm = 1.0 - ems;
      tp = 1.0 + ems;
      t = sb1*sb1/2.0 + 1.0;
      e = tm*(sig*tp - tm - tm);
      sum = sum + u*(y2 - s2*u/2.0)+(sig*tm*(tp*t-e1-e2- ...
                  tm*sb1)*d2 - (tm*(tm*t-e1+e2-tp*sb1) + ...
                  sig*(e1*ems-e2+2.0*sb1*ems))*(d1+d2))/ ...
                  ((sig/dx)^2*e);
   end
end

% Add in the integral from X(IL+1) to X(J) for J =
%   Max(IL+1,IU).

for i = ilp1:max([il,iu-1])
   ip1 = i + 1;
   dx = x(ip1) - x(i);
   if (dx <= 0)
      hi = 0;
      ier = -1;
      return;
   end
   sig = abs(sigma(i));
   if (sig < 1.e-9) 

%   SIG = 0.

      sum = sum + dx*((y(i)+y(ip1))/2.0 - ...
                  dx*(yp(ip1)-yp(i))/12.0);
   elseif (sig <= .50) 

%   0 < SIG <= .5.

      [sm,cm,cmm] = snhcsh(sig);
      e = sig*sm - cmm - cmm;
      sum = sum + dx*(y(i)+y(ip1) - dx*e*(yp(ip1)-yp(i))/ ...
                     (sig*sig*cm))/2.0;
   else

%   SIG > .5.

      ems = exp(-sig);
      sum = sum + dx*(y(i)+y(ip1) - dx*(sig*(1.0+ems)/ ...
                     (1.0-ems)-2.0)*(yp(ip1)-yp(i))/ ...
                     (sig*sig))/2.0;
   end
end

% Add in the integral from X(IU) to XU if IU > IL.

if (il < iu) 
   iup1 = iu + 1;
   dx = x(iup1) - x(iu);
   if (dx <= 0)
      hi = 0;
      ier = -1;
      return;
   end
   u = xu - x(iu);
   if (u == 0) 
      if (a > b), sum = -sum; end
      hi = sum;
      return;
   end
   b2 = u/dx;
   y1 = y(iu);
   s = (y(iup1)-y1)/dx;
   s1 = yp(iu);
   d1 = s - s1;
   d2 = yp(iup1) - s;
   sig = abs(sigma(iu));
   if (sig < 1.e-9) 

%   SIG = 0.

      sum = sum + u*(y1 + u*(6.0*s1 + b2*(4.0*d1 + ...
                    (4.0-3.0*b2)*(d1-d2)))/12.0);
   elseif (sig <= 0.5) 

%   0 < SIG <= .5.

      sb2 = sig*b2;
      [sm,cm,cmm] = snhcsh(sig);
      [sm2,cm2,cmm2] = snhcsh(sb2);
      e = sig*sm - cmm - cmm;
      sum = sum + u*(y1 + s1*u/2.0) + ((cm*cmm2-sm*sm2)* ...
                    (d1+d2) + sig*(cm*sm2-(sm+sig)*cmm2)*d1)/ ...
                    ((sig/dx)^2*e);
   else

%   SIG > .5.

      sb2 = sig*b2;
      sb1 = sig - sb2;
      if (-sb1 > SBIG  ||  -sb2 > SBIG)
         sum = sum + u*(y1 + s*u/2.0);
      else
         e1 = exp(-sb1);
         e2 = exp(-sb2);
         ems = e1*e2;
         tm = 1.0 - ems;
         tp = 1.0 + ems;
         t = sb2*sb2/2.0 + 1.0;
         e = tm*(sig*tp - tm - tm);
         sum = sum + u*(y1 + s1*u/2.0)+(sig*tm*(tp*t-e1-e2- ...
                       tm*sb2)*d1 - (tm*(tm*t-e2+e1-tp*sb2) + ...
                       sig*(e2*ems-e1+2.0*sb2*ems))*(d1+d2))/ ...
                       ((sig/dx)^2*e);
      end
   end
else

% XL and XU are in the same interval (IL = IU), and SUM 
%   contains the integral from XL to X(IL+1).
%   Subtract off the integral from XU to X(IL+1).
%
   dx = x(ilp1) - x(il);
   if (dx <= 0)
      hi = 0;
      ier = -1;
      return;
   end
   y2 = y(ilp1);
   s = (y2-y(il))/dx;
   s2 = yp(ilp1);
   d1 = s - yp(il);
   d2 = s2 - s;
   u = x(ilp1) - xu;
   if (u == 0) 
      if (a > b), sum = -sum; end
      hi = sum;
      return;
   end
   b1 = u/dx;
   sig = abs(sigma(il));
   if (sig < 1.e-9)

%   SIG = 0.

      sum = sum - u*(y2 - u*(6.0*s2 - b1*(4.0*d2 + ...
                    (3.0*b1-4.0)*(d1-d2)))/12.0);
   elseif (sig <= 0.5)

%   0 < SIG <= .5.

      sb1 = sig*b1;
      [sm,cm,cmm] = snhcsh(sig);
      [sm1,cm1,cmm1] = snhcsh(sb1);
      e = sig*sm - cmm - cmm;
      sum = sum - u*(y2 - s2*u/2.0) - ((cm*cmm1-sm*sm1)* ...
                    (d1+d2) + sig*(cm*sm1-(sm+sig)*cmm1)*d2)/ ...
                    ((sig/dx)^2*e);
   else

%   SIG > .5.

      sb1 = sig*b1;
      sb2 = sig - sb1;
      if (-sb1 > SBIG  ||  -sb2 > SBIG) 
         sum = sum - u*(y2 - s*u/2.0);
      else
         e1 = exp(-sb1);
         e2 = exp(-sb2);
         ems = e1*e2;
         tm = 1.0 - ems;
         tp = 1.0 + ems;
         t = sb1*sb1/2.0 + 1.0;
         e = tm*(sig*tp - tm - tm);
         sum = sum - u*(y2 - s2*u/2.0)-(sig*tm*(tp*t-e1-e2- ...
                       tm*sb1)*d2 - (tm*(tm*t-e1+e2-tp*sb1) + ...
                       sig*(e1*ems-e2+2.0*sb1*ems))*(d1+d2))/ ... 
                       ((sig/dx)^2*e);
      end
   end
end

% No errors were encountered.  Adjust the sign of SUM.

if (a > b), sum = -sum; end
hi = sum;
return;

end  % tsintl
