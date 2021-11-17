function [ys,yp] = b2trip(x,y,w,p,d,sd)
% b2trip:  SPD periodic block tridiagonal system solver
%
% USAGE:  [ys,yp] = b2trip(x,y,w,p,d,sd);
%
%   This function solves the order 2(N-1) symmetric posi-
% tive-definite linear system associated with minimizing the
% quadratic functional Q(YS,YP) (described in Function
% SMCRV) with periodic end conditions.  The matrix is block
% tridiagonal except for nonzero blocks in the upper right
% and lower left corners.
%
% On input:
%
%       X = Vector of length N containing a strictly in-
%           creasing sequence of abscissae.  N >= 3.
%
%       Y,W = Vectors of length N-1 containing data values
%             and positive weights, respectively, associated
%             with the first N-1 abscissae.
%
%       P = Positive smoothing parameter defining Q.
%
%       D,SD = Vectors of length N-1 containing positive ma-
%              trix elements.  Letting DX and SIG denote the
%              width and tension factor associated with the
%              interval (X(I),X(I+1)), D(I) = SIG*(SIG*
%              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) =
%              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG)
%              - 2*COSHM(SIG).
%
%   Note that no test is made for a nonpositive value of 
% X(I+1)-X(I), W(I), D(I), or SD(I).
%
% On output:
%
%       YS,YP = Vectors of size(X) containing solution com-
%               ponents:  function and derivative values,
%               respectively, at the abscissae.  YS(N) =
%               YS(1) and YP(N) = YP(1).
%
% Modules required by B2TRIP:  None
%
%***********************************************************

m = size(x);
n = length(x);
nm1 = n - 1;
nm2 = n - 2;
nm3 = n - 3;

% Initialize output arrays.

ys = zeros(m);
yp = zeros(m);

% Work space:

t11 = zeros(nm2,1);
t12 = zeros(nm2,1);
t21 = zeros(nm2,1);
t22 = zeros(nm2,1);
u11 = zeros(nm2,1);
u12 = zeros(nm2,1);
u21 = zeros(nm2,1);
u22 = zeros(nm2,1);

% The forward elimination step consists of scaling a row by
%   the inverse of its diagonal block and eliminating the
%   subdiagonal block for the first N-2 rows.  The super-
%   diagonal is stored in T, the negative of the last column
%   in U, and the right hand side in YS,YP.  For J = 11, 12,
%   and 22, SJI and SJIM1 denote the elements in position J
%   of the superdiagonal block in rows I and I-1, respect-
%   ively.  Similarly, DJI denotes an element in the diago-
%   nal block of row I.

% I = 1:

dx = x(n) - x(nm1);
dnm1 = d(nm1);
s22nm1 = sd(nm1);
s12nm1 = -(dnm1 + s22nm1)/dx;
s11nm1 = 2.0*s12nm1/dx;
dx = x(2) - x(1);
di = d(1);
s22i = sd(1);
s12i = (di + s22i)/dx;
s11i = -2.0*s12i/dx;
r1 = p*w(1);
d11i = r1 - s11nm1 - s11i;
d12i = s12i + s12nm1;
d22i = dnm1 + di;
den = d11i*d22i - d12i*d12i;
t11(1) = (d22i*s11i + d12i*s12i)/den;
t12(1) = (d22i*s12i - d12i*s22i)/den;
t21(1) = -(d12i*s11i + d11i*s12i)/den;
t22(1) = (d11i*s22i - d12i*s12i)/den;
u11(1) = -(d22i*s11nm1 + d12i*s12nm1)/den;
u12(1) = (d12i*s22nm1 - d22i*s12nm1)/den;
u21(1) = (d12i*s11nm1 + d11i*s12nm1)/den;
u22(1) = (d12i*s12nm1 - d11i*s22nm1)/den;
r1 = r1*y(1)/den;
ys(1) = d22i*r1;
yp(1) = -d12i*r1;

% I = 2 to N-2:

for i = 2:nm2
   im1 = i - 1;
   dim1 = di;
   s22im1 = s22i;
   s12im1 = s12i;
   s11im1 = s11i;
   dx = x(i+1) - x(i);
   di = d(i);
   s22i = sd(i);
   s12i = (di + s22i)/dx;
   s11i = -2.0*s12i/dx;
   r1 = p*w(i);
   d11i = r1 - s11im1 - s11i - (s11im1*t11(im1) - ...
          s12im1*t21(im1));
   d12i = s12i - s12im1 - (s11im1*t12(im1) - s12im1*t22(im1));
   d22i = dim1 + di - (s12im1*t12(im1)+s22im1*t22(im1));
   den = d11i*d22i - d12i*d12i;
   t11(i) = (d22i*s11i + d12i*s12i)/den;
   t12(i) = (d22i*s12i - d12i*s22i)/den;
   t21(i) = -(d12i*s11i + d11i*s12i)/den;
   t22(i) = (d11i*s22i - d12i*s12i)/den;
   su11 = s11im1*u11(im1) - s12im1*u21(im1);
   su12 = s11im1*u12(im1) - s12im1*u22(im1);
   su21 = s12im1*u11(im1) + s22im1*u21(im1);
   su22 = s12im1*u12(im1) + s22im1*u22(im1);
   u11(i) = (d12i*su21 - d22i*su11)/den;
   u12(i) = (d12i*su22 - d22i*su12)/den;
   u21(i) = (d12i*su11 - d11i*su21)/den;
   u22(i) = (d12i*su12 - d11i*su22)/den;
   r1 = r1*y(i) - s11im1*ys(im1) + s12im1*yp(im1);
   r2 = -s12im1*ys(im1) - s22im1*yp(im1);
   ys(i) = (d22i*r1 - d12i*r2)/den;
   yp(i) = (d11i*r2 - d12i*r1)/den;
end

% The backward elimination step zeros the first N-3 blocks
%   of the superdiagonal.  For I = N-2,N-3 to 1, T(I) and
%   (YS(I),YP(I)) are overwritten with the negative of the
%   last column and the new right hand side, respectively.

t11(nm2) = u11(nm2) - t11(nm2);
t12(nm2) = u12(nm2) - t12(nm2);
t21(nm2) = u21(nm2) - t21(nm2);
t22(nm2) = u22(nm2) - t22(nm2);
for i = nm3:-1:1
   ip1 = i + 1;
   ys(i) = ys(i) - t11(i)*ys(ip1) - t12(i)*yp(ip1);
   yp(i) = yp(i) - t21(i)*ys(ip1) - t22(i)*yp(ip1);
   t11(i) = u11(i) - t11(i)*t11(ip1) - t12(i)*t21(ip1);
   t12(i) = u12(i) - t11(i)*t12(ip1) - t12(i)*t22(ip1);
   t21(i) = u21(i) - t21(i)*t11(ip1) - t22(i)*t21(ip1);
   t22(i) = u22(i) - t21(i)*t12(ip1) - t22(i)*t22(ip1);
end

% Solve the last equation for YS(N-1),YP(N-1).  SJI = SJNM2
%   and DJI = DJNM1.

r1 = p*w(nm1);
d11i = r1 - s11i - s11nm1 + s11nm1*t11(1) - ...
       s12nm1*t21(1) + s11i*t11(nm2) - s12i*t21(nm2);
d12i = -s12nm1 - s12i + s11nm1*t12(1) - s12nm1*t22(1) + ...
       s11i*t12(nm2) - s12i*t22(nm2);
d22i = di + dnm1 + s12nm1*t12(1) + s22nm1*t22(1) + ...
       s12i*t12(nm2) + s22i*t22(nm2);
den = d11i*d22i - d12i*d12i;
r1 = r1*y(nm1) - s11nm1*ys(1) + s12nm1*yp(1) - ...
     s11i*ys(nm2) + s12i*yp(nm2);
r2 = -s12nm1*ys(1) - s22nm1*yp(1) - s12i*ys(nm2) - s22i*yp(nm2);
ysnm1 = (d22i*r1 - d12i*r2)/den;
ypnm1 = (d11i*r2 - d12i*r1)/den;
ys(nm1) = ysnm1;
yp(nm1) = ypnm1;

% Back substitute for the remainder of the solution
%   components.

for i = 1:nm2
   ys(i) = ys(i) + t11(i)*ysnm1 + t12(i)*ypnm1;
   yp(i) = yp(i) + t21(i)*ysnm1 + t22(i)*ypnm1;
end

% YS(N) = YS(1) and YP(N) = YP(1).

ys(n) = ys(1);
yp(n) = yp(1);
return;

end  % b2trip
