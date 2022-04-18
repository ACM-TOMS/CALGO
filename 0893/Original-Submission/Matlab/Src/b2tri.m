function [ys,yp] = b2tri(x,y,w,p,d,sd)
% b2tri:  SPD block tridiagonal system solver
%
% USAGE:  [ys,yp] = b2tri(x,y,w,p,d,sd);
%
%   This function solves the order 2N symmetric positive-
% definite block tridiagonal linear system associated with
% minimizing the quadratic functional Q(YS,YP) described in
% Function SMCRV.
%
% On input:
%
%       X,Y,W = Vectors of length N containing abscissae,
%               data values, and positive weights, respect-
%               ively.  The abscissae must be strictly in-
%               creasing.
%
%       P = Positive smoothing parameter defining Q.
%
%       D,SD = Vectors of length N-1 containing positive ma-
%              trix entries.  Letting DX and SIG denote the
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
%               respectively, at the abscissae.
%
% Modules required by B2TRI:  None
%
%***********************************************************

m = size(x);
n = length(x);
nm1 = n - 1;

% Initialize output arrays.

ys = zeros(m);
yp = zeros(m);

% Work space:

t11 = zeros(nm1,1);
t12 = zeros(nm1,1);
t21 = zeros(nm1,1);
t22 = zeros(nm1,1);

% The forward elimination step consists of scaling a row by
%   the inverse of its diagonal block and eliminating the
%   subdiagonal block.  The superdiagonal is stored in T and
%   the right hand side in YS,YP.  For J = 11, 12, and 22,
%   SJI and SJIM1 denote the elements in position J of the
%   superdiagonal block in rows I and I-1, respectively.
%   Similarly, DJI denotes an element in the diagonal block
%   of row I.

% Initialize for I = 2.

dx = x(2) - x(1);
dim1 = d(1);
s22im1 = sd(1);
s12im1 = (dim1 + s22im1)/dx;
s11im1 = -2.0*s12im1/dx;
r1 = p*w(1);
d11i = r1 - s11im1;
d12i = s12im1;
d22i = dim1;
den = d11i*d22i - d12i*d12i;
t11(1) = (d22i*s11im1 + d12i*s12im1)/den;
t12(1) = (d22i*s12im1 - d12i*s22im1)/den;
t21(1) = -(d12i*s11im1 + d11i*s12im1)/den;
t22(1) = (d11i*s22im1 - d12i*s12im1)/den;
r1 = r1*y(1)/den;
ys(1) = d22i*r1;
yp(1) = -d12i*r1;

% I = 2 to N-1:

for i = 2:nm1
   im1 = i - 1;
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
   r1 = r1*y(i) - s11im1*ys(im1) + s12im1*yp(im1);
   r2 = -s12im1*ys(im1) - s22im1*yp(im1);
   ys(i) = (d22i*r1 - d12i*r2)/den;
   yp(i) = (d11i*r2 - d12i*r1)/den;
   dim1 = di;
   s22im1 = s22i;
   s12im1 = s12i;
   s11im1 = s11i;
end

% I = N:

r1 = p*w(n);
d11i = r1 - s11im1 - (s11im1*t11(nm1)-s12im1*t21(nm1));
d12i = -s12im1 - (s11im1*t12(nm1) - s12im1*t22(nm1));
d22i = dim1 - (s12im1*t12(nm1) + s22im1*t22(nm1));
den = d11i*d22i - d12i*d12i;
r1 = r1*y(n) - s11im1*ys(nm1) + s12im1*yp(nm1);
r2 = -s12im1*ys(nm1) - s22im1*yp(nm1);
ys(n) = (d22i*r1 - d12i*r2)/den;
yp(n) = (d11i*r2 - d12i*r1)/den;

% Back solve the system.

for i = nm1:-1:1
   ys(i) = ys(i) - (t11(i)*ys(i+1) + t12(i)*yp(i+1));
   yp(i) = yp(i) - (t21(i)*ys(i+1) + t22(i)*yp(i+1));
end
return;

end  % b2tri
