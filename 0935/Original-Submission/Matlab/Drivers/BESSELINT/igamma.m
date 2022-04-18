function fn_val = igamma(alpha, x)
%IGAMMA Incomplete Gamma function.
%   Y = IGAMMA(A,X) evaluates the incomplete gamma function for
%   arbitrary complex A and X. Both A and X should be numbers, not
%   vectors or matrices.
%
%   This implementation is a conversion of the Fortran 77 program
%   cdig by Eric Kostlan & Dmitry Gokhman, described in
%   
%   "KOSTLAN, E., and GOKHMAN, D. A program for calculating the
%   incomplete gamma function. Tech. rep., Dept. of Mathematics,
%   Univ. of California, Berkeley, 1987."
%
%   To speed up performance, the function DNRM in the original Fortran
%   program has been replaced by the built-in Matlab function ABS.

% ----------------------------------------------------------------------
%   Authors: Joris Van Deun & Ronald Cools,
%            Dept. of Computer Science, K.U.Leuven, Belgium.
%
%   Software revision date: September 2, 2005
% ----------------------------------------------------------------------


re = 0.36787944117144232; % exp(-1)
ibuf = 34;

% --- If x is near the negative real axis, then shift to x=1.
if (abs(x) < 1) | (real(x) < 0 & abs(imag(x)) < 1),
  p = 0; q = 0;
  fn_val = re / cdh(alpha, 1);
  ilim = real(x/re);
  for  i = 0 : ibuf - ilim,
    [p,q] = termpq(alpha, x, i, p, q);
    fn_val = fn_val + p * q;
  end
else
  fn_val = exp(-x + alpha*log(x)) / cdh(alpha, x);
end


function fn_val = cdh(alpha, x)
% --- Written By Eric Kostlan & Dmitry Gokhman
% --- March  1986

% --- If Re(alpha-x) is too big, shift alpha.
n = fix(real(alpha-x));
if (n > 0),
  cn = n;
  alpha1 = alpha - cn;
  term = 1 / x;
  sum = term;
  for  i = 1:n - 1,
    cn = n - i;
    term = term * (alpha1 + cn) / x;
    sum = term + sum;
  end
  sum = sum + term * alpha1 / cdhs(alpha1, x);
  fn_val = 1 / sum;
else
  fn_val = cdhs(alpha, x);
end


function fn_val = cdhs(alpha, x)
% --- Written By Eric Kostlan & Dmitry Gokhman
% --- March  1986
% error has been changed from 5e-18 to 10*eps in Matlab-version
tol1 = 1e10; tol2 = 1e-10; error = 10*eps;
ilim = 100000;

q0 = 1;
q1 = 1;
p0 = x;
p1 = x + 1 - alpha;
for i = 1:ilim,
  ci = i;
  if (p0 ~= 0 & q0 ~= 0 & q1 ~= 0),
    r0 = p0 / q0;
    r1 = p1 / q1;
    if (abs(r0-r1) <= abs(r1)*error),
      fn_val = r1;
      return
    end;
% --------- Occasionally renormalize the sequences to avoid over(under)flow.
    if (abs(p0) > tol1 | abs(p0) < tol2 | abs(q0) > tol1  ...
          | abs(q0) < tol2),
      factor = p0 * q0;
      p0 = p0 / factor;
      q0 = q0 / factor;
      p1 = p1 / factor;
      q1 = q1 / factor;
    end
  end
  p0 = x * p1 + ci * p0;
  q0 = x * q1 + ci * q0;
  p1 = p0 + (ci+1-alpha) * p1;
  q1 = q0 + (ci+1-alpha) * q1;
end
% --- If the peripheral routines are written correctly,
% --- the following statements should never be executed.
warning(sprintf('cdhs:  i > %d \n r0,r1= %d,%d', ilim, r0, r1));
fn_val = 0.5 * (r0+r1);


function [p,q] = termpq(alpha, x, i, p, q)
% --- Calculate p*q = -1**i(1 - x**(alpha+i))/(alpha+i)i carefully.

tol = 3e-7; xlim = 39;

if (i == 0), q = 1; end
ci = i;
alphai = alpha + ci;
if (x == 0),
  p = 1 / alphai;
  if (i ~= 0), q = -q / ci; end
  return
end
cdlx = log(x);

%----------------------------------------------------------------------
% We (J+R) believe the following part is wrong in the original version.
%   Example: igamma(-20,0.1) gives
%    4.500507092028393720e+18 (Maple)
%    4.500481621643106e+18 (original igamma)
%    4.500507092028388e+18 (our version with code commented out)
% Note that this is irrelevant for application with besselint,
% because this case will never appear in a call.
% However, we don't know who else will use this file...
%
% 
% --- If (1 - x**alphai) = -x**alphai on the computer,
% --- then change the inductive scheme to avoid overflow.
% if (real(alphai*cdlx) > xlim & i ~= 0),
%   p = p * (alphai - 1) / alphai;
%   q = -q * x / ci;
%   return
% end
%----------------------------------------------------------------------
if (abs(alphai) > tol),
  p = (1 - x^alphai) / alphai;
else
  p = -cdlx * (1 + cdlx*alphai/2);
end
if (i ~= 0), q = -q / ci; end
