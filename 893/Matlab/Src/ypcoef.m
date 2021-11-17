function [d,sd] = ypcoef(sigma,dx)
% ypcoef:  Coefficients for ypc2p and smcrv
%
% USAGE:  [d,sd] = ypcoef(sigma,dx);
%
%   This function computes the coefficients of the deriva-
% tives in the symmetric diagonally dominant tridiagonal
% system associated with the C-2 derivative estimation pro-
% cedure for a Hermite interpolatory tension spline.
%
% On input:
%
%       SIGMA = Nonnegative tension factor (or vector of
%               factors) associated with intervals.
%
%       DX = Positive interval widths in one-to-one corres-
%            pondence with SIGMA.
%
% On output:
%
%       D = Components of the diagonal terms associated with
%           the intervals.  D = SIGMA.*(SIGMA.*COSHM(SIGMA)-
%           SINHM(SIGMA))./(DX.*E), where E = SIGMA.*
%           SINH(SIGMA) - 2*COSHM(SIGMA).
%
%       SD = Subdiagonal (superdiagonal) terms.  SD = 
%            SIGMA.*SINHM(SIGMA)./(DX.*E).
%
% Module required by YPCOEF:  SNHCSH
%
%***********************************************************

k = find(sigma < 1.e-9);

% k = indices for which SIGMA = 0:  cubic interpolant.

d(k) = 4.0./dx(k);
sd(k) = 2.0./dx(k);

% For 0 < SIGMA <= .5, use approximations designed to avoid
%   cancellation error in the hyperbolic functions.

k = find(sigma >= 1.e-9  &  sigma <= 0.5);
[sinhm,coshm,coshmm] = snhcsh(sigma(k));
e = (sigma(k).*sinhm - coshmm - coshmm).*dx(k);
d(k) = sigma(k).*(sigma(k).*coshm-sinhm)./e;
sd(k) = sigma(k).*sinhm./e;

% For SIGMA > .5, scale SINHM and COSHM by 2*EXP(-SIGMA) in 
%   order to avoid overflow when SIGMA is large.

k = find(sigma > 0.5);
ems = exp(-sigma(k));
ssinh = 1.0 - ems.*ems;
ssm = ssinh - 2.0*sigma(k).*ems;
scm = (1.0-ems).*(1.0-ems);
e = (sigma(k).*ssinh - scm - scm).*dx(k);
d(k) = sigma(k).*(sigma(k).*scm-ssm)./e;
sd(k) = sigma(k).*ssm./e;
return;

end  % ypcoef
