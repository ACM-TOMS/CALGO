function t = corr_to_real(s)
%CORR_TO_REAL Maps correlations to reals.
%
%  t = corr_to_real(s) maps the correlation s in [-1,1] to the real number t.
%    Works as the inverse function of real_to_corr, hence
%    corr_to_real(real_to_corr(x)) is (nearly) equal to x.
%
% See also REAL_TO_CORR.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012

  t = zeros(size(s));
  idx = (s>=0);
  t(idx) = 1./(1-s(idx))-1;
  t(~idx) = 1+1./(-1-s(~idx));
end
