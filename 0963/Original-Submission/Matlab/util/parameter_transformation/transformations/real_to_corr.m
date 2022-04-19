function s = real_to_corr(t)
%REAL_TO_CORR Maps reals to correlations.
%
%  s = real_to_corr(t) maps the real number t to the correlation s in [-1,1].
%    Works as the inverse function of corr_to_real, hence
%    corr_to_real(real_to_corr(x)) is (nearly) equal to x.
%
% See also CORR_TO_REAL.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012

  s = zeros(size(t));
  idx = (t>=0);
  s(idx) = 1-1./(t(idx)+1);
  s(~idx) = -1-1./(t(~idx)-1);
end
