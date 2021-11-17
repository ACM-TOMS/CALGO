function s = real_to_pos(t)
%REAL_TO_POS Maps reals to positive numbers.
%
%  s = real_to_pos(t) maps the real number t to the positive number s.
%    Works as the inverse function of pos_to_real, hence
%    pos_to_real(real_to_pos(x)) is nearly equal to x.
%
% See also POS_TO_REAL.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012

  s = zeros(size(t));
  idx = (t>=0);
  s(idx) = t(idx)+1;
  s(~idx) = -1./(t(~idx)-1);
end
