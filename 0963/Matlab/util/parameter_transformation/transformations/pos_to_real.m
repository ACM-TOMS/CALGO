function t = pos_to_real(s)
%POS_TO_REAL Maps positive numbers to reals.
%
%  t = pos_to_real(s) maps the positive number s to the real number t.
%    Works as the inverse function of real_to_pos, hence
%    pos_to_real(real_to_pos(x)) is nearly equal to x.
%
% See also REAL_TO_POS.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012

  t = zeros(size(s));
  idx = (s>=1);
  t(idx) = s(idx)-1;
  t(~idx) = 1-1./s(~idx);
end
