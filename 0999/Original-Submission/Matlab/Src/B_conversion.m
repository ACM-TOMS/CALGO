function ccd = B_conversion(Pd, Ps, ccs)

% Conversion from source to destination B-spline form
%
% INPUT
%   Pd    : destination B-spline patch
%   Ps    : source B-spline patch
%   ccs   : source coefficient vector
%
% OUTPUT
%   ccd   : destination coefficient vector

gg = B_greville(Pd);
Md = B_evaluation_all(Pd, gg);
sss = B_evaluation_spline(Ps, ccs, gg);
ccd = sss / Md;

end
