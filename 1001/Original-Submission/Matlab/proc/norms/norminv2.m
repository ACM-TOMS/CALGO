%% norminv2
% Computes the norm of a vector on the down scaled region of interest as norminv but
% always with q = 2.
%
%% Syntax
%
%   n = norminv2(x,seti)
%
%% See Also
%
% * <norminv.html>

function n = norminv2(x,seti)
n = sqrt(abs(innerinv(x,x,seti)));
end
