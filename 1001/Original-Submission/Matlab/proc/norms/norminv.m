%% norminv
% Computes the norm of a vector on the down scaled region of interest.
%
%% Syntax
%
%   n = norminv(x,seti)
%
%% Description
%
% |n = norminv(x,seti)| does the same as |normroi| for down scaled region of
% interest, i.e. scaling factor |seti.dVinv| instead |seti.dVinv|.
%
%% See Also
%
% * <setGridScale.html>
% * <normroi.html>
%
%% Code
%
function n = norminv(x,seti)
% uses seti.qNorm and seti.dVinv

Q = seti.qNorm;
dV = seti.dVinv;

if Q == 2
    n = sqrt(abs(innerinv(x,x,seti)));
else
    n = norm(x,Q)*dV^(1/Q); % in case Q = 2 the result is the same like above
end

end
