%% normws2
% Computes norm as normws, but always the one, induced by innerhs
% (case pNorm == 2 independent from seti.pNorm input).
%
%% Syntax
%
%   normA = normws2(A,seti)
%
%% See Also
%
% * <normws.html>
%
%% Code
function normA = normws2(A,seti)
% Is normws but uses pNorm = 2.
% Uses seti.dSMeas.
% Norm which is induced by innerhs ($\langle A,B \rangle_\mathrm{dis}$)
normA = sqrt(abs(innerhs(A,A,seti)));
    
end
