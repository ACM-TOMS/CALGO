%% normws
% Computes the weighted Schatten p-norm (WS,P).
%
%% Syntax
%
%   normA = normws(A,seti)
%
%% Description
% |normA = normws(A,seti)| computes the weighted Schatten p-Norm.
% The weight is |seti.dSMeas| and p is |seti.pNorm|.
%
%% Input Arguments
%
% * A :   Complex matrix of size seti.measNb x seti.incNb.
% * seti.pNorm  :   Schatten p-Norm with p = seti.pNorm.
% * seti.dSMeas :   Approximation of the infinitesimal element of a closed contour 
%             with control points.
%             See also <expSetup.html>, <pntsGeometry.html>, and <dS2D.html>.
%
%% More About
%
% *Notation*
% 
% * WS: weighted Schatten norm (case pNorm ~= 2)
% * HS: WS,2 (case pNorm == 2): 
%       weighted Schatten 2-norm is weighted Hilbert-Schmidt norm (WHS)
% * Keep in mind: dSMeas is in HS-Norm:
%   $\|A\|_\mathrm{WHS} = \texttt{ sqrt( trace( A'*diag(seti.dSMeas)*B));}$
%
%
%% See Also
%
% * <innerhs.html>
%
%% Code

function normA = normws(A,seti)

if seti.pNorm == 2
    normA = sqrt(abs(innerhs(A,A,seti)));
    % Norm which is induced by innerhs ($\langle A,B \rangle_\mathrm{dis}$).
else
    S = svd(diag(seti.dSMeas^(1/seti.pNorm))*A);
    normA = norm(S,seti.pNorm);
end
    
end
