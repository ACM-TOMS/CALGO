function [sol, val] = GSC(cc, ff, nn)
%GSC   Greedy submodular cover with submodular cost 
%   [sol, val] = GSC(cc, ff, nn);
%
%   This function implements the greedy algorithm for
%
%   min. cc(S) s.t. ff(S) = ff(true(nn, 1)); 
%
%   where S is a logical vector (called characteristic vector).
%
%   Input:
%       cc, ff: function handle of submodular function
%       nn: dimension of the logical vector. (Cardinarity of the ground set.)
%
%   Output:
%       sol: logical n-dimensional vector
%       val: value of the cost function cc(sol).
%
%   Reference:
%       P.-J. Wan, D.-Z. Du, P. Pardalos, W. Wu.
%       Greedy approximations for minimum submodular cover with submodular cost
%       Computational Optimization and Applications, 45(2):463--474, 2010

S = false(nn, 1);
fS = ff(S);

for ii = 1:nn
    deltamax = -Inf;
    for e = find(~S)'
        Stmp = S;
        Stmp(e) = true;
        fStmp = ff(Stmp);
        delta = (fStmp - fS) / cc(full(sparse(e, 1, true, nn, 1)));
        if delta > deltamax
            Smax = Stmp;
            fSmax = fStmp;
            deltamax = delta;
        end%if
    end%for e
    if deltamax > 0
        S = Smax;
        fS = fSmax;
    else
        break;
    end%if
end%for ii

sol = S;
val = cc(S);

end%function