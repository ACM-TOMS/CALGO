function MP = MDB_patch(PP)

% Construction of an MDB-spline multi-patch from B-spline segments 
%
% INPUT
%   PP    : vector of B-spline patches
%
% OUTPUT
%   MP    : MDB-spline multi-patch

mu = [0 cumsum([PP.n])];
MP = struct('P', PP, 'mu', mu);

end
