function Hl = MDB_extraction_local(MP, H, ip)

% Computation of the local extraction matrix on a patch
%
% INPUT
%   MP    : MDB-spline multi-patch
%   H     : extraction matrix
%   ip    : index of patch
%
% OUTPUT
%   Hl    : local extraction matrix

Hl = H(:,MP.mu(ip)+1:MP.mu(ip+1));

end
