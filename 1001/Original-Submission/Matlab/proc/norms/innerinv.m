%% innerinv
% Computes the inner product on down scaled ROI.
%
%% Syntax
%   i = innerinv(x,y,seti)
%
%% Description
% |innerinv(x,y,seti)| computes the inner product on the down scaled region of interest
% from two complex vectors x, y.
%
%% Input Arguments
%
% * x, y    :   two complex vectors of size seti.nInv^seti.dim
%               (values in region of interest written as a vector).
% * seti.dVinv :    area/volume of the infinitesimal element (pixel/voxel) of down scaled CD,
%                   <see setGrid.html> and especially <setGridScale.html>.
%
%% Output Arguments
%
% i     :   inner product on down scaled ROI, see "More About" in
% <innerroi.html>, but in this case with $h_N^d = \texttt{seti.dVinv}$.
%
%% See Also
%
% * <innerroi.html>
%
%% Code
function i = innerinv(x,y,seti)
i = (x'*y)*seti.dVinv;
end
