%% setGridScale
% Sets functions seti.GD to downscale grid of ROI 
% and seti.GU to upscale grid again.
%
%% Syntax
%
%   seti = setGridScale(seti)
%   seti = setGridScale(seti,dispDepth)
%
%% Description
% |seti = setGridScale(seti)| 
% Sets functions seti.GD to downscale grid of ROI 
% and seti.GU to upscale grid again. The size of the down scaled ROI
% depends on the relation of the number of discretization points of the
% down scaled CD to the upscaled (full) CD.
%
% |seti = setGridScale(seti,dispDepth)| does the same but allows to control
% the depth of displayed messages by dispDepth.
% 
% * If seti.gscale = 1 grid scaling is considered.
%
%% Input Arguments
%
% If the fields was not defined, default values are set.
%
% * seti.gscale     :   1: use a down scaled grid for reconstruction;
%                       0: do not use grid scaling (default)
% * seti.nCDinv     :   number of discretization points
%                       _of the downscaled grid_
%                       for each dimension of CD.
%                       (for nCD see <setGrid.html>).
%                       (default: |floor(seti.nCD/2)|).
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
%
%% Output Arguments
%
% * seti.nInv   :   The analogue to seti.nROI, see <setGrid.html>, for the downscaled grid.
% * seti.hInv   :   The analogue to seti.h, see <setGrid.html>, for the downscaled grid.
% * seti.dVinv  :   The analogoue to seti.dV, see <setGrid.html>, for the downscaled grid.
% * seti.GInv   :   The analogue to seti.G, see <setReshapeVecMat.html>, for the downscaled grid.
%
% * seti.GU     :   function to scale up the grid
%                   (input: down scaled contrast as a vector of size seti.nInv^seti.dim;
%                   output: up scaled contrast as a vector of size seti.nROI^seti.dim)
% * seti.GD     :   function to scale down the grid
%                   (input: up scaled contrast as a vector of size seti.nROI^seti.dim;
%                   output: down scaled contrast as a vector of size seti.nInv^seti.dim)
%
%% More About
%
% *Process in the program*
%
% # seti.GInv or seti.G writes vector as matrix.
% # Scaling by gridUp or gridDown, see <gridUp.html>, <gridDown.html>.
% # seti.iG writes matrix as vector again.
%
% Grid scaling is used in <pda.html>.
%
% *The size of the down scaled grid...*
%
% * Note that ROI is down scaled and not the grid of CD.
% * Therefore the size of the down scaled grid is smaller than 
% seti.nCDinv^seti.dim.
% * So, |seti.nInv = round(seti.nROI/seti.nCD*seti.nCDinv)|.
%
%
%% See Also
%
% * <setGrid.html>
% * <setGeomSim.html>
% * <pda.html>
%
% * <gridUp.html>
% * <gridDown.html>
%
%% Code
function seti = setGridScale(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

%%
% *Check consistency*

seti = checkfield(seti,'gscale',0,dispDepth);

if (seti.gscale == 1) && ~isfield(seti,'nCDinv')
    seti.nCDinv = floor(seti.nCD/2);
    setmessage(seti,'nCDinv',dispDepth);
end

%%
% *Set hInv, dVinv, and nCDinv*
%
if seti.gscale == 1
    % volume of infinitestimal element of CD when nCDinv is used
    seti.hInv = 2*seti.rCD/seti.nCDinv;
    seti.dVinv = seti.hInv^seti.dim;

    seti.nInv = round(seti.nROI/seti.nCD*seti.nCDinv);
else
    seti.nCDinv = seti.nCD;
    seti.nInv = seti.nROI;
    seti.hInv = seti.h;
    seti.dVinv = seti.dV;
end

%%
% *Define seti.GU (scale grid up) and seti.GD (scale grid down)*
%
if seti.gscale == 1
    reshapeVecInv = seti.nInv*ones(1,seti.dim);
    seti.GInv  = @(x) reshape(x,reshapeVecInv); %  G: nInv x 1  -> matrix

    seti.GU = @(qROIdown) seti.iG(gridUp(seti.GInv(qROIdown),seti.nROI,seti.nInv)); % scale grid up
    seti.GD = @(qROI) seti.iG(gridDown(seti.G(qROI),seti.nROI,seti.nInv)); % scale grid down
else % id
    seti.GInv = seti.GROI;

    seti.GU = @(qROI) qROI;
    seti.GD = @(qROI) qROI;
end

end
