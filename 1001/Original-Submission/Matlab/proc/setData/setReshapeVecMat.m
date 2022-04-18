%% setReshapeVecMat
% Defines functions to reshape a vector into a matrix and vice versa.
%
%% Syntax
%
%   seti = setReshapeVecMat(seti)
%
%% Description
% |seti = setReshapeVecMat(seti)| defines functions seti.GROI and seti.GCD
% to reshape vectors into matrices of size seti.nROI
%
% * |seti.GROI| to reshape a vector into a matrix of size 
%
%% Example
%
%   init;
%   seti.dim = 2;
%   seti.nCD = 8;
%   seti.nROI = 3;
%   seti = setReshapeVecMat(seti);
%
%   v = rand(seti.nROI^seti.dim,1) % vector
%   M = seti.GROI(v) % write vector as matrix
%   v = seti.iG(M) % write matrix as vector again
%
%% Input Arguments
%
% * seti.dim    : see <setGrid.html>.
% * seti.nCD    : see <setGrid.html>.
% * seti.nROI   : see <setGrid.html>.
%
%% Output Arguments
%
% * seti.GROI   : reshapes a vector of size seti.nROI^seti.dim 
%                 into a matrix of size seti.nROI x seti.nROI in 2D
%                 and seti.nROI x seti.nROI x seti.nROI in 3D.
% * seti.GCD    : reshapes a vector of size seti.nCD^seti.dim 
%                 into a matrix of size seti.nCD x seti.nCD in 2D
%                 and seti.nCD x seti.nCD x seti.nCD in 3D.
% * seti.G      : same as seti.GROI.
% * seti.iG     : reshapes a matrix into a vector.
%
%% See Also
%
% * <setGeomSim.html>
% * <setGrid.html>
%
%% Code
%
function seti = setReshapeVecMat(seti)

reshapeVecROI = seti.nROI*ones(1,seti.dim);
seti.GROI = @(x) reshape(x,reshapeVecROI); %  G: n x 1  -> matrix

reshapeVecCD = seti.nCD*ones(1,seti.dim);
seti.GCD  = @(x) reshape(x,reshapeVecCD); %  G: n x 1  -> matrix

seti.G = @(q) seti.GROI(q);

seti.iG = @(x) x(:);                  % iG: matrix -> n x 1

end