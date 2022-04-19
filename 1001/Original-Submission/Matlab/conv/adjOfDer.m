%% adjOfDer
% Convenience function to compute the adjoint of the derivative
%
%% Syntax
%
%   [ADFFq,seti] = adjOfDer()
%   [ADFFq,seti] = adjOfDer(seti,qROI,FmeasDelta)
%
%% Description
% |[ADFFq,seti] = adjOfDer()| and |[ADFFq,seti] =
% adjOfDer(seti,qROI,FmeasDelta)| compute the adjoint of the derivative. 
%
% If there are no input arguments, some default values are used.
%
%% Example
%   [ADFFq,seti] = adjOfDer();
%   figure(101); imagesc(real(seti.G(ADFFq))); axis xy; axis square;
%
%% Input Arguments
%
% * seti       :    struct as described in <mimo.html>
% * qROI       :    contrast in ROI (region of interest)
%                   (matrix as vector of size seti.nROI^seti.dim x 1)
% * FmeasDelta :    $\mathcal{F}_\mathrm{meas}^\delta$ measurement data 
%                   (at receivers' positions) with noise level $\delta$
%                   (matrix of size seti.measNb x seti.incNb)
%
%% Output Arguments
%
% * ADFFq : adjoint of derivative
% * seti  : struct with settings
%           (interesting if |adjOfDer()| is called without input arguments)
%
%% Best Practice
%
% This is a convenience function. To avoid unnecessary calls, use
% <mimo.html>.
%
% If you need the scattered field evaluated on receivers positions too, 
% i.e. FFqMeas = uScattRX, then use mimo
%
%   [FFqMeas,ADFFq] = mimo(seti, qROI,seti.FmeasDelta,'adjOfDer');
%
% or (same result)
%
%   [FFqMeas,ADFFq] = mimo(seti, qROI,seti.FmeasDelta);
%
% If you only need uScattRX, use mimo with option 'simo', because it is
% faster.
%
%
%% More About
%
% The adjoint of derivative is
%
% $$\texttt{ADFFq} = [\mathcal{F}'(q)]^\ast [\mathcal{F}(q) -
% F_\mathrm{meas}^\delta]$$
%
% Read also <mimo.html>.
%
%
%% See Also
%
% * <mimo.html>
%
%% Code
function [ADFFq,seti] = adjOfDer(varargin)
%%
% Preparations for mimo

if nargin == 0
    init;
    if ~exist('seti','var')
        seti = struct;
    end
    seti = setGeomSim(seti);
    qROI = seti.qROIexact;
    FmeasDelta = zeros(seti.measNb,seti.incNb); % 0
elseif nargin == 3
    seti = varargin{1}; % seti... see abvoe

    qROI = varargin{2}; % contrast q on ROI as a vector e.g. 8281 x 1
    if size(qROI,1) ~= seti.nROI^seti.dim || size(qROI,2) ~= 1
        error('Make sure that dimension of input qROI is seti.nROI^seti.dim x 1.')
    end
  
    FmeasDelta = varargin{3}; % see formula above...
    if size(FmeasDelta,1) ~= seti.measNb || size(FmeasDelta,2) ~= seti.incNb
        error('Dimension of matrix FmeasDelta must be seti.measNb x seti.incNb.')
    end
else
    error('Function "adjOfDer" needs 0 or 3 input arguments.')
end

%%
% Call of mimo

[~,ADFFq] = mimo(seti,qROI,FmeasDelta,'adjOfDer');
% dimension of ADFFq: seti.nROI^seti.dim x 1 (i.e. qROI as a vector)
end