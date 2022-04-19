%% forward
% Convenience function to compute the result of the forward operator 
% (contrast-to-measurement operator) on the contrast |q|, i.e. 
% the scattered field evaluated at the receivers' positions.
%
%
%% Syntax
%
%   [FFqMeas,FFqROI,seti] = forward()
%   [FFqMeas,FFqROI,seti] = forward(seti,qROI)
%
%
%% Description
% |[FFqMeas,FFqROI,seti] = forward()| and 
% |[FFqMeas,FFqROI,seti] = forward(seti,qROI)| 
% compute the result of the forward operator. 
%
% If there are no input arguments, some default values are used.
%
%% Example
%
% *Example 1*
%
%   [FFqMeas,FFqROI,seti] = forward();
%   figure(101); imagesc(real(seti.G(FFqROI(:,1)))); axis xy; axis square;
%   figure(102); imagesc(real(FFqMeas)); axis xy; axis square;
%
% *Example 2*
%
%   init;
%   seti = struct;
%   seti = setGeomSim(seti);
%   qROI = seti.qROIexact;
%   [FFqMeas,FFqROI,seti] = forward(seti,qROI);
%   figure(101); imagesc(real(seti.G(FFqROI(:,1)))); axis xy; axis square;
%   figure(102); imagesc(real(FFqMeas)); axis xy; axis square;
%
%% Input Arguments
%
% * |seti|  :    struct as described in <mimo.html>
% * |qROI|  :    contrast in ROI (region of interest)
%                (matrix as vector of size seti.nROI^seti.dim x 1)
%
%% Output Arguments
%
% * |FFqMeas| : result of forward operator, i.e. 
%               the scattered field evaluated on receivers' positions (measurements)
% * |FFqROI|  : scattered field evaluated on ROI (region of interest)
%               stored as vector for each transmitter
%               (complex matrix of size seti.nROI^seti.dim x seti.incNb)
% * |seti|    : struct with settings
%               (interesting if |forward()| is called without input arguments)
% 
%% Best Practice
%
% This is a convenience function. To avoid unnecessary calls, use
% <mimo.html>, e.g.
%
%   [FFqMeas,~,FFqROI] = mimo(seti,qROI,'simo')
%
%
%% More About
%
% The forward operator is
%
% $$\texttt{FFqMeas} = \mathcal{F}(q)$$
%
% Read also <mimo.html>.
%
%
%% See Also
%
% * <mimo.html>
%
%% Code
function [FFqMeas,FFqROI,seti] = forward(varargin)
%%
% Preparations for mimo
if nargin == 0
    init;
    if ~exist('seti','var')
        seti = struct;
    end
    seti = setGeomSim(seti);
    qROI = seti.qROIexact;
elseif nargin == 2
    seti = varargin{1}; % seti... see abvoe
    qROI = varargin{2}; % contrast q on ROI as a vector e.g. 8281 x 1
    if size(qROI,1) ~= seti.nROI^seti.dim || size(qROI,2) ~= 1
        error('Make sure that dimension of input qROI is seti.nROI^seti.dim x 1.')
    end
else
    error('Function "forward" needs 0 or 2 input arguments.')
end

%%
% Call of mimo
[FFqMeas,~,FFqROI] = mimo(seti, qROI, 'simo');

end
