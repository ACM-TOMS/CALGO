%% addNoise
% Adds artificial noise to the exact data |FmeasExact|.
%
%% Syntax
%
%   [seti, FmeasDelta] = addNoise(seti,FmeasExact)
%   [seti, FmeasDelta] = addNoise(seti,FmeasExact,dispDepth)
%
%% Description
%
% |[seti, FmeasDelta] = addNoise(seti, FmeasExact)| adds artificial noise
% of noise level |seti.delta| in dependence of the choosen noise type 
% |seti.whichNoise| to the exact data |FmeasExact|. 
% The simulated data with noise is saved in |FmeasDelta|.
%
% |[seti, FmeasDelta] = addNoise(seti, FmeasExact,dispDepth)| does the same
% but displayed messages can be controlled by |dispDepth|.
%
%% Example
%
%   init;
%   seti.measNb = 10;
%   seti.incNb = 5;
%   FmeasExact = rand(seti.measNb,seti.incNb) + 1i*rand(seti.measNb,seti.incNb);
%   seti.delta = 0.01;
%   seti.whichNoise = 'normal';
%   seti.dSMeas = 1; % here set to 1 to have a value, see setMeasPnts.html
%   [seti, FmeasDelta] = addNoise(seti, FmeasExact);
%
%% Input Arguments
%
% * seti        : structure array
% * FmeasExact  : exact scattered field evaluated on receivers' positions 
%                 (complex matrix of size seti.measNb x seti.incNb).
%
% _The following field of the structure array is required:_
%
% * seti.dSMeas         : Approximation of the infinitesimal element of closed 
%                         contour with control points, see <setMeasPnts.html>.
%
% _The following fields of the structure array are optional (otherwise
% default values are set):_
%
% * seti.delta          : Relative (artificial) noise level of data 
%                         (default: 0.01)
% * seti.whichNoise     : Type of used probability density function to 
%                         noise the data:
%                         'laplace', 'uniform', 'normal' (default).
% * seti.seed           : Number to control the random number generator.
%                         Using a non-negative integer (default: 0).
%
% *Optional Input Argument*
%
% * dispDepth           : Depth of displayed messages (0 or greater).
%
%% Output Arguments
%
% * seti    :   structure array
%
% If the fields |delta|, |whichNoise|, and |seed| of |seti| was not set,
% default values are set.
%
% * FmeasDelta  :   FmeasExact with noise (synthetic data)
%                   (complex matrix of size seti.measNb x seti.incNb).
%
%% More About
%
% * |FmeasDelta = FmeasExact + noise|.
% * Noise level is measured in the norm with seti.pNorm.
%
% *In case of |seti.whichNoise = 'normal'|*, see Section 5 in [1]:
%
% $F_\mathrm{meas}^\delta 
%  = F_\mathrm{meas} + \delta \, \frac{
%  \big\| F_\mathrm{meas} \big\|_\mathrm{dis}
%  }{
%  \| N_\mathrm{Re} + \mathrm{i} N_\mathrm{Im} \|_\mathrm{dis}
%  }
%  \big(N_\mathrm{Re} + \mathrm{i} N_\mathrm{Im} \big)$.
%
% Note that 
% $N_{\mathrm{Re}},N_{\mathrm{Im}} \in 
% \bf{R}^{\texttt{seti.measNb} \times \texttt{seti.incNb}}$ 
% are two real matrices sampled from standard, normal distribution.
%
% The considered relative noise level is 
% $\delta = \| F_\mathrm{meas}^\delta - F_\mathrm{meas}\|_\mathrm{dis} / 
% \| F_\mathrm{meas} \|_\mathrm{dis}$.
%
% *Remind*
%
% * If seti.pNorm = 2, then standard, normal distributed noise.
% * if seti.pNorm = 1, then Laplace distributed noise.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <mimo.html>
% * <setData.html>
% * <setMeasPnts.html>
%
%% Code

function [seti, FmeasDelta] = addNoise(seti,FmeasExact,varargin)

if nargin == 3
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

%% 
% *Set noise type and noise level*
%
seti = checkfield(seti,'delta',0.01,dispDepth);
seti = checkfield(seti,'whichNoise','normal',dispDepth); % whichNoise: laplace, normal, uniform
seti = checkfield(seti,'seed',0,dispDepth);

%%
rng(seti.seed); % standard is seed = 0

if strcmp(seti.whichNoise,'laplace')
    % Laplace-distr. noise
    [m,n] = size(FmeasExact);
    noise  = (setLaplDistNoise(m,n,0,1)+1i*setLaplDistNoise(m,n,0,1));
    noise  = noise/normLp(noise,seti)*normLp(FmeasExact,seti)*seti.delta;
    FmeasDelta = FmeasExact + noise;
elseif strcmp(seti.whichNoise,'uniform')
    % uniform noise
    noise  = (rand(size(FmeasExact))+1i*rand(size(FmeasExact)));
    noise  = noise/normLp(noise,seti)*normLp(FmeasExact,seti)*seti.delta;
    FmeasDelta = FmeasExact + noise;
elseif strcmp(seti.whichNoise,'normal')
    % standard, normal distributed noise
    noise  = (randn(size(FmeasExact))+1i*randn(size(FmeasExact)));
    noise  = noise/normws2(noise,seti)*normws2(FmeasExact,seti)*seti.delta;
    FmeasDelta = FmeasExact + noise;
else
    if dispDepth >= 1
        disp('seti.whichNoise not set correctly in addNoise.m - add normally distributed noise')
    end
    noise  = (randn(size(FmeasExact))+1i*randn(size(FmeasExact)));
    noise  = noise/normws2(noise,seti)*normws2(FmeasExact,seti)*seti.delta;
    FmeasDelta = FmeasExact + noise;
end

%%
% *How to compute the delta?*
%
% This can be used, when systematic error is added:
%
%   noise = FmeasDelta-FmeasExact;
%   deltaComp = normws2(noise,seti)/normws2(FmeasDelta,seti);
    
end

%% Code: subfunction: setLaplDistNoise
function y  = setLaplDistNoise(m, n, mu, sigma)
% Generate random numbers drawn from Laplacian distribution
% mean: mu and standard deviation = sigma (Default mu = 0, sigma = 1)
% [m, n]  : the dimension of y 
% Function goes back to laprnd.m by Elvis Chen (bee33@sjtu.edu.cn)

% Check inputs
if nargin < 2
    error('At least two inputs are required');
end

if nargin == 2
    mu = 0; sigma = 1;
end

if nargin == 3
    sigma = 1;
end

% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u).* log(1- 2* abs(u));
end
