%% setMeasKer
%
% Set up measurement kernels
% (such that measurement = k^2 * kernel * solution * voxelVolume).
%
%% Syntax
%   seti = setMeasKer(seti)
%
%% Description
% |seti = setMeasKer(seti)| sets up the measurement kernel
% for each receiver and stores them in struct seti in field |measKer|.
%
%
%% Input Arguments
%
% * seti.measType    :   Type of scattered field: 'nearField' or 'farField'
% * seti.measNb      :   number of receivers
%
% * seti.k      :   wave number, see <setKernel.html>
% * seti.model  :   Model of problem 'helmholtz2D' or 'helmholtz3D',
%                   see <setKernel.html>
%
% * seti.nROI    : discretization points for each dimension
%                  of region of interest (ROI) (in samples),
%                  see <setGrid.html>.
% * seti.gridROI : grid of region of interest (ROI) 
%                  (matrix of size seti.dim x seti.nROI^seti.dim)
%                  see <setGrid.html>.
%
% * seti.measPnts    : coordinates of receivers
%                      (real matrix of size seti.dim x seti.measNb), see <setMeasPnts.html>.
%
%
%% Output Arguments
%
% * |seti.measKer|      :   measurement kernel on region of interest (ROI)
%                           for each receiver
%                           (complex matrix of size seti.measNb x seti.nROI^seti.dim)
%
%% More About
%
% For the measurement kernel see also Section 3.5 in [1]. Here, we give a 
% brief summary.
%
% * Previously definition: 
%   The radiating *fundamental solution* $\Phi(x)$ of the Helmholtz equation is
%
% $$\Phi(x) = \frac{\mathrm{i}}{4} H_0^{(1)}(k |x|) \quad \mathrm{ if\ } x \in \bf{R}^2, \ x
% \neq 0,$$
%
% $$\Phi(x) = \frac{1}{4\pi} \frac{\exp(\mathrm{i} k|x|)}{|x|} \quad \mathrm{ if\ } x \in \bf{R}^3, \ x \neq 0,$$
%
% where $k$ is the wave number and 
% $H_0^{(1)}$ is the Hankel function of first kind and order $0$.
%
%
% * *Near field* (2D and 3D)
%
% The *measurement kernel* in case of *near field* data is
%
% $\texttt{seti.measKer} = 
%   (\Phi(x-y_j))_{x\in \mathrm{ROI},\ j = 1, ..., N_S}$,
%
% where $x$ are points in ROI,
% $y_j$ is the position of the j-th receiver.
%
%
% * *Far field* (2D and 3D)
%
% The *measurement kernel* in case of *far field* data is
%
% $\texttt{seti.measKer} 
%   = \gamma \exp(-\mathrm{i}k\langle x,\theta_j\rangle)_
%     {x\in \mathrm{ROI},\ j = 1, ..., N_S}$,
%
% where $x$ are points in ROI,
% $\Theta_j$ is the direction of the j-th receiver and
% $\gamma = \exp(\mathrm{i} \pi/4) / \sqrt{8 \pi k} 
%  \quad \mathrm{ if\ } x \in \bf{R}^2$,
% $\quad$
% $\gamma = 1/(4\pi)
%  \quad \mathrm{ if\ } x \in \bf{R}^3$.
%
%
% * *Connection between measurement of receivers and measurement kernel*:
%
% *measurement = uScattRX = FmeasDelta =
% seti.k^2*seti.measKer*qROI.*(uIncROI+uScattROI)*seti.dV*
%
% This connection is used in <simo.html>.
%
% With _measurement_ we mean the scattered field at the receivers
% positions (in dependence of the active transmitter).
%
% k             :   wave number, 
% seti.measKer  :   measurement kernel, 
% qROI          :   contrast, 
% seti.dV       :   voxel volume.
%
% In <helmholtz2Dr2data.html> is computed:
% measData = (seti.measKer*fROI)*seti.dV;
%
% In <simo.html> is computed:
% uScattRX = seti.k^2.*helmholtz2Dr2data(fROI, seti);
%
% Note that |fROI| is here |qROI.*(uIncROI+uScattROI)| via function handle
% |QU| defined in <intOpsFuncs.html>.
%
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <expSetup.html>
% * <setKernel.html>
% * <setGrid.html>
% * <setMeasPnts.html>
% * <start.html>
%
%
%% Code: function setMeasKer
%
function seti = setMeasKer(seti)
%
% For future models:
%
% elseif strcmp(seti.model,'helmholtzAniso2D')
% elseif strcmp(seti.model,'helmholtzAniso3D')
% elseif strcmp(seti.model,'maxwell')
%
if strcmp(seti.model,'helmholtz2D') || strcmp(seti.model,'helmholtzHMode2D')
    seti = setMeasKer2D(seti);
elseif strcmp(seti.model,'helmholtz3D')
    seti = setMeasKer3D(seti);
else
    disp('Unknown seti.model in setMeasKer.')
end

end

%% Code: subfunction: setMeasKer2D
% Set measurement kernel in 2D.
%
function seti = setMeasKer2D(seti)

seti.measKer = zeros(seti.measNb,seti.nROI^2);

if strcmp(seti.measType, 'nearField')
    for jj = 1:seti.measNb
        radii = (seti.gridROI(1,:)-seti.measPnts(1,jj)).^2 ...
                + (seti.gridROI(2,:)-seti.measPnts(2,jj)).^2;
        radii = sqrt(radii);
        radii = reshape(radii, 1, seti.nROI^2);
        seti.measKer(jj,:) = (1i/4)*besselh(0,1,seti.k.*radii); % fundamental solution in 2D
    end
    
elseif strcmp(seti.measType, 'farField')
    for jj = 1:seti.measNb
        theta = kron(seti.measPnts(:,jj), ones(1,seti.nROI^2)); % size 2 x nROI^2
        seti.measKer(jj,:) = (exp(1i*pi/4)./sqrt(8*pi*seti.k)) ...
            .* exp(-1i.*seti.k.*(dot(seti.gridROI, theta))).'; % measKer farField 2D
    end
end

end

%% Code: subfunction: setMeasKer3D
% Set measurement kernel in 3D.
%
function seti = setMeasKer3D(seti)
seti.measKer = zeros(seti.measNb,seti.nROI^3);

if strcmp(seti.measType, 'nearField')
    for jj = 1:seti.measNb
        radii = (seti.gridROI(1,:)-seti.measPnts(1,jj)).^2;
        radii = radii + (seti.gridROI(2,:)-seti.measPnts(2,jj)).^2;
        radii = radii + (seti.gridROI(3,:)-seti.measPnts(3,jj)).^2;
        radii = sqrt(radii);
        radii = reshape(radii, 1, seti.nROI^3);
        seti.measKer(jj,:) = exp(1i.*seti.k.*radii)./(4.*pi.*radii); % fundamental solution in 3D

        % Workaround for Inf-Values
        seti.measKer(isinf(seti.measKer)) = 0;
    end
    
elseif strcmp(seti.measType, 'farField')
    for jj = 1:seti.measNb
        theta = kron(seti.measPnts(:,jj), ones(1,seti.nROI^3)); % size 3 x nROI^3
        seti.measKer(jj,:) = (1./(4.*pi)) .* exp(-1i.*seti.k.*(dot(seti.gridROI, theta))).'; % measKer farField 3D
    end
end

end
