%% setIncField
% Set incident fields on region of interest (ROI) for each transmitter.
%
%% Syntax
%
%   seti = setIncField(seti)
%
%% Description
% |seti = setIncField(seti)| evaluates the incident fields on region of
% interest (ROI) for each transmitter and stores them in struct seti in 
% field |incField|.
%
% 
%% Input Arguments
%
% * seti.incType    :   Type of incident field: 'pointSource' or 'planeWave'
% * seti.incNb      :   number of transmitters
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
% * seti.incPnts    : coordinates of transmitters
%                     (real matrix of size 2 x seti.incNb), see <setIncPnts.html>.
%
%
%% Output Arguments
%
% * |seti.incField|   :   incident field on region of interest (ROI)
%                         for each transmitter
%                         (complex matrix of size seti.nROI^seti.dim x seti.incNb)
%
%
%% More About
%
% * In case of experimental data (e.g. seti.expData = 'fresnel'), 
%   the incident fields seti.incField are computed in <loadData.html> 
%   from measurements at receivers positions.
%
% *Theory*
%
% For the incident fields see also Section 3.5 in [1]. Here, we give a 
% brief summary.
%
% *Theory: Plane waves*
%
% The incident field at position $x$ in case of a *plane wave* (in 2D and 3D) 
% in direction $\theta$ is
%
% $$u_\mathrm{plane}^\mathrm{i}(x) = \exp(\mathrm{i} k \langle x, \theta \rangle).$$
%
% The incident fields of $N_i$ plane waves with directions $\theta_j$, 
% $j = 1,...,N_i$, at points $x$ in ROI are stored as a complex matrix
%
% $$\texttt{seti.incField} = 
%   (\exp(\mathrm{i} k \langle x, \theta_j \rangle))_{x \in \mathrm{ROI},\
%   j = 1, ..., N_i}.$$
%
% The complex matrix has size
% seti.nROI^seti.dim x seti.incNb,
% where ROI is a vector instead of a matrix and 
% the number of incident fields is $N_i = \texttt{seti.incNb}$.
%
%
% *Theory: Point sources*
%
% The incident field at position $x$ in case of a *point source* at
% position $y$ is
%
% $$u_\mathrm{point}^\mathrm{i}(x) = \Phi(x-y)$$
%
% with radiating fundamental solution $\Phi(x)$ of the Helmholtz equation,
% points $x$ in ROI and transmitter position $y$.
%
% $$\Phi(x) = \frac{\mathrm{i}}{4} H_0^{(1)}(k |x|) \quad \mathrm{ if\ } x \in \bf{R}^2, \ x
% \neq 0,$$
%
% $$\Phi(x) = \frac{1}{4\pi} \frac{\exp(\mathrm{i} k|x|)}{|x|} \quad \mathrm{ if\ } x \in \bf{R}^3, \ x \neq 0,$$
%
% where $H_0^{(1)}$ is the Hankel function of first kind and order $0$.
%
% Analog to the case of a plane waves the incident fields (at
% points $x$ in ROI) of $N_i$ point sources at positions $y_j$, $j = 1, ...
% N_j$, are stored as a complex matrix
%
% $$\texttt{seti.incField} = (\Phi(x-y))_{x\in\mathrm{ROI},\ j = 1,...,N_i}.$$
%
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <expSetup.html>
% * <setKernel.html>
% * <setGrid.html>
% * <setIncPnts.html>
% * <start.html>
%
%
%% Code
%
% *function: setIncField*
%
function seti = setIncField(seti)

seti.incField = zeros(seti.nROI^2,seti.incNb);

if strcmp(seti.model,'helmholtz2D') || strcmp(seti.model,'helmholtzHMode2D')
    seti = setIncField2D(seti);
elseif strcmp(seti.model,'helmholtz3D')
    seti = setIncField3D(seti);
else
    disp('Unknown seti.model in setIncField.')
end

end

%%
% *subfunction: setIncField2D*
%
function seti = setIncField2D(seti)

if strcmp(seti.incType, 'pointSource')
    for jj = 1:seti.incNb
        x = seti.gridROI(1,:)-seti.incPnts(1,jj);
        y = seti.gridROI(2,:)-seti.incPnts(2,jj);
        [~,radii] = cart2pol(x,y);
        radii = reshape(radii, seti.nROI^2, 1);
        %phi = reshape(phi,seti.nROI^2,1);
        seti.incField(:,jj) = (1i/4)*besselh(0,1,seti.k.*radii);
    end

elseif strcmp(seti.incType, 'planeWave')
    for jj = 1:seti.incNb
        theta = kron(seti.incPnts(:,jj), ones(1,seti.nROI^2)); % size 2 x nROI^2
        seti.incField(:,jj) = exp(1i.*seti.k.*(dot(seti.gridROI, theta)));
    end
end

end

%%
% *subfunction: setIncField3D*
%

function seti = setIncField3D(seti)

seti.incField = zeros(seti.nROI^3,seti.incNb);
if strcmp(seti.incType, 'pointSource')
    for jj=1:seti.incNb
        radii = (seti.gridROI(1,:)-seti.incPnts(1,jj)).^2;
        radii = radii + (seti.gridROI(2,:)-seti.incPnts(2,jj)).^2;
        radii = radii + (seti.gridROI(3,:)-seti.incPnts(3,jj)).^2;
        radii = sqrt(radii);
        radii = reshape(radii, seti.nROI^3, 1);
        seti.incField(:,jj) = exp(1i.*seti.k.*radii)./(4.*pi.*radii);
    end
elseif strcmp(seti.incType, 'planeWave')
    for jj=1:seti.incNb
        theta = kron(seti.incPnts(:,jj), ones(1,seti.nROI^3)); % 3 x nROI^3
        seti.incField(:,jj) = exp(1i.*seti.k.*(dot(seti.gridROI, theta)));
    end
end

end
