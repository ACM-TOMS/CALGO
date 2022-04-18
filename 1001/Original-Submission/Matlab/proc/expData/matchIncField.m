%% matchIncField
% Match an incident field on region (ROI or CD) corresponding to the
% measured data of the incident field on measurements (receivers'
% positions).
%
%% Syntax
%
%   [uIncROI,errC] = matchIncField(uIncRX,seti,region)
%
%% Description
% |[uIncROI,errC] = matchIncField(uIncRX,seti,region)| finds suitable
% incident fields |uIncROI| on a region 
% (region of interest (ROI) or computational domain (CD) depending on
% |region|)
% corresponding to the measured data, |uIncRX|, of each incident field on measurements (receivers positions). 
% The relative error of every matched incident field at receivers positions
% is given too in |errC|.
% Further input parameters are in struct |seti|. 
% Note that this works only in case of a problem in 2 dimensions.
% 
%
%% Example
%
%   init; % Skip it if you are in folder proc/expData/, where matchIncField.m is located.
%
%   seti.dim = 2;   % dimension 2
%   % seti.nROI = 89; % (in case of region = 'ROI')
%   % seti.nCD = 250; % (in case of region = 'CD')
%
%   seti.incNb = 3;     % number of transmitters
%   % seti.measNb = 10;   % number of receivers
%   % seti.radSrc = 4;    % transmitters are on a circle with radius 4
%   % seti.radMeas = 5;   % receivers are on a circle with radius 5
%   seti.k = 250;       % wave number
%
%   % Positions of the transmitters
%   seti.incPnts =  [4.0000   -2.0000   -2.0000;
%                    0    3.4641   -3.4641];
%
%   % Positions of the receivers
%   seti.measPnts = [
%    5.0000    4.0451    1.5451   -1.5451   -4.0451   -5.0000   -4.0451    -1.5451    1.5451    4.0451;
%         0    2.9389    4.7553    4.7553    2.9389    0.0000   -2.9389    -4.7553   -4.7553   -2.9389];
%
%   % Generate the grid for ROI
%   h = 0.0016;
%   a = 0.0704;
%   x1ROI = -a:h:+a;
%   seti.nROI = length(x1ROI); % is 89 in this example as above
%   [X1ROI,X2ROI] = meshgrid(x1ROI,x1ROI);
%   seti.gridROI = [X1ROI(:).'; X2ROI(:).'];
%
%   seti.ampCalc = 1; % method to compute the coefficients
%   seti.nuMax = 3;   % 2*nuMax+1 coefficients are computed for polynomial approximation
%
%   % uIncRX : incident field at receivers positions...
%   %          (complex matrix of size seti.measNb x seti.incNb)
%
%   uIncRX = [
%   -0.0072 + 0.0028i  -0.0022 - 0.0009i  -0.0014 + 0.0012i;
%    0.0125 + 0.0063i   0.0023 + 0.0017i  -0.0017 + 0.0003i;
%    0.0059 + 0.0027i  -0.0002 - 0.0035i  -0.0009 + 0.0057i;
%    0.0002 - 0.0005i  -0.0140 - 0.0085i  -0.0043 - 0.0071i;
%   -0.0021 + 0.0011i  -0.0042 - 0.0011i  -0.0067 - 0.0037i;
%   -0.0008 - 0.0024i  -0.0003 - 0.0039i  -0.0030 + 0.0004i;
%   -0.0004 + 0.0021i   0.0019 - 0.0013i  -0.0034 - 0.0040i;
%    0.0028 + 0.0010i   0.0004 + 0.0013i  -0.0008 + 0.0054i;
%   -0.0020 + 0.0026i   0.0023 - 0.0003i   0.0023 - 0.0022i;
%   -0.0019 - 0.0012i   0.0005 + 0.0028i   0.0029 + 0.0014i];
% 
%   % Compute the incident field uIncROI on ROI
%   [uIncROI,errC] = matchIncField(uIncRX,seti,'ROI');
%
%   % Plots of incident fields from three transmitters
%
%   uIncROI1 = reshape(uIncROI(:,1),[seti.nROI seti.nROI]);
%   uIncROI2 = reshape(uIncROI(:,2),[seti.nROI seti.nROI]);
%   uIncROI3 = reshape(uIncROI(:,3),[seti.nROI seti.nROI]);
% 
%   figure(101); imagesc(real(uIncROI1)); axis xy; axis square;
%   figure(102); imagesc(real(uIncROI2)); axis xy; axis square;
%   figure(103); imagesc(real(uIncROI3)); axis xy; axis square;
% 
%
%% Input Arguments
% 
% * seti.dim        :   dimension of the problem: 2 (3 not available)
% * seti.incNb      :   number of transmitters
% * seti.k          :   wave number
% * seti.incPnts    :   Positions of the transmitters, 
%                       (matrix of size seti.dim x seti.incNb), 
%                       seti.incPnts = [5 -2 3; 0 4 2] 
%                       describes coordinates (5,0), (-2,4), and (3,2).
% * seti.measPnts   :   Positions of the receivers
%                       (matrix of size seti.dim x seti.measNb),
%                       coordinates analogical to seti.incPnts.
% * seti.gridROI    :   grid of region of interest (ROI) 
%                       (in case of region = 'ROI')
%                       (matrix of size seti.dim x seti.nROI^seti.dim)
% * seti.grid       :   grid of computational domain (CD)
%                       (in case of region = 'CD')
%                       (seti.dim x seti.nCD^seti.dim)
%
% * seti.ampCalc    :   method to compute the coefficients c (1, 2 or 3),
%                       we recommend to use method 1 (default)
%                       because it is the most accurate and fast, as shown in [1].
%                       (1: best-approximation to V c = z via MATLAB
%                       function,
%                       2: $c = (\gamma I + V^\ast\ V)^{-1} (V^\ast\
%                       \texttt{uIncRX})$
%                           via linear Tikhonov regularization,
%                       3: Landweber iteration)
% * seti.nuMax      :   match incident field using Hankel functions of first kind and orders
%                       $\nu = -\texttt{nuMax}, ..., -1, 0, 1, ...,
%                       \texttt{nuMax}$ (default: 7)
%
% * uIncRX          :   incident field at receivers' positions for each
%                       transmitter
%                       (complex matrix of size seti.measNb x seti.incNb)
% * region          :   'ROI' or 'CD' (region of interest or computational domain)
%                       (usually we are interested in ROI)
%
%
%% Output Arguments
%
% * uIncROI   :       incident field on ROI (or CD) for each transmitter 
%                     (complex array of size seti.nROI^seti.dim x seti.incNb)
% * errC      :       stores the relative error of 
%                     the matched incident field at receivers positions for
%                     every transmitter
%                     (vector of size 1 x seti.incNb)
%
%
%% More About
%
% *How does it work?* 
%
% # Matching incident field in data |uIncRX| (receivers positions) 
%   with harmonic polynomials for each transmitter by finding corresponding
%   coefficients |c| solving (approximately) $V c = \texttt{uIncRX}$.
% # Computes the resulting incident fields on the grid uIncROI 
%   for each transmitter.
%
% More details of the method are in the subfunctions and in Section 6 of [1].
%
% *Please note*
%
% * The method is for 2 dimensional problems.
% * We assume point sources at transmitters positions.
% * We assume near field data.
% * c are coefficients. (They are kind of amplitudes, but this is not exact.)
% * Avoid same positions for receivers and transmitters. Some entries in
% matrix V will be not a number (NaN).
%
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% Code: matchIncField
%
function [uIncROI,errC] = matchIncField(uIncRX,seti,region)
%%
% |[uIncROI,errC] = matchIncField(uIncRX,seti,region)|
% computes the incident field on grid of region ROI or CD (depends on
% |region|) from given incident field |uIncRX| at receivers postions.
% Relative error of given and matched indicent field |uIncRX| for each 
% transmitter are stored in |errC|.
%   
%   region =  'ROI' or 'CD'
%
nTX = seti.incNb; % nTX: number of transmitters
% nROI = nnz(seti.ROImask); % works for ROI, alternative for ROI and CD
switch region
    case 'ROI'
        nROI = seti.nROI^seti.dim;
    case 'CD'
        nROI = seti.nCD^seti.dim;
end
uIncROI = zeros(nROI,nTX);
errC = zeros(1,nTX);    % store relative error of V c = uIncRX, i.e. every matched incident field at receivers positions
errCAdj = zeros(1,nTX); % store relative error of (\gamma Id + V^\ast V) c = V^\ast uIncRX
% All transmitters are active in experiment (but only one at once)
for activeTX = 1:nTX
    % Matching for single incident field
    [u,errCval,errCAdjval] = matchIncFieldSingleTX(uIncRX,seti,activeTX,region);
    uIncROI(:,activeTX) = u(:);
    errC(1,activeTX) = errCval;
    errCAdj(1,activeTX) = errCAdjval;
end

disp('--')
disp('Check relative error of z = uIncRX because of choosen coefficients:')
fprintf('rel. err. of (gamma Id + V''V)c = V''z    : mean = %g\n',mean(errCAdj));
fprintf('             [min, max] = [%g, %g]\n', min(errCAdj), max(errCAdj));
fprintf('rel. err. of V c = z                    : mean = %g\n',mean(errC));
fprintf('             [min, max] = [%g, %g]\n', min(errC), max(errC));
disp('--')
end

%% Code: matchIncFieldSingle
function [uIncMatched,errC,errCAdj] = matchIncFieldSingleTX(uIncRX,seti,activeTX,region)

%%
% |[uIncMatched,errC,errCAdj] = matchIncFieldSingleTX(uIncRX,seti,activeTX,region)|
% matches the incident field |uIncRX| in data (receivers positions) with 
% harmonic polynomials and then computes the resulting field |uIncMatched| 
% on the grid for one transmitter.
%

% Active receivers depend on the active transmitter
% (not all receivers are active at once).
% Mark active receivers with 1 and inactive with 0: (logical):
activeRX = ~isnan(uIncRX(:,activeTX));
origin = seti.incPnts(:,activeTX); % position of transmitter with number activeTX
%
% seti.measPnts(1,activeRX) : first components of active measPnts
% (depends on choosen active transmitter: activeTX)
% x and y: coordinates away from origin (currently active transmitter)
% so in this coordinate system:
%
% * active transmitter position (0,0)
% * receiver positions (x,y)
%
x = transpose(seti.measPnts(1,activeRX) - origin(1)); % x coordinate
y = transpose(seti.measPnts(2,activeRX) - origin(2)); % y coordinate
z = uIncRX(activeRX,activeTX); % incident field at this receivers (because of active transmitter)
N = seti.nuMax; % degree of matching polynomial (e.g. N = 30)
phi0 = 0; % see function harmonicpolyfit

[c,errC,errCAdj] = harmonicpolyfit(x,y,z,N,seti.k,phi0,seti);

if 0 % store coefficients c (matchIncidentFieldSingleTX is called in for-loop)
    seti.ampCalcNew = 1; % Compute new coefficients (recommended).
                         % Could be stored, but not recommended...
                         % This works only in case of incNb == 1...
    % [pathstr,name,ext] = fileparts(file)
    [~,name,~] = fileparts(seti.fresnelFile);
    freqGHz = seti.fresnelFreq/1E9;
    % seti.ampCalc is not included in filename...
    % so compute the coefficients (and do not use the saved one...)
    fileDataAmp = sprintf('inmat/dataAmp/%s_%01d_GHz_N_%03d.mat',name,freqGHz,N);
    if (seti.ampCalcNew ~= 1 && exist(fileDataAmp, 'file') == 2)
        %file with name exists
        disp('Loading amplitudes of incident field of Fresnel data.');
        sloaded = load(fileDataAmp); % contains amplitudes c
        c = sloaded.c;
        clear sloaded;
    else
        disp('Compute coefficients c')
        [c,errC,errCAdj] = harmonicpolyfit(x,y,z,N,seti.k,phi0,seti);
        fprintf('Store in file: %s\n',fileDataAmp);
        save(fileDataAmp,'c');
    end
end

switch region
    case 'ROI'
        x = transpose(seti.gridROI(1,:) - origin(1));
        y = transpose(seti.gridROI(2,:) - origin(2));
    case 'CD'
        x = transpose(seti.grid(1,:) - origin(1));
        y = transpose(seti.grid(2,:) - origin(2));
end
uIncMatched = harmonicpolyval(c,x,y,seti.k,phi0); % matched field on grid
end

%% Code: Auxiliary Functions

%%
% *harmonicpolyfit*
%
function [c,errC,errCAdj] = harmonicpolyfit(x,y,uIncRX,N,k,phi0,seti)

%%
% |[c,errC,errCAdj] = harmonicpolyfit(x,y,uIncRX,N,k,phi0,seti)| finds the 
% coefficients of a 2D harmonic polynomial of (harmonic) degree |N| 
% that fits the data |uIncRX| at the best in a least-squares sense
% (in case of |seti.ampCalc = 1| -- other methods are available).
% The |c| is a row vector of length 2*N+1 containing the polynomial
% coefficients.
%
% Solution of Helmholtz equation in $\bf{R}^2$ without 0 is
%
% $$\texttt{hmonomial} = 
%     \frac{\mathrm{i}}{4} H_\nu^{(1)} (k r) \exp(\mathrm{i} \nu \varphi), 
%     \quad r > 0, \quad \varphi \in [0,2\pi),$$
%
% where r and phi are the polar coordinates of (x,y) and 
% $H_\nu^{(1)}$ the Hankel function of the first kind and 
% order $\nu \in \bf{Z}$.
%
% Any linear combination of these products are solutions too.
% This linear combination is expressed by matrix V multiplied by 
% vector c, which contains the coefficients.
%
% $$V c = \texttt{uIncRX}$$
%
% See Also HARMONICPOLYVAL.

[phi, r] = cart2pol(x(:),y(:));
V = zeros(length(r),2*N+1);
hmonomial = @(phi,r,n,k) (1i/4)*besselh(n, k*r).*exp(1i*n*(phi+phi0)); % phi+phi0

for n = -N:N
    V(:,n+N+1) = hmonomial(phi,r,n,k);
end

% seti.ampCalc: options 1: c = V\z | 2: using adjoint | 3: landweber
% seti.ampCalc = 1 is recommended
gamma = 1E-4; % used in case 2 and calculation of relative error
seti = checkfield(seti,'ampCalc',1);
switch seti.ampCalc
    case 1
        c = V\uIncRX(:); % 
    case 2
        c = (gamma*eye(size(V'*V)) + V'*V)\(V'*uIncRX(:));
    case 3
        lwN = 2000; % number of landweber iterations
        %lwN = 1E5;
        %lwN = 1E8; %test; very time-consuming
        c = landweber(V,uIncRX,lwN); % find c such that V c \approx z
end

%figure(101);
%plot(real(V*c)); title('V*c approx z = uIncRX, real');

%figure(102);
%plot(imag(V*c)); title('V*c approx z = uIncRX, imag');

% z is incident field uIncRX at receivers positions of the active(!) receivers
errCAdj = norm((gamma*eye(size(V'*V)) + V'*V)*c-V'*uIncRX(:))/norm(V'*uIncRX(:));
errC = norm(V*c-uIncRX(:))/norm(uIncRX(:));

if 0
    disp('Coefficients c = V\z correct? Is V*c-z = 0?')
    % (\varepsilon Id + V^\ast V) c = V^\ast z, so V c = z
    fprintf('rel. err. of (gamma Id + V''V)c = V''z: %g\n',errCAdj);
    fprintf('rel. err. of V c = z: %g\n',errC);
end
%figure(103); plot(real(c));

%figure(104);
%plot(real(c))

%figure(105);
%plot(imag(c))

end

%%
% *harmonicpolyval*
%
function uIncROI = harmonicpolyval(c,x,y,k,phi0)

%%
% |uIncROI = harmonicpolyval(c,x,y,k,phi0)| returns the value of the 2D harmonic 
% polynomial evaluated at points (x,y).
% The c is a vector of length $2\ \xi+1$ ($\xi = \texttt{nuMax}$) 
% whose elements are the coefficients 
% of the polynomial in ascending (!!!) powers.
%
% $$\texttt{uIncROI} = \sum_{\nu = -\xi}^{+\xi}
% c_\nu\ \texttt{hmonomial}(\varphi,r)$$
%
% where $\varphi$ and $r$ are the polar coordinates of (x,y).
%
% See also HARMONICPOLYFIT.
%

[phi, r] = cart2pol(x(:),y(:));
uIncROI = zeros(size(x(:)));
L = length(c);
N = floor( (L-1)/2 );
hmonomial = @(phi,r,n,k) (1i/4)*besselh(n, k*r).*exp(1i*n*(phi+phi0)); % phi+phi_0
for n = -N:(-N+L-1)
    uIncROI = uIncROI + c(n+N+1) * hmonomial(phi,r,n,k);
end;
uIncROI = reshape(uIncROI,size(x));
end

%%
% *landweber*
%
function res = landweber(V,uIncRX,lwN)

%%
% Try to improve fitting of coefficients c by Landweber sheme.
%
% Landweber sheme in general is used to minimize
%
% $\min_x \|Vx-y\|_2^2/2 \quad $ by $\quad x_{k+1} = x_{k} - \omega V^*(V x_k - y)$.
%
% Here we us it with 
%
% * y = uIncRX
% * $x_k = c$ (seeked coefficients c)
%
% Further input:
%
% * lwN   : number of iterations
%
m = size(V,2); % 2*N+1
n = size(uIncRX,2);
xk = zeros(m,n); % zeros(size(c)); % start with c = 0
% $0 < \omega < 2/\sigma_1^2$
s = svd(V);
size(V);
s1 = max(s);
omegaMax = 2/s1^2;
omega = 0.99*omegaMax;
A = V;
AT = transpose(conj(A));
y = uIncRX(:);
for i = 1:lwN
    xkp1 = xk-omega*AT*(A*xk-y);
    xk = xkp1;
    if floor(i/1000) == i/1000
        fprintf('dis = %g\n', norm(A*xk-y));
    end
end
res = xk;
% res is a c (or better c than input)
end
