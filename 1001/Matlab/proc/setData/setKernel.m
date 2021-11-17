%% setKernel
% Computes the Fourier coefficients of the kernel.
%
%% Syntax
%
%   seti = setKernel(seti);
%   seti = setKernel(seti,dispDepth);
%
%% Description
% |seti = setKernel(seti)| computes the Fourier coefficients of the kernel
% and stores them in |seti.kHat|. 
% Note that the factor $k^2$ is _not_ included in Fourier coefficients.
%
% |seti = setKernel(seti,dispDepth)| does the same, but allows to control
% the depth of displayed messages.
%
%% Example
%
%   init;
%   
%   seti.dim = 2;
%   seti.nCD = 256;
%   seti.rCD = 0.2;
%   seti.k = 250;
%   seti.model = 'helmholtz2D';
%   
%   seti = setKernel(seti);
%
% The computed Fourier coefficients of the kernel are in |seti.kHat|.
%
%% Input Arguments
%
% _This arguments was set in rebis before setKernel is called_:
%
% * seti.dim    : see <setGrid.html>.
% * seti.nCD    : see <setGrid.html>.
% * seti.rCD    : see <setGrid.html>.
%
% _Can be set by user (otherwise default is set)_:
%
% * seti.k      : Wave number, default: 250
% * seti.model  : Model of the problem, default: 'helmholtz2D'
%
% *Optional Input Argument*
%
% * dispDepth   : Depth of displayed messages (0: no, 1 or greater: yes).
%
%% Output Arguments
%
% * seti.kHat   : kernel, complex matrix of size nCD x nCD.
%
%% More About
%
% For theory see [1], [2] or Section 3.2 in [3].
%
% Note that the factor $k^2$ is _not_ included in Fourier coefficients in
% this code, but the factor is included in [3].
%
%% References
%
% * [1] Gennadi Vainikko. Fast solvers of the Lippmann-Schwinger equation.
% In Robert P. Gilbert, Joji Kajiwara, and Yongzhi S. Xu, editors, 
% _Direct and Inverse Problems of Mathematical Physics_,
% volume 5 of _International Society for Analysis, Applications and Computation_, pages 423-440.
% Springer, New York, 2000.
% * [2] Thorsten Hohage. 
% On the numerical solution of a three-dimensional inverse medium scattering problem.
% _Inverse Problems_, 17(6):1743, 2001.
% * [3] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <setGeomSim.html>
%
%% Code
%
function seti = setKernel(seti,varargin)

if nargin == 2
    dispDepth = varargin{1};
else
    dispDepth = 0;
end

seti = kernelConsis(seti,dispDepth);

if strcmp(seti.model,'helmholtz2D') || strcmp(seti.model,'helmholtzHMode2D')
    seti.kHat = kernel2D(seti);
elseif strcmp(seti.model,'helmholtz3D')
    seti.kHat = kernel3D(seti);
else
    disp('setKernel.m: Error - Kernel coefficients of chosen model not implemented')
end

end

%% Code: subfunction: kernelConsis
% Checks consistency of kernel settings.

function seti = kernelConsis(seti,dispDepth)

seti = checkfield(seti,'k',250,dispDepth);

if ~isfield(seti,'model')
    if seti.dim == 2
        seti.model = 'helmholtz2D';
    else
        seti.model = 'helmholtz3D';
    end
    setmessage(seti,'model',dispDepth);
end

if strcmp(seti.model,'helmholtz')
    if seti.dim == 2
        seti.model = 'helmholtz2D';
    else
        seti.model = 'helmholtz3D';
    end
    if dispDepth >= 1
        fprintf('   Adapted Helmholtz model to dimension %d.\n', seti.dim);
    end
end

if strcmp(seti.model,'helmholtz2D') && (seti.dim == 3)
    seti.dim = 2;  
    if dispDepth >= 1
        disp('   Parameter "dim"=3 does not fit model "helmholtz2D". Set dim = 2.');
    end
end

if strcmp(seti.model,'helmholtz3D') && (seti.dim == 2)
    seti.dim = 3;  
    if dispDepth >= 1
        disp('   Parameter "dim"=2 does not fit model "helmholtz3D". Set dim = 3.');
    end
end

end

%% Code: subfunction: kernel2D
% Computes the Fourier coefficients of the kernel in 2D case.
%
function y = kernel2D(seti)
    
s = ceil((seti.nCD-1)/2);
j1 = (0:seti.nCD-1)-s;
[j1,j2] = meshgrid(j1,j1);
p = pi*sqrt(j1.^2+j2.^2); % p = pi*|j|

%%
% In this file we use notatio R = seti.rCD 
% (this is the notation in Vainikko, see [1]).
% In [3] we used the notation 2 R = seti.rCD.

kR = seti.k*seti.rCD;

J  = @besselj; %besselj: Bessel function of first kind
% besselh: Bessel function of third kind = Hankel function
H0 = 1i*pi/2*besselh(0,1,kR); % using Hankel function of first kind and order 1
H1 = 1i*pi/2*besselh(1,1,kR); % using Hankel function of first kind and order 0

%%
% The following two lines follow Vainikko (see [1]):
% * Note that we also deal with k != 1, while Vainikko presents only the
% case k = 1 on page 16.
% * Note the missing factor 1/(2R) in Fourier coefficients y = ... because
% Vainikko uses another basis function (with factor 1/(2R)).

y = seti.rCD^2./(p.^2-kR^2).*(1 + H0*p.*J(1,p ) - H1*kR*J(0,p )); % for p != kR
y(p == kR) = seti.rCD^2/2*(   H0*J(0,kR)    + H1*J(1,kR)); % can be proven by l'Hospital from the above formula

s = s*[1 1]; % Do not delete this line!
y = circshift(y,-s);

end

%% Code: subfunction: kernel3D
% Computes the Fourier coefficients of the kernel in 3D case.
%
function y = kernel3D(seti)

kR = seti.k*seti.rCD;

N2 = ceil((seti.nCD-1)/2);
[JX,JY, JZ] = meshgrid((0:seti.nCD-1)-N2);

% pi*|j|
PIJ = pi*sqrt(JX.^2 + JY.^2 + JZ.^2);

% The following computations of K follow [1] and [2], 
% but note the correction in case of PIJ == kR (as in Section 3.2 of [3]).

K = kR^2./(PIJ.^2-kR^2) .* (1 - exp(kR*1.0i) * ( cos(PIJ) - 1.0i*kR./PIJ .* sin(PIJ) ) );
%K(PIJ == kR) = -1.0i/sqrt(2) * kR *(1 - exp(kR*1.0i) / kR * sin(kR)); % wrong
K(PIJ == kR) = 1/2*exp(1i*kR)*(kR*sin(kR) + 1i*(kR*cos(kR) - sin(kR))); % corrected (!) (see Section 3.2 in [3])
K(PIJ == 0) = exp(kR*1.0i) * (1-kR*1.0i) - 1;

% Exclude factor k^2
K =  K/seti.k^2;

% Note that circshift(K,-N2*[1,1,1]) is the same as fftshift(K).
y = circshift(K,-N2*[1,1,1]);

end
