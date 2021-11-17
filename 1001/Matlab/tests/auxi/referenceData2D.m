%% referenceData2D
% Computes scattering data in 2D via series solution for a reference
% contrast q, that is a ball.
%
%% Syntax
%
%   [uScattRXRef, uScattROIRef] = referenceData2D(seti)
%
%% Description
%
% |[uScattRXRef, uScattROIRef] = referenceData2D(seti)| computes in 2D 
% via series solution 
% the reference scattering data |uScattRXRef| at receivers positions and 
% the reference scattering data |uScattROIRef| at grid points in ROI
% for reference contrast q (a homogeneous ball in 2D).
%
% * Currently the reference solution of scatterung is only available for 
% incident plane wave with direction (1,0).
%
% Compute scattered field uScatt outside the ball:
%
% * If type = 'farField' compute far field in direction of receivers positions.
% * If type = 'nearField' compute near field at points RX (receivers positions).
%
%% Example
%
%   init;
%   seti.dim = 2;
%   seti.incType = 'planeWave';
%   seti.incNb = 1;
%   seti.rBall = 0.015;
%   seti.qBall = 0.8;
%   seti.measType = 'farField'; % 'nearField' is possible too...
%   seti = setGeomSim(seti); % set geometry and simulation
%   [uScattRXRef, uScattROIRef] = referenceData2D(seti);
%
%
%% Input Arguments
%
% * seti    : structural array
%
% Several fiels in seti are required, see <setGeomSim.html>.
%
% Important settings:
%
% * seti.dim = 2    : dimension of the problem is 2
% * seti.incType = 'planeWave' because only reference data for plane wave is implemented)
% * seti.incNb = 1 because reference data is ONE plane wave incidence with direction (1,0)) (i.e. seti.incPnts must be (1,0))
% * seti.rBall  :   radius of the ball (the scattering object)
% * seti.qBall  :   contrast value of the ball
%
%% Output Arguments
%
% * uScattRXRef     :   scattered field at receivers positions 
%                       (complex vector of size seti.measNb x 1).
% * uScattROIRef    :   scattered field on ROI written as vector 
%                       (complex vector of size seti.nROI^seti.dim x 1).
%
%% More About
%
% * This function is used in <testMimo.html>.
% * The reference solution of scattering is from an homogeneous ball in 2D.
% * The incident field is $\exp(\mathrm{i} k x_1)$ 
%   (the direction of the wave is hence (1,0)).
%
% For *far field* see "Bemerkung 2.23", p. 47 in [1]:
%
% $u_\infty(\hat{x}) = 4 \gamma_2 \sum_{n \in \bf{Z}} b_n
%  (-\mathrm{i})^{n+1} \exp(\mathrm{i}\, n\, \varphi)$
%
% Note that $\phi$ is $\varphi$ _ $\hat{x}$.
% 
% Note that $\gamma_2$ is missing in last line of "Bemerkung 2.23":
%
% $\gamma_2 = \exp(\mathrm{i} \pi/4)/\sqrt{8 k \pi}$
%
% (Definition of $\gamma_2$, see [1, p. 36]).
%
% For the connection between far field $u_\infty(\hat{x})$ and 
% radiating solution $u(x)$, see [1, Satz 2.16 on p. 36]:
%
% $u(x) = \exp(\mathrm{i} k |x|)/(|x|^{(d-1)/2}) (u_\infty(\hat{x} + \mathcal{O}(|x|^{-1})) \quad \mathrm{for\ } |x| \rightarrow \infty$
%
% with dimension $d = 2$.
%
%
%% References
%
% * [1] Armin Lechleiter: Script to _Zeitharmonische Wellen: Theorie und Anwendungen_.
%     University of Bremen, 2012.
%
%% See Also
%
% * <start.html>
% * <testMimo.html>
% * <referenceData3D.html>
%
%% Code
%
function [uScattRXRef, uScattROIRef] = referenceData2D(seti)

X1 = seti.grid(1,:);
X2 = seti.grid(2,:);

[uScattRXRef,uScattRef] = measurementReference2D(seti.qBall, seti.rBall, seti.k, X1, X2, seti.measPnts, seti.measType);

uScattROIRef = uScattRef(seti.ROImask);
uScattROIRef = uScattROIRef(:);

uScattRXRef = uScattRXRef(:);
end

function [uScattRXRef, uScattROIRef] = measurementReference2D(qBall, rBall, k, X1, X2, RX, type)

nMax = 20;
r = rBall;

nMat = -nMax:nMax;
a = zeros(size(nMat));  
b = zeros(size(nMat));

J = @besselj;
H = @besselh;
c = sqrt(1+qBall); %c is refractive index n because contrast q := n^2-1
kr = k*r;
for m = 1:length(nMat)
  % Solve linear system
  % A*a_n+B*b_n=i^n*U
  % C*a_n+D*b_n=i^n*V
  n = nMat(m);
  T = [J(n,c*kr)                      -H(n,kr);
       c*(J(n-1,c*kr)-J(n+1,c*kr))/2  -(H(n-1,kr)-H(n+1,kr))/2];
  tRight = 1i^n*[J(n,kr); (J(n-1,kr)-J(n+1,kr))/2];
  warning('off', 'MATLAB:nearlySingularMatrix');
  t = T\tRight;
  warning('on', 'MATLAB:nearlySingularMatrix');
  
  a(m) = t(1);
  b(m) = t(2);
end

%%
% Compute scattered field uScatt inside the ball
[phi, z] = cart2pol(X1,X2);
I = (z <= r);
uScattROIRef   = zeros(size(z)); 
for m = 1:length(nMat)
    n = nMat(m);
    uScattROIRef(I) = uScattROIRef(I) + (a(m)*J(n,k*c*z(I)) - 1i^n *J(n,k*z(I))).*exp(1i*n*phi(I));
end

%%
% Compute scattered field uScatt outside the ball
%
[phiRX, rRX] = cart2pol(RX(1,:),RX(2,:));
uScattRXRef = zeros(size(RX(1,:)));
I = ~I; % I inside ball, then ~I outside ball.
for m = 1:length(nMat)
    n = nMat(m);
    uScattROIRef(I) = uScattROIRef(I) + b(m)*H(n,k*z(I)).*exp(1i*n*phi(I));
    if strcmp(type,'nearField')
        uScattRXRef = uScattRXRef + b(m)*H(n,k*rRX).*exp(1i*n*phiRX);
    elseif strcmp(type,'farField')
        gamma2 = exp(1i*pi/4)/sqrt(8*k*pi);
        uScattRXRef = uScattRXRef + 4*gamma2*b(m).*(-1i)^(n+1).*exp(1i*n*phiRX);
    end
end

end
