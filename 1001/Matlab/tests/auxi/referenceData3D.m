%% referenceData3D
% Computes scattering data in 3D via series solution for a reference
% contrast q, that is a ball.
%
%% Syntax
%
%   [uScattRXRef, uScattROIRef] = referenceData3D(seti)
%
%% Description
% |[uScattRXRef, uScattROIRef] = referenceData3D(seti)| computes the same
% in 3D as <referenceData2D.html> in 2D.
%
% See <referenceData2D.html> for details. The differences in 3D are:
%
% * Of course dimension of the problem is seti.dim = 3.
% * Of course the incident plane wave has direction (1,0,0) in 3D.
% * In comparison to 2D (<referenceData2D.html>) no near field reference
% data is available in 3D, because it is not implemented yet.
%
%% Example
%
%   init;
%   seti.dim = 3;
%   seti.incType = 'planeWave';
%   seti.incNb = 1;
%   seti.rBall = 0.015;
%   seti.qBall = 0.8;
%   seti.measType = 'farField'; % Note that 'nearField' is not available
%   seti.nCD = 64; % a high number will take a long computational time
%   seti = setGeomSim(seti); % set geometry and simulation
%   [uScattRXRef, uScattROIRef] = referenceData3D(seti);
%
%% See Also
% * <referenceData2D.html>
%
%% Code
%
function [uScattRXRef, uScattROIRef] = referenceData3D(seti)

X = seti.grid(1,:);
Y = seti.grid(2,:);
Z = seti.grid(3,:);

[uScattRXRef,uScattRef] = measurementReference3D(seti.qBall, seti.rBall, seti.k, X, Y, Z, seti.measPnts, seti.measType);
uScattROIRef = uScattRef(seti.ROImask);

uScattROIRef = uScattROIRef(:);
uScattRXRef = uScattRXRef(:);
end

%%
% *Code: subfunction measurementReference3D*
%
% Computation of reference solution in 3D and evaluation on mesh X, Y, Z.
%
function [uRXScatt, uScatt] = measurementReference3D(q, R, k, X, Y, Z, RX, type)
% k: wave number
% n = q^2 inside B(0,R)
%
% MODEL: \Delta u + k^2 (1+q) u = 0
%

uRXScatt = zeros(size(RX(1,:)));

% normalize RX
[~, NX] = size(RX);
RXdirs = zeros(size(RX)).';

for j=1:NX
    RXdirs(j,:) = RX(:,j) / norm(RX(:,j));
end

% Direction of incoming plane wave
dir = [1 0 0];

nMax=10;

norms = sqrt(X(:).^2+Y(:).^2+Z(:).^2);
directions = [X(:)./norms Y(:)./norms Z(:)./norms];
inside = norms < R;
outside = norms >= R;

zeroIdx = find(norms==0);
uScatt = zeros(size(norms));

for n=0:nMax
    A = [-spherHankel(n,k*R)  spherBessel(n,k*sqrt(1+q)*R)
        -dSpherHankel(n,k*R)  sqrt(1+q)*dSpherBessel(n,k*sqrt(1+q)*R)];
    
    detA = A(1,1)*A(2,2)-A(1,2)*A(2,1);
    AInv = 1/detA * [A(2,2), -A(1,2); -A(2,1), A(1,1)];
    
    rhs = [spherBessel(n,k.*R); dSpherBessel(n,k.*R)];
    t = AInv * rhs;
    
    coeff = zeros(2*n+1,2);
    for m=-n:n
        coeff(m+n+1,:) = (4*pi*1i^n)*conj(spherHarmonic(n, m, dir))*t;
    end
    
    % evaluation outside 
    if ~all(outside==0)
        spherHar = zeros(length(norms(outside)), 1);
        
        for m=-n:n
            spherHar = spherHar + ...
                spherHarmonic(n,m,directions(outside,:))*coeff(m+n+1,1);
        end
        
        uScatt(outside) = uScatt(outside) + spherHankel(n,k*norms(outside)).*spherHar;
    end
    
    % evaluation inside (subtract incident wave)
    if ~all(inside==0)
        spherHarUi = zeros(length(norms(inside)), 1);
        spherHarInside = zeros(length(norms(inside)), 1);
        
        for m=-n:n
            tmp = spherHarmonic(n,m,directions(inside,:));
            
            spherHarInside = spherHarInside + tmp*coeff(m+n+1,2);
            spherHarUi = spherHarUi + tmp*4*pi*1.0i^n*conj(spherHarmonic(n,m,dir));
        end
        
        uScatt(inside) = uScatt(inside) + ...
            spherBessel(n,k*sqrt(1+q)*norms(inside)).*spherHarInside ...
            - spherBessel(n,k*norms(inside)).*spherHarUi;
    end
    
    if strcmp(type,'nearField')
        uRXScatt = 0;
        disp('uRXScatt was set to 0 in case of nearField in referenceData3D.m')
        % no near field implemented...
    elseif strcmp(type,'farField')
        for m=-n:n
            uRXScatt = uRXScatt + 1/k * 1/(1.0i^(n+1)) * spherHarmonic(n,m,RXdirs)' * coeff(m+n+1,1);
        end
    end
end

% not computable at 0, use average
uScatt(zeroIdx) = 0.5*(uScatt(zeroIdx+1)+uScatt(zeroIdx-1));
uScatt = reshape(uScatt, size(X));
end

%%
% *Code: subfunction: spherHarmonic: Spherical harmonics*
%
function y = spherHarmonic(n,m,x)
[phi,theta,~] = cart2sph(x(:,1),x(:,2),x(:,3));

leg = legendre(n,cos(theta+pi/2).');
legM = leg(abs(m)+1,:);

y = sqrt( ((2*n+1)/(4*pi))*(factorial(n-abs(m))/factorial(n+abs(m))) ) ...
    .* legM.' .* exp(1i*m*phi);
end

%%
% *Code: subfunction: spherHankel: 
% Spherical Hankel functions out of normal ones*
%
function u = spherHankel(n,z)
u = sqrt(pi./(2.*z)).*besselh(n+1/2,z);
end

%%
% *Code: subfunction: spherBessel: 
% Spherical Bessel functions out of normal ones*
%
function u = spherBessel(n,z)
u = sqrt(pi./(2.*z)).*besselj(n+1/2,z);
end

%%
% *Code: subfunction: dSpherHankel: Derivative of spherical Hankel functions*
%
function u = dSpherHankel(n,z)
u = -spherHankel(n+1,z) + (n./z).*spherHankel(n,z);
end

%%
% *Code: subfunction: dSpherBessel: Derivative of spherical Bessel functions*
%
function u = dSpherBessel(n,z)
u = -spherBessel(n+1,z) + (n./z).*spherBessel(n,z);
end
