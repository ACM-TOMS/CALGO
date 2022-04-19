%% Single dielectric cylinder (obstacle to Fresnel in 2D: dielTM_dec8f.exp)
% For reference see <incontrastsRef.html>.

%% Position of obstacle
%
% * Theoretical, see [1], position should be left from origin, but Fresnel 
%   data are rotated 270 degrees mathematical positive.
% * Note: rotation seems to be 90 degrees, but take care of y-axis, because
%   it is top negative, when you plot with imagesc
%
%   imagesc(...);
%
% * Solution: Corrected when plotting with:
%
%   axis xy;
%
%% Other parameters of the obstacle
%
% Fresnel, op. 1: Section 4.2 dielectric targets, see [1]:
%
% * Radius a = 15 mm
% * Real part of relative permittivity $\varepsilon_r = 3 \pm 0.3$.
% * Relative permeability $\mu_r$ not given, so guess $\mu_r = 1$.
% * Cylinders are "dielectric", this means imaginary part of contrast is 0.
%
%% References
%
% * [1] Kamal Belkebir and Marc Saillard. Special section on testing inversion algorithms against experimental data. _Inverse Problems_, 17(6):1565-1571, 2001.

function q = fresnel_op1_dielTM(X1,X2,varargin)

finecorrection = 1; % manually correction of obstacles position

qBall = 2 + 1i*0; % dielectricum, so imaginary part is 0 (in theory)

rBall = 15E-3; % 15 mm
dist = 30E-3; % 30 mm distance from origin
r = hypot(X1,X2-dist); % 30 mm top of origin

if finecorrection == 1
    xshift = 1E-3;
    yshift = -3E-3;
    r = hypot(X1-xshift,X2-dist-yshift);
end

q = qBall*( r < rBall );
q = double(q);

end
