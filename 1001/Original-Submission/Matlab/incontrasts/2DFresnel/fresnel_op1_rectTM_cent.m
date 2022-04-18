%% Metallic rectangle (obstacle to Fresnel in 2D: rectTM_cent.exp)
% This file contains an assumption, i.e. is *experimentally*.
%
% For reference see <incontrastsRef.html>.
%
%% Geometry and position of the obstacle
%
% * Geometry of rectangle: 24.5 mm, 12.7 mm.
% * One metallic rectangle centered.
% * Mathematical positive rotation of 270 degrees, see
% <fresnel_op1_dielTM.html>.
% 
%% Other parameters of the obstacles
%
% * The metallic objects in [1] are perfectly conducting, i.e.
%   $\sigma = \infty$, so $\mathrm{imag}(\varepsilon_r) = \infty$
%
% * See, e.g. [2, p. 110] for 
% $\varepsilon_r(x) = \varepsilon_r(x)/\varepsilon_0 + i \sigma(x)/(\omega \varepsilon_0),$
% where $\sigma$ is the conductivity.
%
% * Refractive index n, $n^2 = \varepsilon_r \mu_r$.
%
% * Contrast $q:= n^2-1$ (in other literature you can find $q =
% \sqrt{n^2-1}$).
%
% * Finally, because the obstacle is metal, i.e. perfectly conducting, the
% *imaginary part of the contrast* is $\infty$: $\mathrm{imag}(q) = \infty$.
%
% * We *assume(!)*: q = 0.5 + i 1E6.
%
% * The *real part of contrast* is *assumed* from reconstruction...
% (nothing found in [1])
%
% * Metal: refractive index is highly dependent on wavelength...., i.e.
% our predefined contrast should depend on the choosen frequency 
% so will be a question of choosen frequency
%
% * But: |FmeasData| seems not to be highly dependent from choosen
% qValue...
%
%
%% References
%
% * [1] Kamal Belkebir and Marc Saillard. Special section on testing inversion algorithms against experimental data. _Inverse Problems_, 17(6):1565-1571, 2001.
% * [2] Andreas Kirsch: An integral equation approach and the interior transmission problem for Maxwell's equations. _Inverse Problems and Imaging_,1(1):107-127, 2007.
    
function q = fresnel_op1_rectTM_cent(X1,X2,varargin)

qValue = 0.5 + 1i*1E6; % perfectly conducting, so imaginary part is infinity...
% real part from reconstruction... not sure if correct...

w = 12.7E-3; % width = 12.7 mm
h = 24.5E-3; % height = 24.5 mm

right = 0; % move it right
up    = -3E-3; % move it up (or down with a negative number)
% shift the target 3 mm down... and it does visual better fit to Fresnel data

q = (-w/2+right <= X1) & (X1 <= w/2+right) & (-h/2+up <= X2) & (X2 <= h/2+up);
q = qValue*q;
q = double(q);

end
