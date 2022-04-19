%% inerroi
% Computes the inner product on ROI.
%
%% Syntax
%   i = innerroi(x,y,seti)
%
%% Description
% |innerroi(x,y,seti)| computes the inner product on the region of interest
% from two complex vectors x, y.
%
%% Input Arguments
%
% * x, y    :   two complex vectors of size seti.nROI^seti.dim
%               (values in region of interest written as a vector).
% * seti.dV :   area/volume of the infinitesimal element (pixel/voxel) of CD,
%               <see setGrid.html>.
%
%% Output Arguments
%
% i     :   inner product on ROI, see More About.
%
%% More About
% 
% The inner product is, see [1, Sec. 3.6, eq. (33)]:
%
% $\langle x, y \rangle_\mathrm{roi} = h_N^d \sum_i x_i \overline{y_i}$
%
% with $h_N^d = \texttt{seti.dV}$ 
% for $x,y \in \bf{C}^{N_d}$
% with $N_D = \texttt{seti.nROI}^\texttt{seti.dim}$.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <innerinv.html>
% * <normroi.html>
%
%% Code
function i = innerroi(x,y,seti)
i = (y'*x)*seti.dV;
end
