%% normTVinv1
% Computes the total variation semi-norm.
%
%% Syntax
%
%   n = normTVinv1(gus,seti)
%
%% Description
% |n = normTVinv1(gus,seti)| computes the total variation semi-norm on
% gradient |gus|.
%
%% Input Arguments
%
% * gus         : gradient(u): real matrix of size 2*dim x nInv x nInv (x nInv)
%                 (i.e. u was not on ROI, but down scaled ROI).
% * seti.dim    : dimension of the problem
% * seti.nInv   : see <setGridScale.html>
% * seti.dVinv  : see <setGridScale.html>
%
%% Output Arguments
%
% * n   : result of the TV semi-norm.
%
%% More About
%
% $n = \|\nabla u\|_{\mathrm{tv},\bf{R}} = (\|\nabla u\|_{\mathrm{tv},1})$
%
% where $\nabla u$ was already computed and rewritten as real in |gus| 
% (gradient of u in real, i.e. squared because of identification 
%  $\bf{C} = \bf{R} \times \bf{R}$, see <setIdImagReal.html>).
% 
% In case of 2D (3D is analog) (see [1, Sec. 4.5, eq. (56)]:
%
% $\|a\|_{\mathrm{tv},\bf{R}} := h_N^2 \sum_i
% |(a_i^{(1),\mathrm{Re}}, a_i^{(1),\mathrm{Im}}, 
%  a_i^{(2),\mathrm{Re}}, a_i^{(2),\mathrm{Im}})|$
%
% with $a = \texttt{gus}$.
%
% The semi-norm is used in <setFuncsPda.html> to define TV-penalty |seti.fg|.
%
% * gus is RMD: real matrix down scaled ROI 
%  (size: 2*dim x nInv x nInv in 2D and 2*dim x nInv x nInv x nInv in 3D).
% * guz is CMD: complex matrix down scaled ROI
%  (size: dim x nInv x nInv and dim x nInv x nInv x nInv in 3D).
% 
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.

%% See Also
%
% * <innergrad.html>
%
%% Code
function n = normTVinv1(gus,seti)

absgu = normTVinvAux(gus,seti);

if seti.dim == 2
    n = sum(sum(absgu))*seti.dVinv; % \|grad(u)\|_{TV,1} =  \sum_i,j | (grad(u))_i,j |
elseif seti.dim == 3
    n = sum(sum(sum(absgu)))*seti.dVinv;
end

end
