%% innergrad
% Inner product of gradient evaluations stored as real matrices (arrays).
%
%% Syntax
%
%   res = innergrad(gus,gvs,seti)
%
%% Description
%
% |res = innergrad(gus,gvs,seti)| computes the inner product of the
% gradient evalutations 
% $\langle \mathrm{grad}(u), \mathrm{grad}(v)\rangle|$, which are stored as
% real matrices (arrays).
%
%% Input Arguments
%
% * gus     : gradient(u): real matrix of size 2*dim x nInv x nInv (x nInv)
% * gvs     : gradient(v): real matrix of size 2*dim x nInv x nInv (x nInv)
% * seti.dim    : dimension of the problem
% * seti.dVinv  : see <setGridScale.html>
% * seti.nInv   : see <setGridScale.html>
%
% Note: The name "gus" was choosen because *gradient* from *u* stored as 
% real (*squared*, because $\bf{C} = \bf{R} \times \bf{R}$).
%
%% Output Argument
%
% * res     :   inner product, see More About for definition.
%
%% More About
%
% We consider the inner product of p and q (e.g. p = grad(u), q = grad(v)):
%
% $\langle p,q \rangle = \sum_{i,j} p_{i,j}^1 q_{i,j}^1 + p_{i,j}^2 +
% q_{i,j}^2$ in 2D
%
% and $+ p_{i,j}^3 + q_{i,j}^3$ in 3D.
%
% Note that the exponent is an index(!).
%
% The real and imaginary parts are stored as real values in gus and guv,
% so add all parts (real and imag):
%
% $\langle p,q \rangle = \sum_{i,j} \sum_{l=1}^{2\,\mathrm{dim}} p_{i,j}^l q_{i,j}^l$
%
% Finally, scale it with factor |dVinv|.
%
% *Corresponding notation in [1]*
%
% The function |innergrad| is defined in [1, Sec. 4.5, eq. (59)] as
%
% $\langle x, y \rangle_{\mathrm{tv},\bf{R}} := 
%  \langle x^{(1)}, y^{(1)} \rangle_{\mathrm{roi},\bf{R}} + 
%  \langle x^{(2)}, y^{(2)} \rangle_{\mathrm{roi},\bf{R}}, \quad
%  x, y \in Y_{\mathrm{tv},\bf{R}}$
%
% with the following definitions in [1, Sec. 4.5]:
%
% * $Y_{\mathrm{tv},\bf{R}} = X_{\bf{R}} \times X_{\bf{R}}$.
%
% * $X_{\bf{R}} = \bf{R}^{N_D} \times \bf{R}^{N_D}$.
%
% * So $x = (x^\mathrm{Re},x^\mathrm{Im}) \in X_{\bf{R}}$
% is a rewrited complex vector as real matrix, 
% see also <setIdImagReal.html>.
%
% * $N_D = \texttt{seti.nROI}^\texttt{seti.dim}$.
%
% * $\langle x^{(1)}, y^{(1)} \rangle_{\mathrm{roi},\bf{R}}
%  = h_N^d \sum_i
%  (x_i^\mathrm{Re} y_i^\mathrm{Re} + x_i^\mathrm{Im} y_i^\mathrm{Im})$
% $\quad$ with $x, y \in X_{\bf{R}}$.
%
% * $h_N^d = \texttt{seti.dVinv}$.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <setGridScale.html>
%
%% Code
%
function res = innergrad(gus,gvs,seti)

% guz = seti.T(gus);
% gvz = seti.T(gvs);

add = 0;
for l = 1:2*seti.dim
    add = add + squeeze(gus(l,:,:) .* gvs(l,:,:));
end

if seti.dim == 2
    res = sum(sum(add))*seti.dVinv;
elseif seti.dim == 3
    res = sum(sum(sum(add)))*seti.dVinv;
end

end
