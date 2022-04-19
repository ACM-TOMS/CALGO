%% normroi
% Computes the norm of a vector on the region of interest.
% 
%% Syntax
%
%   n = normroi(x,seti,q)
%   n = normroi(x,seti)
%
%% Description
%
% * |normroi(x,seti,q)| computes the norm of the complex vector x 
% on the region of interest. The norm is a q-Norm with weight |seti.dV|.
%
% * |normroi(x,seti)| computes the norm of the complex vector x 
% on the region of interest. The norm is a seti.qNorm-Norm with weight |seti.dV|.
%
%% Input Arguments
%
% * x, y    :   two complex vectors of size seti.nROI^seti.dim
%               (values in region of interest written as a vector).
% * seti.dV :   area/volume of the infinitesimal element (pixel/voxel) of CD,
%               <see setGrid.html>.
% * seti.qNorm  :   seti.qNorm-Norm
% * q           :   q-Norm 
%                   (third argument is not necessary, if seti.qNorm was defined).
%
%% Output Arguments
%
% * n : evaluated norm
%
%% More About
%
% * In case of seti.qNorm = 2 or q = 2, the norm is the same as the norm
% induced by the inner product <innerroi.html>.
%
% * The norm is defined as, see [1, Sec. 3.6],
%
% $\|x\|_{\mathrm{roi},q} := \left( h_N^d \sum_i |x_i|^q \right)^{1/q}$
% for $x \in \bf{C}^{N_D}$
%
% with $N_D = \texttt{seti.nROI}^\texttt{seti.dim}$ 
% and $h_N^d = \texttt{seti.dV}$.
%
%% References
%
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <innerroi.html>
%
%% Code
function n = normroi(x,seti,q)

if nargin == 2
    Q = seti.qNorm;
elseif nargin == 3
    Q = q; 
end

dV = seti.dV;

if Q == 2
    n = sqrt(abs(innerroi(x,x,seti)));
else
    n = norm(x(:),Q)*dV^(1/Q); % in case Q = 2 the result is the same as above
end

end
