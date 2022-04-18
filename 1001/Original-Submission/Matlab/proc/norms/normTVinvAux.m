%% normTVinvAux
% Auxiliary function called in normTVinv1 and normTVinv2.
%
%% Syntax
%
%   absgu = normTVinvAux(gus,seti)
%
%% Description
% |absgu = normTVinvAux(gus,seti)| computes the auxiliary matrix absgu.
%
%% Input Arguments
%
% The input arguments are described in <normTVinv1.html>:
%
% * gus
% * seti.dim
% * seti.nInv
% * seti.dVinv
%
%% Output Arguments
%
% * absgu : See "More About".
%
%% More About
%
% In case of 2D: 
%
% $\texttt{absgu} = | (\mathrm{grad}(u))_{i,j} | = | (\nabla(u)_{i,j} | = 
% \sqrt{ |(\nabla u)_{i,j}^1|^2 + |(\nabla u)_{i,j}^2|^2 }$,
%
% see [1] or [2, Sec. 4.3, eq. (43)].
%
% Note that the input $\texttt{gus}$ is $\mathrm{grad}(u)$, but rewritten
% as real matrix by identification $\bf{C} = \bf{R} \times \bf{R}$.
%
% The 3D case is analog.
%
%% References
%
% * [1] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <normTVinv1.html>
% * <normTVinv2.html>
%
%% Code

function absgu = normTVinvAux(gus,seti)

guz = seti.T(gus); % transform gus from R x R into C

%sgu = size(gu,2); % 2nd entry of gu is nROI or nInv

if seti.dim == 2 && ~isequal(size(guz),[2 seti.nInv seti.nInv])
    error('normTVinvAux: guz has wrong dimension')
elseif seti.dim == 3 && ~isequal(size(guz),[3 seti.nInv seti.nInv seti.nInv])
    error('normTVinvAux: guz has wrong dimension')
end

if 0 % test
    disp('inside normTVinvAux')
    size(guz)
    gus = seti.S(guz);
    size(gus)
    max(abs(squeeze(real(guz(1,:,:))-gus(1,:,:)))) % OK
    max(abs(squeeze(real(guz(2,:,:))-gus(2,:,:)))) % OK
    max(abs(squeeze(imag(guz(1,:,:))-gus(3,:,:)))) % OK
    max(abs(squeeze(imag(guz(2,:,:))-gus(4,:,:)))) % OK

    error('stop')
end

if seti.dim == 2
    %absgu = zeros(seti.nInv,seti.nInv);
    absgu = squeeze( sqrt( abs(guz(1,:,:)).^2 + abs(guz(2,:,:)).^2 ) ); % absgu = | (grad(u))_i,j |
elseif seti.dim == 3
    %absgu = zeros(seti.nInv,seti.nInv,seti.nInv);
    absgu = squeeze( sqrt( abs(guz(1,:,:)).^2 + abs(guz(2,:,:)).^2 + abs(guz(3,:,:)).^2) );
end

end
