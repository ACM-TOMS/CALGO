function [A, b]  = make_system( eps, w, r, h );

%
% Generate test problem. The matrix is the FDM-discretisation of the (constant coefficient) 
% Convection-diffusion-reaction equation with zero-order reaction term.
% The right-hand-side is determined by the pre-defined solution u = xyz(1-x)(1-y)(1-z).
%
% The routine is called as folows:
%     [A, b]  = make_system( eps, w, r, h );
%
%     Input parameters:
%        eps: diffusion coefficient (scalar)
%        w:   convections coefficients (3-dimensional array, entries correspond to w_x, w_y and w_z)
%        r:   reaction parameter
%        h:   gridsize
%
% make-system is called in test_idrs.m.
%
%   Martin van Gijzen
%   Version August 31, 2010
%
%   This software is distributed under the
%   ACM Software Copyright and License Agreement.
%

% Show parameters:
disp('FDM discretisation of a 3D convection-diffusion-reaction problem on a unit cube');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

disp(' ');
disp('The parameters of the problem are :');
disp(['Gridsize h = ',num2str(h),';']);
disp(['Diffusion parameter eps = ', num2str(eps),';']);
disp(['Convection parameters w = (',num2str(w(1)),',',num2str(w(2)),',',num2str(w(3)),');']);
disp(['Reaction parameter r = ',num2str(r),' (Note: linear reaction term, gives negative shift to matrix);']);
disp(' ');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(' ');

% Generate matrix
m = round(1/h)-1;
n = m*m*m;
e = ones(m,1);
Sx = spdiags([(-eps/h^2-w(1)/(2*h))*e (2*eps/h^2)*e (-eps/h^2+w(1)/(2*h))*e], -1:1, m, m );
Sy = spdiags([(-eps/h^2-w(2)/(2*h))*e (2*eps/h^2)*e (-eps/h^2+w(2)/(2*h))*e], -1:1, m, m );
Sz = spdiags([(-eps/h^2-w(3)/(2*h))*e (2*eps/h^2)*e (-eps/h^2+w(3)/(2*h))*e], -1:1, m, m );
Is = speye(m,m);
I = speye(n,n);
A = kron(kron(Is,Is),Sx) + kron(kron(Is,Sy),Is)+ kron(kron(Sz,Is),Is) -r*I;

x = linspace(h,1-h,m);
sol = kron(kron(x.*(1-x),x.*(1-x)),x.*(1-x))';
b = A*sol;

return
