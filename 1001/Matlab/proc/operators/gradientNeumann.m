%% gradientNeumann
% Discretization of the gradient using standard finite differences with
% Neumann boundary conditions.
%
%% Syntax
%
%   res = gradientNeumann(xnRVD,h,N,G,seti)
%
%% Description
% |res = gradientNeumann(xnRVD,h,N,G,seti)| computes the discretized
% gradient |res| of |xnRVD| (matrix on down scaled ROI stored as a vector)
% with length of the infinitesimal element |h| and 
% discretization points for each dimension |N|.
%
%% Example
%
% From <setFuncsPda.html>:
%
% gradientNeumann(xnRVD,seti.hInv,seti.nInv,seti.GInv,seti)
%
%% Input Arguments
%
% * |xnRVD| : matrix on down scaled ROI stored as a real vector
%            (i.e. in the code: firstly we identify $\bf{R} \times \bf{R}$ with $\bf{C}$
%             and secondly reshape the vector to a matrix in the down scaled ROI).
% * h       : length of the infinitesimal element, see <setGrid.html>
% * N       : discretization points for each dimension, see <setGrid.html>
% * G       : function to reshape a vector to a matrix, see
% <setReshapeVecMat.html>
% * seti    : structural array, in this file we need the transformation 
%             operators |seti.T| and |seti.S| to identify 
%             $\bf{C}$ with $\bf{R} \times \bf{R}$.
%
%% Output Arguments
%
% * res     : gradient of xnRVD
%             (real array of size 2*2 x N x N in case of 2D and
%              size 2*3 x N x N x N in case of 3D).
%
%% More About
%
% Discretization of $\nabla: X \to Y$ using standard 
% finite differences with Neumann boundary conditions, 
% see [1, Sec. 6.1] or [2, Sec. 4.3].
%
% We present the 2D case (this can be extended to 3D).
%
% $(\nabla u)_{i,j} = ( (\nabla u)_{i,j}^{(1)}; (\nabla u)_{i,j}^{(2)} )$
%
% with 
%
% * $(\nabla u)_{i,j}^{(1)} = (u_{i+1,j}-u_{i,j})/h \quad \mathrm{if\ } i < N$.
% * $(\nabla u)_{i,j}^{(1)} = 0 \quad \mathrm{if\ } i = N$
%
% and
%
% * $(\nabla u)_{i,j}^{(2)} = (u_{i,j+1}-u_{i,j})/h \quad \mathrm{if\ } j < N$.
% * $(\nabla u)_{i,j}^{(2)} = 0 \quad \mathrm{if\ } j = N$
% 
%% References
%
% * [1] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
% * [2] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
%
% * <gradientNeumannAdj.html>
%
%% Code

function res = gradientNeumann(xnRVD,h,N,G,seti)

% xnRVD: real, stored as vector, down scaled grid
xnCVD = seti.T(xnRVD); % seti.T, because we use complex values in this function

u = G(xnCVD); % u: matrix: nInv^dim
% u = xnCMD

% h = seti.hInv
% N = seti.nInv
% G = seti.G (reshape vector to matrix)

% test this code:
% 1. run start (to set some vectors and functions)
% 2. gradientNeumann(seti.S(seti.qROIexact),seti.hInv,seti.nInv,seti.G,seti)

%%
% *Compare slow and fast version of computation of gradient u*

if 0
    disp('test gradient: slow and fast version')
    for dim = 2:3
        
        fprintf('dimension %g\n',dim)
        N = 100;
        if dim == 2
            u = rand(N,N);
        elseif dim == 3
            u = rand(N,N,N);
        end
        h = 0.1;

        tic
        disp('  gradient u calculation (slow version)')
        gu = guCalcSlow(u,h,N);
        toc

        tic
        disp('  gradient u calculation (fast version)')
        gu2 = guCalc(u,h,N);
        toc

        fprintf('  max difference: %g \n',max(unique(abs(gu-gu2))))
    end
    error('stop: test in gradient neumann end')
end

%%
% *Computation of the gradient*

gu = guCalc(u,h,N); % gu is a COMPLEX array of size 
                    % 2 x N x N in case of 2D 
                    % and size 3 x N x N x N in case of 3D.

%%
% *Save as real array*

res = seti.S(gu);   % identify complex C with real R x R, then res is real 
                    % of size 2*2 x N x N in case of 2D 
                    % and size 2*3 x N x N x N in case of 3D.
end

%%
% *Code: subfunction: guCalc (fast computation)*

function res = guCalc(u,h,N)
dim = length(size(u));

if dim == 2
    gu = zeros(dim,N,N);
    
    % gu(1,:,:) = grad(u)_{i,j}^1 = 1/h*(u(i+1,j)-u(i,j)) if i < N (0 else)
    uRowShift = circshift(u,[-1 0]);
    uRowShift(end,:) = uRowShift(end-1,:); % to get 0 entry if i = N
    gu(1,:,:) = 1/h*(uRowShift-u);

    % gu(2,:,:) = grad(u)_{i,j}^2 = 1/h*(u(i,j+1)-u(i,j)) if j < N (0 else)
    uColShift = circshift(u,[0 -1]);
    uColShift(:,end) = uColShift(:,end-1); % to get 0 entry if j = N
    gu(2,:,:) = 1/h*(uColShift-u);

    res = gu;
elseif dim == 3
    gu = zeros(dim,N,N,N);
    
    % gu(1,:,:,:) = grad(u)_{i,j,k}^1 = 1/h*(u(i+1,j,k)-u(i,j,k)) if i < N (0 else)
    u1Shift = circshift(u,[-1 0 0]);
    u1Shift(end,:,:) = u1Shift(end-1,:,:); % to get 0 entry if i = N
    gu(1,:,:,:) = 1/h*(u1Shift-u);

    % gu(2,:,:) = grad(u)_{i,j,k}^2 = 1/h*(u(i,j+1,k)-u(i,j,k)) if j < N (0 else)
    u2Shift = circshift(u,[0 -1 0]);
    u2Shift(:,end,:) = u2Shift(:,end-1,:); % to get 0 entry if j = N
    gu(2,:,:,:) = 1/h*(u2Shift-u);

    % gu(3,:,:) = grad(u)_{i,j,k}^3 = 1/h*(u(i,j,k+1)-u(i,j,k)) if k < N (0 else)
    u3Shift = circshift(u,[0 0 -1]);
    u3Shift(:,:,end) = u3Shift(:,:,end-1); % to get 0 entry if k = N
    gu(3,:,:,:) = 1/h*(u3Shift-u);
    
    res = gu;
else
    error('gradientNeumann.m: dimension must be 2D or 3D.')
end
end

%%
% *Code: subfunction: guCalcSlow (slow computation)*
%
% Is not used any more.

function res = guCalcSlow(u,h,N)

dim = length(size(u));
if dim == 2
    gu = zeros(dim,N,N); % grad(u) stored in components gu(1,:,:), gu(2,:,:), (3D: gu(3,:,:))
    % 2D
    for l = 1:dim
        gu(l,:,:) = gradulCalc2D(l,u,h,N); % (grad u)_{i,j}^l
    end
    res = gu;
elseif dim == 3
    gu = zeros(dim,N,N,N); 
    % 3D
    for l = 1:dim
        gu(l,:,:,:) = gradulCalc3D(l,u,h,N); % (grad u)_{i,j,k}^l
    end
    res = gu;
else
    error('gradientNeumann.m: dimension of u must be 2D or 3D.')
end
end

%%
% *Code: subfunction: gradulCalc2D*
%
% gradulCalc2D is required in guCalcSlow.

function res = gradulCalc2D(l,u,h,N)
% compute (grad u)_{i,j}^l
% l = 1 and 2 in 2D
% l = 1, 2, 3 in 3D
gradul = zeros(N,N);
for i = 1:N
    for j = 1:N
        if (l == 1 && i < N) || (l == 2 && j < N)
            switch l
                case 1
                    gradul(i,j) = (u(i+1,j)-u(i,j))/h; % first component
                case 2
                    gradul(i,j) = (u(i,j+1)-u(i,j))/h; % second component
            end
        else %case i = N
            gradul(i,j) = 0;
        end
    end
end
res = gradul;
end

%%
% *Code: subfunction: gradulCalc3D*
%
% gradulCalc3D is required in guCalcSlow

function res = gradulCalc3D(l,u,h,N)
% compute (grad u)_{i,j}^l
% l = 1 and 2 in 2D
% l = 1, 2, 3 in 3D
gradul = zeros(N,N,N);
for i = 1:N
    for j = 1:N
        for k = 1:N
            if (l == 1 && i < N) || (l == 2 && j < N) || (l == 3 && k < N)
                switch l
                    case 1
                        gradul(i,j,k) = (u(i+1,j,k)-u(i,j,k))/h; % first component
                    case 2
                        gradul(i,j,k) = (u(i,j+1,k)-u(i,j,k))/h; % second component
                    case 3
                        gradul(i,j,k) = (u(i,j,k+1)-u(i,j,k))/h; % second component
                end
            else %case i = N
                gradul(i,j,k) = 0;
            end
        end
    end
end
res = gradul;
end
