%% gradientNeumannAdj
% Computes the discretized adjoint of the gradient with Neumann boundary
% conditions.
%
%% Syntax
%
%   res = gradientNeumannAdj(yg,seti)
%
%% Description
% |res = gradientNeumannAdj(yg,seti)| computes the adjoint of the 
% gradient using standard finite differences with Neumann boundary conditions.
%
% This function corresponds to <gradientNeumann.html>.
%
%% Input Arguments
% 
% * |yg|    is in 2D: 2*dim x nInv x nInv, see <gradientNeumann.html>.
% * |yg|    is in 3D: 2*dim x nInv x nInv x nInv, see <gradientNeumann.html>.
% * seti    : structural array.
% * seti.T, seti.S  : identify complex with real, see <setIdImagReal.html>.
% * seti.hInv, seti.nInv : see <setGridScale.html>.
% * seti.iG : write matrix as vector, see <setReshapeVecMat.html>.
%
%% Output Arguments
%
% * res : resulting adjoint of gradient
%
%% More About
%
% *Analytical*
%
% Analytical we have the relation
%
% $\langle \nabla u, p\rangle_Y = -\langle u, \mathrm{div} p\rangle_X$,
%
% i.e. $\nabla^\ast = -\mathrm{div}$, see [1, Sec. 6.1]$.
%
% *Discretization*
%
% For discretization we use explicitly  in case of 2D:
%
% $(\nabla y)^{(1)} = \texttt{gradAdjComp(1,:,:)} =$
%
% * (case i = 1):     $\quad$ $-y_{1,j}/h$
% * (case 1 < i < n): $\quad$ $-(y_{i,j}-y_{i-1,j})/h$
% * (case i = n):     $\quad$ $y_{n-1,j}/h$
%
% The case $(\nabla y)^{(2)}$ is analog.
%
% An extension to 3D is straightforward.
%
%% References
%
% * [1] Antonin Chambolle and Thomas Pock. A first-order primal-dual algorithm for convex problems with applications to imaging. _Journal of Mathematical Imaging and Vision_, 40(1):120-145, 2011.
%
%% See Also
% * <gradientNeumann.html>
%
%% Code

function res = gradientNeumannAdj(yg,seti)

ygz = seti.T(yg); % ygz: dim x nInv x nInv (x nInv): so GD and complex
% if seti.dim == 2
%     x = divergence(squeeze(ygz(1,:,:)),squeeze(ygz(2,:,:)));
% elseif seti.dim == 3
%     x = divergence(squeeze(ygz(1,:,:,:)),squeeze(ygz(2,:,:,:)),squeeze(ygz(3,:,:,:)));
% end
% % x is a complex matrix of dimension 2 or 3
% res = seti.iG(x); % write matrix as vector
% res = seti.S(res); % result: matrix in R^2
% %res is xnRVD

h = seti.hInv;
N = seti.nInv;
% G = seti.G (vector -> matrix)

gAdj = gradAdjCalc(ygz,h,N);
res = seti.iG(gAdj); % write matrix as vector
res = seti.S(res); % result: matrix in R^2: res is xnRVD

end

%% 
% *Code: subfunction: guAdjCalc (fast)*
function res = gradAdjCalc(y,h,N)
% N = nInv
dim = size(y,1); % y is ygz, so dim x nInv x nInv (x nInv)

if dim == 2
    gradAdjComp = zeros(dim,N,N);
    
    % gradAdjComp(1,:,:) = gradAdj(y)_{i,j}^1
    % = (case i = 1):     -y(1,j)/h
    %   (case 1 < i < n): -(y(i,j)-y(i-1,j))/h
    %   (case i = n):     y(n-1,j)/h
    yRowShift = circshift(y,[0 +1 0]);
    gradAdjComp(1,:,:) = squeeze(-(y(1,:,:)-yRowShift(1,:,:))/h); % case 1 < i < n
    gradAdjComp(1,1,:) = squeeze(-y(1,1,:)/h); % case i = 1
    gradAdjComp(1,end,:) = squeeze(y(1,end-1,:)/h); % case i = n

    % gradAdjComp(2,:,:) = gradAdj(y)_{i,j}^2
    yColShift = circshift(y,[0 0 +1]);
    gradAdjComp(2,:,:) = squeeze(-(y(2,:,:)-yColShift(2,:,:))/h); % case 1 < j < n
    gradAdjComp(2,:,1) = squeeze(-y(2,:,1)/h); % case j = 1
    gradAdjComp(2,:,end) = squeeze(y(2,:,end-1)/h); % case j = n
    
    % gradAdj = gradAdj(y)_{i,j}^1 + gradAdj(y)_{i,j}^2
    gradAdj = gradAdjComp(1,:,:) + gradAdjComp(2,:,:); % size N x N
    res = gradAdj;
elseif dim == 3
    gradAdjComp = zeros(dim,N,N,N);
    
    % gradAdjComp(1,:,:,:) = gradAdj(y)_{i,j,k}^1
    y1Shift = circshift(y,[0 +1 0 0]);
    gradAdjComp(1,:,:,:) = squeeze(-(y(1,:,:,:)-y1Shift(1,:,:,:))/h); % case 1 < i < n
    gradAdjComp(1,1,:,:) = squeeze(-y(1,1,:,:)/h); % case i = 1
    gradAdjComp(1,end,:,:) = squeeze(y(1,end-1,:,:)/h); % case i = n

    % gradAdjComp(2,:,:) = gradAdj(y)_{i,j,k}^2
    y2Shift = circshift(y,[0 0 +1 0]);
    gradAdjComp(2,:,:,:) = squeeze(-(y(2,:,:,:)-y2Shift(2,:,:,:))/h); % case 1 < j < n
    gradAdjComp(2,:,1,:) = squeeze(-y(2,:,1,:)/h); % case j = 1
    gradAdjComp(2,:,end,:) = squeeze(y(2,:,end-1,:)/h); % case j = n

    % gradAdjComp(3,:,:) = gradAdj(y)_{i,j,k}^3
    y3Shift = circshift(y,[0 0 0 +1]);
    gradAdjComp(3,:,:,:) = squeeze(-(y(3,:,:,:)-y3Shift(3,:,:,:))/h); % case 1 < k < n
    gradAdjComp(3,:,:,1) = squeeze(-y(3,:,:,1)/h); % case k = 1
    gradAdjComp(3,:,:,end) = squeeze(y(3,:,:,end-1)/h); % case k = n

    % gradAdj = gradAdj(y)_{i,j,k}^1 + gradAdj(y)_{i,j,k}^2 + gradAdj(y)_{i,j,k}^3
    gradAdj = gradAdjComp(1,:,:) + gradAdjComp(2,:,:) + gradAdjComp(3,:,:); % size N x N x N
    res = gradAdj;

else
    error('gradientNeumannAdj.m: dimension must be 2D or 3D.')
end
end

