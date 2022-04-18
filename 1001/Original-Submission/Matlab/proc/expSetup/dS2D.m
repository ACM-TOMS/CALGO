%% dS2D
% Approximation of the infinitesimal element of a closed contour with control points.
%
%% Syntax
%
%   dSp = dS2D(p)
%
%% Description
% |dSp = dS2D(p)| computes approximation of the infinitesimal element of a
% closed contour with control points.
% p is a real matrix of size Nb x 2. 
% dSp is a vector of size Nb x 1.
%
% * If all elements of dSp are almost the same dSp is a value.
% * If input is only one point, it is set dSp = 1 (a closed contour does
% not make sense).
%
% Note that p is a *matrix of size Nb x 2 and not 2 x Nb*.
% This is why in pntsGeometry is used:
%
%   dS = dS2D(transpose(Pnts));
%
%
%% Examples
%
% *Example 1*
%
% Point coordinates: (0,0), (1,0), (1,1)
%
%   p = [0 0; 1 0; 1 1];
%   dS2D(p)
%
% Result:
%
%   1.2071
%   1.0000
%   1.2071
%
% *Example 2*
%
% Point coordinates: (0,0), (1,0), (1,1), (0,1) (a square)
%
%   p = [0 0; 1 0; 1 1; 0 1];
%   dS2D(p)
%
% Result:
%
%   1
%
% Result is a value and not a vector because if all elements of dSp are 
% almost the same size, it is reduced to a single number.
%
%
%% Input Arguments
%
% * p     :   coordinates of points (real matrix of size Nb x 2)
%
%% Output Arguments
%
% * dSp   :   approximation of the infinitesimal element of a closed contour 
%             with control points. dSp is a vector of size Nb x 1. 
%             If all elements of dSp are almost the same dSp is a value.
%             If input is only one point, it is set dSp = 1 (a closed contour does
% not make sense).
%
%
%% More About
%
% Please note that |dSp| does not contain the distances between the points,
% because then it would depend on the orientation of the parametrization. 
% To omit this problem we use an average value (similar to the central 
% difference in comparison to the forward or backward difference). It is
% explained in the following.
% 
% *First, but wrong idea*
%
% The first (but kind of wrong) idea to discretize a closed contour is to 
% compute the distances between the points, e.g. 
% $dS(p_n) = \|p_{n+1} - p_n\|_2$ for points $p_n$.
%
% Then, the distances for the point coordinates $p_1 = (0,0)$, $p_2 = (1,0)$
% and $p_3 = (1,1)$ are 1, $\sqrt{2}$ and 1.
%
% The problem is that the infinitesimal element then depends on the
% orientation of the parametrization. We consider the parametrization 
% $[p_1, p_2, p_3, (p_1)]$ and the (equivalent) one $[p_3, p_2, p_1, (p_3)]$. 
% Then the length of the infinitesimal element of $p_2$ changes: 
% In the first case it is $dS(p_2) = \|p_3 - p_2\|_2 = \sqrt{2}$ and 
% in the second $dS(p_2) = \|p_1 - p_2\|_2 = 1$.
%
% *Second, orientation independent idea*
%
% The second idea avoids the orientation's dependence using the forward and
% backward definition of the distances between the points, i.e. 
% $dS(p_n) = \|p_{n+1} - p_n\|_2$ and $dS(p_n) = \|p_n - p_{n-1}\|_2$ 
% and compute the average, i.e.
%
% $$ dS(p_n) = \frac{1}{2} \left( \|p_{n+1} - p_n\|_2 + \|p_n - p_{n-1}\|_2 \right).$$
%
% This results in the infinitesimal elements are $\quad$
% $1+(\sqrt(2)-1)/2 \approx 1.2071$, $\quad$
% $1$ $\quad$ and $\quad$
% $1+(\sqrt(2)-1)/2 \approx 1.2071$.
% 
% The definition is chosen such that the point $p_n$ is in the center of the 
% closed contour's element belonging to $dS(p_n)$. 
% This is kind of a 'staggered definition' between $dS$ and $p$ 
% (similar to the central difference in comparison to the forward or backward difference).
%
%
%% See Also
% * <pntsGeometry.html>
% * <dS3D.html>
%
%% Code
%
function dSp = dS2D(p)

if numel(p) == 2 % i.e. ONE coordinate (a distance does not make sense)
    dSp = 1;
elseif numel(p) > 2
    pm1 = circshift(p, 1); % p minus 1
    pp1 = circshift(p,-1); % p plus 1
    n = @(x) sqrt(sum( x.^2 ,2));
    dSp  = (  n(p-pm1) + n(p-pp1)  )/2;
    % if all elements of dSp almost the same size,
    % reduce dSp to a single number
    if var(dSp) < 10E-5*mean(dSp)
        dSp = mean(dSp);
    end
else
    disp('Error in function dSp: size of input array too small')
end
end
