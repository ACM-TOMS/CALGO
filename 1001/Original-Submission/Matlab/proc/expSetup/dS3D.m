%% dS3D
% Approximation of the infinitesimal element of a closed contour with control points.
%
%% Syntax
%
%   dS = dS3D(Nb, rad)
%
%% Description
%
% |dSp = dS3D(Nb, rad)| 
% computes approximation of the infinitesimal element on a
% surface containing points in three dimensional space.
%
% Currently, it is only written for spheres with radius rad on which Nb
% points are arranged on.
%
% *For simplicity, uniform points on surface are assumed.*
%
%% Input Arguments
%
% * Nb      : Number of points.
% * rad     : Radius of sphere transmitters are arranged on.
%
%% Output Arguments
%
% * dS   :   Approximation of the infinitesimal element of a closed contour 
%            with control points. 
%            Because assumption of uniform points dS is a value. 
%            If input is only one point, it is set dS = 1 (a closed contour does
% not make sense).
%
%
%% See Also
% * <dS2D.html>
% * <expSetup.html>
%
%% Code
%
function dS = dS3D(Nb, rad)

if Nb == 1
    dS = 1;
elseif Nb > 1
    % For simplicity, assume uniform points on surface
    dS = 4*pi*rad^2/Nb;
else
    disp('Error in function dS3D.m: size of input array too small.')
end

end
