
function [XD,YD,XM,YM] = slt_particle(A,U,V,dt)
% Semi-Lagrangian transport particle tracking 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% A is the extended slt_grid
% A.X,A.Y are the grid values on the interior produced by meshgrid
% A.Xe,A.Ye are the extended grid values produced by meshgrid
% U = u*cos(lat) is the spherical velocity over the time interval
% V = v*cos(lat) is the spherical velocity over the time interval
% dt is the time interval from level n to n+1 
% Output: 
% XD,YD is the resulting departure points at time level n
% XM,YM is the resulting mid points at time level n+0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the departure points by integrating the position equation
% using the extrapolated half time velocities
% find the midpoint M position and half velocity at midpoint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  -------------------------------------
%  |                                   |
%  |       C                D          |
%  |  +++++++++++++++N++++++++++++++   +1.0
%  |  +                            +   |
%  |  +                            +   |
%  |  +            interior        + periodic
%  |  +                            +   |
%  |  +                            +   |
%  |  +++++++++++++++S++++++++++++++   -1.0 (= sin lat)
%  | 0.0   A                B      2*pi
%  |                                   |
%  |                                   |
%  -------------------------------------
%     +nxpt+1                      +nxpt+nx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UM = U';  %transpose for correct meshgrid/interp matrix
VM = V';
%  Particle tracking must be done in Cartesian space to avoid pole problems.
%  Cartesian points are then injected back into the interior extended grid
%  for interpolation.  The tensor product 2-D interpolation is too nice to give 
%  up for 3-D interpolation on the manifold.
%   ( XC,YC,ZC are Cartesian particle arrays.)
%   ( UC,VC,WC are Cartesian velocity arrays.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [XC,YC,ZC] = S2Cscalar(A.X,A.Y);
%XC pause  % debug
%YC pause  % debug
% ZC pause % debug
for iteration = 1:3  % only a few iterations needed for midpoint
  [UC,VC,WC] = S2Cvector(A.X,A.Y,UM,VM);
%UC pause  % debug
%VC pause  % debug
%WC pause  % debug
  XCM = XC - 0.5*dt*UC;   % note half dt step for midpoint
  YCM = YC - 0.5*dt*VC;
  ZCM = ZC - 0.5*dt*WC;
% Project back to interpolation domain [0:2*pi]x[-1,1]
  [XM,YM] = C2Sscalar(XCM,YCM,ZCM);
%XM pause  % debug
%YM pause  % debug
% Interpolate velocities at the midpoint on the extended slt_grid
  Ue  = slt_extend(A,U);
  UM=interp2(A.Xe,A.Ye,Ue,XM,YM,'linear');  % linear interpolation of velocity
%UM pause  % debug
  Ve  = slt_extend(A,V);
  VM=interp2(A.Xe,A.Ye,Ve,XM,YM,'linear');
%VM pause  % debug
end
% find the departure point taking the full dt
  [UC,VC,WC] = S2Cvector(A.X,A.Y,UM,VM);
  XCM = XC - dt*UC;
  YCM = YC - dt*VC;
  ZCM = ZC - dt*WC;
% Project back to interpolation domain [0:2*pi]x[-1,1]
  [XD,YD] = C2Sscalar(XCM,YCM,ZCM);
%XD pause  % debug
%YD pause  % debug
%  end function slt_particle


function [XC,YC,ZC] = S2Cscalar(X,Y);
%Usage: [XC,YC,ZC] = S2Cscalar(X,Y);
%Purpose: Transform the Spherical points in X,Y to Cartesian points on the unit sphere
% Input:  X  	 longitude in [0,2*pi] 
%         Y  	 sin(latitude) in [-1,1] 
% Output: XC,YC,ZC Cartesian coordinates of points
% Uses MATLAB function sph2cart
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = ones(size(X));
[XC,YC,ZC] = sph2cart(X,asin(Y),R);
%end function S2Cscalar

function  [UC,VC,WC] = S2Cvector(X,Y,U,V);
%Usage: [UC,VC,WC] = S2Cvector(X,Y,U,V);
%Purpose: Transform the Spherical vectors in U,V to Cartesian vectors on the unit sphere
% Input:  X  	 longitude in [0,2*pi] 
%         Y  	 sin(latitude) in [-1,1] 
%         U  	 u*cos(lat) spherical velocity
%         V  	 v*cos(lat) spherical velocity
% Output: UC,VC,WC Cartesian velocity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref:  Paul Swarztrauber,
%       Vector functions and their Derivatives on the Sphere, 
%       SIAM J. Numer Anal,  Vol 18, No 2, p 195, (1981)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transform matrix Q
W = zeros(size(V));
UC = zeros(size(V));
VC = zeros(size(V));
WC = zeros(size(V));
for i = 1: numel(X)
lat = asin(Y(i));
lon = X(i);
Q = [ -sin(lon)  -sin(lat)*cos(lon)  cos(lat)*cos(lon) ; 
       cos(lon)  -sin(lat)*sin(lon)  cos(lat)*sin(lon) ;
       0.0            cos(lat)          sin(lat)       ];
u = U(i)/cos(lat); v = V(i)/cos(lat); w = W(i);
vel = Q * [u v w]';
UC(i) = vel(1);
VC(i) = vel(2);
WC(i) = vel(3);
end
%end function S2Cvector

function  [X,Y] = C2Sscalar(XC,YC,ZC);
%Usage: [X,Y] = C2Sscalar(XC,YC,ZC);
%Purpose: Transform the Cartesian points on unit sphere to Spherical points in X,Y
% Input: XC,YC,ZC Cartesian coordinates of points
% Output: X  	 longitude in [0,2*pi] 
%         Y  	 sin(latitude) in [-1,1] 
% Uses MATLAB function cart2sph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%A
[X,Y,R] = cart2sph(XC,YC,ZC);
%ditch R and return only lon and sin(lat) coordinates
X = mod(2.0*pi+X,2.0*pi);  %make the interval [0,2pi]
Y = sin(Y);                % \mu = \sin(lat)
%end function C2Sscalar
