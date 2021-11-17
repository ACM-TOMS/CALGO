% 
%  Test program for the slt class.
%  The semi-lagrangian transport class defines and extended grid
%  on which to perform interpolation of fields, does particle tracking
%  along trajectories determined from a velocity input, and interplates
%  field values at the departure points of the trajectories.
%
%   This test program exercises the slt_grid, slt_particle and slt_interp
%   methods.  The private function slt_extend is utilized to fill halo regions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Uniform grid initialization
time = 0.0 % model time in days
dt=0.1;    % time step (also in days)
nx = 64;
ny = 32;
delx = 2*pi/nx; dely = 2.0/(ny+1);  % x is longitude, y is sin(latitude)
x = 0:delx:2*pi-delx;
y = -1.0+dely:dely:1.0-dely;
%  Set the slt grid
A = slt_grid(x,y);
%  Pull the lon and lats from the slt_grid.  Should be a get method.
Y = get(A,'Y');
X = get(A,'X');
Y = asin(Y);
%  Define a velocity field on grid interior as in the SW Advection Test Case 1
alpha = pi/3  % define rotation angle for test case 1
[U,V,F] = init_case1(alpha,nx,ny,X,Y);  % note conversion to lat from mu in Y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
surf(F);
 pause % debug
%  Compute departure points
[XD,YD,XM,YM] = slt_particle(A,U,V,dt);
%  Interpolate field to departure points
Fadv = slt_interp(A,F,XD,YD);
time = time + dt
surf(Fadv)
% Play it again.
while (time+dt<12.0)
 Fadv = slt_interp(A,Fadv,XD,YD);
 time = time + dt
 surf(Fadv)
 pause
end
 surf(Fadv-F)
% end test




