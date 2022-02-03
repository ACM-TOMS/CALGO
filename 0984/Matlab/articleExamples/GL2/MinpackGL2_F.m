function f=MinpackGL2_F(x,Prob)
% MinpackGL2_F: Computes objective for the 2D Ginzburg Landau (GL2) problem
%                  from the MINPACK-2 collection.
% USE:
%               f = MinpackGL2_F(x,Prob)
% where
% x    : vector of length Prob.user.n
% Prob : structure created by MinpackGL2_Prob with components
%    Prob.user.nx=nx;         % no. of gridpoints in x-direction
%    Prob.user.ny=ny;         % no. of gridpoints in y-direction
%    Prob.user.vornum=vornum; % no. of vortices
%    Prob.user.tkappa=tkappa; % Ginzburg-Landau constant
%    Prob.user.hx=hx;         % mesh spacing in x-direction
%    Prob.user.hy=hy;         % mesh spacing in x-direction
%    Prob.user.sqn=sqn;       % sqn=nx*ny
%
% OUTPUT:
%     f: scalar objective function
%
% See also MinpackGL2_Prob, MinpackGL2_FG, , MinpackGL2_Hv, , MinpackGL2_Sp

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 7/01/09
% Copyright 2009-2009: S.A. Forth, Cranfield University 
% REVISIONS:
% DATE  WHO   WHAT

% Original Fortran Header comments follow
%
%      subroutine dgl2fg(nx,ny,x,f,fgrad,task,w,vornum)
%      character*(*) task
%      integer nx, ny, vornum
%      double precision f
%      double precision x(4*nx*ny), fgrad(4*nx*ny), w(4*(nx+1)*(ny+1))
% **********
%
% Subroutine dgl2fg
%
% This subroutine computes the function and gradient of the
% Ginzburg-Landau (2-dimensional) superconductivity problem.
%
% The subroutine statement is
%
%   subroutine dgl2fg(nx,ny,x,f,fgrad,task,w,vornum)
%
% where
%
%   nx is an integer variable.
%     On entry nx is the number of grid points in the first
%        coordinate direction.
%     On exit nx is unchanged.
%
%   ny is an integer variable.
%     On entry ny is the number of grid points in the second
%        coordinate direction.
%     On exit ny is unchanged.
%
%   x is a double precision array of dimension 4*nx*ny.
%     On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
%        Otherwise x need not be specified.
%     On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
%        x is set according to task.
%
%   f is a double precision variable.
%     On entry f need not be specified.
%     On exit f is set to the function evaluated at x if task = 'F'
%        or 'FG'.
%
%   fgrad is a double precision array of dimension 4*nx*ny.
%     On entry fgrad need not be specified.
%     On exit fgrad contains the gradient evaluated at x if
%        task = 'G' or 'FG'.
%
%   task is a character variable.
%     On entry task specifies the action of the subroutine:
%
%        task               action
%        ----               ------
%         'F'     Evaluate the function at x.
%         'G'     Evaluate the gradient at x.
%         'FG'    Evaluate the function and the gradient at x.
%         'XS'    Set x to the standard starting point xs.
%
%     On exit task is unchanged.
%
%   w is a double precision work array of dimension 4*(nx+1)(ny+1).
%
%   vornum is an integer variable.
%     On entry vornum specifies the number of vortices.
%     On exit vornum is unchanged.
%
% MINPACK-2 Project. November 1993.
% Argonne National Laboratory and University of Minnesota.
% Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright.
%
% **********

% if nargin<2
%     error('MINPack:GL2:MinpackGL2_FG:nargin',...
%         'Arguments x and Prob must be supplied')
% end

n=length(x);
nx=Prob.user.nx;
ny=Prob.user.ny;

% % check arguments
if (n ~= 4*nx*ny)
    error ('MINPack:GL2:MinpackGL2_FG:n',...
        'MINPack:GL2:MinpackGL2_FG:length(x): x not of length 4*nx*ny')
end

% extract parameters from Prob structure
vornum =Prob.user.vornum;
tkappa =Prob.user.tkappa;
hx =Prob.user.hx;
hy =Prob.user.hy;
sqn =Prob.user.sqn;

% w=reshape(x,nx,ny,4);
% w(nx+1,ny+1,4)=0;
% x=w(:,:,1);
% y=w(:,:,2);
% vpotx=w(:,:,3);
% vpoty=w(:,:,4);
Z = reshape(x,nx*ny,4);


znx = zeros(nx,1); znyp1 = zeros(1,ny+1);
x     = [reshape(Z(:,1),nx,ny) znx; znyp1];
y     = [reshape(Z(:,2),nx,ny) znx; znyp1];
vpotx = [reshape(Z(:,3),nx,ny) znx; znyp1];
vpoty = [reshape(Z(:,4),nx,ny) znx; znyp1];

arg = 2*pi*vornum*(0:ny)/(ny);
x(nx+1,1:ny+1) = x(1,1:ny+1).*cos(arg) - y(1,1:ny+1).*sin(arg);
y(nx+1,1:ny+1) = x(1,1:ny+1).*sin(arg) + y(1,1:ny+1).*cos(arg);
vpotx(nx+1,1:ny+1) = vpotx(1,1:ny+1);
vpoty(nx+1,1:ny+1) = vpoty(1,1:ny+1) + 2*pi*vornum/(ny*hy);

% Top face for order parameter and vector potential.
x(1:nx+1,ny+1) = x(1:nx+1,1);
y(1:nx+1,ny+1) = y(1:nx+1,1);
vpotx(1:nx+1,ny+1) = vpotx(1:nx+1,1);
vpoty(1:nx+1,ny+1) = vpoty(1:nx+1,1);

%    Compute the Condensation Energy Density
delsq = x(1:nx,1:ny).^2 + y(1:nx,1:ny).^2;
fcond =  sum(sum(- delsq + (delsq.^2)./2))/sqn;

%    Compute the Kinetic Energy Density.
x1a = sum(sum((x(2:nx+1,1:ny) - x(1:nx,1:ny).*cos(hx*vpotx(1:nx,1:ny)) +...
    y(1:nx,1:ny).*sin(hx*vpotx(1:nx,1:ny))).^2));
x2a = sum(sum((y(2:nx+1,1:ny) - y(1:nx,1:ny).*cos(hx*vpotx(1:nx,1:ny)) -...
    x(1:nx,1:ny).*sin(hx*vpotx(1:nx,1:ny))).^2));

x1b = sum(sum((x(1:nx,2:ny+1) - x(1:nx,1:ny).*cos(hy*vpoty(1:nx,1:ny)) +...
    y(1:nx,1:ny).*sin(hy*vpoty(1:nx,1:ny))).^2));
x2b = sum(sum((y(1:nx,2:ny+1) - y(1:nx,1:ny).*cos(hy*vpoty(1:nx,1:ny)) -...
    x(1:nx,1:ny).*sin(hy*vpoty(1:nx,1:ny))).^2));

fkin = ((x1a+x2a)/(hx^2) + (x1b+x2b)/(hy^2))/sqn;

%    Compute the Magnetic Field Energy Density.
ffield = sum(sum(((vpoty(2:nx+1,1:ny)-vpoty(1:nx,1:ny))/hx -...
    (vpotx(1:nx,2:ny+1)-vpotx(1:nx,1:ny))/hy).^2))*(tkappa^2)/sqn;

% sum contributions
f = fcond + fkin + ffield;