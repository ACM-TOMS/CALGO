function [Prob,nuse]=MinpackGL2_Prob(varargin)
% MinpackGL2_Prob: Minpack 2D Ginzburg Landau (GL2) problem data.
%   
%    Computes start point and assigns physical constants etc for the 
%    2D Ginzburg Landau (GL2) problem from the MINPACK-2 collection.
%
% USE:
%       [Prob,nuse] = MinpackGL2_Prob(nx,ny,vornum)
%       [Prob,nuse] = MinpackGL2_Prob(nx,ny)
%       [Prob,nuse] = MinpackGL2_Prob(n)
% where
%   nx     : number of gridpoints in x-direction [10].
%   ny     : number of gridpoints in y-direction [nx].
%   vornum : number of vortices [8].
%   n      : problem size then nx=floor(sqrt(n/4)); ny=nx; vornum=8
%
% Returns
%    Prob.x_0         : starting point
%    Prob.user.nx     : no. of gridpoints in x-direction
%    Prob.user.ny     : no. of gridpoints in y-direction
%    Prob.user.vornum : no. of vortices
%    Prob.user.tkappa : Ginzburg-Landau constant
%    Prob.user.hx     : mesh spacing in x-direction
%    Prob.user.hy     : mesh spacing in x-direction
%    Prob.user.sqn    : sqn=nx*ny
%    nuse             : nuse=nx*ny
%
% See also MinpackGL2_F, MinpackGL2_FG, MinpackGL2_G, MinpackGL2_Hv,
% MinpackGL2_Sp, MinpackGL2_Fnovec, MinpackGL2_FGnovec

% AUTHOR: S.A.Forth & K. Lenton
% DATE: 7/01/09
% Copyright 2009-2009: S.A. Forth, Cranfield University 
% REVISIONS:
% DATE  WHO   WHAT
% 13/1/08 SAF removed input x
% 24/6/09 SAF single argument case for problem size and nuse returned.
%             header info standardised

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

% set any missing inputs
ny=[];vornum=[];
if nargin==0
    error('MADMinpack:MinpackGL2_Prob:nargin',...
        'MADMinpack:MinpackGL2_Prob:nargin at least one argument must be set')
elseif nargin==1
    n=varargin{1};
    nx=floor(sqrt(n/4));
    ny=nx;
else 
    nx=varargin{1};
    if nargin>=2
        ny=varargin{2};
    end
    if nargin>=3
        vornum=varargin{3};
    end
end
if isempty(nx)
    nx=10;
end
if isempty(ny)
    ny=nx;
end
if isempty(vornum)
    vornum=8;
end

% Initialization.
tkappa = 5;
hx = sqrt(vornum/2)*3/nx;
hy = sqrt(vornum/2)*3*sqrt(3)/ny;
sqn = nx*ny;
bave = 2*pi*vornum*tkappa/(sqn*hx*hy);
sqrtv = sqrt(vornum)*pi;

%    Initial Order Parameter.
ypt=(0:ny-1)*hy;
xpt=(0:nx-1)'*hx;
xp = 1 - (sin(sqrtv*xpt/(2*3)) * sin(sqrtv*ypt/(2*sqrt(3)*3))) .^2;

%    Initial Vector Potential. 
vpoty=bave*(0:nx-1)'*hx/tkappa*ones(1,ny);

x=[xp zeros(nx,ny*2) vpoty];

% Compute the standard starting point
Prob.x_0 = x(:);

% store required values in Prob.user
Prob.user.nx=nx;
Prob.user.ny=ny;
Prob.user.vornum=vornum;
Prob.user.tkappa=tkappa;
Prob.user.hx=hx;
Prob.user.hy=hy;
Prob.user.sqn=sqn;
nuse=4*nx*ny;