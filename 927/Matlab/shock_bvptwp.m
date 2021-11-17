function shock_bvptwp(solver)
%shock_bvptwp  The solution has a shock layer near x = 0
%   This is an example used in U. Ascher, R. Mattheij, and R. Russell,
%   Numerical Solution of Boundary Value Problems for Ordinary Differential
%   Equations, SIAM, Philadelphia, PA, 1995,  to illustrate the mesh
%   selection strategy of COLSYS. 
%
%   For 0 < e << 1, the solution of 
%
%       e*y'' + x*y' = -e*pi^2*cos(pi*x) - pi*x*sin(pi*x)
%
%   on the interval [-1,1] with boundary conditions y(-1) = -2 and y(1) = 0
%   has a rapid transition layer at x = 0.
%
%   For this problem,
%   analytical partial derivatives are easy to derive and the solver benefits
%   from using them.  
%
%   By default, this example uses the 'twpbvpc_l' solver. Use syntax 
%   SHOCK_BVPTWP(solver) to solve this problem with the another solver
%     available solvers are: 'twpbvp_m', 'twpbvpc_m', 'twpbvp_l', 'twpbvpc_l', 
%                            'acdc', 'acdcc'
%
%   See also bvptwp, bvptwpset, bvptwpget, bvpinit, function_handle.
%   THIS MFILE IS ADAPDET FORM THE SHOCKBVP OF
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2007 The MathWorks, Inc.
%   $Revision: 1.10.4.3 $  $Date: 2007/05/23 18:53:57 $
%
%
%  
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universit\`a di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%

 if nargin < 1
   solver = 'twpbvpc_l';
 end

  
% The differential equations written as a first order system and the
% boundary conditions are coded in shockODE and shockBC, respectively. Their
% partial derivatives are coded in shockJac and shockBCJac and passed to the
% solver via the options. The option 'Vectorized' instructs the solver that
% the differential equation function has been vectorized, i.e.
% shockODE([x1 x2 ...],[y1 y2 ...]) returns [shockODE(x1,y1) shockODE(x2,y2) ...]. 
% Such coding improves the solver performance. 

options = bvptwpset('FJacobian',@shockJac,'BCJacobian',@shockBCJac,'Solver',solver,'reltol',[1e-4;0],'Linear','on');

% A guess for the initial mesh and the solution
%sol = bvpinit([-1 -0.5 0 0.5 1],[1 0]);
solinit =bvpinit(linspace(-1,1,9),zeros(1,2)); 
e = 1e-6;   
if  strcmp(solver,'acdc')|| strcmp(solver,'acdcc')
     options=bvptwpset(options,'Lambdamin',e);
end

sol = bvptwp(@shockODE,@shockBC,solinit,options);

% The final solution 
figure;
plot(sol.x,sol.y(1,:));
axis([-1 1 -2.2 2.2]);
title(['There is a shock at x = 0 when \epsilon =' sprintf('%.e',e) '.']); 
xlabel('x');
ylabel('solution y');

  % -----------------------------------------------------------------------
  % Nested functions -- e is shared with the outer function.
  %
  
  function dydx = shockODE(x,y,ExtraArgs)
  %SHOCKODE  Evaluate the ODE function (vectorized)
  if nargin==3
     e=ExtraArgs;
   end
    pix = pi*x;
    dydx = [                 y(2,:)
             (-x.*y(2,:) - e*pi^2*cos(pix) - pix.*sin(pix))/e ];   
  end
  % -----------------------------------------------------------------------

  function res = shockBC(ya,yb,ExtraArgs)
  %SHOCKBC  Evaluate the residual in the boundary conditions
  if nargin==3
     e=ExtraArgs;
   end
    res = [ ya(1)+2
            yb(1)  ];
  end
  % -----------------------------------------------------------------------

  function jac = shockJac(x,y,ExtraArgs)
  %SHOCKJAC  Evaluate the Jacobian of the ODE function
  %  x and y are required arguments.
        if nargin==3
            e=ExtraArgs;
        end
    jac = [ 0   1
            0 -x/e ];
  end
  % -----------------------------------------------------------------------

  function [dBCdya,dBCdyb] = shockBCJac(ya,yb,ExtraArgs)
  %SHOCKBCJAC  Evaluate the partial derivatives of the boundary conditions
  %  ya and yb are required arguments.
        if nargin==3
            e=ExtraArgs;
        end
    dBCdya = [ 1 0
               0 0 ];

    dBCdyb = [ 0 0
               1 0 ];
  end
end  % shock_bvptwp  

