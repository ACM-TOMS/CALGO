function [x,fval,output] = aditest_GL2_solve(nx)
% function [x,output] = aditest_GL2_solve(nx)
% This function solves the GL2 minimization problem using ADiGator and
% fminunc
% inputs:
% nx - scalar that defines problem dimension (default is 8)
% outputs:
% x      - optimal value of decision vector
% fval   - objective at optimal solution
% output - output of fminunc
%
% solves in approx 4sec with default value of 8, will increase as nx is
% increased.

if nargin == 0
  nx = 8;
end
n = 4*nx^2;
Prob = MinpackGL2_Prob(n);

setup.numvar = n;
setup.auxdata = Prob;
setup.objective = 'MinpackGL2_F';
setup.order = 2;

% Generate ADiGator Files

funcs = adigatorGenFiles4Fminunc(setup);
xx.f = Prob.x_0;
xx.dx = ones(length(xx.f),1);
ff = MinpackGL2_F_ADiGatorHes(xx,Prob);
Hpat = sparse(ff.dxdx_location(:,1),ff.dxdx_location(:,2), ones(size(ff.dxdx)),n,n);

solveflag = exist('fminunc','file');
if solveflag
  options = optimoptions(@fminunc,'GradObj','on','Hessian','on',...
    'Algorithm','trust-region','HessPattern',Hpat,'TolFun',sqrt(eps),'Display','Iter');
  problem.objective = funcs.hessian;
  problem.x0 = Prob.x_0;
  problem.solver = 'fminunc';
  problem.options = options;
  [x,fval,~,output] = fminunc(problem);
else
  fprintf('Problem requires optimization toolbox\n');
end
end