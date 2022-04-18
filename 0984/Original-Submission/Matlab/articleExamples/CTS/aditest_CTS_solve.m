% Solves the CTS problem using lsqnonlin, supplies derivatives with
% ADiGator, also solves without supplying derivatives.

clc
Prob = MinpackCTS_Prob();
n = length(Prob.x_0);

if ~exist('MinpackCTS_F_ADiGatorJac.m','file')
  adisetup.numvar = n;
  adisetup.function = 'MinpackCTS_F';
  adisetup.auxdata = Prob;
  adigatorGenFiles4Fsolve(adisetup);
end

yy.f = Prob.x_0;
yy.dx = ones(n,1);
f1 = MinpackCTS_F_ADiGatorJac(yy,Prob);
Jpat = sparse(f1.dx_location(:,1),f1.dx_location(:,2),ones(size(f1.dx)),f1.dx_size(1),n);


solveflag = exist('lsqnonlin','file');
if solveflag
  funtol = 10^(-6);
  xtol = 10^(-6);
  opts = optimset('Jacobian','on','JacobPattern',Jpat,'TolFun',funtol,'TolX',xtol,'Display','Iter');
  
  problem.objective = @(x)MinpackCTS_F_Jac(x,Prob);
  problem.x0 = Prob.x_0;
  problem.lb = [];
  problem.ub = [];
  problem.solver = 'lsqnonlin';
  problem.options = opts;
  
  tic
  [x,resnorm,residual,exitflag,output,lambda,jacobian] = lsqnonlin(problem);
  adisolve = toc
  
  opts = optimset('JacobPattern',Jpat,'TolFun',funtol,'TolX',xtol,'Display','Iter');
  problem.objective = @(x)MinpackCTS_F(x,Prob);
  problem.x0 = Prob.x_0;
  problem.lb = [];
  problem.ub = [];
  problem.solver = 'lsqnonlin';
  problem.options = opts;
  
  tic
  [x2,resnorm2,residual2,exitflag2,output2,lambda2,jacobian2] = lsqnonlin(problem);
  fdsolve = toc
else
  fprintf('Problem requires optimization toolbox\n');
end