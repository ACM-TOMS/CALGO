% Solves the CPF problem using fsolve, supplies derivatives with ADiGator,
% also solves without supplying derivatives.
clc
Prob = MinpackCPF_Prob();
n = length(Prob.x_0);


if ~exist('MinpackCPF_F_ADiGatorJac.m','file')
  adisetup.numvar = n;
  adisetup.function = 'MinpackCPF_F';
  adisetup.auxdata = Prob;
  adigatorGenFiles4Fsolve(adisetup);
end

yy.f = Prob.x_0;
yy.dx = ones(n,1);
f1 = MinpackCPF_F_ADiGatorJac(yy,Prob);
Jpat = sparse(f1.dx_location(:,1),f1.dx_location(:,2),ones(size(f1.dx)),f1.dx_size(1),n);

solveflag = exist('fsolve','file');
if solveflag
  opts = optimset('Jacobian','on','JacobPattern',Jpat);
  problem.objective = @(x)MinpackCPF_F_Jac(x,Prob);
  problem.x0 = Prob.x_0;
  problem.lb = Prob.x_L;
  problem.ub = [];
  problem.solver = 'fsolve';
  problem.options = opts;
  
  tic
  [x,fval,exitflag,output,jacobian] = fsolve(problem);
  adisolve = toc
  
  opts = optimset('JacobPattern',Jpat);
  problem.objective = @(x)MinpackCPF_F(x,Prob);
  problem.x0 = Prob.x_0;
  problem.lb = Prob.x_L;
  problem.ub = [];
  problem.solver = 'fsolve';
  problem.options = opts;
  
  tic
  [x2,fval2,exitflag2,output2,jacobian2] = fsolve(problem);
  fdsolve = toc
else
  fprintf('Problem requires optimization toolbox\n');
end