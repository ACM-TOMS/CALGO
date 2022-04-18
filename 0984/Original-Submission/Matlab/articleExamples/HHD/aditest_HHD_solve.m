% Solves the HHD problem using fsolve, supplies derivatives with ADiGator,
% also solves without supplying derivatives.
clc
Prob = MinpackHHD_Prob();
n = length(Prob.x_0);


if ~exist('MinpackHHD_F_ADiGatorJac.m','file')
  adisetup.numvar = n;
  adisetup.function = 'MinpackHHD_F';
  adisetup.auxdata = Prob;
  adigatorGenFiles4Fsolve(adisetup);
end

yy.f = Prob.x_0;
yy.dx = ones(n,1);
f1 = MinpackHHD_F_ADiGatorJac(yy,Prob);
Jpat = sparse(f1.dx_location(:,1),f1.dx_location(:,2),ones(size(f1.dx)),f1.dx_size(1),n);

solveflag = exist('fsolve','file');
if solveflag
funtol = 10^(-6);
xtol = 10^(-6);
opts = optimset('Jacobian','on','JacobPattern',Jpat,'TolFun',funtol,'TolX',xtol,'Display','Iter');
problem.objective = @(x)MinpackHHD_F_Jac(x,Prob);
problem.x0 = Prob.x_0;
problem.lb = Prob.x_L;
problem.ub = Prob.x_U;
problem.solver = 'fsolve';
problem.options = opts;

tic
[x,fval,exitflag,output,jacobian] = fsolve(problem);
adisolve = toc

opts = optimset('JacobPattern',Jpat,'TolFun',funtol,'TolX',xtol,'Display','Iter');
problem.objective = @(x)MinpackHHD_F(x,Prob);
problem.x0 = Prob.x_0;
problem.lb = Prob.x_L;
problem.ub = Prob.x_U;
problem.solver = 'fsolve';
problem.options = opts;

tic
[x2,fval2,exitflag2,output2,jacobian2] = fsolve(problem);
fdsolve = toc
else
  fprintf('Problem requires optimization toolbox\n');
end