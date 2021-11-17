function ratios = aditest_HHD(numeval)
% function ratios = aditest_HHD(numeval)
% Gets AD ratios for ADiGator/ADiMat/INTLAB/MAD tools
% input: numeval - number of evaluations (default 100)
% output: ratios - AD ratios
if nargin < 1
  numeval = 100;
end

Prob = MinpackHHD_Prob();
n = length(Prob.x_0);

setup.numvar = n;
setup.function = 'MinpackHHD_F';
setup.auxdata = Prob;

adimatFlag = exist('admColorSeed','file');
madFlag    = exist('fmad','file');
intlabFlag = exist('gradientinit','file');

tic
adifuncs = adigatorGenFiles4Fsolve(setup);
gentime = toc;

yy.f = Prob.x_0;
yy.dx = ones(n,1);
f1 = MinpackHHD_F_ADiGatorJac(yy,Prob);
Jpat = sparse(f1.dx_location(:,1),f1.dx_location(:,2),ones(size(f1.dx)),f1.dx_size(1),n);

if adimatFlag
  [S,c] = admColorSeed(Jpat);
  adimat_derivclass('scalar_directderivs');
else
  [c,S] = adigatorColor(Jpat);
end

for k = 0:numeval
  if k < 2
    time.adigator = 0;
    time.intlab   = 0;
    time.mad      = 0;
    time.func     = 0;
    time.hand     = 0;
    time.adimat   = 0;
  end
  x = rand(n,1);
  
  %S = sparse(eye(n));
  
  if intlabFlag
    tic
    gx = gradientinit(x);
    gf = MinpackHHD_F(gx,Prob);
    time.intlab = time.intlab+toc;
  end
  
  tic
  [F,J2] = MinpackHHD_F_Jac(x,Prob);
  time.adigator = time.adigator + toc;
  
  if madFlag
    tic
    mx = fmad(x,S);
    mf = MinpackHHD_F(mx,Prob);
    time.mad = time.mad+toc;
  end
  
  tic
  f = MinpackHHD_F(x,Prob);
  time.func = time.func+toc;
  
  tic
  [F,J] = MinpackHHD_FJ(x,Prob);
  time.hand = time.hand+toc;
  
  if adimatFlag
    tic
    J3 = zeros(size(Jpat,1),size(S,2));
    for i = 1:size(S,2)
      J3(:,i) = g_MinpackHHD_F(S(:,i),x,Prob);
    end
    time.adimat = time.adimat+toc;
  end
end
time.adigator = time.adigator./numeval;
time.intlab   = time.intlab./numeval;
time.mad      = time.mad./numeval;
time.func     = time.func./numeval;
time.hand     = time.hand./numeval;
time.adimat     = time.adimat./numeval;

ratios.adigator = time.adigator./time.func;
ratios.intlab   = time.intlab./time.func;
ratios.mad      = time.mad./time.func;
ratios.hand     = time.hand./time.func;
ratios.adimat     = time.adimat./time.func;
ratios.sfd = size(S,2)+1;