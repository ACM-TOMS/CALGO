function ratios = burgersRatios_loop(nn,numeval)
% function ratios = burgersRatios_loop(nn,numeval)
% This file computes Burgers ODE Jacobian numerous times using ADiGator,
% ADiMat, INTLAB, and MAD using the ORIGINAL burgers ode function (with
% loops). User must have these packages installed in order to run this
% file.
%
% inputs: nn - problem dimension, can be vector of dimensions (dafault
%              2.^5)
%         numeval - number of Jac evaluations to average over (default 100)
%         
% output: ratios - Jac/fun CPU ratios for each tool.
if ~exist('numeval','var')
  numeval = 100;
end
if ~exist('nn','var')
  nn = 2.^5;
end

adimatFlag = exist('admColorSeed','file');
madFlag    = exist('fmad','file');
intlabFlag = exist('gradientinit','file');

if adimatFlag
  adimat_derivclass('scalar_directderivs');
end

count = 0;
for n = nn
  count = count+1;
N = n/2;
h = 1/(N+1);
xinit = h*(1:N);
% u(x,0) at grid points
ainit = sin(2*pi*xinit) + 0.5*sin(pi*xinit);

y0 = [ainit xinit];
y0 = y0(:);
f0 = burgersfun(0,y0,N);
ay = adigatorCreateDerivInput(size(y0),'y');
output = adigatorGenJacFile('burgersfun',{1,ay,N},adigatorOptions('overwrite',1));
Jpat = output.JacobianStructure;

ay = [];

if adimatFlag
  [S,c] = admColorSeed(Jpat);
else
  [c,S] = adigatorColor(Jpat);
end

ti.n = n;
for i = 0:numeval
  % Add some noise
  y = y0.*(1+.1*rand(n,1)-.1*rand(n,1));
  if rem(i,10) == 0
    fprintf('noloop - n = %1.0d -- i = %1.0f\n',n,i)
  end
  if i < 2
    ti.adigator = 0;
    ti.intlab = 0;
    ti.adimat = 0;
    ti.mad    = 0;
    ti.fun    = 0;
  end
  
  tic
  f = burgersfun(1,y,N);
  ti.fun = ti.fun+toc;
  
  if intlabFlag
    tic;
    gy = gradientinit(y);
    gf = burgersfun(1,gy,N);
    ti.intlab = ti.intlab+toc;
  end
  
  tic;
  ay.f = y;
  ay.dy = ones(n,1);
  f1 = burgersfun_ADiGatorJac(1,ay,N);
  ti.adigator = ti.adigator+toc;
  
  if madFlag
    tic;
    my = fmad(y,S);
    mf = burgersfun(1,my,N);
    ti.mad = ti.mad + toc;
  end
  
  if adimatFlag
    tic;
    JJ = zeros(size(Jpat,1),size(S,2));
    for j = 1:size(S,2)
      g_y = S(:,j);
      JJ(:,j) = g_burgersfun(1,g_y,y,N);
    end
    ti.adimat = ti.adimat + toc;
  end
end
ti.adigator = ti.adigator/numeval;
ti.adimat   = ti.adimat/numeval;
ti.mad      = ti.mad/numeval;
ti.intlab   = ti.intlab/numeval;
ti.fun      = ti.fun/numeval;
if count == 1
  time = ti;
else
  time(count) = ti;
end
end
ratios = struct('adimat',[],'adigator',[],'mad',[],'intlab',[]);
for i = 1:length(time)
  ratios.adigator(i) = time(i).adigator/time(i).fun;
  ratios.adimat(i) = time(i).adimat/time(i).fun;
  ratios.mad(i) = time(i).mad/time(i).fun;
  ratios.intlab(i) = time(i).intlab/time(i).fun;
end
