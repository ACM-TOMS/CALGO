function rats = aditest_GL2(numx,numeval)
% rats = aditest_GL2(numx,numeval)
% This function will compare results of ADiGator/ADiMat/INTLAB/MAD for GL2
% problem.
% inputs-
% numx :  scalar which defines problem dimension, if numel(numx) > 1,
%         results will be obtained for each element. (default is 8)
% numeval: number of evaluations to compare against. (default is 100)
% output-
% rats : structure containing ratios
%
% run time is approx 35sec with defaults, will increase as numx and numeval
% are increased.

if nargin < 1
  numx = 8;
end
if nargin < 2
  numeval = 100;
end
numx = numx(:).';

adimatFlag = exist('admColorSeed','file');
madFlag    = exist('fmad','file');
intlabFlag = exist('gradientinit','file');

rats = [];
% if adimatFlag
%   admTransform(@MinpackGL2_G,admOptions('i',1));
%   admTransform(@MinpackGL2_F,admOptions('i',1,'m','r'));
% end
count = 0;

for nx = numx
count = count+1;
n = 4*nx^2;
Prob = MinpackGL2_Prob(n);

setup.numvar = n;
setup.auxdata = Prob;
setup.objective = 'MinpackGL2_F';
setup.order = 2;

% Generate ADiGator Files
ax = adigatorCreateDerivInput([n 1],'x');
ax2.f = ax;
ax2.dx = ones(n,1);
t = tic;
adigator('MinpackGL2_F',{ax,Prob},'MinpackGL2_Fdx',adigatorOptions('overwrite',1));
gentimef1 = toc(t);

t = tic;
adigator('MinpackGL2_Fdx',{ax2,Prob},'MinpackGL2_Fdxdx',adigatorOptions('overwrite',1));
gentimef2 = toc(t);

t = tic;
adigator('MinpackGL2_G',{ax,Prob},'MinpackGL2_Gdx',adigatorOptions('overwrite',1));
gentimeg1 = toc(t);



t1 = cputime;
yy.f = rand(n,1);
yy.dx = ones(n,1);
ff1 = MinpackGL2_Fdxdx(yy,Prob);
Hpat = sparse(ff1.dxdx_location(:,1),ff1.dxdx_location(:,2),ones(size(ff1.dxdx)),n,n);

ff2 = MinpackGL2_Gdx(yy,Prob);
ax= [];

if adimatFlag
  [S,c] = admColorSeed(Hpat);
else
  [c,S] = adigatorColor(Hpat);
end

%disp(['2nd deriv gen time + solve time: ',num2str(gentime+solvetime)]);
jac = [];
grd = [];
hes = [];


for i = 0:numeval
  if i < 2
    grd.adigator = 0;
    grd.adimat = 0;
    grd.intlab = 0;
    grd.mad = 0;
    grd.hand = 0;
    
    jac.adigator = 0;
    jac.adimat = 0;
    jac.intlab = 0;
    jac.mad = 0;
    
    hes.adigator = 0;
    hes.adimat = 0;
    hes.intlab = 0;
    hes.mad = 0;
    
    fun = 0;
  end
  
  if rem(i,10) == 0
    fprintf('n = %1.0d -- i = %1.0f\n',n,i);
  end
  x = rand(n,1);
  
  t = tic;
  f = MinpackGL2_F(x,Prob);
  fun = fun+toc(t);
  
  
  % ------------------------- GRADIENT COMPUTATIONS --------------------- %
  t = tic;
  ax.f = x;
  ax.dx = ones(n,1);
  ff = MinpackGL2_Fdx(ax,Prob);
  grd.adigator = grd.adigator+toc(t);
  
  if intlabFlag
    t = tic;
    gx = gradientinit(x);
    gf = MinpackGL2_F(gx,Prob);
    grd.intlab = grd.intlab+toc(t);
  end
  
  if madFlag
    SF1 = sparse(eye(n));
    t = tic;
    gx = fmad(x,SF1);
    gf = MinpackGL2_F(gx,Prob);
    grd.mad = grd.mad+toc(t);
  end
  
  t = tic;
  G = MinpackGL2_G(x,Prob);
  grd.hand = grd.hand+toc(t);
  %adimat_derivclass('opt_derivclass');
  %adimat_derivclass('arrderivclass');
  if adimatFlag
    adimat_derivclass('scalar_directderivs');
    adimat_adjoint('scalar')
    pause(.1);
    rehash
    t = tic;
    a_y = 1;
    [a_f,fff] = a_MinpackGL2_F(x,Prob,a_y);
    grd.adimat = grd.adimat+toc(t);
  end

  % ----------------------- JACOBIAN COMPUTATIONS --------------------- %
  t = tic;
  ax.f = x;
  ax.dx = ones(n,1);
  ff = MinpackGL2_Gdx(ax,Prob);
  jac.adigator = jac.adigator+toc(t);
  
  if intlabFlag
    t = tic;
    gx = gradientinit(x);
    gf = MinpackGL2_G(gx,Prob);
    jac.intlab = jac.intlab+toc(t);
  end
  
  if madFlag
    t = tic;
    gx = fmad(x,S);
    gf = MinpackGL2_G(gx,Prob);
    jac.mad = jac.mad+toc(t);
    adimat_derivclass('scalar_directderivs');
  end
  
  if adimatFlag
    t = tic;
    JJ = zeros(n,size(S,2));
    for j = 1:size(S,2)
      JJ(:,j) = g_MinpackGL2_G(S(:,j),x,Prob);
    end
    jac.adimat = jac.adimat+toc(t);
  end
  
  t = tic;
  ax.f = x;
  ax.dx = ones(n,1);
  ff = MinpackGL2_Fdxdx(ax,Prob);
  hes.adigator = hes.adigator+toc(t);
  
  if intlabFlag && nx < 2^6
  t = tic;
  gx = hessianinit(x);
  gf = MinpackGL2_F(gx,Prob);
  hes.intlab = hes.intlab+toc(t);
  else
    hes.intlab = inf;
  end
  
  SF1 = sparse(eye(n));
  if madFlag && nx < 2^6
  t = tic;
  gx = fmad(x,SF1);
  gx = fmad(gx,S);
  gf = MinpackGL2_F(gx,Prob);
  hes.mad = hes.mad+toc(t);
  else
    hes.mad = inf;
  end
  
  if adimatFlag
    t = tic;
    H = admHessian(@MinpackGL2_F,S,x,Prob,admOptions('i',1,'functionResults',{ff.f},'hessianStrategy','t1rev'));
    hes.adimat = hes.adimat+toc(t);
  end
  
end
rats.grd.adigator(count) = grd.adigator/fun;
rats.grd.adimat(count) = grd.adimat/fun;
rats.grd.intlab(count) = grd.intlab/fun;
rats.grd.mad(count) = grd.mad/fun;
rats.grd.hand(count) = grd.hand/fun;
rats.grd.sfd(count) = n+1;

rats.jac.adigator(count) = jac.adigator/fun;
rats.jac.adimat(count) = jac.adimat/fun;
rats.jac.intlab(count) = jac.intlab/fun;
rats.jac.mad(count) = jac.mad/fun;
rats.jac.sfd(count) = rats.grd.hand(count)*(size(S,2)+1);

rats.hes.adigator(count) = hes.adigator/fun;
rats.hes.adimat(count) = hes.adimat/fun;
rats.hes.intlab(count) = hes.intlab/fun;
rats.hes.mad(count) = hes.mad/fun;
rats.hes.sfd(count) = 2*n;

rats.n(count) = n;
rats.colors(count) = size(S,2);
end