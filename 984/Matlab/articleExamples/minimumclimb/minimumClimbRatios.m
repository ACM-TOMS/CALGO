function ratios = minimumClimbRatios(NN,numeval)
% function ratios = minimumClimbRatios(NN,numeval)
% This function will compute ratios for ADiGator/ADiMat/INTLAB/MAD tools
% for first and second derivatives of minimumClimb problem.
%
% inputs: NN - number of collocation points, can be a vector of
%              dimensions, will run example for each element of NN 
%              (default 2^5)
%         numeval - number of evaluations to average over
% outputs: ratios - ratios for each of the tested tools.
%
% test values of N = 2^4,2^5,...,2^12
if ~exist('NN','var')
  NN = 2.^5;
end
if ~exist('numeval','var')
  numeval = 100;
end

adimatFlag = exist('admColorSeed','file');
madFlag    = exist('fmad','file');
intlabFlag = exist('gradientinit','file');

% Problem bounds to sample from
feettometer = .3048;
hscale = 1000;
tscale = 200;
hmin =  0*feettometer       /hscale;
hmax =  69000*feettometer   /hscale;
vmin =  1*feettometer       /hscale*tscale;
vmax =  2000*feettometer    /hscale*tscale;
fpamin = -90/180*pi;
fpamax = -fpamin;
umin = -10;
umax =  10;

xMin = [hmin vmin fpamin umin];
xMax = [hmax vmax fpamax umax];

% Set different flags for different tools..
nx = 4;

auxdata = minimumClimbAuxdata();
% Differentiate Minimum Climb Dynamics
ax = adigatorCreateDerivInput([Inf nx],'X'); % 3 states 1 control
opts = adigatorOptions('overwrite',1,'comments',0);
axdot_X = adigator('minimumClimbDynamics',{ax,auxdata},'minimumClimbDynamics_X',opts);
ifx = axdot_X{1}.deriv.nzlocs(:,1);
jfx = axdot_X{1}.deriv.nzlocs(:,2);
ax2.f = ax; ax2.dX = adigatorCreateAuxInput([Inf nx]);
axdot_XX = adigator('minimumClimbDynamics_X',{ax2,auxdata},'minimumClimbDynamics_XX',opts);
ifxx = ifx(axdot_XX{1}.dX.deriv.nzlocs(:,1));
jfxx = jfx(axdot_XX{1}.dX.deriv.nzlocs(:,1));
kfxx = axdot_XX{1}.dX.deriv.nzlocs(:,2);



jti = [];
fti = [];
hti = [];
fun = struct('time',[],'N',[]);
jac = struct('adigator',[],'adigator_nonvec',[],'adimat',[],'intlab',[],'mad',[],'N',[]);
hes = struct('adigator',[],'adigator_nonvec',[],'adimat',[],'intlab',[],'mad',[],'N',[]);

for j = 1:length(NN)
  N = NN(j);
  
  % Colors
  c = repmat(1:4,[N 1]);
  c = c(:).';
  
  % Seed Matrix
  SI = speye(N*nx);
  I = ones(N,1);
  II = repmat({I},1,nx);
  S = blkdiag(II{:});
  
  % Adigator vectorized projections
  [IFX,JFX] = adigatorProjectVectLocs(N,ifx,jfx);
  [IFXX,JFXX,KFXX] = adigatorProjectVectLocs(N,ifxx,jfxx,kfxx);
  JFX = c(JFX);
  JKFXX = sub2ind([nx*N,nx],JFXX,c(KFXX));
  
  sFXX = sparse(IFXX,JKFXX,ones(size(IFXX)),N*3,N*nx*nx);
  [iFXX,jFXX] = find(sFXX);
  
  % GENERATE ADIGATOR NONVECTORIZED DERIV FILES
  nz = N*(nx);
  derivinfo.vodname = 'z';
  derivinfo.vodsize = [nz 1];
  derivinfo.nzlocs = [(1:nz).' (1:nz).'];
  ax = adigatorCreateDerivInput([N nx],derivinfo);
  axdot_z = adigator('minimumClimbDynamics',{ax,auxdata},'minimumClimbDynamics_z',opts);
  IFXnv = axdot_z{1}.deriv.nzlocs(:,1);
  JFXnv = c(axdot_z{1}.deriv.nzlocs(:,2));


  ax2.f = ax; ax2.dz = ones(nz,1);
  axdot_zz = adigator('minimumClimbDynamics_z',{ax2,auxdata},'minimumClimbDynamics_zz',opts);
  
  
  for i = 0:numeval
    if rem(i,10) == 0
      display(['N = ',num2str(N),' -- i = ',num2str(i)]);
    end
    if i < 2
      jti.adigator = 0;
      jti.adigator_nonvec = 0;
      jti.adimat   = 0;
      jti.intlab   = 0;
      jti.mad      = 0;
      
      hti.adigator = 0;
      hti.adigator_nonvec = 0;
      hti.adimat   = 0;
      hti.intlab   = 0;
      hti.mad      = 0;
      
      fti.time     = 0;
    end
    
    nz = size(S,1);
    
    % Generate Random Value of X
    x = zeros(N,nx);
    for ii = 1:nx
      x(:,ii) = rand(N,1).*(xMax(ii)-xMin(ii))+xMin(ii);
    end
    % Randome value of lambda
    lambda = rand(N*3,1);
    
    
    tic
    % FUNCTION
    out1 = minimumClimbDynamics(x,auxdata);
    fti.time = fti.time+toc;
    
    
    % ------------------ Jacobians -------------------------------------- %
    if adimatFlag
      adimat_derivclass('scalar_directderivs');
      tic
      % ADIMAT
      J_adimat = zeros(N*3,nx);
      for k = 1:nx
        Sk = S(:,k);
        g_x = reshape(Sk,N,nx);
        kout = g_minimumClimbDynamics(g_x,x,auxdata);
        J_adimat(:,k) = kout(:);
      end
      jti.adimat = jti.adimat + toc;
    end
    
    % ADIGATOR
    tic
    adix = struct('f',x,'dX',ones(N,nx));
    out2 = minimumClimbDynamics_X(adix,auxdata);
    J_adigator = sparse(IFX,JFX,out2.dX,N*3,nx);
    jti.adigator = jti.adigator+toc;
    
    % ADIGATOR NONVECTORIZED
    tic
    adix = struct('f',x,'dz',ones(N*nx,1));
    out3 = minimumClimbDynamics_z(adix,auxdata);
    J_adigator_nv = sparse(IFXnv,JFXnv,out3.dz,N*3,nx);
    jti.adigator_nonvec = jti.adigator_nonvec+toc;
    
    % INTLAB
    if intlabFlag
      tic
      gx = gradientinit(x(:));
      gx = reshape(gx,[N nx]);
      G_out = minimumClimbDynamics(gx,auxdata);
      J_intlab = G_out(:).dx;
      jti.intlab = jti.intlab + toc;
    end
    
    %MAD
    if madFlag
      tic
      gx = fmad(x(:),S);
      gx = reshape(gx,[N nx]);
      out = minimumClimbDynamics(gx,auxdata);
      J_mad = getinternalderivs(out);
      jti.mad = jti.mad + toc;
    end

    % ------------------- hessians -------------------------------------- %
    % ADIGATOR
    tic
    adix = struct('f',x,'dX',ones(N,nx));
    out2 = minimumClimbDynamics_XX(adix,auxdata);
    FXX = sparse(IFXX,JKFXX,out2.dXdX,N*3,N*nx*nx);
    hti.adigator = hti.adigator+toc;
    
    tic
    H_adigator = reshape(lambda.'*FXX,N*nx,nx);
    sumtime = toc;
    
    tic
    adix = struct('f',x,'dz',ones(N*nx,1));
    out3 = minimumClimbDynamics_zz(adix,auxdata);
    FXX  = sparse(iFXX,jFXX,out3.dzdz,N*3,N*nx*nx);
    hti.adigator_nonvec = hti.adigator_nonvec+toc;
    
    L = lambda.'*out3.f(:);

    % INTLAB
    if intlabFlag
      tic
      gx = hessianinit(x(:));
      gx = reshape(gx,[N nx]);
      out = minimumClimbDynamics(gx,auxdata);
      hti.intlab = hti.intlab + toc;
    end
    
    %MAD
    if madFlag
      tic
      gz1 = fmad(x(:),S);
      gz = fmad(gz1,S);
      gz = reshape(gz,[N nx]);
      out = minimumClimbDynamics(gz,auxdata);
      FXX_mad = getinternalderivs(getinternalderivs(out));
      hti.mad = hti.mad + toc;
    end

    % ADIMAT
    if adimatFlag
      adimat_adjoint('scalar');
      tic
      H_adimat = admHessian(@minimumClimbSum,S,x,lambda,auxdata,admOptions('i',1,'functionResults',{L},'hessianStrategy','t1rev'));
      hti.adimat = hti.adimat + toc;
    end
    
    % Add the sum time to all forward over forward mode methods
    hti.adigator = hti.adigator+sumtime;
    hti.adigator_nonvec = hti.adigator_nonvec+sumtime;
    if intlabFlag
      hti.intlab  = hti.intlab + sumtime;
    end
    if madFlag
      hti.mad     = hti.mad + sumtime;
    end
  end
  fti.time    = fti.time/numeval;
  fti.N       = N;
  fun(j) = fti;
  
  
  jti.adimat   = jti.adimat/numeval;
  jti.adigator = jti.adigator/numeval;
  jti.adigator_nonvec = jti.adigator_nonvec/numeval;
  jti.intlab  = jti.intlab/numeval;
  jti.mad     = jti.mad/numeval;
  jti.N       = N;
  
  jac(j) = jti;
  
  hti.adimat = hti.adimat/numeval;
  hti.adigator = hti.adigator/numeval;
  hti.adigator_nonvec = hti.adigator_nonvec/numeval;
  hti.intlab  = hti.intlab/numeval;
  hti.mad     = hti.mad/numeval;
  hti.N       = N;
  
  hes(j) = hti;
end

ratios = [];
for j = 1:length(NN)
  ratios.hes.adimat(j) = (hes(j).adimat/fun(j).time);
  ratios.hes.adigator(j) = (hes(j).adigator/fun(j).time);
  ratios.hes.adigator_nonvec(j) = (hes(j).adigator_nonvec/fun(j).time);
  ratios.hes.intlab(j) = (hes(j).intlab/fun(j).time);
  ratios.hes.mad(j) = (hes(j).mad/fun(j).time);
  
  ratios.jac.adimat(j) = (jac(j).adimat/fun(j).time);
  ratios.jac.adigator(j) = (jac(j).adigator/fun(j).time);
  ratios.jac.adigator_nonvec(j) = (jac(j).adigator_nonvec/fun(j).time);
  ratios.jac.intlab(j) = (jac(j).intlab/fun(j).time);
  ratios.jac.mad(j) = (jac(j).mad/fun(j).time);
end
