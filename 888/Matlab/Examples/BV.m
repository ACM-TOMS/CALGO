 
function BV              
%  Barotropic Vorticity Equation simulation on a Sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Program BV using a Spectral-SLT method
%
%   This program uses a pseudo spectral method 
%   to integrate the initial value problem (with periodic initial data)
%   of the shallow water equations
%    (1)         d/dt (eta)              =  0.0
%   Variables:  
%      eta	- potential vorticity, (= xi + f)
%      f	- Coriolis term ( = 2 \Omega \sin \theta)
%
%   The inversion relation with velocity is
%                L(psi) = eta -f 
%                v = k x grad(psi)
%   Above eta = k*(curl v), delta = div(v) = 0.0
%
%   The program uses the leapfrog (centered in time) integration routine
%   A two time level, semi-Lagrangian transport scheme is used for the advective terms.

%  Algorithm:
%   Particle tracking uses Hortal (1999) varient
%     dX/dt =  V :  X_D^n = X_A^(n+1) - dt * V_M^(n+0.5)
%              where V_M^(n+0.5) = 0.5*( (2V_D^n - V_D^(n-1)) + V_A^n )
%    The underscore A,D and M refer to arrival, departure and midpoints of the 
%    particle trajectories.
%     (D1)       eta_A^(n+1) = eta_D^n - (eta*delta)_M^(n+0.5)
%
%   Input:  
%   User enters the grid parameter for an NxN mesh and the number of
%   sample points of the initial data and the time T at which the
%   the solution is to be viewed. Solution at time T and initial data
%   are plotted together.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Remarks: 
%     Storage:  The spectral_field is not overloaded for operations of 
%       addition and multiplication, so it is necessary to do these as 
%       matrix operations at the model level.  This also makes it clear
%       what calcuations are going on and all are in grid point space.
%     Spectral calculations:  All the operations that require transforms and
%       manipulations in spectral space are done with methods of the 
%       spectral_field class.  Thus it is necessary to set values.  
%     Time history is not held in copies of the spectral_field, since there
%       are no ODE integrations taking place with the spectral coefficients.
%     Transport calculations:  The slt_grid class provides the interpolation
%       and particle tracking on an extended grid that is hidden from the
%       model level.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Begin method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
string = 'Begin'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control initialization  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nj = input('enter the value of N   ')
%ndays = input('enter the number of days to integrate ')
nj=16;     % number of latitudes
ndays = 3;  % Number of days to integrate the model
T = 86400*ndays   % time interval of the solution in seconds
%nsteps=5;  % number of steps to take
%dt = T/nsteps;   % model time step
dt = 3600   % model time step (20 minutes)
nsteps = T/dt
nout=3;    % frequency of model output in number of steps
t = 0:dt:T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Grid initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Domain Spherical grid (ni,nj)
G = gauss_grid('T10',nj);
aradius = 6.37122e6; % units - m (earth radius)
G = set(G,'radius',aradius);  % set the Earths radius for the grid operators
ni = get(G,'ni');
xg = get(G,'xg');
yg = get(G,'yg');
%setup the Slt extended grid
A = slt_grid(xg,yg);
X = get(A,'X');
Y = get(A,'Y');
Y = asin(Y);  % adjust to longitude for initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set Initial Fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Define the spectral_fields on the spectral_grid G
%Prognostic:
ETA= spectral_field('Potential Vorticity',G);
%Diagnostic
PSI= spectral_field('Streamfunction',G);
XI= spectral_field('Vorticity',G);
ff= spectral_field('Scratch',G);
DELTA= spectral_field('Divergence',G);
Z= zeros(ni,nj);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 2 initializiation: Steady zonal flow
%alpha = 0.0;
%alpha = pi/4;
%[Ucos,Vcos,Xi,gp,gp1,Fcor] = init_case2(alpha,ni,nj,X,Y); % gp is PSI, gp1 is PHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case 6 initialization: Rossby Wave (hmm)
[Ucos,Vcos,Xi,gp,gp1,Fcor] = init_case6(ni,nj,X,Y); % gp is PSI, gp1 is PHI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PSI=set(PSI,'gp',gp);
FCOR = set(XI,'gp',Fcor);
%XI=set(XI,'gp',Xi);
%gp = get(XI,'gp');
%PSI= del2inv(XI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve Laplaces equation for vorticity
XI = del2(PSI);
DELTA=set(DELTA,'gp',Z);
%Velocities (by inverting the vorticity and divergence relationship)
% But U = u*cos(lat), V = v*cos(lat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U,V] = UVinv(XI,DELTA);
Un = get(U,'gp');
Vn = get(V,'gp');
% Adjust for potential vorticity (need Coriolis term)
gp1 = get(XI,'gp');
gp = gp1 + Fcor;  % eta = xi + fcor
ETA=set(XI,'gp',gp);
%%%%%%%%%%%debug -
plot(U,' U ');
pause
ff = set(ff,'gp',Ucos);
plot(ff,' Ucos ');
pause
Xi=Ucos - Un ;
ff = set(ff,'gp',Xi);
plot(ff,'error Un ');
pause
plot(V,' V ');
pause
ff = set(ff,'gp',Vcos);
plot(ff,' Vcos ');
pause
Xi=Vcos - Vn;
ff = set(ff,'gp',Xi);
plot(ff,'Error Vcos ');
pause
Unm1=Un;
Vnm1=Vn;    % - debug
% end initialize method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(FCOR,'FCOR - Coriolis forcing term');
pause;
plot(ETA,'Initial ETA - Potential Vorticity');
pause;
plot(PSI,'Initial PSI - Stream Function');
pause;
plot(U,'Initial U - velocity');
pause;
plot(V,'Initial V - velocity');
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the old time level values needed for extrapolations, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Unm1 = get(U,'gp');
Vnm1 = get(V,'gp');
R1nm1 = Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Begin time stepping  -- Run method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
string = 'Run'
for nstep = 1:nsteps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute the SLT departure points from known velocities at time level n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Un = get(U,'gp');
Vn = get(V,'gp');
% extapolate the velocity to the half time level 
%   (and on unit sphere for the slt)
Uhalf = ((1.5)*Un - (0.5)*Unm1)/aradius; 
Vhalf = ((1.5)*Vn - (0.5)*Vnm1)/aradius;
[XD,YD,XM,YM] = slt_particle(A,Uhalf,Vhalf,dt);
%XD  %debug
%YD
%pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Explicit terms: Calculate right hand sides of equations 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Interpolation to the departure point of level n values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Etan = get(ETA,'gp');
  EtaD   = slt_interp(A,Etan,XD,YD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Calculate prognostic potential energy tendencies
%  KEM = 0.5*(Un.^2 + Vn.^2 );
%  ff = set(ff,'gp',KEM);
%   Extrapolate right hand sides to the half time level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  R1n = Z;    % eta right hand side
  R1M = 1.5*R1n - 0.5*R1nm1;
%   Interpolate right hand sides to midpoint XM
  R1M = slt_interp(A,R1M,XM,YM);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the explicit tendency term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  EtaA = EtaD + dt* R1M;                              % (Eqn D1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Prognostics are at the new time level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ETA = set(ETA,'gp',EtaA);
% Store time level n temporaries to nm1
  Unm1 = Un;
  Vnm1 = Vn;
  R1nm1 = R1n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update diagnostic relations:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Fcor = get(FCOR,'gp');
  gp = EtaA - Fcor;
  XI = set(XI,'gp',gp);
%Calculate velocities by inverting the UV relation
  [U,V] = UVinv(XI,DELTA);
% Solve for stream function 
  PSI = del2inv(XI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot every nout timesteps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (mod(nstep,nout) == 0 ) 
  nstep
  day = nstep*dt/86400
  plot(ETA,'ETA - potential vorticity');
  pause;
%  plot(XI); 
%  title('XI - vorticity');
%  pause;
%  plot(PSI); 
%  title('PSI - stream function');
%  pause;
%  plot(U); 
%  title('U - velocity');
%  pause;
%  plot(V); 
%  title('V - velocity');
%  pause;
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end   % time stepping loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% end run method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the final solution  - Finalize method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
string = 'Finalize'
day = nstep*dt/86400
plot(ETA,'Final ETA - vorticity');
pause;
plot(PSI,'Final PSI - stream function');
pause;
plot(U,'Final U - velocity');
pause;
plot(V,'Final V - velocity');
pause;
string = 'Done'
% end finalize method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end BV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

