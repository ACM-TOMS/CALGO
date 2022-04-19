function [z, F,probinfo] = setupMinimumClimb(numLGR,numIntervals,sol)
% This function sets up the minimum climb optimal control problem - builds
% initial guess/bounds and setups up collocation


h0 = 0;
hf = 20;
v0 = 26;
vf = 59;
fpa0 = 0;
fpaf = 0;

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
t0min = 0;
t0max = 0;
tfmin = 100 /tscale;
tfmax = 350 /tscale;

bounds.state   = [ [hmin ; hmax] [vmin; vmax] [fpamin; fpamax] ];
bounds.control = [umin ; umax];
bounds.tf      = [tfmin ; tfmax];
guess.state    = bounds.state;
guess.control  = [0 ; 0];
guess.tf       = tfmax;


% --------------------- Call Direct Setup ------------------------------- %
n = size(guess.state,2);
m = size(guess.control,2);
probinfo.nstates   = n;
probinfo.ncontrols = m;

variablecount   = 0;

% ---------------------- Get LGR Stuff ---------------------------------- %
meshPoints = linspace(-1,1,numIntervals+1).';
polyDegrees = numLGR*ones(numIntervals,1);
[tau,w,D] = lgrPS(meshPoints,polyDegrees);
nLGR = length(w);

probinfo.LGR.nLGR    = nLGR;
probinfo.LGR.Dmatrix = D;
probinfo.LGR.weights = w;
probinfo.LGR.tau     = tau;

% ------------- Get Reference Indices off of Decision Vector ------------ %
xInds = variablecount+(1:n*(nLGR));
xInds = reshape(xInds,nLGR,n);
variablecount = xInds(end);
probinfo.map.state = xInds;

uInds = variablecount+(1:m*nLGR);
uInds = reshape(uInds,nLGR,m);
variablecount = uInds(end);
probinfo.map.control = uInds;

xNInds = variablecount+(1:n);
variablecount = xNInds(end);
probinfo.map.stateNp1 = xNInds;

tfInd = variablecount+(1);
variablecount = tfInd(end);
probinfo.map.tf = tfInd;

% -- Create Dummy Variables and Get Defects/Path Constraint Locations --- %
defectInds  = 1:n*nLGR;
probinfo.map.defect = defectInds;

probinfo.nconstraints = n*nLGR;
probinfo.nvariables   = variablecount;

% ---------------- Set Bounds on Decision Variable ---------------------- %
zmin   = zeros(variablecount,1);
zmax   = zeros(variablecount,1);
zguess = zeros(variablecount,1);
for i = 1:n
  xiInds = xInds(:,i);
  zmin(xiInds) = bounds.state(1,i);
  zmax(xiInds) = bounds.state(2,i);
  if nargin == 3
    zguess(xiInds) = interp1(sol.tau,[sol.state(:,i);sol.stateNp1(i)],tau(1:end-1));
  else
    zguess(xiInds) = linspace(guess.state(1,i),guess.state(2,i),nLGR);
  end
end

for i = 1:m
  uiInds = uInds(:,i);
  zmin(uiInds) = bounds.control(1,i);
  zmax(uiInds) = bounds.control(2,i);
  if nargin == 3
    zguess(uiInds) = interp1(sol.tau(1:end-1),sol.control(:,i),tau(1:end-1),'linear','extrap');
  else
    zguess(uiInds) = linspace(guess.control(1,i),guess.control(2,i),nLGR);
  end
end

for i = 1:n
  xiInds = xNInds(i);
  zmin(xiInds) = bounds.state(1,i);
  zmax(xiInds) = bounds.state(2,i);
  if nargin == 3
    zguess(xiInds) = sol.stateNp1(i);
  else
    zguess(xiInds) = bounds.state(2,i);
  end
end

zmin(tfInd) = bounds.tf(1);
zmax(tfInd) = bounds.tf(2);
if nargin == 3
  zguess(tfInd) = sol.tf;
else
  zguess(tfInd) = guess.tf;
end

z.guess = zguess;
z.min   = zmin;
z.max   = zmax;
z.mul   = zeros(probinfo.nvariables,1);
z.state = z.mul;
% -------------------- Set Bounds on Constraints ------------------------ %
Fmin = zeros(probinfo.nconstraints,1);
Fmax = zeros(probinfo.nconstraints,1);

F.min   = Fmin;
F.max   = Fmax;
F.mul   = zeros(probinfo.nconstraints,1);
F.state = F.mul;

xInds = probinfo.map.state;
z.min(xInds(1,:)) = [h0 v0 fpa0];
z.max(xInds(1,:)) = [h0 v0 fpa0];
z.min(xInds(end,:)) = [hf vf fpaf];
z.max(xInds(end,:)) = [hf vf fpaf];

end

function [x,w,P]=lgrnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lgrnodes.m
%
% Computes the Legendre-Gauss-Radau nodes, weights and the LGR Vandermonde
% matrix. The LGR nodes are the zeros of P_N(x)+P_{N+1}(x).
%
% References on LGR nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
%   F. B. Hildebrand , "Introduction to Numerical Analysis," Section 8.11
%   Dover 1987
%
% Written by Greg von Winckel - 05/02/2004
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N1=N+1;

% Use Chebyshev-Gauss-Radau nodes as initial guess for LGR nodes
x=-cos(2*pi*(0:N)/(2*N+1))';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1+1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.

xold=2;

% Free abscissae
free=2:N1;

while max(abs(x-xold))>eps
  
  xold=x;
  
  P(1,:)=(-1).^(0:N1);
  
  P(free,1)=1;    P(free,2)=x(free);
  
  for k=2:N1
    P(free,k+1)=( (2*k-1)*x(free).*P(free,k)-(k-1)*P(free,k-1) )/k;
  end
  
  x(free)=xold(free)-((1-xold(free))/N1).*(P(free,N1)+P(free,N1+1))...
    ./(P(free,N1)-P(free,N1+1));
end

% The Legendre-Gauss-Radau Vandermonde
P=P(1:N1,1:N1);

% Compute the weights
w=zeros(N1,1);
w(1)=2/N1^2;
w(free)=(1-x(free))./(N1*P(free,N1)).^2;
end

function [tau,w,D] = lgrPS(s,N);

%--------------------------------------------------------%
% This function computes the points, weights, and        %
% differentiation matrix for use with a multiple-segment %
% Radau pseudospectral method.                           %
% INPUTS                                                 %
%  s:  a vector of length K+1 that contains mesh points  %
%      on the closed interval [-1,+1]. s must be         %
%      monotonically increasing on [-1,+1].              %
%  N:  a vector of length K that contains the degree of  %
%      the polynomial within each mesh interval          %
% OUTPUTS                                                %
%  tau: a vector of LGR points on [-1,+1]                %
%  w:   a vector of LGR weights                          %
%  D:   the LGR differentiation matrix                   %
%--------------------------------------------------------%

% s = [-1; s];
Ntotal = sum(N);
tau = zeros(Ntotal+1,1);
w   = zeros(Ntotal,1);
D   = sparse(Ntotal,Ntotal+1);
rowshift = 0;
colshift = 0;
for i=1:length(N)
  [xcurr,wcurr,Pcurr] = lgrnodes(N(i)-1);
  rows = rowshift+1:rowshift+N(i);
  cols = colshift+1:colshift+N(i)+1;
  scurr = (s(i+1)-s(i))*xcurr/2+(s(i+1)+s(i))/2;
  tau(rows) = scurr;
  w(rows) = (s(i+1)-s(i))*wcurr/2;
  scurr = [scurr; s(i+1)];
  Dcurr = collocD(scurr);
  Dcurr = Dcurr(1:end-1,:);
  D(rows,cols) = Dcurr;
  rowshift = rowshift+N(i);
  colshift = colshift+N(i);
end;
tau(end)=+1;
end

function D=collocD(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% collocD.m
%
% Computes the pseudospectral/collocation differentiation matrix for the
% arbitrary nodes stored in the vector x. Uses the lagrange polynomial
% formulation.
%
% Reference:
%
% Jean-Paul Berrut & Lloyd N. Trefethen, "Barycentric Lagrange Interpolation"
% http://web.comlab.ox.ac.uk/oucl/work/nick.trefethen/berrut.ps.gz
%
% Written by: Greg von Winckel       07/18/04
% Contact:    gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make x a column vector if it isn't already and order it
% and get the number of nodes
x=sort(x(:));                       N=length(x); N1=N+1; N2=N*N;

% Compute the barycentric weights
X=repmat(x,1,N);                    Xdiff=X-X'+eye(N);
W=repmat(1./prod(Xdiff,2),1,N);     D=W./(W'.*Xdiff);
D(1:N1:N2)=1-sum(D);                D=-D';
end