function [x,phi,psi]=cascade(wavelet,q)
%CASCADE  Compute Daubechies scaling function and wavelet at dyadic rationals
%
%[x,phi,psi]=cascade(D,q)
%  Input 
%    wavelet:   Wavelet genus - optional, default D=4.
%         N.B.: If D is a vector it will be interpreted as a vector of 
%             low pass filters (lp).
%    q:   Desired dyadic resolution, default q=7.
%
%   Output
%     x:   Abscissae: (0, 1/2^q, ... , (D-1)*2^q)
%     phi: Vector containing corresponding values of scaling function
%     psi: Vector containing corresponding values of wavelet
%
%    Dependencies
%    none
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Parameter check%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
  q = 7;
  if nargin < 1
    wavelet = 4;
  end
end;    
len = length(wavelet);  %Determine the nature of the input argument #1:
if len > 1
  D  = len;             % Use actual filter length as wavelet genus
  lp = wavelet;         % First argument is a vector of low pass filters 
else
  D  = wavelet;         % First argument is the wavelet genus
  [lp,hp] = wfilters(['db' num2str(D/2)],'r');     % Low and high pass filters  
  lp=lp';hp=hp';
end;    
lp = lp(:);             % Force column vector

if isempty(lp)
  error('ERROR (cascade.m): invalid wavelet genus\n');
  fprintf('ERROR (cascade.m): invalid wavelet genus\n');
  return;
end  

if q < 0
  error('ERROR (cascade.m): resolution q must be non negative\n');
  fprintf('ERROR (cascade.m): resolution q must be non negative\n');
  return;
end  

if rem(D,2) > 0
  error('ERROR (cascade.m): Filter must have even length\n');
  fprintf('ERROR (cascade.m): Filter must have even length\n');
  return;
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shift = 1;              % Matlab hack to start vectors in 0

stride  = 2^q;          % Distance between rationals at which phi and psi are calculated
stride2 = stride/2;

L    = (D-1)*stride;    %  Support interval [0:L]
x    = ((0:L)/stride)'; % Represent points from 0 to D-1 (both included)
phi  = zeros(size(x));  % Corresponding phi - last point never touched


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Build matrixes A0 and A1 %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0 = zeros(D-1,D-1); % Dilation matrix for eigenvalues and "down" computation
A1 = zeros(D-1,D-1); % Dilation matrix for "up" computation
E = lp(1:2:D-1);     % Even numbered coefficients a0,a2,a4,...
O = lp(2:2:D);       % Odd numbered coefficients a1,a3,a5,...
Firstrow = 1;                    % Matrix A0
Lastrow  = D/2;
A0(Firstrow:Lastrow,1)     = E; 
for j=2:2:D-1
  Firstrow = Firstrow+1;
  Lastrow  = Lastrow+1;  
  A0(Firstrow:Lastrow,j)   = O; 
  A0(Firstrow:Lastrow,j+1) = E;   
end;
A0 = sqrt(2)*A0;
Firstrow = 1;                    % Matrix A1
Lastrow  = D/2;
for j=1:2:D-3
  A1(Firstrow:Lastrow,j)   = O; 
  A1(Firstrow:Lastrow,j+1) = E;   
  Firstrow = Firstrow+1;
  Lastrow = Lastrow+1;  
end;
A1(Firstrow:Lastrow,D-1)   = O; 
A1 = sqrt(2)*A1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Find initial phi using eigenvalue method%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,Dia]=eig(A0);      % Get eigen-vectors and -values  
[maxeig,index]=max(diag(Dia));  % Find largest eigenvalue (=1), eigen values of A0 are of the form 2^-m
v = V(:,index);                 % Pick corresponding eigen vector i.e. we are finding \Phi_{0}
v = v/sum(v);                   % Normalize s.t. sum(v) = 1, (sum(\phi_{k})=1)
for k=0:D-2 
  phi(k*stride + shift) = v(k+shift);   % Initialize phi at integers
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compute phi at successively finer levels%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if q > 0
  %%%%%%% Compute phi at half integers%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  start = stride2;                      
  phi(start+shift:stride:L+shift) = A1*phi(0+shift:stride:L-1+shift);
  %%%%%%% Do succesive levels: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  while start > 1
    step  = start;           % Distance between subvectors
    start = start/2;         % Leftmost starting point (1/2, 1/4, 1/8,...)

    for k=start:step:stride2
      k2 = 2*k;              %Starting point for source vector
      ks = k + stride2;      %Starting point for high bit (A1) target vector

      phi(k+shift :stride:L+shift) = A0*phi(k2+shift:stride:L+shift);
      phi(ks+shift:stride:L+shift) = A1*phi(k2+shift:stride:L+shift);
    end;   
  end;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compute psi from phi%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = zeros(size(phi));
   srt=sqrt(2);
if q > 0 
  L2   = L/2;
  phi2 = phi(1:2:L);        % Phi at even numerators 

  first = 1;
 
  for k=1:D
 psi(first:first+L2-1) = psi(first:first+L2-1) + srt*hp(k)*phi2;
    first = first + stride2;
  end;
else                        % q = 0  (a special case)
  for m=0:D-1
    for k=0:D-1
      index = 2*m-k;
      if (index >= 0) & (index <= D-1) 
     psi(m+shift) = psi(m+shift) + srt*hp(k+shift)*phi(index+shift);
      end;	
    end;  
  end;
end;




