function [x,phi_d,psi_d]=cascade_der(wavelet,q,d)
fprintf('d=%d',d);
%function [x,phi_d,psi_d]=cascade_der(wavelet,q,d)
% General for calculating any order of derivative...............
%*********************************************************************
%CASCADE  Compute Daubechies scaling function and wavelet at dyadic rationals
%[x,phi_d,psi_d]=cascade1(wavelet,q,d)
%
%  Input 
%    wavelet:   Wavelet genus - optional, default D=4.
%         N.B.: if D is a vector it will be interpreted as a vector of 
%             low pass filters (lp).
%    q:   Desired dyadic resolution, default q=7.
%
%    d: d is the order of differentiation required
%
%
%   Output
%     x:   Abscissae: (0, 1/2^q, ... , (D-1)*2^q)
%     phi_d: Vector containing corresponding values of derivatives of scaling function.
%     psi_d: Vector containing corresponding values of derivatives of wavelet function.
%
%  See also daubfilt, low2hi
%Dependencies 
% none

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Parameter check%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
  q = 7;
  if nargin < 1
    wavelet = 4;
  end
end;    

len = length(wavelet);  %Determine the nature of argument #1:
if len > 1
  D  = len;             % Use actual filter length
  lp = wavelet;         % First argument is a vector of low pass filter coefficients 
else
  D  = wavelet;         % First argument is the wavelet genus
   [lp,hp] = wfilters(['db' num2str(D/2)],'r');     % Low and high pass filters    
  lp=lp';hp=hp';
  lp = sqrt(2)*lp;     % Low pass filter  
end;    
lp = lp(:);             % Force column vector

if isempty(lp)
  error('ERROR (cascade.m): invalid wavelet genus\n');
  fprintf('ERROR (cascade.m): invalid wavelet genus\n');
  x=[]; phi_d = []; psi_d=[];
  return;
end  

if q < 0
  error('ERROR (cascade.m): resolution q must be non negative\n');
  fprintf('ERROR (cascade.m): resolution q must be non negative\n');
  x=[]; phi_d = []; psi_d=[];
  return;
end  

if rem(D,2) > 0
  error('ERROR (cascade.m): Filter must have even length\n');
  fprintf('ERROR (cascade.m): Filter must have even length\n');
  x=[]; phi_d = []; psi_d=[];
  return;
end  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shift = 1;              % Matlab hack to start vectors in 0

stride  = 2^q;          % Distance between integers in phi and psi
stride2 = stride/2;

L    = (D-1)*stride;    %Support interval [0:L]
x    = ((0:L)/stride)'; %Represent points from 0 to D-1 (both included)
phi_d  = zeros(size(x));  %Corresponding phi - last point never touched


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Build matrixes A0 and A1%%%%%%%%%%%%%%%%%
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
  A0(Firstrow:Lastrow,j)   = O; % O is odd numbered low pass filter coefficients
  A0(Firstrow:Lastrow,j+1) = E;  % E is even numbered low pass filter coefficients
end;
A0 = A0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% note that we have multiplied our vector containig low pass filter coefficients by sqrt(2),
%%%% so we dont need to multiply the matrix A0 by sqrt(2).. Or we can do the other way round..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Firstrow = 1;                    % Matrix A1
Lastrow  = D/2;
for j=1:2:D-2
  A1(Firstrow:Lastrow,j)   = O; 
  A1(Firstrow:Lastrow,j+1) = E;   
  Firstrow = Firstrow+1;
  Lastrow = Lastrow+1;  
end;
A1(Firstrow:Lastrow,D-1)   = O; 
A1 =(2^(d))*A1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Find initial phi using eigenvalue method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,Dia]=eig(A0);      % Get eigen-vectors and -values  
if(d==0)
[maxeig,index]=max(abs(diag(Dia)));  % Find largest eigenvalue (=1)
else
eig1=diag(Dia);
for i=1:size(eig1)
%    if(2^(-d)==eig1(i))
if(abs(eig1(i)-(2^(-d)))<=10^(-6))
            index=i;
        eig2=eig1(i);
        break;
    end
end
end
v = V(:,index);        % Pick corresponding eigen vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%normalization of the eigenvector obtained%%%%%%%%%%%%%%%%%%%%%%%%
%sum2=0;
%for k=1:size(v)
%sum2=sum2+((-k)^(d))*v(k);                   
%end
%v=(factorial(d)*v)/sum2;
if d == 0
  Mrow = ones (1,D-1);
else
  Mom  = moments (lp, d);
  
  Mrow = zeros (1,D-1);
  for c=1:D-1,
    l       = c-1;              
    TMom    = tmoments (Mom,l);
    Mrow(c) = TMom(d);
  end       
end;
constant = (-1)^d*factorial(d)/((v)'*Mrow');
v=v*constant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=0:D-2 
  phi_d(k*stride + shift) = v(k+shift);   % Initialize phi at integers
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compute phi at successively finer levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if q > 0
  %%%%%%% Compute phi at half integers 
  start = stride2;                      
  phi_d(start+shift:stride:L+shift) = A1*phi_d(0+shift:stride:L-1+shift);

  %%%%%%% Do succesive levels: 
  while start > 1
    step  = start;           % Distance between subvectors
    start = start/2;         % Leftmost starting point (1/2, 1/4, 1/8,...)

    for k=start:step:stride2
      k2 = 2*k;              %Starting point for source vector
      ks = k + stride2;      %Starting point for high bit (A1) target vector
phi_d(k+shift :stride:L+shift) = (2^d)*A0*phi_d(k2+shift:stride:L+shift);
phi_d(ks+shift:stride:L+shift) = A1*phi_d(k2+shift:stride:L+shift);
    end;   
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compute psi_d from phi_d%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

psi_d = zeros(size(phi_d));
  srt=sqrt(2);
if q > 0 
  L2   = L/2;
  phi_d2 = phi_d(1:2:L);        % Phi_d at even numerators 

  first = 1;
  for k=1:D
    psi_d(first:first+L2-1) = psi_d(first:first+L2-1) + (2^(d))*srt*hp(k)*phi_d2;
    first = first + stride2;
  end;
else                        % q = 0  (a special case)
  for m=0:D-1
    for k=0:D-1
      index = 2*m-k;
      if (index >= 0) & (index <= D-1) 
        psi_d(m+shift) = psi_d(m+shift) + (2^(d))*srt*hp(k+shift)*phi_d(index+shift);
      end;	
    end;  
  end;
end;

