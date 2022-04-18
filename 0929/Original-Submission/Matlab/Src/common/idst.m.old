function f=idst(c, Wavelet, N)
%IDST Inverse Discrete Scaling-function Transform.
%
%  f = idst(c, D, N)
%
%  Input
%    c:  Matrix, where columns are scaling function coefficients
%    D:  Wavelet genus - optional, default D=4.
%    N:  Desired resolution of f - optional, default N=length(c).
%
%  Output
%    f:  Matrix, where columns are the transformed vectors.
%
%
%  See also dst, dstmat.

shift=1;
%%%%%%%%%%%%%%% Parameter check%%%%%%%%%%%%%%%%
M=size(c,1)  %% number of rows of matrix c
if M == 0
  fprintf('ERROR (idst.m): Null vector c\n');
  f = [];
  return;
end;  
[rows, cols] = size(c);
if rows == 1
  fprintf('ERROR (idst.m): c is a row vector - no transform made\n');
  f = [];  
  return;
end;  
% Factor lengths into powers of two and odd integers
Lc = M;              % Factor M=Lc*2^j
while rem(Lc,2)==0,  
  Lc=Lc/2;
end;  
j = log2(M/Lc);  % Get resolution of c

if nargin < 3
  r = j;               % Default resolution of f is that of c
else
  Lf = N;              % Factor N=Lf*2^r
  while rem(Lf,2)==0,  
    Lf=Lf/2;
  end;  
  r = log2(N/Lf);  % Get resolution of f  
  
  if Lf ~= Lc
    fprintf('ERROR (idst.m): Resolutions of c and f are incompatible,\n');
    fprintf('                N must be a multiple of length(c).\n');    
    f = [];
    return
  end;
end;    
  
if r < j,
  fprintf('Warning (idst.m):\n   Resolution (r) of f is less than resolution (j) of c. Using r=j.\n');
  r = j;
end;

if nargin < 2
  Wavelet = wfilters('db2','r'); % Default filters 
end;
  
len = length(Wavelet);  % Determine the nature of argument #2:
if len > 1
  D  = len;
  hk = Wavelet;     % First argument is the wavelet filter (from daubfilt.m)
else
  D  = Wavelet;
  hk = wfilters(['db' num2str(D/2)],'r'); % First argument is the wavelet genus
end;    

%%%%%%%%%%%%%%% End parameter check%%%%%%%%%%%%%%%%%%
  
q = r-j;                 % Sufficient resolution for phi.
T = dstmat(hk,r,j,q,Lc);

f=T*c;


