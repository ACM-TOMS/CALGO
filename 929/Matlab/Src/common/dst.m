
function c=dst(f, Wavelet)
%DST   Discrete Scaling-function Transform.
%
%  c = dst(f, D)
%
%  Input
%    f:  Matrix, where columns are the vectors to be transformed.
%    D:  Wavelet genus - optional, default D=4.
%
%  Output
%    c:  Matrix, where columns are the scaling function coefficients
%        corresponding to f
%
% Dependencies
% dstmat
%
%
%
%  See also idst, dstmat.



%%%%%%%%%%%%%%% Parameter check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
  Wavelet = wfilters('db2','r'); % Default filters (low pass filter coefficients)
end;
  
len = length(Wavelet);  % Determine the nature of argument #2:
if len > 1
  D  = len;
  hk = Wavelet;         % First argument is a vector of low pass filter coefficients
else
  D  = Wavelet;
  hk = wfilters(['db' num2str(D/2)],'r');     % First argument is the wavelet genus
end;    
N = size(f,1);       % number of rows of f
if N == 0
  fprintf('ERROR (dst.m): Null vector f\n');
  c = [];
  return;
end;  
[rows, cols] = size(f);
if rows == 1
  fprintf('ERROR (dst.m): f is a row vector - no transform made\n');
  c = [];  
  return;
end;  

%%%%%%%%%%%%%%% End parameter check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = N;                       % Factor N=K*2^m
while rem(L,2)==0,  
  L = L/2;
end;  
r = log2(N/L);              % Get resolution of f

j=r;    %Get as many coefficients as samples.
% j=input('dyadic resolution for scaling function coefficients i.e.  j');
q=r-j;  %Get sufficient resolution for phi
%q=input('enter the value of q(dyadic resolution of scaling function) such that j+q-r >= 0')
% while(j+q-r < 0)
%    q=input('enter the value of q(dyadic resolution of scaling function) such that j+q-r >= 0');
% end
T = dstmat(D,r,j,q,L);
c = T\f;               %Solve f=T*c for c. 
