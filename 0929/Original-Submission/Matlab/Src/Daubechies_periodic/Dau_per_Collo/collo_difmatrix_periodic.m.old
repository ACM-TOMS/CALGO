function T=collo_difmatrix_periodic(wavelet,r,j,q,K,d)
%.................................................................
% Get matrix T with derivative of scaling function values.
% It is also the differential matrix for collocation method
%.................................................................
%   T=diff_collo_daub(D,r,j,q,K,d_n)
%
%  Input
%    wavelet:  Wavelet genus - optional, default D=4.
%    r:  Dyadic resolution of function.
%    j:  Dyadic resolution of scaling function coefficients.
%    q:  Dyadic resolution of scaling function.
%    K:  Integer.
%    d: order of derivative required
%  Output
%    T:  Matrix transforming between a function f in V_r
%        and a vector of scaling function coefficients c.
%   
%  T will have K*2^r rows and K*2^j columns.
%
%  The following must hold:  j+q-r >= 0 
%
%  Example:
%
%    T = collo_difmatrix_periodic(4,3,3,0,1,0)
% Dependencies
%cascade_der.m
%

if nargin < 5
  K = 1;      % Factor which are not a positive pwr of two (integer)
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Parameter check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = length(wavelet);  %Determine the nature of argument #1:
if len > 1
  D  = len;             % Use actual filter length
  lp = wavelet;         % First argument is the wavelet filter (from filters.m)
else
  D  = wavelet;         % First argument is the wavelet genus
   [hk,gk] = wfilters(['db' num2str(D/2)],'r');      % Low and high pass filters
   hk=hk';gk=gk';
end;    

shift=1;

if j+q-r < 0 
  error('Resolution of phi is too coarse. Adjust j, q, or p so j+q-r >= 0');
end;

% if (2^j< D-1)
% error('Adjust j and D so that  2^j >= D-1');
% end
[x,phi_d] = cascade_der(D,q,d); % Get phi_d:		
phi_d;

N=K*2^r;    % Number of samples
P=K*2^j;    % Number of coefficients

Q=2^q;      % The resolution of phi
            % (D-1)*Q is the number of dyadic points where phi is known.

%T = zeros(N,P);
T = spalloc (N,P,P*(D-1)*2^(r-j)); % There will be (D-1)*2^(r-j) entries 
                                 % in each column

S = P*Q/N;  % step between values needed in phi is 2^(j-r+q) 

phid_s = phi_d(1:S:length(phi_d));
len   = length(phid_s);

if len < N                % phi_s fits into a column
  for l=0:P-1
    firstrow = l*N/P;
    lastrow  = firstrow + len - 1;
    wrap = lastrow - (N-1);
    if wrap > 0, 
      lastrow = lastrow - wrap;
      T(0+shift:wrap-1+shift, l+shift) = phid_s(len-wrap+1:len);     
    else 
      wrap=0; 
    end;
    T(firstrow+shift:lastrow+shift, l+shift) = phid_s(1:len-wrap);
  end;  
else     % RARE (If vectors are small): 
         % phid_s does not fit into a column (use modulus wrapping)
  for l=0:P-1   % For each column
    firstrow = rem(l*N/P,N);
    for k=0:len-1                % For element in phid_s
      row = rem(firstrow + k, N);
      T(row+shift,l+shift) = T(row+shift, l+shift) + phid_s(k+shift);     
    end;
  end;
end;      
T=(2^((d*j)+(j/2)))*T;

