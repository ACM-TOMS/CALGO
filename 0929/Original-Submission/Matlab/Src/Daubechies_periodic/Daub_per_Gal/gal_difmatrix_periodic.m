function DM=gal_difmatrix_periodic(d,N,L,D,flag)
%DIFMATRIX  differentiation matrix wrt scaling functions.
%
%  DM=difmatrix(d,N,L,D)
%
%  Input
%    d: Order of differentiation
%    N: Size  of matrix
%    L: Interval length
%    D: Genus of wavelet
%
%  Output
%    DM: Resulting matrix (banded, circulant)
%
% See also conn.


% Get connection coefficients
Gamma = conn(d,D);            % Get <phi'',phi_k> 
Gamma = (N/L)^d*Gamma;        % Scale according to x
 
  
% Setup differentiation matrix:

NN = N;
while NN <= D-2 %Protect the modulus function (rem) from negative arguments
  NN = 2*NN;  
end;         

DM = spalloc(N,N,N*(2*D-3));  	    
for l=1:N
  for d = -(D-2):D-2
    k = rem(NN+l-d-1,N)+1;
    DM(l,k) = DM(l,k) + Gamma(d+D-1);
  end
end





