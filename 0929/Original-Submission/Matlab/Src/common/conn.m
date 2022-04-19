function Gamma = conn (d,D)
%CONN  2-term connection coefficients Gamma(d,D)
%
%  Gamma = conn(d,D)
%
%  Compute 2-term connection coefficients 
%
%  Input:
%         d: Order of differentiation.
%            (Restriction: d < D/2) 
%         D: Genus (default 6) 
%            (Restriction D >= 4)         
%
%         The case d=2, D=4 is not valid.
%
%  Output:
%         Gamma: The connection coefficients. 


makeeigtest = 0;       % Find max d for which 2^(-d) is in eig(A)
tol         = 1.0e-1;  % related to eigtest

if nargin < 2
  D = 6;
end

if D < 4
  error('undefined for D<4')
end  
if D == 4 & d == 2
  error('undefined for D=4 and d=2')
end; 
if d >= D/2
error('d>=D/2')
end  

  
M  = 2*(D-2)+1;              % Number of coefficients   
[hk,gk] = wfilters(['db' num2str(D/2)],'r');     % Low and high pass filters    
hk=hk';gk=gk';


% Setup homogeneous system

A=zeros(M, M);
for r=1:M,
  for c=1:M,
    l = r-(D-1); 
    q = c-(D-1); 

    A(r,c) = 0; 
    for i = 1:D;
      p=i-1;
      qq = p+q-2*l;
      if (qq >= 0) & (qq <= D-1)
        A(r,c) = A(r,c) + hk(p+1)*hk(qq+1);  
      end
    end
  end
end  

if makeeigtest
  %--------------------------------------------------------
  % Make a numerical test of the maximal d for which 2^(-d) 
  % is an eigenvalue.
  %
  dvec = (0:M-1)';
  x = 2.^(-dvec)
  eigval = sort(eig(A))
  i = 0;
  j = 0;
  while i < M & j <= M
    maxind = i;
    i = i+1;
    j = search(eigval,x(i),tol);
  end;
  if maxind > 0 
    fprintf('Maximal d = %d, D = %d\n', dvec(maxind), D);    
    j = search(eigval,x(maxind),tol);
    fprintf('Eigenvalue = %e\n', eigval(j));
    fprintf('2^(-%5d) = %e\n', maxind-1, x(maxind));  
  else
    fprintf('No admissible d, D = %d\n', D);    
  end
end;

A = A - eye(M)/2^(d);


% Setup heterogeneous moment equation:

if d == 0
  Mrow = ones (1,M);
else
  Mom  = moments (hk, d);
  
  Mrow = zeros (1,M);
  for c=1:M,
    l       = c-(D-1);              % Translate into true coordinates  
    TMom    = tmoments (Mom,l);
    Mrow(c) = TMom(d);
  end       
end;


%--------------------
% Solve system
%
[U,S,V] = svd(A);
c = (-1)^d * factorial(d) / (Mrow * V(:,M));
Gamma = c * V(:,M);



