function resid=TP_Lobattoresid(xprolwprol, c)
global exact_integral   psicoeffsymm  ;
     global xproltemp;

N=length(xprolwprol) + 1;

	 
resid = exact_integral;

if mod(N,2)==1
M=length(xprolwprol)/2 - 1;   % One grid point is always
                                                      % at x=1 for Lobatto grid
						      % M is the number of non-negative, non-unit grid points

xprol=xprolwprol(1:M);   wprol=xprolwprol(M+1:2*M+2);
else
   M = (N-2)/2;   % N is even
   xprol=xprolwprol(1:M);     wprol=xprolwprol(M+1:2*M+1);
end % if

% First step: compute grid point values of the Legendre polynomials.
[NEIG,psicoeffcoll]= size(psicoeffsymm);
if mod(N,2)==1
	xprolfull = [0 xprol 1];   % 
else  % N even
	  %   disp('N is even !!!!!!!!!!!!!!!')
      xprolfull=[xprol   1]  ;     % xprolfull includes the fixed grid pt., x=1
  end % if
  
PN2dquad=PN_Legendre_vectN(xprolfull,2*NEIG);  % evaluate normalized Legendre
                      % at points xprolfull;  P+N2d(ii,j)=P_{j-1}(xprolfull(ii)).
								   
PN2dquadsymm=PN2dquad(:,1:2:2*NEIG);   % grid point values of symetric functions only 
% compute grid point values of the prolate functions at the quadrature points

psi_j = PN2dquadsymm * psicoeffsymm;  %psi_j =   psi( x_i, j )

if mod(N,2)==0
	approx_integral=2* psi_j' * wprol'  ;  
else
	% if N is odd, then xprol(1)=0, and we must not double-count it.
%approx_integral=2*   (psi_j(:,2:N)  )' * (wprol (2:N))'  + psi_j(:,1)*wprol(1);
approx_integral=zeros(N-1,1);

for j=1:N-1
	approx_integral(j)=psi_j(1,j)*wprol(1);
	%  disp([' j=',int2str(j),'  ',num2str(approx_integral(j)),' psi(j,1)=',num2str(psi_j(1,j)), 'wprol(1)=',num2str(wprol(1))   ])     
	for k=2:length(wprol)
	approx_integral(j)= approx_integral(j) + 2*psi_j(k,j) * wprol(k);
 % disp([num2str(k),'  ',num2str(approx_integral(j)),num2str(2*psi_j(k,j)*wprol(k)),' psi(j,k)=',num2str(psi_j(k,j)), 'wprol(k)=',num2str(wprol(k))   ])     
end % k
end % j

end % if

resid= exact_integral' -  approx_integral;
