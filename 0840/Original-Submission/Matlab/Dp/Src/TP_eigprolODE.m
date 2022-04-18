function [chi_array,B_array]=TP_eigprolODE(c,nmax);
% Computes the prolate spheroidal functions of order zero.
% These are normalized so that int_{-1}^{1} [\psi_{j}(x)]**2 dx=1
%     Input: c=bandwidth
%            nmax=degree of highest prolate function desired.
%            [number of prolate functions computed is (nmax+1)]

% output: chi_array contains the eigenvalues of the prolate 
%                   differential equation chi, sorted so the 
%                   smallest is first
%          B_array contains the NORMALIZED Legendre coefficients of the 
%                   NORMALIZED eigenfunctions. The j-th column is the 
%                   eigenfunction corresponding the chi_array(j), i. e.,
%                   psi_{j-1}(x; c).

neig=max(2*nmax+2,30);
N=nmax+1;
% matrix element from Eq.(62), pg. 814, of Xiao, Rokhlin and Yarvin
% Inverse Problems (2002)

  maindiag=ones(neig,1);    sup=ones(neig,1);   sub=ones(neig,1);
for j=1:neig
%	k=2*j - 2; % for symmetric modes
	k=j-1;   % for general case
maindiag(j)=k*(k+1) + c*c*( 2*k*(k+1)-1) /   (   (2*k+3)*(2*k-1)   ) ;
sup(j)= c*c*(k+1)*(k+2)/ ( (2*k+3)*sqrt( (2*k+1)*(2*k+5)) );
sub(j)=c*c*(k-1)*(k)/ ( (2*k-1)*sqrt( (2*k-3)*(2*k+1)) );

end % 
AA=zeros(neig,neig);

for j=3:neig-2
	AA(j,j)=maindiag(j);   AA(j,j+2)=sup(j);    AA(j,j-2)=sub(j);
end % j
AA(1,1)=maindiag(1);  AA(1,3)=sup(1);
AA(2,2)=maindiag(2);  AA(2,4)=sup(2);
AA(neig,neig)=maindiag(neig);   AA(neig-2,neig)=sub(neig);
AA(neig-1,neig-1)=maindiag(neig-1);   AA(neig-3,neig-1)=sub(neig-1);

[VV,DD]=eig(full(AA));
[chi_arrayfull,iindex]=sort(diag(DD));
chi_array=chi_arrayfull(1:N);
B_array=zeros(neig,N);

for j=1:N
	B_array_unnorm(:,j)=VV(:,iindex(j));
end % j

% The coefficients are in terms of normalized Legendre functions,
% but we now need to normalize the eigenfunctions so as to have 
% unit L2 norm on [-1, 1]. 
for j=1:N
	psi_squared(j) = B_array_unnorm(:,j)' * B_array_unnorm(:,j);
end % j
psinormfact = ones(1,N) ./ sqrt(psi_squared);

for j=1:N
	B_array(:,j) = B_array_unnorm(:,j) * psinormfact(j);
end % j

