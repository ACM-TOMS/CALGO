function [xprol,wprol,Deriv,Deriv2d]=TP_LobptswtsPERT(N,c)
% is called by Main program; 

% input:      N = number of grid points
%                    c = bandwidth parameter c
% output:   xprol = vector of length(1,N) containing the prolate-Lobatto grid points
%                    wprol = vector of length(1,N) containing the prolate-Lobatto quadrature weights
%                    Deriv(i,j) = d C_{j}/dx(x_i)
%                   Deriv2(i,j) = d**2 C_{j}/dx**2 (x_i)


cstar=  (pi/2) * (( N-1) + 0.5)
if c  >  cstar
	disp('WARNING! c > cstar!   Function returns NULL result')
else

%  Approximate Legendre-Lobatto points (i. e., prolate points and weights for c=0)
	ta=pi*( 0:(1/(N-1)):1);
xguess = - cos(ta);
wguess = ( pi/N) * (1- xguess.^2).^(1/2);
wguess(1)= 2/(  (N-1)*N);
wguess(N)=2/(  (N-1)*N);

 if N <= 12
	% direct call
% [xprol,wprol,Deriv,Deriv2d]=	TP_Lobptwtder2dGUESS(N,c,xguess,wguess);
[xprol,wprol]=TP_LobattoQuad(N,c,xguess,wguess);
else
			 [xL,wL]=TP_LobattoQuad(N,0,xguess,wguess);
[xprol1,wprol1]=TP_LobattoQuad(N,cstar/10,xguess,wguess);
			 xguess= xL + 100*(xprol1 - xL) * (c/cstar)^2;
		wguess=wL +  100*(wprol1 - wL) * (c/cstar)^2;
		[xprol,wprol]=TP_LobattoQuad(N,c,xguess,wguess);
end % if

end % if  on whether c is too large or not





% *******************************************************************
 % Next compute derivative matrix.
 
 [chi_array,B_array]=TP_eigprolODE(c,N-1);
 %          B_array contains the NORMALIZED Legendre coefficients of the 
%                   NORMALIZED eigenfunctions. The j-th column is the 
%                   eigenfunction corresponding the chi_array(j), i. e.,
%                   psi_{j-1}(x; c).
[psicoll,dummy]=size(B_array) ; 
PN2d=PN_Legendre_vectN(xprol,psicoll - 1) ;  % evaluate normalized Legendre
                      % at points xproll;  PN2d(ii,j)=P_{j-1}(xprol(ii)).
PSI_Gram= PN2d * B_array;  % PSI_Gram(i, j) =   psi( x_i, j ); matrix Psi in notes.

PNx2d=PNx_Legendre_vectN(xprol,psicoll-1);
PNxx2d=PNxx_Legendre_vectN(xprol,psicoll-1);

PSI_x= PNx2d * B_array;   % matrix Psi^{(x)} in notes
                                                % PSI_x(i,j) = d psi/dx(x_i, j)
												
PSI_xx= PNxx2d * B_array;   % matrix Psi^{(xx)} in notes
                                                % PSI_xx(i,j) = d**2 psi/dx**2(x_i, j)
																								
C_matrix=inv(PSI_Gram);   % Matrix   C in the noted.
Deriv = PSI_x * C_matrix;   % First derivatives of prolate CARDINAL basis.
Deriv2d = PSI_xx * C_matrix;   % Second derivatives of prolate CARDINAL basis.
