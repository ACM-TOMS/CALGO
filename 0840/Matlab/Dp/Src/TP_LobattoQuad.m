function [xprol,wprol]=TP_LobattoQuad(N,c,xguess,wguess);
global exact_integral  psicoeffsymm;
% TP_LobattoQuad computes the prolate-Lobatto-Gaussian quadrature points and weights
%      input: N = number of collocation points
%                 c = bandwidth parameter
%                xguess --- vector of length(1,N) containing first guess for the Lobatto-grid points
%                wguess -- vector of length(1,N) containing first guess for Lobatto-prolate quadrature weights.

%
%  output:   xprol is vector of  size(1,N) with the final prolate-Lobatto points
%                 wprol is vector of size(1,N), the corresponding quadrature points

M = (N/2) - 1 - (1/2)*rem(N,2);  %if mod(N,2)==0  % npts is even  %M= N/2 - 1
                                                 %else   % npts is ODD   %M= (N-3)/2;  %end % if

% M is the number of non-negative quadrature abscissas

%                 xandw0 = vector of length(N-1) containing a first guess for the unknowns in 
%                                 the Newton iteration. The first M entries are the grid points on the
%                                 interval   0 < x_k < 1. The remaining entries are the quadrature weights
%                                 associated with the nonnegative collocation points, i. e., weights for 
%                                 x=1 and, if N is odd, x=0 are included

if mod(N,2)==0
	% N is even
xnonneg=xguess(M+2:2*M+1);       wnonneg=wguess(M+2:N);
else
	xnonneg=xguess(M+3:N-1);   	wnonneg=wguess(M+2:N);
end % if
xandw0 = [xnonneg wnonneg];


% step one: compute the exact integral of each prolate function 
% of degree 0, 1, ... (2*npts-1)

% Because of symmetry, we only use the integrals of EVEN prolate functions.
       nmax=2*N;   % degree of higest prolate function that will be integrated exactly.

[chi_array,B_array]=TP_eigprolODE(c,nmax);
	  [psirow,psicol]=size(B_array) ;
	  
psicoeffsymm = B_array(1:2:psirow,1:2:(2*(N-1)-1)  ) ;

% Because the integral of overline(P_{j}) on [-1,1] is the inner product of 
%  P_{j} with itself and the Legendre are orthogonal, the integral 
% is just the constant in the prolate normalized Legendre series
exact_integral = sqrt(2) * psicoeffsymm(1,1:N-1);

% apply a Newtoniteration until the quadrature is exact
NewtMax=20;
nalpha=4;  %    nalpha=0;
JacMod=1;

%  [xandw,residratio]=TP_Newton(N-1,nalpha,NewtMax,JacMod,xandw0,lambda,'TP_Lobattoresid');
xandw=TP_NewtonJac(xandw0,c,'TP_Lobattoresid');

% process roots: xandw is a vector of dimension (N-1) contains only nonnegative x
%                        and quadrature weights

if mod(N,2)==0
	% N is even
xprol=[-1   -xandw(M:-1:1) xandw(1:M)  1] ;
wprol=[xandw(N-1:-1:M+1)  xandw(M+1:N-1)]  ;

else
	% N is ODD
xprol=[-1   -xandw(M:-1:1)  0     xandw(1:M)  1];
wprol=[xandw(2*M+2:-1:M+2)   xandw(M+1:2*M+2)];
end % if

% xprol, wprol are each N-dimensional vectors containing 
% all grid points and weights, respectively.
