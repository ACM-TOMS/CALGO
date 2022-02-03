function PNxa=PNx_Legendre_vectN(x,n)
%This evaluates the  FIRST DERIVATIVE of the 
%  NORMALIZED, ORTHOGONORMAL Legendre polynomials P_{n)(x)
% for all degrees up to and including N.
% x may be either a scalar or a vector

%         Input: x= scalar or vector of grid points where Legendre polynomials
%                             are to be evaluated.
%                       n=degree of highest Legendre polynomial needed.
%
%         Output: PNxa is a  size(x) times (n+1) array

% Example:   let x = [0 0.3 0.9], n=3. Then the output will be the 3 x 4 array
%          PNxa= |  dP_0/dx(x(1)=0)     dP_{1}/dx(x(1))   dP_{2}/dx(x(1))  dP_{3}/dx(x(1))  |
%                    |  dP_0/dx(x(2)=0.3)  dP_{1}/dx(x(2))   dP_{2}/dx(x(2))  dP_{3}/dx(x(2))  | 
%                    |  dP_0/dx(x(3)=0.9)  dP_{1}/dx(x(3))   dP_{2}/dx(x(3))  dP_{3}/dx(x(3))  |


% Checked by LegendreFirstDeriv.maple;
%  Recurrence is that for Gegenbauer polynomials of order 3/2

PNxa=zeros(length(x),n);

	PNxa(:,1)=0;

if n > 0
	PNxa(:,2)=1 ;
end % if

if n > 1
	PNxa(:,3)=3*x';
end % if

if n > 2
for j=1:(n-2)
PNxa(:,j+3) =  (1/(j+1))   *  (   2*(j+3/2)*x' .*PNxa(:,j+2) -  (j+2)*PNxa(:,j+1)   )   ;
end % j	
end % if	



for j=1:(n+1)
	PNxa(:,j)=PNxa(:,j) * sqrt( j-1/2);
end % j
