function PNa=PN_Legendre_vectN(x,n)
%This evaluates the NORMALIZED, ORTHOGONORMAL Legendre polynomials P_{n)(x)
% for all degrees up to and including n.
% x may be either a scalar or a vector
%         Input: x= scalar or vector of grid points where Legendre polynomials
%                             are to be evaluated.
%                       n=degree of highest Legendre polynomial needed.
%
%         Output: PNa is a  size(x) times (n+1) array

% Example:   let x = [0 0.3 0.9], n=3. Then the output will be the 3 x 4 array
%          PNa= |  P_0(x(1)=0)  P_{1}(x(1))  P_{2}(x(1))  P_{3}(x(1))  |
%                    |  P_0(x(2)=0.3)  P_{1}(x(2))  P_{2}(x(2))  P_{3}(x(2))  | 
%                    |  P_0(x(3)=0.9)  P_{1}(x(3))  P_{2}(x(3))  P_{3}(x(3))  |

PNa=zeros(length(x),n+1);

	PNa(:,1)=1;
	if n > 0
	PNa(:,2)=x';
end % if
if  n > 1
for j=1:(n-1)
PNa(:,j+2) = (1/(j+1))*( (2*j+1)*x' .*PNa(:,j+1) - j*PNa(:,j));
end % j
end % if

for j=1:(n+1)
	PNa(:,j)=PNa(:,j) * sqrt( j  -   1/2); 
end % j
