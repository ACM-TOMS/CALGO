function PNxxa=PNxx_Legendre_vectN(x,n)
%This evaluates the  SECOND DERIVATIVE of the 
%  NORMALIZED, ORTHOGONORMAL Legendre polynomials P_{n)(x)
% for all degrees up to and including N.
% x may be either a scalar or a vector

%         Input: x= scalar or vector of grid points where Legendre polynomials
%                             are to be evaluated.
%                       n=degree of highest Legendre polynomial needed.
%
%         Output: PNxxa is a  size(x) times (n+1) array

% Example:   let x = [0.1 0.3 ], n=3. Then the output will be the 2 x 4 array
%      PNxxa= |  P_{0,xx}(x(1)=0.1)  P_{1,xx}(x(1))  P_{2,xx}(x(1))  P_{3,xx}(x(1))  |
%                    |  P_{0,xx}(x(2)=0.3)  P_{1,xx}(x(2))  P_{2,xx}(x(2))  P_{3,xx}(x(2))  | 


% P"_{n+3} = (1/(n+1)) * { 2 (n+5/2)*x*P"_{n+2} - (n+4) P"_{n+1}   }
% This is the recurrence for Gegenbauer of degree (n+1) and order (5/2)
%    on pg. 503 of my book. 

% Checked by  LegendreSecondDeriv.maple;

PNxxa=zeros(length(x),n);

	PNxxa(:,1)=0;

if n > 0,	   PNxxa(:,2)=0;           end % if

if n > 1,    PNxxa(:,3)=3;            end % if

if n > 2,    PNxxa(:,4)=15*x' ;   end % if

if n > 3
for j=1:(n-3)
PNxxa(:,j+4) =  (1/(j+1))   *  (   2*(j+5/2)*x' .*PNxxa(:,j+3) -  (j+4)*PNxxa(:,j+2)   )  ; 
end % j	
end % if	

for j=1:(n+1),   PNxxa(:,j)=PNxxa(:,j) * sqrt( j-1/2);       end % j
