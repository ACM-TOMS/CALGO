function newzero=nth_zero(zero2,a,b,rho,tau,type,h1_or_h2)
% Find the n-th zeros of h_1 or h_2 to define the subintervals to be
% used subsequently as endpoints of the integrals used in acceleration 
% algorithms. We use this approximation as a more precise zero-finding
% algorithm has insignificant impact in the overall efficiency of the
% toolbox. 

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input:
% zero2 = second zero derived from first_2_zeros.m
% func = function handle
% a,b = a,b from IIPBF.m
% rho,tau = rho,tau from IIPBF.m
% type = 'JJ', 'JY' or 'YY'
% h1_or_h2 = 1 if the zero derived is for integrating h_1, 2 for �
% integrating h_2.
% 
% Output:
% newzero = the next zero as approximated, assuming the functions behave
% like the sin function.
% 
% NOTE: Prior versions was some variation of Newton's method, which
% involved approximating derivatives and consequently required additional
% function evaluations (which is our yardstick for measuring algorithm
% efficiency). Empirical tests show there is no significant differences in
% accuracy of the quadratures if a more accurate zero is chosen.
% 
% Author: Jung Hun Kim
% Last Revision Date: 30/08/2012

newzero=0;
k = 0;t = 2;
if h1_or_h2 == 2
    newzero=zero2+pi/(abs(rho-tau));
    return;
end 
while (zero2 >= newzero)
    switch type
        case 'JJ'
            newzero = (k*pi + ((a+b+1)*t*pi))/(rho+tau); 
        case 'JY'
            newzero = (k*pi + ((a+b+1)*t*pi) + pi/2)/(rho+tau);
        case 'YY'
            newzero = (k*pi + ((a+b+1)*t*pi) )/(rho+tau);
        otherwise 
            error('No such type')
    end
    k = k+1;
end
