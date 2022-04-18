function case13_func = case13(u,rho,tau)
% Helper function in test case 13 of reference paper below
%
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J.,
%              and Lucas, S. K. 2011. Algorithm: IIPBF, a MATLAB toolbox
%              for infinite integral of products of two Bessel functions.
%              Pending review from ACM Transactions on Mathematical
%              Software
%
% The function is used to verify the accuracy of the estimation of the case
% 13 with the same parameters, and is used in TestCases_1to23.m. It has no
% role in the actual estimation.
%
% Input:
% u,rho,tau: parameter for f(x)=exp(-u*x)*besselj(0,rho*x)*besselj(0,tau*x)
% 
% Output:
% case13func: function handle.
% 
% NOTE: The formula needs to be revised as there are conflicting versions
% from various sources (??)
%
% Author: Sirong Zhang & Jung Hun Kim
% Revision Date: 30/08/2012

    function Fy = case13_ans(y)
        
        u2 = power(u,2);
        rho2 = power(rho,2);
        tau2 = power(tau,2);
        y2 = power(y,2);
        temp = power(u2 + rho2 - y2,2) + 4*u2.*y2;
        num1 = sqrt(temp) +u2 +rho2 -y2;
        den1= sqrt(temp);
        Fy =-(2/pi).*sqrt(num1/2)./den1./sqrt(y2 - tau2);
        
    end
case13_func = @case13_ans;
end