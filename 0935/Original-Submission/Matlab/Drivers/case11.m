function case11_func = case11(u,rho,tau)
% Helper function in test case 11 of reference paper below
%
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J.,
%              and Lucas, S. K. 2011. Algorithm: IIPBF, a MATLAB toolbox
%              for infinite integral of products of two Bessel functions.
%              Pending review from ACM Transactions on Mathematical
%              Software
%
% The function is used to verify the accuracy of the estimation of the case
% 11 with the same parameters, and is used in TestCases_1to23.m. It has no
% role in the actual estimation.
%
% Input:
% u,rho,tau: parameter for f(x)=exp(-u*x)*bessely(0,rho*x)*bessely(0,tau*x)
% 
% Output:
% case11func: function handle.
% 
% NOTE: The formula needs to be revised as there are conflicting versions
% from various sources (??)
%
% Author: Sirong Zhang & Jung Hun Kim
% Revision Date: 30/08/2012


    function y = case11_ans(x)
        u2 = power(u,2);
        s2 = power(rho,2);
        t2 = power(tau,2);
        x2 = power(x,2);
        alpha1 = sqrt(sqrt(power((u2 + t2 - x2),2) + 4.*u2.*x2) + ...
            (u2 + t2 - x2)) ./ sqrt(2);
        alpha2 = sqrt(sqrt(power((u2 + t2 - x2),2) + 4.*u2.*x2) - ...
            (u2 + t2 - x2)) ./ sqrt(2);
        alpha1_2 = power(alpha1,2);
        alpha2_2 = power(alpha2,2);
        alpha2x = alpha2 + x;
        alpha1u = alpha1 + u;
        alpha2x_2 = power(alpha2x,2);
        alpha1u_2 = power(alpha1u,2);
        num = sqrt(alpha1u_2 + alpha2x_2);
        gx = -(2 ./ (pi.*(alpha1_2 + alpha2_2))) .* ...
            ((alpha1 .* log(num ./ tau) ) + (alpha2 .* atan(alpha2x./alpha1u)));
        y = -(2/pi).* gx./sqrt(x2 - s2);
    end
case11_func = @case11_ans;
end

