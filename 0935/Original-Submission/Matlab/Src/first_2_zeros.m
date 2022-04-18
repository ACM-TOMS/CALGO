function [zero1,zero2,its]=first_2_zeros(func,h_type,nu1,nu2,rho,tau)
% Calculates the first zero of the h_1 and h_2 after Y_a(x) or Y_b(x) no
% longer dominates. We then use a step function to get a more refined zero.
% For details, refer to Lucas [1995] in paper.

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input:
% func= h_1 or h_2
% h_type= 1 for h_1 or 2 for h_2, not to be confused with type from IIPBF
% nu1= equivalent to a in IIPBF.m
% nu2= equivalent to b in IIPBF.m
% rho,tau = equivalent to rho,tau in IIPBF.m
%
% Output:
% zero1= first zero after Y_a(x) or Y_b(x) no longer dominates
% zero2= second zero
% its= iterations needed to derive the zeros.
% 
% NOTE: A if-then condition was added for cases of extreme values of rho
% and tau, as this would lead zero1 and zero2 to be extremely large and
% likely lead to inaccurate and/or non-convergent quadrature.
% 
% Author: Sirong Zhang
% Last update: Tilak Ratnanather 
% Revision Date: 25/09/2012

if (nu1==0)   %for Y_a zeros
    yazero=0.8936/rho;
else
    yazero=nu1+0.9315768*nu1^(1/3)+0.260351*nu1^(-1/3)+0.01198*nu1^(-1)...
        -0.0060*nu1^(-5/3)-0.001*nu1^(-7/3);
    yazero=yazero/rho;
end

if (nu2==0)  %for Y_b zeros
    ybzero=0.8936/tau;
else
    ybzero=nu2+0.9315768*nu2^(1/3)+0.260351*nu2^(-1/3)+0.01198*nu2^(-1)...
        -0.0060*nu2^(-5/3)-0.001*nu2^(-7/3);
    ybzero=ybzero/tau;
end

% choose the zero
% In cases where |a-b|>5, choosing the larger zero may improve results, 
% since simple oscillations may occur after a larger interval in these cases
ymint = min(yazero,ybzero);
ymaxt = max(yazero,ybzero);
if(abs(ymaxt/ymint)>100)
   yzero=ymint;
else
   yzero=ymaxt;
end

% use small step estimate to find the next zeros.
[zero1,its1]=step_zero(func,h_type,yzero,rho,tau); % first zero;
[zero2,its2]=step_zero(func,h_type,zero1,rho,tau); % second zero;
its=its1+its2;
end

function [newzero,its]=step_zero(func,h_type,oldzero,rho,tau)
% compute the new zeros by  step method.
its=1;
if strcmp(h_type,'1')
    step=pi/(4*(rho+tau));
else
    step=pi/(4*abs(rho-tau));
end
oldzero=oldzero+step;
x1=oldzero+step;
fx0=feval(func,oldzero);
fx1=feval(func,x1);
while (fx0*fx1>=0)
    x1=x1+step;
    its=its+1;
    fx0=fx1;
    fx1=feval(func,x1);
end

newzero=0.5*(x1+x1-step);
end
