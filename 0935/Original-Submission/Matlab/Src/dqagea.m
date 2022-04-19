function [result,err,neval]=dqagea(func,a,b,abserr,relerr)
% Deals with definite integrals. Partions [a,b] into subintervals which are
% evaluated separately using quadgk. The subinterval with highest estimated
% error is partitioned again and reevaluated separately with quadgk.
% Iteration stops when required error bound is met. Note the similarity
% with dqagiea.m (for indefinite integrals)
%
% Essentially, this is a simplified version of SLATEC's dqage.m due to the
% recent development of quadgk.m, and no need to call a separate error
% sorting routine.

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input
% func = the function to be integrated
% a, b = endpoints of integration, NOT a,b from IIPBF.m
% abserr = required absolute error
% relerr = required relative error
%
% Output
% result = estimate of \int_a^b func
% err = estimated error
% neval = number of function evaluations
%
% Dependencies:
% quadgk.m
%
% Author: Sirong Zhang
% Last modified: 10/01/2010 by Jung Hun Kim
% - error sorting
% - corrected count of iterations (was originally being under-counted)
% - incorporate quadgk.m instead of dqk15.m
% initialize

% maximum number of partitions
maxlist=5000;

% vector containing the left endpoints of the integrals.
lends=zeros(maxlist,1);

% vector containing the right endpoints of the integrals
rends=zeros(maxlist,1);

% vector with the estimates from each intervals
sums=zeros(maxlist,1);

% vector with the error estimates from each intervals
errs=zeros(maxlist,1);

%first estimate
last=1;
if (a>b)
    t=b; b=a; a=t;
end
lends(1)=a;rends(1)=b;
[sums(last),errs(last)]=quadgk(func,a,b ,'RelTol',0);
if isnan(sums(last))
    errsum=realmax;
else
    errsum=errs(last);
end
i=1;neval=1;maxerr=1;
area=sums(last);
errbnd=max(abserr,relerr*area);


%main loop
while ((last < maxlist && errsum>errbnd) || isnan(errsum))
    
    % bisect the interval with max error in to [a1,b1][b1,c1]
    a1=lends(maxerr);
    b1=0.5*(lends(maxerr)+rends(maxerr));
    c1=rends(maxerr);
    
    % re-estimate each of the new partitions with quadgk.
    [area1,err1]=quadgk(func,a1,b1,'RelTol',0);
    [area2,err2]=quadgk(func,b1,c1,'RelTol',0);
    
    neval=neval+1;
    narea=area1+area2;
    nerr=err1+err2;
    
    % subtract the estimates from the old partition.
    errsum=errsum+nerr-errs(maxerr);
    area=area+narea-sums(maxerr);
    errbnd=max(abserr,relerr*area);
    last=last+1;
    
    % insert the new estimates
    if (err2>err1)
        lends(maxerr)=b1;
        lends(last)=a1;
        rends(last)=b1;
        sums(maxerr)=area2;
        sums(last)=area1;
        errs(maxerr)=err2;
        errs(last)=err1;
    else
        rends(maxerr)=b1;
        lends(last)=b1;
        rends(last)=c1;
        sums(last)=area2;
        sums(maxerr)=area1;
        errs(maxerr)=err1;
        errs(last)=err2;
    end
    % find which interval has the largest error
    [~,maxerr]=max(errs);
    i=i+1;
    
end
if (last==maxlist)
    error('No convergence in the finite integration');
end
result=sum(sums);
err=sum(errs);
