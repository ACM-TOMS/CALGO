function [result,err,neval]=dqagiea(func,bound,abserr,relerr)
% Deals with improper integrals. Simplified version of SLATEC's dqagie.m
% Essentially combines the bisection/quadrature method with the epsilon
% algorithm for faster convergence.

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input:
% func = function to be integrated
% bound = lower bound of integrand
% abserr = required absolute error
% relerr = required relative error
% 
% Output
% result = estimate of \int_0_\infty func
% err = estimated error
% neval = number of evaluations
% 
% Dependencies
% quadgk.m - note this is the same quadgk.m used for definite integrals
% dqext.m - SLATEC's epsilon algorithm.
% 
% Author: Sirong Zhang
% Last modified: 10/01/2010 by Jung Hun Kim
% - incorporate quadgk.m
% - error sorting
% - correct call of epsilon algorithm
% Last modified: 26/08/2012 by Jung Hun Kim
% - added absolute and relative tolerance in the first estimate in [0,inf]

maxlist=500;
lends=zeros(maxlist,1);
rends=zeros(maxlist,1);
sums=zeros(maxlist,1);
errs=zeros(maxlist,1);

%first estimate in [0,infinity];
last=1;
lends(1)=bound;rends(1)=1;
[sums(last),errs(last)]=quadgk(func,bound,Inf,'AbsTol',abserr,'RelTol',relerr);

if (isnan(sums(last)) || isnan(errs(last)))
    errs(last)=100; sums(last)=100;
end

errsum=errs(last);area=sums(last);
neval=1;maxerr=1;
errbnd=max(abserr,relerr*abs(area));

%set up the epsilon table and last three results to prepare for eps algo.
maxtab=50;
etab=zeros(1,maxtab);
etab(1)=area;
nres=0;ntab=1;
err=realmax;old3r=zeros(1,3); 

%main loop in case quadgk above fails or not accurate enough
while (last < maxlist && errsum>errbnd)
    %bisect the interval with max error in to [a1,b1][b1,c1]
    a1=lends(maxerr);
    b1=0.5*(lends(maxerr)+rends(maxerr));
    c1=rends(maxerr);
    %get new estimates;
    [area1,err1]=quadgk(func,a1,b1);
    [area2,err2]=quadgk(func,b1,c1);

    %update the results and calculate new errorbound.    
    narea=area1+area2; nerr=err1+err2;
    errsum=errsum+nerr-errs(maxerr);
    area=area+narea-sums(maxerr);
    errbnd=max(abserr,relerr*abs(area));
    
    neval=neval+2;
    %delete the last maxerr element and save the old error;
    last=last+1;
    errlast=errs(maxerr);
    %add the new interval;
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
    
    %sort the error and save the position of each interval into pos;
    [~,pos]=sort(errs(1:last));
    maxerr=pos(last);
    
    %check errorbound conditions exit main loop if estimated error is small
    if (errsum <errbnd)
        continue;
    end
    if (last==2)
        %initial large intervals estimate;
        small=0.375;
        errlarge=errsum;
        errtest=errbnd;
        etab(2)=area;
    end
    
    errlarge=errlarge-errlast;
    if (abs(b1-a1) >small)
        errlarge=errlarge+nerr;
    end
    if (abs(rends(maxerr)-lends(maxerr)) >small)
        continue
    end
    
    % Epsilon algorithm if necessary
    if (errlarge < errtest)
        ntab=ntab+1; etab(ntab)=area;
        nres=nres+1;
        [ntab,etab,result,eperr,old3r,nres]=dqelga(ntab,etab,old3r,nres);
        if (eperr<err)
            err=eperr;
            errtest=max(abserr,relerr*result);
            if (err<errtest)
                errbnd=realmin;
                continue;
            end
        end
        
        if result<eps
            error('epsi error');
        end
        
        %prepare for bisecting the smallest interval;
        maxerr=pos(last);
        small=small*0.5;
        errlarge=errsum;    
    else
        %find the new maximum error of of all large interval.
        nmax=2;
        find=0;
        while (nmax<last && find==0)
            maxerr=pos(last+1-nmax);
%             errmax=errs(maxerr);
            if (abs(rends(maxerr)-lends(maxerr)) <= small) 
                nmax= nmax+1;
            else
                find=1;
            end
        end    
    end
end

if (last==maxlist)
    error('Minimum step size reached. Convergent error for infinite integral.');
end

%set the final result
if (err>errsum || err/abs(result)>(errsum/abs(area)))
    result=sum(sums);
end

if (result<eps);
        result=0;
end
err=errsum;
