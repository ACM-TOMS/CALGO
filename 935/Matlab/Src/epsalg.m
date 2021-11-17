function [result,err,neval]=epsalg(zero1,zero2,func,abserr,relerr,rho,tau)
% The wrapper for the epsilon algorithm. We use this routine to generate
% the sequence to be input in the epsilon algorithm. The routine works in 
% the following way:
% 1. Obtain an initial estimate of \int_0^Inf f(x) h_2 dx using quadgk
% 2. If the initial estimate has high\isnan then
%    a. Obtain the next zero of f(x) h_2 
%    b. quadgk( f(x) h_2, zero(i), zero(i+1))  
%    c. Run epsilon algorithm with new sequence
%    d. Back to a if sum of error is larger than error bound. 

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input:
% zero1,zero2: first two zeros as the beginning point and to get later zeros.
% func: function handle;
% abserr,relerr: error bounds.
% 
% Output:
% result: estimate of the integral;
% error:  estimate of the error based on the algorithm;
% neval:  the number of function evaluation.
% 
% Author: Sirong Zhang.
% Last Modified: 10/02/2010 Jung Hun Kim
% - inclusion of quadgk
% - eliminate redundancy
% - resolved function name conflicts
% - incorrect counting of function evaluations
% Dependencies: 
% dqelga - SLATEC's version of epsilon algorithm

%initial values;
tabmax=500;
parsum=0;      
errorbound=0;
diff=realmax;
left=zero1;
right=zero2;
nres=0;  
count=0;
persistent ntab;
ntab=0;  
neval=0;

%the table for epsilon-algorithm;
persistent tab;
tab=zeros(1,tabmax); 

old3r=zeros(1,3);
%main loop
while (count<tabmax && diff>errorbound);
    area=quadgk(func,left,right,'AbsTol',abserr);
    if isnan(area)
        area=0;
    end
    newzero=right+pi/(abs(rho-tau));its=0;
    left=right;
    right=newzero;
    ntab=ntab+1;
    parsum=parsum+area;
    tab(ntab)=parsum;
    if (ntab>=3)
        [ntab,tab,result,diff,old3r,nres]=dqelga(ntab,tab,old3r,nres);
        errorbound=max(abserr,relerr*abs(result));
    end
    count=count+1;
    neval=neval+its+1;

end
if(count==tabmax && diff>errorbound)
    error('Warning: no convergence for extrapolation in epsilon-algorithm');
end
err=diff;
neval=neval+nres;        
end   

    
    
 

    
    
    


