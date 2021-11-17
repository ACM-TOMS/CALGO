function [result,ierror,neval]=mWtrans(zero1,zero2,func,abserr,relerr,a,b,rho,tau,type)
% mW-Transform implementation to estimate infinite oscillatory integrals.

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input:
% zero1 = the first zero derived from first2zeros.m
% zero2 = second zero derived from first2zeros.m
% func = function handle
% abserr = required absolute error for integration
% relerr = required relative error for integration
% a,b = a,b from IIPBF.m
% rho,tau = rho,tau from IIPBF.m
% type = 'JJ', 'JY' or 'YY'
%
% Output:
% result = estimate via the extrapolation
% ierror = estimate error
% neval = number of iterations
%
% NOTE: Vectorization of the mWtrans procedure proved to be difficult due
% to the complication arising from vectorizing quadgk. Subsequent
% calculation of the estimate has been vectorized however. The assumption
% being all reasonable 'func' with the given parameters will converge
% within 100 iterations.
% 
% Author: Sirong Zhang
% Last Modified: 10/01/2010 by Jung Hun Kim
% - incorporate quadgk.com
% - criteria for termination
%
% Dependencies:
% quadgk.m - replaces dqk15.m
% nth_zero.m - for finding successive n-th zeros.
%

%% Initialize
% max number of iterations
datamax=100;

% tables for m-W reaformation
M = zeros(datamax);
N = zeros(datamax);

% estimates of the integral
W=zeros(datamax);

% matrix to be used in updating W
u = [1 -1; -1 1]; temp=tril(repmat(u,datamax/2,datamax/2));

% endepoints of integrals
x=zeros(datamax,1);

% make sure we integrate at least three times;
W(1)=1;W(2)=1;

left=zero1;
right=zero2;

% First estimate
fxs=quadgk(func,left,right,'AbsTol',abserr);

% Shift interval with the new zero
newzero=nth_zero(right,a,b,rho,tau,type,1);
left=right; right=newzero; x(1)=right;

% Note 'neval' = number of evaluations.
neval=1;count=1;
DIFF=realmax; errorbound=0;



%% Main loop
while (count<datamax-2 && DIFF>errorbound)
    newzero=nth_zero(right,a,b,rho,tau,type,1);
    [phixs,~]=quadgk(func,left,right,'AbsTol',abserr);
    if (phixs == 0 || abs(phixs) < 1e-200)
        DIFF = 0; break
    end
    if (phixs ~= eps)
        left=right;right=newzero;x(count+1)=right;
        neval=neval+1;
        
        % update the partial sum
        M(1,count)=fxs/phixs;N(1,count)=1/phixs;
        fxs=fxs+phixs;
        
        if count == 2
        % Hard coded when count=2 since W has to be updated at least once.
            denom = [1;1./(1./x(1) - 1/x(2))];
            mVec = [M(1,2);M(1)];
            M([datamax+1 2]) = ([1 0; -1 1] * mVec) .* denom;
            nVec = [N(1,2);N(1)];
            N([datamax+1 2]) = ([1 0; -1 1] * nVec) .* denom;
        elseif count>2
            idx = count-1:-1:1;
            
            % indices of M and N used in updating W
            idx_mat = sub2ind(size(M),1:count-1, idx)'; 
            
            % Note that elements of denom has to be multiplied recursively
            % and cumulatively to select elments of M and N.
            denom =[1;1./(1./x(idx) - 1/x(count))];
            temp_a=(repmat(cumprod(denom),1,count));
            temp_b=(repmat(cumprod(denom)',count,1));
            temp_ab=tril([cumprod(denom) temp_a(:,2:end)./temp_b(:,1:end-1)]).*temp(1:count,1:count);
            
            mVec = [M(1,count);M(idx_mat)];
            M(sub2ind(size(M),1:count,count:-1:1)) = temp_ab * mVec ;
            nVec = [N(1,count);N(idx_mat)];
            N(sub2ind(size(N),1:count,count:-1:1)) = temp_ab * nVec ;
        end
                
        % derive the new estimate and calculate the error bound
        W(count)=M(count,1)/N(count,1);
        errorbound=max(abserr,relerr*abs(W(count)));
        
        % compare with the previous 2 estimates.
        if (count>2)
            DIFF=max(abs(W(count)-W(count-1)),abs(W(count)-W(count-2)));
        end
        if isnan(DIFF)
            error('Warning: Divergent partial sums for mW tranformation.');
        end
        count=count+1;
    else
        if count > 1
            W(count) = W(count-1);
            DIFF=max(abs(W(count)-W(count-1)),abs(W(count)-W(count)));
        else
            DIFF = 0;
        end
        break
    end
end

if ((count==datamax-2) && (DIFF > errorbound) )
    error('Warning: Possible non-convergence of mW transformation');
end

if count > 1
    result=W(count-1);
else
    result=W(count);
    if count ==1
        result=0;
    end
end
ierror = DIFF;

