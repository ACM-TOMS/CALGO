% COMPVS  COMPVS(COEFF,X) evaluates a polynomial curve represented in 
%            the Bernstein basis with the coefficients COEFF at the points 
%            in vector X by the compensated VS algorithm
% Input arguments:
%        - coeff is the vector of coefficients of the polynomial to be 
%          evaluated respect to the VS basis
%        - x is the vector of points where the polynomial must be evaluated
% Output arguments:
%        - y provides the values obtained for the evaluation of obj at the
%           points in x 
%        - errBound provides upper bounds for the relative errors of the
%           corresponding evaluations. If for some point in x, the 
%           corresponding element in errBound is NaN it means that algoritm 
%           has not been able to obtain a reliable bound  
function [y,errBound] = CompVs(coeff,x)
    n = length(coeff)-1;
    for i=0:n
        coeffAbs(i+1)=abs(coeff(i+1));
    end
    u=eps/2;

    [rowsNo,columnsNo] = size(x);
    if rowsNo~=1 && columnsNo~=1
        disp('x must be a row or a column vector');
        return;
    elseif rowsNo==1
        dim = 2;
        y = ones(1,columnsNo);
        bound = ones(1,columnsNo);
        alpha = zeros(n,columnsNo);
        p = zeros(n,columnsNo);
    else
        dim = 1;
        y = ones(rowsNo,1);
        bound = ones(rowsNo,1);
        alpha = zeros(rowsNo,n);
        p = zeros(rowsNo,n);
    end

    % Compute r=1-x
    [r,rho] = TwoSum(1,-x);

    % Initialization for points < 0.5
    ind1 = find(x<0.5);
    if ~isempty(ind1)
        x1 = x(ind1);
        r1 = r(ind1);
        [coc1,beta1] = DivRem(x1,r1);
        y1 = y(ind1);
        rho1 = rho(ind1);
        pi1 = zeros(size(y1));
        sigma1 = zeros(size(y1));
        if dim==1
            p1 = p(ind1,:);
        else
            p1 = p(:,ind1);
        end
        bound1 = bound(ind1);
    end

    % Initialization for points >= 0.5
    ind2 = find(x>=0.5);
    if ~isempty(ind2)
        x2 = x(ind2);
        r2 = r(ind2);
        [coc2,beta2] = DivRem(r2,x2);
        y2 = y(ind2);
        rho2 = rho(ind2);
        pi2 = zeros(size(y2));
        sigma2 = zeros(size(y2));
        if dim==1
            p2 = p(ind2,:);
        else
            p2 = p(:,ind2);
        end
        bound2 = bound(ind2);
    end

    % Evaluation at points <0.5
    if ~isempty(ind1)
        [s,pi1] = TwoProduct(coc1,coeff(n+1));
        [y1,sigma1] = TwoSum(s,coeff(n));
        if dim==1
            p1(:,n) = ((beta1-rho1.*coc1)./r1) * coeff(n+1) + pi1 + sigma1;
        else
            p1(n,:) = ((beta1-rho1.*coc1)./r1) * coeff(n+1) + pi1 + sigma1;
        end
        bound1 = coc1*coeffAbs(n+1)+coeffAbs(n);
        for i=1:n-1
            [s,pi1] = TwoProduct(coc1,y1);
            aux_y1 = y1;
            [y1,sigma1] = TwoSum(s,coeff(n-i));
            if dim==1
                p1(:,n-i) = ((beta1-rho1.*coc1)./r1) .* aux_y1 + pi1 + sigma1;
            else
                p1(n-i,:) = ((beta1-rho1.*coc1)./r1) .* aux_y1 + pi1 + sigma1;
            end
            bound1 = coc1.*bound1 + coeffAbs(n-i);
        end
       
        switch dim
            case 1
                for i=1:n
                    [y1,alpha(ind1,i)] = TwoProduct(y1,r1);                    
                end
                vs1 = zeros(length(ind1),1);
                horner1 = zeros(length(ind1),1);
                for i=1:length(ind1)
                    vs1(i) = Vs(p1(i,:),x1(i));
                    horner1(i) = Horner(fliplr(alpha(ind1(i),:)),n-1,r1(i));                
                end                  
            case 2
                for i=1:n
                    [y1,alpha(i,ind1)] = TwoProduct(y1,r1);                    
                end
                vs1 = zeros(1,length(ind1));
                horner1 = zeros(1,length(ind1));               
                for i=1:length(ind1)
                    vs1(i) = Vs(p1(:,i)',x1(i));
                    horner1(i) = Horner((flipud(alpha(:,ind1(i))))',n-1,r1(i));
                end
        end           
        y(ind1) = y1+r1.*vs1+horner1;       
        bound1 = bound1.*r1.^n;
        bound(ind1) = 2*u*abs(y(ind1)) + 64*n^2*u^2*bound1;
    end
    % Evaluation at points >=0.5   
    if ~isempty(ind2)
        [s,pi2] = TwoProduct(coc2,coeff(1));
        [y2,sigma2] = TwoSum(s,coeff(2));
        if dim==1
            p2(:,1) = ((rho2+beta2)./x2) * coeff(1) + pi2 + sigma2;
        else
            p2(1,:) = ((rho2+beta2)./x2) * coeff(1) + pi2 + sigma2;
        end
        bound2 = coc2*coeffAbs(1)+coeffAbs(2);
        for i=2:n
            [s,pi2] = TwoProduct(coc2,y2);
            aux_y2 = y2;
            [y2,sigma2] = TwoSum(s,coeff(i+1));
            if dim==1
                p2(:,i) = ((rho2+beta2)./x2) .* aux_y2 + pi2 + sigma2;
            else
                p2(i,:) = ((rho2+beta2)./x2) .* aux_y2 + pi2 + sigma2;
            end
            bound2 = coc2.*bound2 + coeffAbs(i+1);
        end

        switch dim
            case 1
                for i=1:n
                    [y2,alpha2(ind2,n+1-i)] = TwoProduct(y2,x2);
                end
                vs2 = zeros(length(ind2),1);
                horner2 = zeros(length(ind2),1);
                for i=1:length(ind2)                  
                    vs2(i) = Vs(p2(i,:),x2(i)); 
                    horner2(i) = Horner(alpha(ind2(i),:),n-1,x2(i));                
                end                
            case 2
                for i=1:n
                    [y2,alpha2(n+1-i,ind2)] = TwoProduct(y2,x2);
                end
                vs2 = zeros(1,length(ind2));
                horner2 = zeros(1,length(ind2));
                for i=1:length(ind2)
                    vs2(i) = Vs((p2(:,i))',x2(i)); 
                    horner2(i) = Horner(alpha(:,ind2(i))',n-1,x2(i));
                end                
        end
        y(ind2) = y2+x2.*vs2+horner2;
        bound2 = bound2.*x2.^n;
        bound(ind2) = 2*u*abs(y(ind2)) + 64*n^2*u^2*bound2;
    end
    ind = find(abs(y)<=bound);
    errBound(ind) = NaN;
    ind = find(abs(y)>bound);
    errBound(ind) = bound(ind) ./ abs(y(ind)); 
end 
