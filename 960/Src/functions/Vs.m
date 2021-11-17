% VS  [Y,ERRBOUND]=VS(COEFF,X) evaluates a polynomial curve represented in 
%            the Bernstein basis with the coefficients COEFF at the points 
%            in vector X by the VS algorithm
% Input arguments:
%        - coeff is the vector of coefficients of the polynomial to be 
%          evaluated respect to the Bernstein basis
%        - x is the vector of points where the polynomial must be evaluated
% Output arguments:
%        - y provides the values obtained for the evaluation of obj at the
%           points in x
%        - errBound provides upper bounds for the relative errors of the
%           corresponding evaluations. If for some point in x, the 
%           corresponding element in errBound is NaN it means that algoritm 
%           has not been able to obtain a reliable bound  
function [y,errBound] = Vs(coeff,x)
    n = length(coeff)-1;
    
    [rowsNo,columnsNo] = size(x);
    if rowsNo~=1 && columnsNo~=1
        disp('x must be a row or a column vector');
        return;
    else
        y = zeros(size(x));
        errBound = zeros(size(x));
    end

    ind1 = find(x<0.5);
    if ~isempty(ind1)
        x1 = x(ind1);
        y1 = y(ind1);
        bound1 = errBound(ind1);
        t1 = x1 ./ (1-x1);
    end

    ind2 = find(x>=0.5);
    if ~isempty(ind2)
        x2 = x(ind2);
        y2 = y(ind2);
        bound2 = errBound(ind2);
        t2 = (1-x2) ./ x2;
    end

    % Evaluation at points <0.5
    if ~isempty(ind1)
        y1 = coeff(n+1)*t1+coeff(n);
        bound1 = 2*abs(coeff(n+1))*abs(t1)+3*abs(y1);
        for i=1:n-1
            y1 = y1.*t1+coeff(n-i);
            bound1 = bound1.*abs(t1)+3*abs(y1);
        end
        [pow1,k1] = powerDC(1-x1,n);
        %pow1= (1-x1).^n;k1=n-1;
        aux1 = y1;
        y1 = y1 .* pow1;
        bound1 = (pow1.*(bound1-2*abs(aux1))+k1*abs(aux1).*abs(pow1)+abs(y1))*(eps/2);
    end
    % Evaluation at points >=0.5
    if ~isempty(ind2)
        y2 = coeff(1)*t2+coeff(2);
        bound2 = 2*abs(coeff(1))*abs(t2)+3*abs(y2);
        for i=2:n
            y2 = y2.*t2+coeff(i+1);
            bound2 = bound2.*abs(t2)+3*abs(y2);
        end
        [pow2,k2] = powerDC(x2,n);
        %pow2=x2.^n;k2=n-1;        
        aux2 = y2;
        y2 = y2 .* pow2;
        bound2 = (pow2.*(bound2-2*abs(aux2))+k2*abs(aux2).*abs(pow2)+abs(y2))*(eps/2);
    end

    if ~isempty(ind1)
        y(ind1) = y1;
        errBound(ind1) = bound1;
    end
    if ~isempty(ind2)
        y(ind2) = y2;
        errBound(ind2) = bound2;
    end

    ind1 = find(abs(y)<=errBound);
    ind2 = find(abs(y)>errBound);
    errBound(ind1) = NaN;
    errBound(ind2) = errBound(ind2)./abs(y(ind2));
end
