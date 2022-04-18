% CASTELJAU  CASTELJAU(COEFF,X) evaluates a polynomial curve represented in 
%            the Bernstein basis with the coefficients COEFF at the points 
%            in vector X by the de Casteljau algorithm
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
function [y,errBound] = Casteljau(coeff,x)
    n = length(coeff)-1;

    [rowsNo,columnsNo] = size(x);
    if rowsNo~=1 && columnsNo~=1
        disp('x must be a row or a column vector');
        return;
    elseif rowsNo==1
        dim = 2;
        eval = zeros(n,columnsNo);
        bound = zeros(n+1,columnsNo);
    else
        dim = 1;
        eval = zeros(rowsNo,n);
        bound = zeros(rowsNo,n+1);
    end



    switch dim
        case 1
            % First step of the de Casteljau algorithm (r=1)
            for i=1:n
                eval(:,i) = (1-x)*coeff(i) + x*coeff(i+1);                
                bound(:,i) = (1-x).*bound(:,i)+x.*bound(:,i+1)+(1-x).*abs(coeff(i))+x.*abs(coeff(i+1))+abs(eval(:,i));
            end
            % Steps 2 to n of de Casteljau algorithm (r=2:n)
            for r=2:n
                for i=1:n+1-r
                    aux = eval(:,i);
                    eval(:,i) = (1-x).*eval(:,i) + x.*eval(:,i+1);
                    bound(:,i) = (1-x).*bound(:,i) + x.*bound(:,i+1) + (1-x).*abs(aux)+x.*abs(eval(:,i+1))+abs(eval(:,i));
                end
            end
            y = eval(:,1);
            errBound = bound(:,1)*eps/2;
        case 2
            % First step of the de Casteljau algorithm (r=1)
            for i=1:n
                eval(i,:) = (1-x)*coeff(i) + x*coeff(i+1);
                bound(i,:) = (1-x).*bound(i,:)+x.*bound(i+1,:)+(1-x).*abs(coeff(i))+x.*abs(coeff(i+1))+abs(eval(i,:));
            end
            % Steps 2 to n of de Casteljau algorithm (r=2:n)
            for r=2:n
                for i=1:n+1-r
                    aux = eval(i,:);
                    eval(i,:) = (1-x).*eval(i,:) + x.*eval(i+1,:);
                    bound(i,:) = (1-x).*bound(i,:) + x.*bound(i+1,:) + (1-x).*abs(aux)+x.*abs(eval(i+1,:))+abs(eval(i,:));
                end
            end
            y = eval(1,:);
            errBound = bound(1,:)*eps/2;
    end

    ind1 = find(abs(y)<=errBound);
    ind2 = find(abs(y)>errBound);
    errBound(ind1) = NaN;
    errBound(ind2) = errBound(ind2) ./ abs(y(ind2));
end

