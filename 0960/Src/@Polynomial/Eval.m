% EVAL  Y=EVAL(OBJ,X,PREC) evaluates polynomial OBJ at the points in x with
%       a pretended relative precision given by prec.
% Input arguments:
%        - obj is the polynomial to be evaluated
%        - x is the vector of points where the polynomial must be evaluated
%        - prec is an upper bound goal for the relative error of the
%           evaluations. In the case that prec is not provided the
%           algorithm takes prec as 1e-12
% Output argument:
%        - y is a matrix providing the values obtained for the evaluation 
%          of obj at the points in x, upper bounds of the relative errors 
%          for each one of these evaluations (if for some point in x, the 
%          corresponding element in errBound is NaN it means that algoritm 
%          has not been able to obtain a reliable bound) and a sequence of 
%          flags (flag=1 means that the relative error of the corresponding 
%          evaluation is lower than prec and a flag=0 means that precision 
%          requirement either is not satisfied or is not known if it is 
%          satisfied). If the number of elements in x is n then y is a 3 by 
%          n matrix (y=[evaluation; error bound; flags]) in the case that x
%          is a row vector, whereas y is a n by 3 matrix (y=[evaluation,
%          error bound, flags]) in the case that x is a column vector or an
%          scalar
function y = Eval(obj,x,varargin)
    % If argument tol it is no provided it is fixed to 1e-12
    if nargin < 3
        prec = 1e-12;
    else
        prec = varargin{1};
    end
    
    sizeX = size(x);
    if sizeX(1)~=1 && sizeX(2)~=1
        error('Polynomial.Eval: x must be a vector');
    end
    res = zeros(sizeX);
    errBound = zeros(sizeX);
    flags = zeros(sizeX);

    coeff = obj.getCoeff();
    n = obj.getDegree();
    coeffVS = zeros(size(coeff));
    for i=1:n+1
        coeffVS(i) = nchoosek(n,i-1)*coeff(i);
    end
    
    % Now we evaluate the polynomial by the VS algorithm
    [res,errBound] = Vs(coeffVS,x);
    flags = errBound<prec;
    ind = find(flags==0);
    if length(ind)>0
        % At the points where evaluation is not accurate enough
        % the polynomial is evaluated by other more accurate
        % algorithms                
        if n<=32
            [res(ind),errBound(ind)] = Casteljau(coeff,x(ind));
            flags(ind) = errBound(ind)<prec;
            ind = find(flags==0);
            if length(ind)>0
                % At the points where evaluation is not accurate
                % enough the polynomial is evaluated by the
                % compensated VS algorithm
                [res(ind),errBound(ind)] = CompVs(coeffVS,x(ind));
                flags(ind) = errBound(ind)<prec;
            end
        else
            [res(ind),errBound(ind)] = Vs(coeffVS,x(ind));
            flags(ind) = errBound(ind)<prec;
            ind = find(flags==0);
            if length(ind)>0
                % At the points where evaluation is not accurate
                % enough the polynomial is evaluated by the
                % compensated VS algorithm
                [res(ind),errBound(ind)] = CompVS(coeffVS,x(ind));
                flags(ind) = errBound(ind)<prec;
            end
        end
    end
    
    if sizeX(2) == 1
        y=[res,errBound,flags];
    else
        y=[res;errBound;flags];
    end
end
