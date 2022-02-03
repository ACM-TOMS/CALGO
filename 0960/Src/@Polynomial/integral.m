% INTEGRAL  INTEGRAL(OBJ) provides the definite integral of the
%           polynomial obj in the interval [0,1]
function q = integral(obj)
    coeff = obj.getCoeff();
    n = obj.getDegree();
    q = (1/(n+1))*sum(coeff);
end

