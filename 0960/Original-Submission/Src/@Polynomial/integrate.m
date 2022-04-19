% INTEGRATE  INTEGRATE(OBJ) provides the indefinite integral of the
%            polynomial obj like an object of the class Polynomial
function r = integrate(obj)
    coeff = obj.getCoeff();
    n = obj.getDegree();
    coeffInt = zeros(1,n+2);
    for i=1:n+1
        coeffInt(i+1) = (1/(n+1))*sum(coeff(1:i));
    end
    r = Polynomial(coeffInt,'c');
end % integrate

