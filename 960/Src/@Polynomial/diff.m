% DIFF  DIFF(OBJ) is the derivative of the polynomial obj.
function r = diff(obj)
    c = obj.getCoeff();
    n = obj.getDegree();
    r = Polynomial(n*(c(2:end)-c(1:end-1)),'c');
end % diff

