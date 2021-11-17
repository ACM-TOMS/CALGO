% MTIMES   MTIMES(OBJ1,OBJ2) implements obj1 * obj2 for Polynomials
function r = mtimes(obj1,obj2)
    obj1 = Polynomial(obj1);
    obj2 = Polynomial(obj2);
    m = obj1.getDegree();
    n = obj2.getDegree();
    for i=1:m+1
        obj1.BernsCoeff(i) = obj1.BernsCoeff(i) * nchoosek(m,i-1);
    end
    for i=1:n+1
        obj2.BernsCoeff(i) = obj2.BernsCoeff(i) * nchoosek(n,i-1);
    end
    % Product is performed by using convolution
    coefProduct = conv(obj1.BernsCoeff,obj2.BernsCoeff);
    for i=1:n+m+1
        coefProduct(i) = coefProduct(i) / nchoosek(n+m,i-1);
    end
    r = Polynomial(coefProduct,'c');
end % mtimes
