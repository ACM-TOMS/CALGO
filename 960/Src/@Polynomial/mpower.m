% MPOWER  MPOWER(OBJ,K) implements obj^k
function r = mpower(obj,k)
    coef = obj.BernsCoeff;
    n = obj.getDegree();
    for i=1:n+1
        coef(i) = coef(i) * nchoosek(n,i-1);
    end
    coefPower = coef;
    for i=1:k-1
        coefPower = conv(coefPower,coef);
    end
    for i=1:k*n+1
        coefPower(i) = coefPower(i) / nchoosek(n*k,i-1);
    end
    r = Polynomial(coefPower,'c');
end

