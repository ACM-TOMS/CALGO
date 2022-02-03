% PLUS    PLUS(OBJ1,OBJ2) implements obj1+obj2 for Poylnomials
function r = plus(obj1,obj2)
    % PLUS  implements obj1 + obj2 for Polynomials
    obj1 = Polynomial(obj1);
    obj2 = Polynomial(obj2);
    k = obj2.getDegree() - obj1.getDegree();
    % Increase degree
    if k>0
        obj1 = obj1.degreeElevation(k);
    elseif k<0
        obj2 = obj2.degreeElevation(-k);
    end
    r = Polynomial(obj1.BernsCoeff+obj2.BernsCoeff,'c');
end % plus

