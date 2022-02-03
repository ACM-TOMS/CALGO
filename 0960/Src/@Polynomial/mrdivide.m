%MRDIVIDE  MRDIVIDE(OBJ1,OBJ2) implements obj1/obj2
function [q,r] = mrdivide(obj1,obj2)
    poly1 = obj1.degreeReduction();
    poly2 = obj2.degreeReduction();
    m = poly1.getDegree();
    n = poly2.getDegree();
    if m<n
        q = Polynomial();
        r = poly1;
    else
        A = zeros(m+1);
        coeff1 = poly1.getCoeff();
        coeff2 = poly2.getCoeff();
        for i=0:m
            for j=max(0,i-n):min(m-n,i)
                A(i+1,j+1) = ((nchoosek(m-n,j)*nchoosek(n,i-j))/nchoosek(m,i))*coeff2(i-j+1);
            end
            for j=max(0,i-m+n-1):min(n-1,i)
                A(i+1,m-n+j+2) = ((nchoosek(m-n+1,i-j)*nchoosek(n-1,j))/nchoosek(m,i));
            end
        end
        sol = A\(coeff1');
        q = Polynomial(sol(1:m-n+1),'c');
        r = Polynomial(sol(m-n+2:end),'c');
    end
end

