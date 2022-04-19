% DEGREEREDUCTION     DEGREEREDUCTION(OBJ) return a polynomial object
% represented in the Bernstein basis of the true degree of the
% polynomial obj (it is no obvious the true degree of a polynomial
% represented in the Bernstein basis)
function poly = degreeReduction(obj)
    n = obj.getDegree();
    coeff = obj.getCoeff();
    differences = coeff;
    % Compute finite differences
    for r=1:n
        differences(r+1:end) = differences(r+1:end)-differences(r:end-1);
    end
    % Compute true degree
    degree = n;
    while abs(differences(degree+1))<30*eps
        degree=degree-1;
    end
    % Reduction of the degree
    if degree==n
        poly = obj;
    else
        coeffRed = zeros(1,degree+1);
        r=n-degree;
        for k=1:degree+1
            for j=1:k
                coeffRed(k) = coeffRed(k)+(-1)^(k-j)*nchoosek(k-j+r-1,r-1)*nchoosek(n,j-1)*coeff(j);
            end
            coeffRed(k) = coeffRed(k)/nchoosek(degree,k-1);
        end
        poly = Polynomial(coeffRed,'c');
    end
end

