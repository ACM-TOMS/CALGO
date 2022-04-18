% DEGREEELEVATION   DEGREEELEVATION(OBJ,K) implements k degrees elevation 
% of a polynomial represented in the Bernstein basis
function poly = degreeElevation(obj,k)
    n = obj.getDegree();
    coef = obj.BernsCoeff;
    for i=1:n+1
        coef(i) = coef(i)*nchoosek(n,i-1);
    end
    % Degree elevation is performed by using convolution
    coefUnity = zeros(1,k+1);
    for i=1:k+1
        coefUnity(i) = nchoosek(k,i-1);
    end
    coefDE = conv(coefUnity,coef);
    for i=1:n+k+1
        coefDE(i) = coefDE(i) / nchoosek(n+k,i-1);
    end
    poly = Polynomial(coefDE,'c');
end
