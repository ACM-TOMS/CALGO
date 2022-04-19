function [exponents] = generate_homogeneous_exponents(vertices,degree)
% [exponents] = generate_homogeneous_exponents(vertices,degree)
% Recursive algorithm to generate the exponents of a homogeneous polynomial
exponents = [];
if (vertices == 1)
    exponents = degree;
    return
end
if (degree == 0)
    exponents = zeros(1,vertices);
    return
end
for cont=degree:-1:0
    [expaux] = generate_homogeneous_exponents(vertices-1,degree-cont);
    exponents = [exponents; cont*ones(size(expaux,1),1) expaux];
end
return