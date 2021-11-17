function [resul] = product_poly(poly1,poly2)
%Algorithm to perform a product between two polynomials.
%The present algorithm considers only homogenous polynomials.
%poly1.label   -> string with the name of the variable (eg. 'P')
%poly1.data(n) -> data correspondent to the n-th monomy
%poly1.data(n).exponent -> values of the exponents of the monomy (eg. [0 2])
%poly1.data(n).value    -> value of the monomy

%Second version: using hash tables
vertices = length(poly1.data(1).exponent);
degree = sum(poly1.data(1).exponent)+sum(poly2.data(1).exponent);
[tabexponents] = generate_homogeneous_exponents(vertices,degree);
if (degree>1)
    base = (degree*ones(1,vertices)).^(vertices-1:-1:0);
else
    base = (2*ones(1,vertices)).^(vertices-1:-1:0);
end
for cont=1:size(tabexponents,1)
    indhash(cont) = sum(base.*tabexponents(cont,:));
    resul.data(cont).value = 0;
end
resul.label = strcat(poly1.label,'*');
resul.label = strcat(resul.label,poly2.label);
for cont1=1:length(poly1.data)
    for cont2=1:length(poly2.data)
        exponent = poly1.data(cont1).exponent + poly2.data(cont2).exponent;
        
        hash = sum(base.*exponent);
        indresul = find(indhash==hash);
        
        resul.data(indresul).exponent = exponent;
        resul.data(indresul).value = resul.data(indresul).value + (poly1.data(cont1).value*poly2.data(cont2).value);
    end
end
resul.vertices = max(poly1.vertices,poly2.vertices);

return