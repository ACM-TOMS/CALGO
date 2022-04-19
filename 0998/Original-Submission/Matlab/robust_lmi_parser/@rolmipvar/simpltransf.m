function out = simpltransf(poly,simplexin,simplexout,T)

% Compute the polynomial resultant from a linear transformation on simplex
% variables.
%
% [out] = simpltransf(poly,simplexin,simplexout,T) returns the polynomial
% out when the input polynomial poly, dependent on simplex variables of the
% simplexin domain, is transformed to be dependent on simplex variables of
% the simplexout domain, being such transformation given by matrix T.
%
% For instance, let poly = a1*P1 + a2*P2, being a1 and a2 variables of
% simplex 1, and suppose that such variables are transformed to another
% simplex 2, with variables b1, b2 and b3, through the relation
% a1 = b1 + b2 - 3*b3, a2 = 2*b1 + 5*b2 - 2*b3; then 
% [out] = simpltransf(poly, 1, 2, [1 1 -3; 2 5 -2]).

if ~isa(poly, 'rolmipvar')
    error('The first argument must be a rolmipvar variable');
end
if (nargin ~= 4)
    error('Wrong input; please refer to help simpltransf');
end

for cont=1:size(T,1)
    %Create each polynomial
    a{cont} = rolmipvar(T(cont,:),'aux',size(T,2),1);
    a{cont} = fork(a{cont},'aux',1,simplexout);
end

deg = sum(poly.data(1).exponent{simplexin});
outini = poly;
[outini,outdeg] = removesimplex(outini,simplexin);
%out = resetpoly(out);
aux = rolmipvar(ones(1,size(T,2)),'aux',size(T,2),1);
um = aux;
for cont=1:deg-1
    aux = aux*um;
end
aux = fork(aux,'aux',1,simplexout);
out = 0;
%out = out*aux;


for cont=1:length(poly.data)
    expaux = poly.data(cont).exponent{simplexin};
    polaux = 1;
    for contexp=1:length(expaux)
        temp = expaux(contexp);
        while (temp >= 1)
            polaux = polaux*a{contexp};
            temp = temp - 1;
        end
    end


    %Isolates a monomial
    for contmon=1:length(outini.data(cont).exponent)
        monomi{1}{contmon} = outini.data(cont).exponent{contmon};
    end
    monomi{1}{contmon+1} = outini.data(cont).value;
    monom = rolmipvar(monomi,'aux',outini.vertices,outdeg);
    
    %Multiplies the monomial by the transformated simplex
    out = out + monom*polaux;

    
end

out.label = poly.label;



return

function [out,deg] = removesimplex(poly,simplex)
% Internal function to set all the exponents of the simplex variables
% ''simplex'' to zero, which removes the simplex.

for cont=1:length(poly.data)
    poly.data(cont).exponent{simplex} = zeros(size(poly.data(cont).exponent{simplex}));
end
out = poly;

for cont=1:length(out.data(1).exponent)
    deg(cont) = sum(out.data(1).exponent{cont});
end
return
