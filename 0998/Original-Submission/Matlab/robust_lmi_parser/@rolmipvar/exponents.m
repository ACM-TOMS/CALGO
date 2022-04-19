function c = exponents(X)
%  Retrieve the exponents of the polynomial struct
%  to get more flexibility.
%
    if isa(X, 'rolmipvar')
        for i = 1:length(X.data)
            c{i} = X.data(i).exponent;
        end
    else
        c = {X};
    end
end