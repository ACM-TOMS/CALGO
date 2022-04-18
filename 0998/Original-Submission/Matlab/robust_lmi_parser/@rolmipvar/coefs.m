function c = coefs(X)
%  Retrieve the coefficients of the polynomial struct
%  to get more flexibility.
%
if isa(X, 'rolmipvar')
    if (length(X.data) == 1)
        c = X.data(1).value;
    else
        for i = 1:length(X.data)
            c{i} = X.data(i).value;
        end
    end
else
    c = {X};
end
return
