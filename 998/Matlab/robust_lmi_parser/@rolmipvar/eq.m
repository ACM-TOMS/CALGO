function out = eq(A,B)
% Sets equality conditions on rolmipvar variables

out = [];
if ((isa(A,'rolmipvar')) && (isa(B,'rolmipvar')))
    if (length(A.data) ~= length(B.data))
        error('Both sides of the equality constraint must have the same number of monomials');
    else
        for cont=1:length(A.data)
            out = [out, constraint(A.data(cont).value,'==',B.data(cont).value)];
        end
    end
elseif (isa(A,'rolmipvar'))
    for cont=1:length(A.data)
        out = [out, constraint(A.data(cont).value,'==',B)];
    end
elseif (isa(B,'rolmipvar'))
    for cont=1:length(B.data)
        out = [out, constraint(A,'==',B.data(cont).value)];
    end
end
return