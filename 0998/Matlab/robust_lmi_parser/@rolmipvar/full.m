function out = full(X)

out = X;
for cont=1:length(X.data)
    out.data(cont).value = full(X.data(cont).value);
end

return