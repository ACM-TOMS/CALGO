function s = power(s, n)

[row, col] = size(s);
if row~=1 && col~=1
    error('First argument should be an array.')
end

len = length(s);

for i=1:len
    s(i) = s(i)^n;
end

end