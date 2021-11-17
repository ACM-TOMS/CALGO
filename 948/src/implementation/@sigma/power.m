function s = power(s, n)

num = numel(s);

for i=1:num
    s(i) = s(i)^n;
end

end