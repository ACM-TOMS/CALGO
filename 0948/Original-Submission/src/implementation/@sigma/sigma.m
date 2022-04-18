function  s = sigma(a)


if a > 0
    s.size = a;
    s.vector(1:a) = -Inf;
    s = class(s, 'sigma');
else
    error('No sigma object can be constructed. Input size must be positive integer!');
end

end