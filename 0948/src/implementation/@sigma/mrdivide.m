function s = mrdivide(t1, t2)

[m1, n1] = size(t1);
if ~isscalar(t2)
    error('Wrong size of second argument!');
end

if isa(t1, 'numeric') && isa(t2, 'sigma') % t1 is sigma object
    s(1:m1,1:n1) = t2;
elseif isa(t1, 'sigma') && isa(t2, 'numeric') % t2 is sigma object
    s = t1;
elseif isa(t1, 'sigma') && isa(t2, 'sigma') % both are sigma object
    s = t1;
    num = m1*n1;
    for i=1:num
        if s(i).size~=t2.size, error('Size mismatch!');end
        s(i).vector = max(t2.vector, s(i).vector);
    end
else
    error('Wrong class of objects.');
end

end