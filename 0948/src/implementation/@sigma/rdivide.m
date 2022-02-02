function s = rdivide(t1, t2)

[m1, n1] = size(t1);
[m2, n2] = size(t2);

if m1==m2 && n1==n2 % both matrices
    if isa(t1, 'numeric') && isa(t2, 'sigma') % t2 is sigma object
        s(1:m2,1:n2) = t2;
        
    elseif isa(t1, 'sigma') && isa(t2, 'numeric') % t1 is sigma object
        s(1:m1,1:n1) = t1;
        
    elseif isa(t1, 'sigma') && isa(t2, 'sigma') % both are sigma object
        s(1:m1,1:n1) = t1;
        num = m1*n1;
        for i=1:num
            if s(i).size~=t2(i).size, error('Size mismatch!');end
            s(i).vector = max(t2(i).vector, s(i).vector);
        end
    else
        error('Wrong class of objects.');
    end
    
elseif isscalar(t1) % t1 is scalar, t2 is matrix
    if isa(t1, 'numeric') && isa(t2, 'sigma') % t2 is sigma object
        s = t2;
    elseif isa(t1, 'sigma') && isa(t2, 'numeric') % t1 is sigma object
        s(1:m2,1:n2) = t1;
    elseif isa(t1, 'sigma') && isa(t2, 'sigma') % both are sigma object
        s = t2;
        num = m2*n2;
        for i=1:num
            if s(i).size~=t1.size, error('Size mismatch!');end
            s(i).vector = max(t1.vector, s(i).vector);
        end
    else
        error('Wrong class of objects.');
    end
    
elseif isscalar(t2) % t1 is matrix, t2 is scalar
    if isa(t1, 'numeric') && isa(t2, 'sigma') % t2 is sigma object
        s(1:m1,1:n1) = t2;
    elseif isa(t1, 'sigma') && isa(t2, 'numeric') % t1 is sigma object
        s = t1;
    elseif isa(t1, 'sigma') && isa(t2, 'sigma') % both are sigma object
        s = t1;
        num = m1*n1;
        for i=1:num
            if s(i).size~=t2.size, error('Size mismatch!');end
            s(i).vector = max(t2.vector, s(i).vector);
        end
    else
        error('Matrix dimensions must agree.');
    end
else
    error('Matrix dimensions must agree.');
end
end