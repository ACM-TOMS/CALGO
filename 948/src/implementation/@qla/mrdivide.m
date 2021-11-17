function s = mrdivide(t1, t2)

[m1, n1] = size(t1);
if ~isscalar(t2)
    error('Wrong size of second argument!');
end

if isa(t1, 'numeric') && isa(t2, 'qla') % t1 is qla object
    s(1:m1,1:n1) = t2;
elseif isa(t1, 'qla') && isa(t2, 'numeric') % t2 is qla object
    s = t1;
elseif isa(t1, 'qla') && isa(t2, 'qla') % both are qla object
    s = t1;
    num = m1*n1;
    for i=1:num
        if s(i).size~=t2.size, error('Size mismatch!');end
        s(i).offset = min(s(i).offset, t2.offset);
        temp = max([s(i).type; t2.type]);
        indx = s(i).type==2 & t2.type==2;
        temp(indx) = 3;
        s(i).type = temp;
    end
else
    error('Wrong class of objects.');
end

end