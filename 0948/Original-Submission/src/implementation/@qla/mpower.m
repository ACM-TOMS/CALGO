function s = mpower(s, n)

num = numel(s);
if n==1
    return;
elseif n~=0
    for i=1:num
        index = s(i).type==2;
        s(i).type(index) = 3;
    end
end
end