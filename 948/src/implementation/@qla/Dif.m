function t = Dif(t, d)

if nargin == 1, d = 1; end

t.offset = t.offset - d;
if ~all(t.offset)
    index = t.type==1;
    t.type(index) = 2;
end
end
