function  s = qla(a)

if a > 0
    s.size = a;
    s.offset = inf(1,a);
    s.type = zeros(1, a);
    
    s = class(s, 'qla');
else
    error('No qla object can be constructed. Input size must be positive integer!');
end
end