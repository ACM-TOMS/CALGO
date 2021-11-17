function n = getSize(s)

if isa(s, 'qla')
    n = s.size;
else
    error('Input for getType must be a qla object!');
end

end