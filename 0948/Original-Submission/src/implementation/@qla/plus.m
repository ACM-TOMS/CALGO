function s = plus(s, t)
if isa(s, 'numeric') && isa(t, 'qla') % t is qla object
    s = t;
elseif isa(s, 'qla') && isa(t, 'qla') % both are qla objects
    s.offset = min(s.offset, t.offset);
    s.type = max(s.type,t.type);
end
end