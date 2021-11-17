function s = plus(s, t)
if isa(s, 'numeric') && isa(t, 'sigma') % t is qla object
    s = t;
elseif isa(s, 'sigma') && isa(t, 'sigma') % both are qla objects
    s.vector = max(s.vector, t.vector);
end
end