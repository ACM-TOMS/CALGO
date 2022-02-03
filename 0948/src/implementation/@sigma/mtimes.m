function s = mtimes(s, t)

if isa(s, 'numeric') && isa(t, 'sigma') % t is sigma object
    s = t;
elseif isa(s, 'sigma') && isa(t, 'sigma') % both are sigma object
    s.vector = max(s.vector, t.vector);
end
end