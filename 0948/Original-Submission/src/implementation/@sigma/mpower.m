function s = mpower(s,t)
if isa(t, 'numeric') && isscalar(t)
else
    error('The second argument must be a numeric scalar!');
end
end