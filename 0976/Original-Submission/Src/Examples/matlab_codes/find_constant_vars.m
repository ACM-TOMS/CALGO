% find any variables which have zero value uniformly.  should be
% improved to find constant variables.
function ind = find_constant_vars(data)
%
zerothresh = 1e-7;


biggies = abs(max(data)-min(data));
ind = find(biggies>zerothresh);

end

