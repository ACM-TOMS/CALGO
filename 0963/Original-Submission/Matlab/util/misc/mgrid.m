function g = mgrid(margin, d)
% MGRID Creates the d-dimensional grid with margin as 1-dimensional margin.
%
%  g = mgrid(margin, d) creates a d-dimensional grid by 'folding' margin
%    in d dimensions. The result is a d x length(margin)^d matrix representing
%    all grid points (as column vectors) of the set margin^d.
%
%    INPUT margin: margin of the grid to be created
%               d: dimension of the grid
%
% created by Benedikt Rudolph
% DATE: 05-Sep-2012
    
    m = length(margin);
    g = zeros(d, m^d); % pre-init grid
    for i=1:d % iterate over all dimensions
        g(i,:) = reshape( repmat(margin, m^(i-1), m^(d-i)), 1, m^d );
    end
end
