function [varargout] = addsimplex(varargin)
% Procedure to insert a new multisimplex variable in an already existent
% polynomial variable
%
% [polyout] = addsimplex(polyin,vertices)
% Inserts, in the variable poly, the multisimplex with the given number of
% vertices, which corresponds to a row vector.

if (nargin ~= 2)
    error('Input error. Type ''help addsimplex'' for more details');
    return
else
    varargout{1} = insert_simplex(varargin{1},varargin{2});
end

return