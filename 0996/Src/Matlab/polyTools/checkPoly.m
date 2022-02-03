function flag = checkPoly( poly )
%%CHECKPOLY    Check if the input is simplified sparsePOP format.
%   This function returns a logical value.

flag = true;

if ~isfield(poly,'supports') || ~isfield(poly,'coef')
    flag = false;
    warning('checkPoly: lack of fields');
    %return; 
end

if ~iscolumn(poly.coef)
    flag = false;
    warning('checkPoly: .coef must be a column vector');
    %return;
end

if size(poly.supports,1) ~= length(poly.coef)
    flag = false; 
    warning('checkPoly: dimensions must agree');
    %return; 
end

end

