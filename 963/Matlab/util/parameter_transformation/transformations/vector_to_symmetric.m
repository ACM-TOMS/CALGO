% This program has been developed in a research project at Risklab, Toronto.
% (c) Benedikt Rudolph, 2012

function M = vector_to_symmetric(v, d)
    % Transforms the vector v into the symmetric d x d matrix M
    % by interpreting v as the stacked columns of the lower triangular.
    %
    % It is expected that 0.5*d*(d+1) == length(v),
    % otherwise an error will occur.
    %
    % M = vector_to_symmetric(v, d)
    
    M = zeros(d); % init M
    M(logical(tril(ones(d)))) = v; % set lower triangular to v
    M = M + tril(M,-1)'; % add upper triangular symmetrically
end
