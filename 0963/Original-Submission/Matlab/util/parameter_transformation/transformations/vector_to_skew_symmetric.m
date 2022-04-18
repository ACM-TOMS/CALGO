function M = vector_to_skew_symmetric(v, d)
%VECTOR_TO_SKEW_SYMMETRIC Transforms vector v into skew symmetric matrix M.
%
%  M = vector_to_symmetric(v, d) Transforms the vector v into the skew symmetric
%    d x d matrix M by interpreting v as the stacked columns of the lower triangular.
%    It is expected that 0.5*d*(d+1) == length(v), otherwise an error will occur.
%
% created by Benedikt Rudolph
% DATE: 20-Aug-2012
  
  M = zeros(d); % init M
  M(logical(tril(ones(d),-1))) = v; % set lower triangular to v
  M = M - tril(M,-1)'; % subtract upper triangular symmetrically
  M(logical(eye(d))) = 0; % set diagonal elements to zero
end
