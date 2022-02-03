% ADD_MU_DER  Add derivatives with respect to mu
%
%   [Fd1, Fd2...] = ADD_MU_DER(r, Fd1, Fd2...) modifies the derivatives Fd1,
%   Fd2,... from being w.r.t. theta to being w.r.t. theta and mu, by 
%   increasing the 3rd dimension of each Fdj by r and filling the new elements
%   with zeros. An Fdi may also be a cell array, and then each cell is modified.

function [varargout] = add_mu_deriv(r, varargin)
  for j = 1:length(varargin)
    if iscell(varargin{j})
      for i=1:length(varargin{j})
        varargout{j}{i} = add_mu_deriv(r, varargin{j}{i}); 
      end
    else
      [m,n,k] = size(varargin{j});
      varargout{j} = cat(3, varargin{j}, zeros(m,n,r));
    end
  end
end
