function A = tensor_as_matrix(varargin)
%TENSOR_AS_MATRIX Constructor for matrix representation of a tensor
%
%   TENSOR_AS_MATRIX(T, RDIMS) creates a matrix representation of a
%   tensor T by rearranging and reshaping T.  The dimensions (or modes)
%   specified in RDIMS map to the rows of the matrix, and the
%   remaining dimensions (in ascending order) map to the columns.
%
%   TENSOR_AS_MATRIX(T, RDIMS, CDIMS) creates a matrix representation
%   of tensor T.  The dimensions specified in RDIMS map to the rows of
%   the matrix, and the dimensions specified in CDIMS map to the
%   columns, in the order given.
%
%   TENSOR_AS_MATRIX(T, RDIM, STR) creates the same matrix
%   representation as above, except only one dimension in RDIM maps to
%   the rows of the matrix, and the remaining dimensions span the
%   columns in an order specified by the string argument STR as
%   follows:
%
%     'fc' - Forward cyclic.  Order the remaining dimensions in the
%            columns by [RDIM+1:ndims(T), 1:RDIM-1].  This is the
%            ordering defined by Kiers.
%
%     'bc' - Backward cyclic.  Order the remaining dimensions in the
%            columns by [RDIM-1:-1:1, ndims(T):-1:RDIM+1].  This is the
%            ordering defined by De Lathauwer, De Moor, and Vandewalle.
%
%   TENSOR_AS_MATRIX(A, RDIMS, CDIMS, TSIZE) creates a TENSOR_AS_MATRIX
%   from a matrix A along with the mappings of the row (RDIMS) and column
%   indices (CDIMS) and the size of the original tensor (TSIZE). 
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR.

%Brett W. Bader and Tamara G. Kolda, Released under SAND2004-5189,
%Sandia National Laboratories, 2004.  Please address questions or
%comments to: tgkolda@sandia.gov.  Terms of use: You are free to copy,
%distribute, display, and use this work, under the following
%conditions. (1) You must give the original authors credit. (2) You may
%not use or redistribute this work for commercial purposes. (3) You may
%not alter, transform, or build upon this work. (4) For any reuse or
%distribution, you must make clear to others the license terms of this
%work. (5) Any of these conditions can be waived if you get permission
%from the authors.

if (nargin < 2)  ||  (nargin > 4)
  error('Incorrect number of arguments.');
end

T = varargin{1};
rdims = varargin{2};

% Case I: Called to convert a matrix to a tensor_as_matrix
if (nargin == 4)

  data = T;
  if ~isa(data,'numeric') || (ndims(data) ~= 2)
    error('A must be a matrix.');
  end
  cdims = varargin{3};
  tsize = varargin{4};
  
  % Error check
  n = numel(tsize);
  if ~isempty(rdims) && ((max(rdims) > n) || (min(rdims) < 1))
      error('One or more of the specified dimensions is out of range');
  elseif ~isempty(cdims) && ((max(cdims) > n) || (min(cdims) < 1))
      error('One or more of the specified dimensions is out of range');
  elseif (length([rdims cdims]) ~= n) || ~isequal(1:n, sort([rdims cdims]))
      error('Incorrect specification of dimensions');
  elseif (prod(tsize(rdims)) ~= size(data,1))
    error('SIZE(A,1) does not match size specified by RDIMS and SIZE.');
  elseif (prod(tsize(cdims)) ~= size(data,2))
    error('SIZE(A,2) does not match size specified by CDIMS and SIZE.');
  end
  
  % Save class variables
  A.tsize = tsize;
  A.rindices = rdims;
  A.cindices = cdims;
  A.data = data;
  A = class(A, 'tensor_as_matrix');
  return;
  
end

% Case II: Called to convert an MDA to a tensor_as_matrix
if isa(T,'double')

  A = tensor_as_matrix(tensor(T),varargin{2:nargin});
  return;

end
% Case III: Convert a tensor to a tensor_as_matrix

% Save the size of T and the number of dimensions
tsize = size(T);
n = ndims(T);

% Compute the default set of column indices (cdims might get changed later)
cdims = setdiff(1:n, rdims);

% Process optional third argument
if (nargin == 3)
    arg3 = varargin{3};
    if isa(arg3,'char')
        if (numel(rdims) ~= 1)
            error('Only one row dimension when third argument is a string.');
        end
        if strcmp(arg3,'fc')        % Forward cyclic.
            cdims = [rdims+1:n, 1:rdims-1];
        elseif strcmp(arg3,'bc')    % Backward cyclic.
            cdims = [rdims-1:-1:1, n:-1:rdims+1];
        else
            error('Unrecognized option');
        end
    elseif isa(arg3,'numeric')
        cdims = arg3;
    end
end

% Error check
if ~isempty(rdims) && ((max(rdims) > n) || (min(rdims) < 1))
    error('One or more of the specified dimensions is out of range');
elseif ~isempty(cdims) && ((max(cdims) > n) || (min(cdims) < 1))
    error('One or more of the specified dimensions is out of range');
elseif (length([rdims cdims]) ~= n) || ~isequal(1:n, sort([rdims cdims]))
    error('Incorrect specification of dimensions');
end

% Permute T so that the dimensions specified by RDIMS come first
T = permute(T,[rdims cdims]);

% Convert T to a matrix
data = reshape(T.data, prod(tsize(rdims)), prod(tsize(cdims)));


% Save class variables
A.tsize = tsize;
A.rindices = rdims;
A.cindices = cdims;
A.data = data;
A = class(A, 'tensor_as_matrix');
