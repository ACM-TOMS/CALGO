function t = tensor(varargin)
%TENSOR Tensor object constructor.
% 
%   T = TENSOR creates an empty, dense tensor object.
%
%   T = TENSOR(Z) creates a tensor from the multidimensional array Z.
%
%   T = TENSOR(Z,DIMS) creates a tensor from the multidimensional
%   array Z. The DIMS argument is used to specify any trailing
%   singleton dimensions or if there is only one dimension. If DIMS is
%   empty, then the result is a tensor that is the size of Z.
%
%   T = TENSOR(S) copies a tensor S.
%
%   T = TENSOR(A) converts a cp_tensor, tucker_tensor, or
%   tensor_as_matrix object to a tensor.
%
%   Examples
%   Create a random third order tensor of size 3 x 4 x 2
%      T = tensor(rand(3,4,2));
%   Create a random first order tensor of size 3
%      T = tensor(rand(3),3);
%   Create a random fourth order tensor of size 3 x 4 x 2 x 1
%      T = tensor(rand(3,4,2,1)); <-- Doesn't work!
%      T = tensor(rand(3,4,2,1),[3 4 2 1]); <-- Works.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR_AS_MATRIX, CP_TENSOR, TUCKER_TENSOR.


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


% EMPTY/DEFAULT CONSTRUCTOR
if size(varargin) == 0

    t.data = []; 
    t.size = 0;
    t = class(t, 'tensor');
    return;
    
end

% COPY CONSTRUCTOR    
if isa(varargin{1}, 'tensor')
    
    t.data = varargin{1}.data;
    t.size = varargin{1}.size;
    t = class(t, 'tensor');
    return;

end

% CONVERT CP_TENSOR OR TUCKER_TENSOR
if isa(varargin{1}, 'cp_tensor') || isa(varargin{1}, 'tucker_tensor') 
    
    t = full(varargin{1});
    return;

end

% CONVERT TENSOR_AS_MATRIX
if isa(varargin{1},'tensor_as_matrix')

    a = varargin{1};
    idx = [a.rindices a.cindices];
    if isempty(idx)
      idx = 1;
    end
    sz = a.tsize;
    if ~isempty(sz)
        sz = sz(idx);
    end
    t.data = reshape(a.data, [sz 1 1]);
    t.size = sz;
    t = class(t, 'tensor');
    [sidx indx] = sort(idx);
    t = permute(t,indx);
    return;

end

% CONVERT A MULTIDIMENSIONAL ARRAY
if (nargin == 1) || (nargin == 2) 
    
    if ~isa(varargin{1},'numeric') && ~isa(varargin{1},'logical')
        error('Z must be a multidimensional array.')
    end
    
    t.data = varargin{1};

    t.size = [];
    if nargin == 1
        t.size = size(t.data);
    else
        t.size = varargin{2};

        if (length(t.size) > 2) && (size(t.size,1) ~= 1)
            error('DIMS must be a row vector.');
        end
	
        % -- Error Check --
        % First, check that the matching dimensions do indeed match
        sz = size(t.data);
        j = min( [length(sz),length(t.size)] );
        for i = 1 : j
            if (sz(i) ~= t.size(i))
                error('Specified size is incorrect.');
            end
        end
	
        % Second, check that the remaining dimensions are okay
        for i = j + 1 : length(sz)
            if (sz(i) ~= 1)
                error('Specified size is incorrect.');
            end
        end

        for i = j + 1 : length(t.size)
            if (t.size(i) ~= 1)
                error('Specified size is incorrect.');
            end
        end
    end
    
    t = class(t, 'tensor');
    return;

end


error('Unsupported use of function TENSOR.');


