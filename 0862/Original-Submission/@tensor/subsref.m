function a = subsref(t,s)
%SUBSREF Subscripted reference for tensors.
%
%   Examples
%      T = tensor(rand(3,4,2,1),[3 4 2 1]);
%      T(1,1,1,1) <-- produces a scalar
%      T(1,1,1,:) <-- produces a tensor of order 1 and size 1
%      T(:,1,1,:) <-- produces a tensor of size 3 x 1
%      T(1:2,[2 4],1,:) <-- produces a tensor of size 2 x 2 x 1
%      T(1:2,[2 4],1,1) <-- produces a tensor of size 2 x 2
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

switch s(1).type    
    case '.'
        switch s(1).subs
            case 'data'
		if numel(s) == 1
		    a = t.data;
		else
		    a = subsref(t.data, s(2:end));
		end
	    case 'size'
		if numel(s) == 1
		    a = t.size;
		else
		    a = subsref(t.size, s(2:end));
		end
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    
    case '()'

	% Extract the data
	a = t.data(s.subs{:});

	% Determine the subscripts
	subs = s(1).subs;
	asiz = [];
        kpdims = [];
        rmdims = [];
        for i = 1:length(subs)
          if ischar(subs{i}) && (subs{i} == ':')
            asiz = [asiz size(t,i)];
            kpdims = [kpdims i];
          elseif numel(subs{i}) > 1
            asiz = [asiz numel(subs{i})];
            kpdims = [kpdims i];
          else
            rmdims = [rmdims i];
          end
        end

	% Do we return a scalar or not?
	if ~isempty(asiz)
	    a = tensor(permute(a, [kpdims rmdims]),asiz);
	end
	
    case '{}'
	error('Subscript cell reference cannot be used for dense tensors.')
    
    otherwise
        error('Incorrect indexing into tensor.')
end
