function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for a tensor.
%
%   Examples
%      T = tensor(rand(3,4,2,));
%      T(1:2,1:2,1) = ones(2,2); <-- Calls SUBSASGN
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR, TENSOR/SUBSREF.

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


switch s.type    
    case '.'
	error(['Cannot change field ', s.subs, ' directly.']);
    case '()'
	data = t.data;
        if isa(b,'tensor')
	  data(s.subs{:}) = b.data;
        else
          data(s.subs{:}) = b;
        end
	t = tensor(data, t.size);
    case '{}'      
	error('Subscript cell reference not supported for tensor.');
    otherwise
        error('Incorrect indexing into tensor.')
end


