function a = subsref(t,s)
%SUBSREF Subscripted reference for a cp_tensor.
%
%   Examples
%      A = CP_TENSOR([3; 2], rand{4,2), rand(5,2), rand(3,2));
%      A.lambda returns the lambda array ([3;2]).
%      A.U returns a cell array of 3 matrices.
%      A.U{1} returns the matrix corresponding to the first mode.
%      A(2,3,1) calculates and returns that single element of A.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also CP_TENSOR.


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
            case 'lambda'
		if length(s) == 1
		    a = t.lambda;
		else
		    a = subsref(t.lambda, s(2:end));
		end
            case {'u','U'}
		if length(s) == 1
		    a = t.u;
		else
		    a = subsref(t.u, s(2:end));
		end
            otherwise
                error(['No such field: ', s(1).subs]);
        end
    case '()'
	a = 0;
	for k = 1 : length(t.lambda)
	    b = 1;
	    for i = 1 : length(s.subs)
		b = b * t.u{i}(s.subs{i},k);
	    end
	    a  = a + t.lambda(k) * b;
	end          
    case '{}'
	a = subsref(t.u,s);
   otherwise
        error('Invalid subsref.');
end
