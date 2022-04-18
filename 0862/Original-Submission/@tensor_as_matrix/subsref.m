function a = subsref(t,s)
%SUBSREF Subscripted reference for tensor_as_matrix.
%
%   Copyright 2005, Tamara Kolda and Brett Bader, Sandia National Labs
%
%   See also TENSOR_AS_MATRIX.

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
		if length(s) == 1
		    a = t.data;
		else
		    a = subsref(t.data, s(2:end));
		end
	    case 'tsize'
		if length(s) == 1
		    a = t.tsize;
		else
		    a = subsref(t.tsize, s(2:end));
		end
	    case 'rindices'
		if length(s) == 1
		    a = t.rindices;
		else
		    a = subsref(t.rindices, s(2:end));
		end
	    case 'cindices'
		if length(s) == 1
		    a = t.cindices;
		else
		    a = subsref(t.cindices, s(2:end));
		end
            otherwise
                error(['No such field: ', s.subs]);
        end
    case '()'
	a = t.data(s.subs{:});
	
    otherwise
        error('Invalid subsref into tensor_as_matrix.')
end
