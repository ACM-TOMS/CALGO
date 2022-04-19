function ret = isvectorm(v);
% ret = isvectorm(v);
% returns true if v has exactly one dimension whose size is equal 1.
%
% E.g.: isvectorm(ones(1,1,3))
%       isvectorm(ones(0))
%       isvector(ones(1,1,3))
%       isvector(ones(0))
%
% Note: isvectorm([]) returns false, since the empty set is no vector space
%
% See also: isvector
%
% Written by: tommsch, 2018
val=sizem(v);
if(nnz(val>1)==1 || isequal(val,1))
    ret=true;
else
    ret=false;
end
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   