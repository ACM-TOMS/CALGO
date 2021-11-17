function B = squeezem(A)
% [ B ] = squeezem( A )
% Consistent behaviour of squeeze, with regards to multi-dimensional applications.
% Returns an array B with the same elements as A but with all the singleton dimensions removed.  
% A singleton is a dimension such that size(A,dim)==1.
% Row-vectors become column vectors.
%
% E.g.: squeezem([2 3 4])
%       squeezem([2 1 3])
%
% See also: squeeze
%
% Written by: tommsch, 2018
% Uses code from: matlabs squeeze

  siz = size(A);
  siz(siz==1) = []; % Remove singleton dimensions.
  siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
  B = reshape(A,siz);
  
end

function dummy; end %#ok<DEFNU>  %Generates an error, if the 'end' of a function is missing.