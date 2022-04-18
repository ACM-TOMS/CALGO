function [sub] = ind2subm(sze, idx)
% sub = ind2subm( sze, idx )
% Multiple subscripts from linear indices.
% Similar to matlabs built-in function, but returns an array of subscripts. Each column are the subscripts for one index.
% XX Function is very slow and should be vectorized
%
% Input:
%   sze         Size of the array
%   idx         linear indices to be converted
%
% Output:
%   sub         array, each column are the subscripts for the corresponding indices in idx
%
% E.g.: ind2subm([2 1 3],[1 6])
%
% Written by tommsch, 2019

sub = cell(1,numel(sze));
[sub{:}] = ind2sub(sze,idx);
sub = vertcat(sub{:});
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.


