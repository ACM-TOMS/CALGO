function [sortedarray, I] = fastSort(array)
%%FASTSORT    Sort values in array
%   [sortedarray, I] = fastSort(array);
%   
%   This function treats each row of array as a single entity
%   and returns the rows of array that are sorted by sum(array, 2). 

rng(1);

deg = sum(array, 2);
val = array * rand(size(array, 2), 1);
val = val + deg * 2^ceil(log2(max(val)));
[~, I] = sort(val);
sortedarray = array(I, :);

end