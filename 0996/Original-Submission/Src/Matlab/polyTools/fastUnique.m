function [uniquearray, ia1, ic1] = fastUnique(array)
%%FASTUNIQUE    Find unique values in array
%   [uniquearray, ia1, ic1] = fastUnique(array);
%   
%   This function treats each row of array as a single entity
%   and returns the unique rows of array (uniquearray). 
%
%   The rows of uniquearray are NOT in sorted order.
%   But sum(uniquearray,2) is in sorted order.
%
%   Caution: This function may fail with negligibly small probability.

rng(1);
maxit = 5;

iter = 0;
isSuccess = false;
while ~isSuccess && iter < maxit
    deg = sum(array, 2);
    val = array * sort(rand(size(array, 2), 1),'ascend');
    %val = array * rand(size(array, 2), 1);
    val = val + deg * 2^ceil(log2(max(val)));
    [~, ia1, ic1] = unique(val);
    uniquearray = array(ia1, :);
    
    isSuccess = isequal(array, uniquearray(ic1, :));
    iter = iter + 1;
end

if ~isSuccess 
    warning('fastUnique fails -> grevlexUnique'); 
    [uniquearray, ia1, ic1] = grevlexUnique(array);
end

end