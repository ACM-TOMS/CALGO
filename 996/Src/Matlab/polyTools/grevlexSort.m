function [array,index] = grevlexSort(array)
%GREVLEXSORT   Sort rows in grevlex order
%

tmp = [sum(array,2),-array];
if issorted(tmp,'rows')
    index = (1:size(array,1))';
else
    [tmp,index] = sortrows(tmp); %unique(tmp,'rows');
end
array = -tmp(:,2:end);

end