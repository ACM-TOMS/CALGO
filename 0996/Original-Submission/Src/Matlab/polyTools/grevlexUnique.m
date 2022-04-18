function [uniquearray, ia1, ic1] = grevlexUnique(array)
%GREVLEXUNIQUE   Apply unique() and sort rows in grevlex order

% if 0 
tmp = [sum(array, 2), -array];
[tmp, ia1, ic1] = unique(tmp, 'rows');
uniquearray = -tmp(:, 2:end);
% else
%     [dim, I] = sort(full(sum(array, 2)));
%     [~, revI] = sort(I);
%     
%     num = accumarray(dim+1, 1);
%     cumnum = cumsum([0; num]);
%     tmp = mat2cell(-array(I, :), num, size(array, 2));
%     ic1 = cell(length(tmp), 1);
%     ia1 = cell(length(tmp), 1);
%     uniquearray = cell(length(tmp), 1);
% 
%     iasum = 0;
%     for ii = 1:length(tmp)
%         [uniquearray{ii}, ia1{ii}, ic1{ii}] = unique(tmp{ii}, 'rows');
%         ic1{ii} = ic1{ii} + iasum;
%         ia1{ii} = ia1{ii} + cumnum(ii);
%         iasum = iasum + length(ia1{ii});
%     end
%     ic1 = cell2mat(ic1); ic1 = ic1(revI);
%     ia1 = cell2mat(ia1); ia1 = I(ia1);
%     uniquearray = -cell2mat(uniquearray);
% end

end

