function [idx1,idx2] = myMeshGrid(vec1,vec2)
%%MYMESHGRID a simple meshgrid function with vector format output
% [idx1,idx2] = myMeshGrid(vec1,vec2)

if ~isrow(vec1); vec1 = vec1'; end
if ~isrow(vec2); vec2 = vec2'; end

idx1 = ones(length(vec2),1) * vec1;
idx2 = vec2' * ones(1,length(vec1));

idx1 = idx1(:);
idx2 = idx2(:);

end

