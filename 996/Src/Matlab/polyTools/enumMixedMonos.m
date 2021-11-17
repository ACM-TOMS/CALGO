function monos = enumMixedMonos(numVar, relaxOrder, Ibin)
%ENUMMIXEDMONOS   Generate indices of moment matrix.
%
%   INPUT:
%      numVar     :
%      relaxOrder : relaxation order of DNN hierarchy.
%      Ibin       : logical vector which specifies binary variables.
%
%   OUTPUT:
%      monos: two dimensional array
%                   whose row represents the supports of ,,,

Ibin = logical(Ibin);
if(numVar ~= length(Ibin)); error('dimensions must agree'); end
Ibox = ~Ibin;
numBoxVar = sum(Ibox);
numBinVar = sum(Ibin);

tmp = enumMonosUptoDeg(numBoxVar+1, relaxOrder);
trashIdx = tmp(:, 1) > numBinVar;
degBin = tmp(~trashIdx, 1);
boxMonos = tmp(~trashIdx, 2:end);

binMonos = enum01MonosUptoDeg(numBinVar, max(degBin));
sumBin = sum(binMonos, 2);

ddmax = max(degBin);
boxIdx = cell(ddmax+1, 1);
binIdx = cell(ddmax+1, 1);
for dd = 0:ddmax
    vec1 = find(degBin == dd);
    vec2 = find(sumBin == dd);
    [ddBoxIdx, ddBinIdx] = myMeshGrid(vec1, vec2);
    boxIdx{dd+1} = ddBoxIdx;
    binIdx{dd+1} = ddBinIdx;
end
boxIdx = cell2mat(boxIdx);
binIdx = cell2mat(binIdx);

mixedSupports = sparse(length(boxIdx),numVar);
mixedSupports(:, Ibox) = boxMonos(boxIdx, :);
mixedSupports(:, Ibin) = binMonos(binIdx, :);
monos = grevlexSort(mixedSupports);

end