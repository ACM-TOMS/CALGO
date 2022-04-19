function X = projOntoPolyConeMex(Z, polyCone)
%PROJONTOPOLYCONE   Projection onto a polyhedral cone which has a chain structure.
%   
%   Usage:
%      X = projOntoPolyCone( Z, polyCone );


BB = polyCone.varStructure.BMat;
%chain = polyCone.chain;
%if max(sum(chain, 2)) > 1; error('Chains are not disjoint'); end

cardi = polyCone.varStructure.cardinality;
aveZ = (Z(:)' * BB)' ./ cardi;

yy = zeros(size(BB, 2), 1);

% Singleton (chain with length 1)
singleton = polyCone.singleton;
yy(singleton) = aveZ(singleton); % ./ cardi(singleton);

% Non singleton
% Use a pre-permuted matrix since permuting aveZ at each function call is too time consuming
%w
cardiPerm = polyCone.cardinalityPerm;
%b
aveZPerm = (Z(:)' * polyCone.BMatPerm ./ cardiPerm)'; 
if ~isempty(aveZPerm)
    yytmp = mexProxMonotonicNew(aveZPerm, cardiPerm', polyCone.I');
    yy(polyCone.permIdx) = yytmp;
end

projlbd = yy<polyCone.lbd;
yy(projlbd) = polyCone.lbd(projlbd);
projubd = yy>polyCone.ubd;
yy(projubd) = polyCone.ubd(projubd);
%
% if polyCone.nonneg
%     yy(yy<0) = 0;
% end

%
X = BB * yy;
X = X'; %for BP
end

%chainCell = polyCone.chainCell;
% Slow examples

% % Non singleton
% % each cell contains an input or output vector of each isotonic problem
% %input w
% cardiChainCell = polyCone.cardiChainCell;
% %input b
% aveZcell = cellfun(@(ch) aveZ(ch), chainCell, 'UniformOutput', false);
% %aveZcell = polyCone.computeAveZCell(Z); % mat2cell is slow
% %output yy
% yycell = cellfun(@mexProxMonotonicRev, aveZcell, cardiChainCell, 'UniformOutput', false);
% for pp = 1:length(yycell)
%     yy(chainCell{pp}) = yycell{pp};
% end

% % Non singleton
% yy2 = yy;
% aveZcell = cell(length(polyCone.cardinalityOfChainsCell), 1);
% aveZcell = cellfun(@(BB,ca) (full(Z(:)' * BB) ./ ca)', ...
%     polyCone.BMatOfChainsCell, polyCone.cardinalityOfChainsCell,...
%     'UniformOutput', false);
% yycell = cellfun(@mexProxMonotonicRev, aveZcell, polyCone.cardiChainCell, ...
%     'UniformOutput', false);
% for pp = 1:length(yycell)
%     yy2(chainCell{pp}) = yycell{pp};
% end
% if ~isequal(yy,yy2); error('error!'); end