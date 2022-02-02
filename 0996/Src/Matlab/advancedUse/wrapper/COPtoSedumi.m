function [ A, b, c, K ] =...
    COPtoSedumi( Q0vec, H0vec, H1vec, PSDcone, polyCone, params )
%COPTOSEDUMI   Convert cone optimization problem into sedumi format.
%   Usage:
%      [A,b,c,K] = COPtoSedumi(Q0vec,H0vec,PSDcone,polyCone,params);
%
%   The original problem is
%   min. Q0vec'*x
%   s.t. H0vec'*x = 1
%        H1vec'*x = 0
%        x\in PSDcone \cup polyCone
%



    
% use the basis representation:
%
% x = BMat*y \in PSDcone
%
% and the dual representation:
%
% max. -(Q0vec'*BMat)*y
% s.t. 1 - (H0vec'*BMat)*y = 0
%      0 - (H1vec'*BMat)*y = 0
%      0 - EQ*y = 0                (complementarity constr)
%      0 - (-I)*y \geq 0           (DNN constr)
%      0 - (-E)*y \geq 0           (chain constr)
%      0 - (-BMat)*y \in PSDcone   (varStructure, PSDcone)

if isfield(polyCone.varStructure, 'BMat')
    BMat = polyCone.varStructure.BMat;
else
    BMat = genBasisMat(polyCone.varStructure);
end
nn = size(BMat, 2);


b = -Q0vec'*BMat;

if norm(H1vec) == 0
    Ateq = H0vec' * BMat;
    K.f = 1;
else
    Ateq = [H0vec' * BMat; H1vec' * BMat];
    K.f = 2;
end

% eq0 (complementarity)
eq0Var = find(polyCone.eq0Var);
numeq0 = length(eq0Var);
Ateq0 = sparse(1:numeq0, eq0Var, 1, numeq0, nn);
K.f = K.f + numeq0;

% DNN
K.l = 0;
if polyCone.nonneg
    AtDNN = -speye(nn); %-I
    K.l = K.l + nn;
else
    AtDNN = sparse(0, nn);
end

% chain
numIneq = 0;
numchain = size(polyCone.chain, 2);
row = cell(1, numchain);
col = cell(1, numchain);
pointer = 0;
for ii = 1:numchain
    chainii = find(polyCone.chain(:, ii));
    numIneqii = length(chainii) - 1;
    numIneq = numIneq + numIneqii;

    if numIneqii
        row{ii} = pointer + [1:numIneqii; 1:numIneqii];
        col{ii} = [chainii(1:end-1)'; chainii(2:end)'];
    end
    pointer = pointer + numIneqii;
end
row = cell2mat(row);
col = cell2mat(col);
val = [ones(1, numIneq); -ones(1, numIneq)];
Atchain = -sparse(row(:), col(:), val(:), numIneq, nn); % -E
% sparse([],[],[],0,nn) works.
K.l = K.l + numIneq;

% PSD constr
AtPSD = -BMat;
K.s = PSDcone;

A = [Ateq; Ateq0; AtDNN; Atchain; AtPSD]';

c = sparse(1, 1, 1, size(A, 2), 1);
    
end

%% archive (Null representation)
% varStructure = polyCone.varStructure;
% computeMomentIdx = varStructure.computeMomentIdx;
% momentic = varStructure.momentic;
% function [momentIdx,val] = findVar(varidx)
%     where = momentic == varidx;
%     blkIdx = varStructure.blkIdx(where);
%     rowIdx = varStructure.rowIdx(where);
%     colIdx = varStructure.colIdx(where);
%     momentIdx = [computeMomentIdx(blkIdx, rowIdx, colIdx);...
%         computeMomentIdx(blkIdx, colIdx, rowIdx)];
%     val = (1 / length(momentIdx)) * ones(length(momentIdx), 1);
% end
% %   The major difference lies in a way to represent varStructure.
% switch params.COPrep
%     case 'null'
%         % introduce the slack variables u,v\geq 0
%         %
%         % and use the primal representation
%         %
%         % min. [0,0,Q0vec']*[u;v;x]
%         % s.t. [0,0,H0vec']*[u;v;x] = 1
%         %      [0,0,H1vec']*[u;v;x] = 0 
%         %      [0,0,NMat]*[u;v;x] = 0    (varStructure)
%         %      [0,-I,I]*[u;v;x] = 0  (DNN constr)
%         %      [-I,0,CH]*[u;v;x] = 0 (chain constr)
%         %      [0,0,EQ]*[u;v;x] = 0 (eq0Var(complementarity) constr)
%         %      u,v\geq 0;
%         %      x\in PSDcone
%         NMat = genNullMat(polyCone.varStructure);
%         [mm,nn] = size(NMat);
%         
%         %Achain
%         numIneq = 0;
%         numchain = size(polyCone.chain,2);
%         for ii = 1:numchain
%             chainii = find(polyCone.chain(:,ii));
%             numIneqii = length(chainii) - 1;
%             numIneq = numIneq + numIneqii;
%         end
%         pointer = 0;
%         row = cell(numIneq,1);
%         col = cell(numIneq,1);
%         val = cell(numIneq,1);
%         for ii = 1:numchain
%             chainii = find( polyCone.chain(:,ii) );
%             numIneqii = length(chainii) - 1;
%             [leftmomentIdx,leftcoef] = findVar(chainii(1));
%             for jj = 2:(numIneqii+1)
%                 pointer = pointer + 1;
%                 [rightmomentIdx,rightcoef] = findVar(chainii(jj));
%                 col{pointer} = [leftmomentIdx;rightmomentIdx];
%                 val{pointer} = [leftcoef;-rightcoef];
%                 row{pointer} = pointer * ones(length(col{pointer}),1);
%                 
%                 leftmomentIdx = rightmomentIdx;
%                 leftcoef = rightcoef;
%             end
%         end
%         row = cell2mat(row);
%         col = cell2mat(col);
%         val = cell2mat(val);
%         CH = sparse(row,col,val,numIneq,nn); 
%         Achain = [-speye(numIneq),sparse(numIneq,nn),CH];
%         
%         % verStructure
%         Avar = [sparse(mm,numIneq+nn),NMat];
%         
%         % DNN
%         if polyCone.nonneg
%             ADNN = [sparse(nn,numIneq),-speye(nn),speye(nn)];
%         else
%             ADNN = sparse(0,numIneq+nn+nn);
%         end
%         
%         % eq0var
%         eq0Var = find(polyCone.eq0Var);
%         numeq0 = length(eq0Var);
%         row = cell(numeq0,1);
%         col = cell(numeq0,1);
%         val = cell(numeq0,1);
%         pointer = 0;
%         for ii = 1:numeq0
%             pointer = pointer + 1;
%             [eq0momentIdx,eq0coef] = findVar(eq0Var(ii));
%             col{pointer} = eq0momentIdx;
%             val{pointer} = eq0coef;
%             row{pointer} = pointer * ones(length(col{pointer}),1);
%         end
%         row = cell2mat(row);
%         col = cell2mat(col);
%         val = cell2mat(val);
%         EQ = sparse(row,col,val,numeq0,nn);
%         Aeq0 = [sparse(numeq0,numIneq+nn),EQ];
%         
%         K.l = numIneq + nn;
%         K.s = PSDcone;
%         if norm(H1vec) == 0
%             Aeq = [sparse(1,numIneq+nn), H0vec'];
%         else
%             Aeq = [sparse(1,numIneq+nn), H0vec'; sparse(1,numIneq+nn), H1vec'];
%         end
%         A = [Aeq; Avar; ADNN; Achain; Aeq0];
% 
%         c = [sparse(numIneq + nn, 1); Q0vec];
%         b = sparse(1, 1, 1, size(A, 1), 1);
% 
% 
%     case 'basis'
%         % use the basis representation:
%         %
%         % x = BMat*y \in PSDcone
%         %
%         % and the dual representation:
%         %
%         % max. -(Q0vec'*BMat)*y
%         % s.t. 1 - (H0vec'*BMat)*y = 0
%         %      0 - (H1vec'*BMat)*y = 0
%         %      0 - EQ*y = 0                (complementarity constr)
%         %      0 - (-I)*y \geq 0           (DNN constr)
%         %      0 - (-E)*y \geq 0           (chain constr)
%         %      0 - (-BMat)*y \in PSDcone   (varStructure, PSDcone)
%         
%         if isfield(polyCone.varStructure, 'BMat')
%             BMat = polyCone.varStructure.BMat;
%         else
%             BMat = genBasisMat(polyCone.varStructure);
%         end
%         nn = size(BMat, 2);
% 
% 
%         b = -Q0vec'*BMat;
%         
%         if norm(H1vec) == 0
%             Ateq = H0vec' * BMat;
%             K.f = 1;
%         else
%             Ateq = [H0vec' * BMat; H1vec' * BMat];
%             K.f = 2;
%         end
%         
%         % eq0 (complementarity)
%         eq0Var = find(polyCone.eq0Var);
%         numeq0 = length(eq0Var);
%         Ateq0 = sparse(1:numeq0, eq0Var, 1, numeq0, nn);
%         K.f = K.f + numeq0;
%         
%         % DNN
%         K.l = 0;
%         if polyCone.nonneg
%             AtDNN = -speye(nn); %-I
%             K.l = K.l + nn;
%         else
%             AtDNN = sparse(0, nn);
%         end
%         
%         % chain
%         numIneq = 0;
%         numchain = size(polyCone.chain, 2);
%         row = cell(1, numchain);
%         col = cell(1, numchain);
%         pointer = 0;
%         for ii = 1:numchain
%             chainii = find(polyCone.chain(:, ii));
%             numIneqii = length(chainii) - 1;
%             numIneq = numIneq + numIneqii;
% 
%             if numIneqii
%                 row{ii} = pointer + [1:numIneqii; 1:numIneqii];
%                 col{ii} = [chainii(1:end-1)'; chainii(2:end)'];
%             end
%             pointer = pointer + numIneqii;
%         end
%         row = cell2mat(row);
%         col = cell2mat(col);
%         val = [ones(1, numIneq); -ones(1, numIneq)];
%         Atchain = -sparse(row(:), col(:), val(:), numIneq, nn); % -E
%         % sparse([],[],[],0,nn) works.
%         K.l = K.l + numIneq;
%         
%         % PSD constr
%         AtPSD = -BMat;
%         K.s = PSDcone;
% 
%         A = [Ateq; Ateq0; AtDNN; Atchain; AtPSD]';
% 
%         c = sparse(1, 1, 1, size(A, 2), 1);
%         
%     otherwise
%         error('params.COPrep must be ''basis'' or ''null''.')
% end