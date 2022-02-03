function [x, y, info] = COPsdpnalplus(Q0vec, H0vec, H1vec, PSDcone, polyCone, params)
%COPSDPNAL Solve COP below by SDPNAL+
%   min. Q0vec'*x
%   s.t. H0vec'*x = 1
%        H1vec'*x = 0
%        x \in PSDcone
%        x \in polyCone
%
%   x \in polyCone is converted to x >= 0, A*x = b, and B*x >= 0(chain).
    
    BMat = polyCone.varStructure.BMat;
    [mm, nn] = size(BMat);
    scaledBMat = BMat * spdiags(1 ./ sum(BMat, 1)', 0, nn, nn);
    
    NMat = genNullMat(polyCone.varStructure);
    
    K.s = PSDcone;
    c = Q0vec;

    if norm(H1vec) == 0
        A = [H0vec'; ...
            BMat(:, polyCone.eq0Var)';...
            NMat];
    else
        A = [H0vec'; ...
            H1vec';
            BMat(:, polyCone.eq0Var)';...
            NMat];
    end
    
    % chain
    numIneq = max(0, sum(sum(polyCone.chain) - 1));
    numchain = size(polyCone.chain, 2);
    B = sparse(numIneq, mm);
    idx = 1;
    for ii = 1:numchain
        chainii = find(polyCone.chain(:, ii));
        endidx = idx + length(chainii) - 2;
        B(idx:endidx, :) = (scaledBMat(:, chainii(1:end-1)) - scaledBMat(:, chainii(2:end)))';
        idx = endidx + 1;
    end
    
    AB = [A; B];
    dummyb = [1; sparse(size(AB, 1) - 1, 1)];
    [blknal, ABtnal, Cnal, dummybnal, perm]= read_sedumi(AB, dummyb, c, K, 1);
    
    aa = size(A, 1);
    Atnal = cellfun(@(x) x(:, 1:aa), ABtnal, 'UniformOutput', false);
    Btnal = cellfun(@(x) x(:, aa+1:end), ABtnal, 'UniformOutput', false);
    bnal = dummybnal(1:aa);

    %OPTIONS.tol = 1e-8;
    OPTIONS.stopoption = 2;
    OPTIONS.maxtime = params.maxtime;
    [obj, X, s, y, S, Z, y2, v, info, runhist] ...
        = sdpnalplus(blknal,Atnal,Cnal,bnal,0,[],Btnal,0,[],OPTIONS); 

    if isfield(params, 'UbdIX')
        K1.s = PSDcone;
        z = cell2mat(cellfun(@(x) x(:), Z, 'UniformOutput', false));
        solEig = blkEigBP(c - A'*y - B'*v - z, K1);
        mu = min(min(solEig), 0);
        info.LBv = y(1) + mu * params.UbdIX;
    end
    
    info.cputime = info.timeSSN + info.timeADM;
    info.objval = obj;
    x = cell2mat(cellfun(@(x) x(:), X, 'UniformOutput', false));
    y = (x' * scaledBMat)';
        
end

% function [x, y, info] = COPsdpnalbasis(Q0vec, H0vec, H1vec, PSDcone, polyCone)
% %COPSDPNAL Summary of this function goes here
% %   Detailed explanation goes here
%     
%     BMat = polyCone.varStructure.BMat(:, ~polyCone.eq0Var);
%     [mm, nn] = size(BMat);
%     scaledBMat = BMat * spdiags(1 ./ sum(BMat, 1)', 0, nn, nn);
%     
%     K.l = nn;
%     K.s = PSDcone;
%     c = [BMat' * Q0vec / 2; Q0vec / 2];
% 
%     if norm(H1vec) == 0
%         A = [sparse(1, nn), H0vec'; ...
%             BMat, -speye(mm)...
%             ];
%     else
%         A = [sparse(1, nn), H0vec'; ...
%             sparse(1, nn), H1vec';...
%             BMat, -speye(mm)
%             ];
%     end
%     
%     % chain
%     chain = polyCone.chain(~polyCone.eq0Var, :);
%     numIneq = max(0, sum(sum(chain) - 1));
%     numchain = size(chain, 2);
%     B = sparse(numIneq, nn + mm);
%     topidx = 1;
%     for ii = 1:numchain
%         chainii = find(chain(:, ii));
%         endidx = topidx + length(chainii) - 2;
%         rowidx = repmat(1:length(chainii)-1, 2, 1);
%         colidx = [chainii(1:end-1); chainii(2:end)];
%         val = [ones(1, length(chainii)-1); -ones(1, length(chainii)-1)];
%         B(topidx:endidx, :) = ...
%             [sparse(rowidx(:), colidx(:), val(:), length(chainii)-1, nn), ...
%             (scaledBMat(:, chainii(1:end-1)) - scaledBMat(:, chainii(2:end)))'];
%         topidx = endidx + 1;
%     end
%     
%     AB = [A; B];
%     dummyb = [1; sparse(size(AB, 1) - 1, 1)];
%     [blknal, ABtnal, Cnal, dummybnal, perm]= read_sedumi(AB, dummyb, c, K, 1);
%     
%     aa = size(A, 1);
%     Atnal = cellfun(@(x) x(:, 1:aa), ABtnal, 'UniformOutput', false);
%     Btnal = cellfun(@(x) x(:, aa+1:end), ABtnal, 'UniformOutput', false);
%     bnal = dummybnal(1:aa);
% 
%     OPTIONS.tol = 1e-8;
%     [obj, X, s, y, S, Z, y2, v, info, runhist] ...
%         = sdpnalplus(blknal,Atnal,Cnal,bnal,[],[],Btnal,0,[],OPTIONS); 
% 
%     info.objval = obj;
%     X = cell2mat(cellfun(@(x) x(:), X, 'UniformOutput', false));
%     x = X(nn+1:end);
%     y = X(1:nn);
% end
% function [x, y, info] = COPsdpnalpre(Q0vec, H0vec, H1vec, PSDcone, polyCone)
% %COPSDPNAL Summary of this function goes here
% %   Detailed explanation goes here
%     
%     K.s = PSDcone;
%     c = Q0vec;
%     
%     BMat = polyCone.varStructure.BMat;
%     [mm, nn] = size(BMat);
%     scaledBMat = BMat * spdiags(1 ./ sum(BMat, 1)', 0, nn, nn);
%     
%     if norm(H1vec) == 0
%         A = [H0vec'; scaledBMat(:, polyCone.eq0Var)'];
%     else
%         A = [H0vec'; H1vec'; scaledBMat(:, polyCone.eq0Var)'];
%     end
%     
%     % chain
%     numIneq = max(0, sum(sum(polyCone.chain) - 1));
%     numchain = size(polyCone.chain, 2);
%     B = sparse(numIneq, mm);
%     idx = 1;
%     for ii = 1:numchain
%         chainii = find(polyCone.chain(:, ii));
%         endidx = idx + length(chainii) - 2;
%         B(idx:endidx, :) = (scaledBMat(:, chainii(1:end-1)) - scaledBMat(:, chainii(2:end)))';
%         idx = endidx + 1;
%     end
%     AB = [A; B];
%     dummyb = [1; sparse(size(AB, 1) - 1, 1)];
%     [blknal,ABtnal,Cnal,dummybnal,perm]= read_sedumi(AB,dummyb,c,K);
%     
%     aa = size(A, 1);
%     Atnal = cellfun(@(x) x(:, 1:aa), ABtnal, 'UniformOutput', false);
%     Btnal = cellfun(@(x) x(:, aa+1:end), ABtnal, 'UniformOutput', false);
%     bnal = dummybnal(1:aa);
%     
% 
%     [obj, X, s, y, S, Z, y2, v, info, runhist] ...
%         = sdpnalplus(blknal,Atnal,Cnal,bnal,0,[],Btnal,0,[]); 
% 
%     dummy = sparse(size(BMat, 2), 1);
%     [~,scaledBMatnal,~,~,~]= read_sedumi(scaledBMat',dummy,c,K);
%     
%     info.objval = obj;
%     x = X;
%     
%     
%     
% end
% 
% 
% 
% function [x, y, info] = COPsdpnalprepre(Q0vec, H0vec, H1vec, PSDcone, polyCone)
% %COPSDPNAL Summary of this function goes here
% %   Detailed explanation goes here
%     
%     K.s = PSDcone;
%         BMat = polyCone.varStructure.BMat;
%         nn = size(BMat, 2);
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
%         A = [Ateq; Ateq0; Atchain; AtPSD]';
% 
%         c = sparse(1, 1, 1, size(A, 2), 1);
%         
%     [blknal,Atnal,Cnal,bnal,perm]= read_sedumi(A,b,c,K);
%     
%     [obj, X, s, y, S, Z, y2, v, info, runhist] ...
%         = sdpnalplus(blknal,Atnal,Cnal,bnal,0,[],[],[],[]); 
% 
%     info.objval = obj;
%     x = X;
%     
%     
%     
% end
