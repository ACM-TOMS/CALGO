function solveQAP(printFileName,instance, solver, lambda)
%SOLVEQAP Solve QAP
% Usage
%   solveQAP(printFileName,instance,solver,lambda);
% Input
%   printFileName : a file name for output, for exampe, 'temp', then output is written in 'temp.csv'  
%   instance : the instance name in the directory instances/QAP/qapdata
%              for example, instance = 'chr12a'; 
%   solver   : 'BP, 'sdpnalplus', 'sdpnalAW'
%   lambda   : the Lagrangian multiplier, for example lambda = 1.0e5
% All instances in the directory instances/QAP/qapdata are from QAPLIB: 
% P. Hahn and M. Anjos. QAPLIB: a quadratic assignment problem library. 
% http://www.seas.upenn.edu/qaplib.
%
% For the DNN relaxation, see 
% S. Kim, M. Kojima, and K.-C. Toh.
% "A Lagrangian-DNN relaxation: a fast method for computing tight lower bounds 
% for a class of quadratic optimization problems," 
% Math. Program., 156 (2016) 161-187.
% 
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% "BBCPOP: A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems 
% with Binary, Box and Complementarity Constraints," arXiv:1804.00761, 
% to appear in ACM Transaction on Mathematical Software.
%
% Example:  
% >> printFileName='temp'; instance='chr12a'; solver='BP'; lambda=100000; 
% >> solveQAP(printFileName,instance,solver,lambda);
%

%% read data of QAP for our DNN formulation
% filename = [instance '.dat'];

if nargin == 2
    solver = 'BP';
    lambda = 100000;
elseif nargin == 3
    lambda = 100000;
end

switch solver
    case 'BP'
        %% solve our DNN by BP
        % Lagrangian relaxation
        % lambda = 1e5;
        [objPoly, I01, Icomp, relaxOrder, params] = qapreadBP(instance,lambda);

        % Lag-DNN
        n = size(objPoly.supports, 2);
        %         params = [];
        %         params.sparseSW = 0;
        %         relaxOrder = 1;
        [Q0vec, H0vec, PSDcone, polyCone, UbdIX] ...
            = BBCPOPtoDNN(objPoly, I01, Icomp, relaxOrder, params);
        UbdIX = 1 + sqrt(n);

        % solve
        %         params.delta = 0.1;
        %         params.delta1 = 1e-5;
        %         params.delta2 = 1e-6;
        %         params.UbdObjVal = optvalQAP(instance) * 1.5;
        params.solver = 'BP';
        
        [sol, info] = solveDNN(Q0vec, H0vec, PSDcone, polyCone, UbdIX, params);
        %         info.iterBP = infoBP.iter;
        %         info.timeBP = infoBP.timeBisection;
        %         info.iterAPGR = infoBP.APGiter;
        %         info.termcodeBP = infoBP.break_yes;
        eTime =info.timeBP; 
        iter =info.iterAPGR;
        termcode = info.termcodeBP;
        iterBP = info.iterBP;
        optVal = params.optVal;
        %         info.iter = info.iterBP;
        %         info.breakyes = info.termcodeBP;
        %         info.time = info.timeBP;
    case {'sdpnal','sdpnalplus'}

        %% solve our DNN by sdpnal+ without Lagrangian relaxation

        [objPoly, I01, Icomp, penaltyPoly, C, dd] = qapread(instance);

        % DNN
        n = size(objPoly.supports, 2);
        params = [];
        params.relaxOrder = 1;
        [Q0vec, H0vec, H1vec, PSDcone, polyCone, UbdIX] ...
            = BBPOPtoCOP(objPoly, {penaltyPoly}, ~I01, Icomp, params);
        params.UbdIX = 1 + sqrt(n);

        % solve
        params.delta = 0.1;
        params.UbdObjVal = optvalQAP(instance) * 1.5;
        params.solver = 'sdpnalplus';
        [sol, info] = solveCOP(Q0vec, H0vec, H1vec, PSDcone, polyCone, params);
        eTime =info.totaltime; 
        iter =info.iter;
        termcode = info.termcode;
        iterBP = 0; 
        optVal = [];
    case 'sdpnalAW'
        %% solve AW+(standard DNN formulation for sdpnal+) by sdpnal+
        clear blk At C b L U Bt l u

        [AA,BB] = qapreadnal(instance); 
        r = size(AA, 1);

        % DNN
        [blk,At,C,b,AAscale,BBscale] = qapAW(AA,BB,2);                
        L = 0; U = []; %% X >= 0
        Bt = []; l = []; u = [];

        % solve
        OPTIONS.tol = 1e-6;     
        OPTIONS.stopoption = -2; 
        OPTIONS.maxtime = 1e5;
        [obj,X,s,y,Z1,Z2,y2,v,info,runhist] = ...
            sdpnalplus(blk,At,C,b,L,U,Bt,l,u,OPTIONS);
        
        c = C{1}(:);
        Amat = smat(repmat(blk,size(At{1},2),1),mat2cell(At{1},size(At{1},1),ones(size(At{1},2),1)),1)';
        Avec = cell2mat(cellfun(@(x) x(:), Amat, 'UniformOutput', false));
        z = Z2{1}(:);
        K1.s = sqrt(size(Avec,1)); 
        mu = min(min(blkEigBP(c-Avec*y-z, K1)), 0);
        sol.LBv = (y'*b + mu * r) * AAscale * BBscale;
%         sol.LBv = [];
%         sol.UB = [];
        obj2 = obj * AAscale * BBscale;
        sol.LB = min(obj2);
        sol.UB = max(obj2);
        eTime = info.totaltime;
        iter = info.iter; 
        termcode = info.termcode;
        iterBP = 0;
        lambda = [];
        optVal = [];
        
end
%%
fid = fopen(printFileName,'a');
%fprintf(fid, 'instance,optVal,method,lambda,LB,LBv,UB,time,iter(APGR),termcode,iterBP\n');
fprintf(fid, '%s,%d,%s,%1.2e,%1.8e,%1.8e,%1.8e,%1.2e,%d,%d,%d\n',...
    instance,optVal,solver,lambda,sol.LB,sol.LBv,sol.UB,eTime,iter,termcode,iterBP);
fclose(fid);

end