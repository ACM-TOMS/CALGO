function solvePOPchordal(printFileName,degree,nDim,radiorange,isBin,addComplement,solver)
%SOLVEPOPCHORDAL Solve a POP with chordal sparsity.
% Usage
%   solvePOPchordal(printFileName,degree,nDim,radiorange,isBin,addComplement,solver); 
%
% Input
%   printFileName : a file name for output, for exampe, 'temp', then output is written in 'temp.csv' 
%   degree  : the degree of the polynomial objective function
%   nDim    : the number of variables
%   radiorange : radiorange which controls the sparsity of the chordal graph 
%                to be generated. 
%                for example, radiorange = 0.1
%   isBin   : 1 for the binary constraint, 0 for the box constraint
%   addComplement   : 1 for adding random complementarity constraints 
%   solver   : 'BP, 'sdpnalplus', 'sedumi'
%
% For the DNN relaxation, see 
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% "BBCPOP: A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems 
% with Binary, Box and Complementarity Constraints," arXiv:1804.00761, 
% to appear in ACM Transaction on Mathematical Software.
% 
% Example:  
% >> printFileName='temp'; degree=3; nDim=100; radiorange=0.1; isBin=1; addComplement=0; solver='BP'; 
% >> solvePOPchordal(printFileName,degree,nDim,radiorange,isBin,addComplement,solver);
% 

[objPoly, I01, Icomp] = genPOPchordal(degree,nDim,radiorange,isBin,addComplement);
    
%     Icomp = [];
%     if addComplement
%         Icomp = complementarityConstraint(2, 2, clique, randSeed);
%     end
    
    if nargin == 4
        solver = 'BP'; 
    end
    
    if strcmp(solver, 'newBP')
        params.solver = 'BP';
        params.newBPswitch = true;
    else
        params.solver = solver;
        params.newBPswitch = false;
    end
    
    %% to DNN
    eqPolySys = [];
    params.relaxOrder = ceil(degree / 2);
    [Q0vec, H0vec, H1vec, PSDcone, polyCone, UbdIX] ...
            = BBPOPtoCOP(objPoly, eqPolySys, I01, Icomp, params);

    %% solve
    params.UbdIX = UbdIX;
    params.UbdObjVal = evalPoly(objPoly, zeros(nDim, 1));
    params.delta1 = 1.0e-4;
    params.delta2 = 0;
    params.delta = 0;
    [sol, info] = solveCOP(Q0vec, H0vec, H1vec, PSDcone, polyCone, params);
    
    
    %% print
    fid = fopen(printFileName, 'a');
    % fprintf(fid,'degree,n,rRange,isBin,addComp,rOder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');    
    if strcmp(solver, 'BP')
        %         info.iterBP = infoBP.iter;
        %         info.timeBP = infoBP.timeBisection;
        %         info.iterAPGR = infoBP.APGiter;
        %         info.termcodeBP = infoBP.break_yes;
        fprintf(fid,'%d,%d,%f,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d,%d\n',...
            degree, nDim, radiorange, isBin, addComplement, params.relaxOrder, solver,...
            sol.LBv, info.timeBP, info.iterAPGR, info.termcodeBP, info.iterBP); 
    %     elseif strcmp(solver,'newBP')
    %         fprintf(fid,'%d,%d,%f,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d,%d\n',...
    %             degree, nDim, radiorange, isBin, addComplement, params.relaxOrder, solver,...
    %             sol.LBv, info.cputime, info.APGiter, info.break_yes, info.iter);
    elseif strcmp(solver, 'sdpnalplus') || strcmp(solver, 'sdpnal') || strcmp(solver, 'sedumi')
        fprintf(fid,'%d,%d,%f,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d\n',...
            degree, nDim, radiorange, isBin, addComplement, params.relaxOrder, solver,...
            sol.LBv, info.cputime, info.iter, info.termcode);
    end
    fclose(fid);

end

