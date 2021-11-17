function solvePOPdense(printFileName,degree,nDim,isBin,addComplement,solver)
%SOLVEPOPDENSE Solve a dense POP
% Usage
%   solvePOPdense(printFileName,degree,nDim,isBin,addComplement,solver); 
% Input
%   printFileName : a file name for output, for exampe, 'temp', then output is written in 'temp.csv' 
%   degree  : the degree of the polynomial objective function
%   nDim    : the dimension of the polynomial objective function
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
% >> printFileName='temp'; degree=3; nDim=10; isBin=1; addComplement=0; solver='BP'; 
% >> solvePOPdense(printFileName,degree,nDim,isBin,addComplement,solver); 
%

    if nargin == 5
        solver = 'BP'; 
    end
    
    [objPoly,I01,Icomp] = genPOPdense(degree,nDim,isBin,addComplement); 
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
    params.delta1 = 1.0e-4;
    params.delta2 = 0;
    params.delta = 0;
    params.UbdIX = UbdIX;
    params.UbdObjVal = evalPoly(objPoly, zeros(nDim, 1));
    %     params.deltaRelBP = 1.0e-4;
    %     params.deltaAbsBP = 0;
    %     params.deltaDifRelBP = 0;
    [sol, info] = solveCOP(Q0vec, H0vec, H1vec, PSDcone, polyCone, params);
    
    
    %% print
    fid = fopen(printFileName, 'a');
    % fprintf(fid,'degree,n,isBin,addComp,rOrder,solver,LBv,sec,iter,termcode\n');
    if strcmp(solver, 'BP')
        %         info.iterBP = infoBP.iter;
        %         info.timeBP = infoBP.timeBisection;
        %         info.iterAPGR = infoBP.APGiter;
        %         info.termcodeBP = infoBP.break_yes;
        fprintf(fid,'%d,%d,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d,%d\n',...
            degree, nDim, isBin, addComplement, params.relaxOrder, solver,...
            sol.LBv, info.timeBP, info.iterAPGR, info.termcodeBP, info.iterBP);
    %     elseif strcmp(solver,'newBP')
    %         fprintf(fid,'%d,%d,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d,%d\n',...
    %             degree, nDim, isBin, addComplement, params.relaxOrder, solver,...
    %             sol.LBv, info.cputime, info.APGiter, info.break_yes, info.iter);
    elseif strcmp(solver, 'sdpnalplus')  || strcmp(solver, 'sdpnal') || strcmp(solver, 'sedumi')
        fprintf(fid,'%d,%d,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d\n',...
            degree, nDim, isBin, addComplement, params.relaxOrder, solver,...
            sol.LBv, info.cputime, info.iter, info.termcode);
    end
    fclose(fid);

end

