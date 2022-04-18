function solvePOParrow(printFileName,degree,a,b,c,l,isBin,addComplement,solver)
%SOLVEPOPARROW Solve a POP having arrow-type sparsity
% Usage
%   solvePOParrow(printFileName,degree,a,b,c,l,isBin,addComplement,solver); %This add the result to solvePOParrow.csv.
% Input
%   printFileName : a file name for output, for exampe, 'temp', then output is written in 'temp.csv'  
%   degree  : the degree of the polynomial objective function
%   a       : the block size
%   b       : the number of conseccutive variables 
%   c       : the number of common variables over all the blocks
%   l       : blockSize, the number of blocks
%   isBin   : 1 for the binary constraint, 0 for the box constraint
%   addComplement   : 1 for adding random complementarity constraints 
%   solver   : 'BP, 'sdpnalplus', 'sedumi'
%
% a, b, c, l determines the Hessian matrix of the objective function;
% the size of the Hessian matrix = the number of variabls is given by 
%     n = (a - b) * l + b + c;
% for example, if a = 3, l = 2, b = 1 and c = 2, then the Hessian matrix 
% if 7 x 7 matorix of the form 
%  ***00**
%  ***00**
%  *******
%  00*****
%  00*****
%  *******
%  *******
% a, b, c and l must satisfy 
% 0 <= b < a, 0 <= c, 2 <= l
%
% block-diagonal-type
% if a >=1, b = 0, c = 0 and l >= 2. 
% simple arrow-type
% if a = 1, b = 0, c >= 1 and l >= 2. 
% block-bandwidth-type
% if 2 <= a, 1 <= b < a, c = 0 and l >= 2. 
%
% For the DNN relaxation, see 
% N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
% "BBCPOP: A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization Problems 
% with Binary, Box and Complementarity Constraints," arXiv:1804.00761, 
% to appear in ACM Transaction on Mathematical Software.
% 
% Example:  
% >> printFileName='temp'; degree=3; a=10; b=2; c=2; l=10; isBin=1; addComplement=0; solver='BP'; 
% >> solvePOParrow(printFileName,degree,a,b,c,l,isBin,addComplement,solver);
%

%     randSeed = 2017;
      n = (a - b) * l + b + c;
%     I01 = isBin .* true(1, n);
%     probData.nDim = n;
%     probData.degree = d;
%     probData.randSeed = randSeed;
%     probData.blockSize = a;
%     probData.consecVarSize = b;
%     probData.commonVarSize = c;
%     probData.noBlocks = l;
%     [objPoly, clique] = genPOParrow(probData, I01);
%     
%     Icomp = [];
%     if addComplement
%         Icomp = complementarityConstraint(2, 2, clique, randSeed);
%     end

    if nargin == 8
        solver = 'BP'; 
    end

    % [objPoly, I01, Icomp]=genPOParrow2(d,a,b,c,l,isBin,addComplement)
    [objPoly, I01, Icomp] = genPOParrow(degree,a,b,c,l,isBin,addComplement); 
    
    params.solver = solver;   
    params.newBPswitch = false;
    
    %% to DNN
    eqPolySys = [];
    params.relaxOrder = ceil(degree / 2);
    [Q0vec, H0vec, H1vec, PSDcone, polyCone, UbdIX] ...
            = BBPOPtoCOP(objPoly, eqPolySys, I01, Icomp, params);

    %% solve
    params.UbdIX = UbdIX;
    params.UbdObjVal = evalPoly(objPoly, zeros(n, 1));
    params.delta1 = 1.0e-4;
    params.delta2 = 0;
    params.delta = 0;
    [sol, info] = solveCOP(Q0vec, H0vec, H1vec, PSDcone, polyCone, params);
        
    %% print
    filename = printFileName;
    fid = fopen(filename, 'a');
    % fprintf(fid,'degree,n,a,b,c,l,isBin,addComp,rOrder,solver,LBv,sec,iter(APGR),termcode,iterBP\n');
    if strcmp(solver, 'BP')
        %         info.iterBP = infoBP.iter;
        %         info.timeBP = infoBP.timeBisection;
        %         info.iterAPGR = infoBP.APGiter;
        %         info.termcodeBP = infoBP.break_yes;
        fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d,%d\n',...
            degree, n, a, b, c, l, isBin, addComplement, params.relaxOrder, solver,...
            sol.LBv, info.timeBP, info.iterAPGR, info.termcodeBP, info.iterBP);
    elseif strcmp(solver, 'sdpnalplus') || strcmp(solver, 'sdpnal') || strcmp(solver, 'sedumi')
        fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%1.8e,%1.2e,%d,%d,\n',...
            degree, n, a, b, c, l, isBin, addComplement, params.relaxOrder, solver,...
            sol.LBv, info.cputime, info.iter, info.termcode);
    end
    fclose(fid);

end

