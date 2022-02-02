function solveBIQ(printFileName,instance, solver, lambda)
%SOLVEBIQ Solve Binary QOP
% Usage
%   solveBIQ(printFileName,instance,solver,lambda);
% Input
%   printFileName : a file name for output, for exampe, 'temp', then output is written in 'temp.csv'  
%   instance : the instance name in the directory instances/BIQ/BIQdata
%              for example, instance = 'chr12a'; 
%   solver   : 'BP, 'sdpnalplus', 'sdpnalAW'
%   lambda   : the Lagrangian multiplier, for example lambda = 1.0e5
% All instances in the directory instances/BIQ/BIQdata are from BIQMAQ: 
%
% For the DNN relaxation, see 
% S. Kim, M. Kojima, and K.-C. Toh.
% "A Lagrangian-DNN relaxation: a fast method for computing tight lower bounds 
% for a class of quadratic optimization problems," 
% Math. Program., 156 (2016) 161-187.
%
% Example:  
% >> printFileName='temp'; instance='bqp100-1'; solver='BP'; lambda=100000; 
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
isBinary = 1;

switch solver
    case 'BP'
        %% solve our DNN by BP
        % Lagrangian relaxation
        % lambda = 1e5;
        [objPoly, I01, Icomp, relaxOrder, params] = biqreadBP(instance,lambda,isBinary);
        [Q0vec, H0vec, PSDcone, polyCone, UbdIX] ...
            = BBCPOPtoDNN(objPoly, I01, Icomp, relaxOrder, params);
        params.solver = 'BP';        
        [sol, info] = solveDNN(Q0vec, H0vec, PSDcone, polyCone, UbdIX, params);
        eTime =info.timeBP; 
        iter = info.iterAPGR;
        termcode = info.termcodeBP;
        iterBP = info.iterBP;
        optVal = params.optVal;
    case {'sdpnal','sdpnalplus'}
        %% solve our DNN by sdpnal+ without Lagrangian relaxation
        [objPoly, I01, Icomp, relaxOrder, params] = biqreadBP(instance,lambda,isBinary);
        [Q0vec, H0vec, PSDcone, polyCone, UbdIX] ...
            = BBCPOPtoDNN(objPoly, I01, Icomp, relaxOrder, params);
        params.solver = 'sdpnalplus';
        [sol, info] = solveDNN(Q0vec, H0vec, PSDcone, polyCone, UbdIX, params);
        eTime =info.totaltime; 
        iter =info.iter;
        termcode = info.termcode;
        iterBP = 0; 
        optVal = [];
end
%%
fid = fopen(printFileName, 'a');
% fprintf(fid, 'instance,optVal,method,lambda,LB,LBv,UB,time,iter(APGR),termcode,iterBP\n');
fprintf(fid, '%s,%1.8e,%s,%1.2e,%1.8e,%1.8e,%1.8e,%1.2e,%d,%d,%d\n',...
    instance,optVal,solver,lambda,sol.LB,sol.LBv,sol.UB,eTime,iter,termcode,iterBP);
fclose(fid);

end