(1)
BBCPOP.m can be used with any of general SDP solvers such as 
SDPNAL+, SeDuMi, SDPT3, SDPA, and MOSEK by specifying 

params.solver = 'sdpnalplus', 'sedumi', 'sdpt3', 'sdpa' and 'mosek’.
 
(The SDP solver should be installed.) 

For example:  

>>[objPoly, I01, Icomp] = example1;
>>relaxOrder = 2; params.solver='sdpnalplus';
>>[sol, info] = BCCPOP(objPoly,I01,Icomp,relaxOrder,params);

(2)
To solve randomly generated POP instances in the directory ``instances/POPrandom”,
binary QOP instances (from BIQ mac library) in the directory ``instances/BIQ”, 
and QAP instances (from QAPLIB) in the directory ``instances/QAP/qapdata”, 
the user can use the following functions:

advancedUse/experiment/POPrandom/
solvePOPdense(printFileName,degree,nDim,isBin,addComplement,solver);
% where params.solver can be one of ’BP','sdpnalplus' and 'sedumi’.

advancedUse/experiment/POPrandom/
solvePOParrow(printFileName,degree,a,b,c,l,isBin,addComplement,solver);
% where params.solver can be one of ’BP','sdpnalplus' and 'sedumi’.

advancedUse/experiment/POPrandom/
solvePOPchordal(printFileName,degree,nDim,radiorange,isBin,addComplement,solver); 
% where params.solver can be one of ’BP','sdpnalplus' and 'sedumi’.

advancedUse/experiment/BIQ
solveBIQ(printFileName,instance,solver,lambda)
% where params.solver can be one of ’BP' and ’sdpnalplus’.

advancedUse/experiment/QAP
solveQAP(printFileName,instance,solver,lambda); 
% where params.solver can be one of ’BP','sdpnalplus’and ’sdpnalAW’.

(3)
For numerical experiment of randomly generated POP instances:  
advancedUse/experiment/POPrandom/expPOP

(4) 
For numerical experiment of binary QOP instances: 
advancedUse/experiment/BIQ/expBIQ


(5)
For numerical experiment of QAP instances: 
advancedUse/experiment/QAP/expQAP
