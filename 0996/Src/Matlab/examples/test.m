% Test script for BBCPOP

reltol = 1e-1;
%% Example1
[objPoly, I01, Icomp] = example1;
relaxOrder = 2;
params = [];
[sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

if (abs(((-2.4701) - sol.LBv) / sol.LBv) > reltol)
    error('Large numerical error at Example 1.')
end
%% Dense POP example
degree=3; nDim=5; isBin=true; addComplement=false; 
[objPoly,I01,Icomp,relaxOrder,params] = genPOPdense(degree,nDim,isBin,addComplement); 
[sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);
if abs(((-5.4614) - sol.LBv) / sol.LBv) > reltol
    error('Large numerical error at Dense POP example.')
end
%% Arrow type sparse POP example 
degree=4; a=10; b=2; c=2; l=3; isBin=0; addComplement=true; 
[objPoly, I01, Icomp, relaxOrder, params] = genPOParrow(degree,a,b,c,l,isBin,addComplement);
[sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);
if abs(((-17.9045) - sol.LBv) / sol.LBv) > reltol
    error('Large numerical error at Arrow type sparse POP example.')
end
%% Chordal graph type sparse POP example: 
degree=3; nDim=100; radiorange=0.1; isBin=1; addComplement=false;
[objPoly, I01, Icomp, relaxOrder, params] = genPOPchordal(degree,nDim,radiorange,isBin,addComplement);
[sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);
if abs(((-54.5594) - sol.LBv) / sol.LBv) > reltol
    error('Large numerical error at Chordal graph type sparse POP example.')
end
%% QAP example:
instance='chr12a'; lambda=10000; 
[objPoly, I01, Icomp, relaxOrder, params] = qapreadBP(instance,lambda);
[sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);
if abs(((9551.2) - sol.LBv) / sol.LBv) > reltol
    error('Large numerical error at QAP example.')
end
%% BIQ example:
instance='bqp100-1'; lambda=10000;
[objPoly, I01, Icomp, relaxOrder, params] = biqreadBP(instance,lambda);
[sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);
if abs(((-8041.3) - sol.LBv) / sol.LBv) > reltol
    error('Large numerical error at BIQ example.')
end
%%
fprintf('\n');
fprintf('All tests were passed successfully!\n');
fprintf('Now you are ready for BBCPOP.\n');
