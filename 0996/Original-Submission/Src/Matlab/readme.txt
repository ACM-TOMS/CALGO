** To install **

1) Download BBCPOPv1.0.tar.gz ad unpack it. 

2) Set MATLAB-path to BBCPOPv1.0 and its subdirectories.  
(Or execute startup). 

3) Executable mex functions have been prepared for 64bit Windows, OS X, and Linux.
If they don’t work, recompile with recompile.m.

4) Now you are ready to use BBCPOP. Test it by typing
>> [objPoly, I01, Icomp, relaxOrder, params] = example1;
>> [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

To test BBCPOP on more examples listed below, type 
>> test

** Some other examples. 

Dense POP example:  
>> degree=3; nDim=5; isBin=true; addComplement=false; 
>> [objPoly,I01,Icomp,relaxOrder,params] = genPOPdense(degree,nDim,isBin,addComplement); 
>> [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

Arrow type sparse POP example: 
>> degree=4; a=10; b=2; c=2; l=3; isBin=0; addComplement=true; 
>> [objPoly, I01, Icomp, relaxOrder, params] = genPOParrow(degree,a,b,c,l,isBin,addComplement);
>> [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

Chordal graph type sparse POP example: 
>> degree=3; nDim=100; radiorange=0.1; isBin=1; addComplement=false;
>> [objPoly, I01, Icomp, relaxOrder, params] = genPOPchordal(degree,nDim,radiorange,isBin,addComplement);
>> [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

QAP example:
>> instance='chr12a'; lambda=100000; 
>> [objPoly, I01, Icomp, relaxOrder, params] = qapreadBP(instance,lambda);
>> [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

BIQ example:
>> instance='bqp100-1'; lambda=10000;
>> [objPoly, I01, Icomp, relaxOrder, params] = biqreadBP(instance,lambda);
>> [sol, info] = BBCPOP(objPoly, I01, Icomp, relaxOrder, params);

Help messages are available for most functions. 
Try, for example, 
>> help BBCPOP

References

N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
"BBCPOP: A Sparse Doubly Nonnegative Relaxation of Polynomial Optimization 
Problems with Binary, Box and Complementarity Constraints," arXiv:1804.00761, 
https://arxiv.org/abs/1804.00761.

N. Ito, S. Kim, M. Kojima, A. Takeda, and K.-C. Toh.
User Manual for BBCPOP:  A Sparse Doubly Nonnegative Relaxation
of Polynomial Optimization Problems with Binary, Box and Complementarity 
Constraints. https://sites.google.com/site/bbcpop1/download.