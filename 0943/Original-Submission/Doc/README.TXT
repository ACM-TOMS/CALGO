The main driver for this code is "demo.m"

To run: type "demo" at the MATLAB command line

demo.m calls the following functions:

mss.m  -- the algorithm proposed in this paper (More'-Sorensen sequential method)
ms.m   -- an implementation of the More'-Sorensen method for L-BFGS matrices

These methods make use of the following routines:
unrolling.m  
two_loop.m
recursionSolve_qN.m  (only used by mss.m)
qN_B.m  (only used by ms.m)  
qN_multiply.m (only used when "show=1" in mss.m and ms.m)