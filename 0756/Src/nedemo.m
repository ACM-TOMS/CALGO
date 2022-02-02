%
% NEDEMO demonstrates the Nonlinear Equations sovler NESOLVE.
%
% Written by Richard T. Behrens, July 1988.
%
format compact
echo on
clc
%
% Nonlinear Equations DEMO:  Demonstration of the use of NESOLVE.
%
%
% First, we need a set of simultaneous nonlinear equations whose solution
% we wish to compute numerically.  Let's sovle the 'helical valley function'
% defined as:
%
%                       0 = 10 * (x3 - 10*theta)
%                       0 = 10 * (sqrt(x1^2 + x2^2) - 1)
%                       0 = x3
%
%                       where:  theta = (1/(2*pi)) * atan(x2/x1)
%
%
% We need to write a function file to define the above set of equations.
% Given a 3-vector X, it should evaluate the right-hand sides and return
% a 3-vector F, which will be zero when a solution is found.  Both X and
% F need to be column-vectors (i.e. n x 1, NOT 1 x n).
%
pause
clc
%
% Here is the function file:
%
%
type netestf1
pause
clc
%
% Next, we need an ititial guess of the solution we are looking for.
% Nonlinear equations can have multiple solutions, so the starting point
% should be chosen as close as possible to the desired solution (it is
% also possible to have no solutions at all).  Since we have no idea where
% the solution is, we choose an arbitrary starting point. But note that
% certain starting points like [0;0;0] are unsuitable in this case, because
% the 'helical valley function' is undefined there.
%
% Here is our starting point:
%
echo off
x0 = [10;10;10]
pause
clc
echo on
% The function file and the starting point are the only required inputs
% to NESOLVE.  We will specify additional inputs (DETAILS) which will
% cause the intermediate results to be printed and tell the package
% that the function is cheap (quick) to evaluate.
echo off
details = zeros(16,1);
details(1) = 1;
details(3) = 1
disp(' ')
disp('Press any key to begin finding a solution . . .')
pause
%
clc
disp('Please stand by ... loading & compiling the functions may take a while.')
flops(0);
time1 = clock;
[xf,termcode,path] = nesolve('netestf1',x0,details);
time2 = clock;
clc
f = flops;
disp('The solution is:  ')
xf
disp('The number of flops required to find it was:  ')
f
disp('The time it took was:  ')
etime(time2,time1)
disp(' ')
disp('Did you notice that the third component converged in only')
disp('one iteration?  That is because that equation was linear.')
disp(' ')
disp('Press any key to see a plot of the sequence of iterates . . .')
pause
clc
plot(path)
title('Convergence Path')
pause

