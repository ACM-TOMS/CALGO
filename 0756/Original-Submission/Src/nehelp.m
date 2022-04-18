%
% Help file for NESOLVE.M and its system of subroutines.
% To access the help information, execute this file by typing NEHELP
% from within MATLAB.
%
clc
echo on
%
%                    NONLINEAR EQUATION SOLVER
%
% INTRODUCTION
%
% NESOLVE and its system of subordinate functions is a software package
% designed to solve systems of nonlinear equations.  Everything must be
% REAL; the algorithms are not designed for complex numbers (complex
% functions could be handled by breaking each equation and each variable
% into two parts, real and imaginary, and solving a system twice as large).
% Newton's method, with a few modifications, is used to find the solution.
% A line search has been included to ensure global convergence (convergence
% from poor initial guesses).  The necessary derivatives of the function are
% computed using finite differences unless the user supplies a function
% which computes them from analytically obtained equations.
%
% Detailed descriptions of the algorithmic modules as well as the theory
% behind them may be found in "Numerical Methods for Unconstrained
% Optimization and Nonlinear Equations" by J. E. Dennis, Jr. and
% R. B. Schnabel, 1983.  The user is refered to that book if such
% information is needed.
%
pause
clc
%
% HOW TO USE NESOLVE
%
% To use the Nonlinear Equations package, you must first create a function
% file xxxxxx.M which evaluates the simultaneous nonlinear equations you
% wish to solve.  You may name this function anything you like, just supply
% the name (in single quotes) as the FVEC argument to NESOLVE.  To develop
% the function file, first write each equation as
%
%                      'expression' = 0.
%
% Next write a function file which evaluates each of these 'expressions' and
% returns the results as an n-vector F(X), where the input X is also an
% n-vector.  Note that the number of variables and the number of equations
% must be the same, but that dummy variables or trivial equations (0=0)
% could be added to accomplish this.  It is possible to write some functions
% in such a way that the number of variables, n, is determined at run time
% by the length of the vector X0, the starting point (see netestf4).
%
% A solution is any vector X for which F(X) = 0.
%
pause
clc
%                    NESOLVE INPUT ARGUMENTS
%
% REQUIRED INPUTS
%
% FVEC, the first input argument, must be the name of the function file
% defining the system of nonlinear equations (as described on the previous
% screen).  The name must be placed in single quotation marks or contained
% in a variable as a string.  The '.m' should not be included.
%
% X0, the second input augument, must be a column-vector containing an
% initial guess of the solution.  The size of X0 must be compatible with
% the function named in FVEC.
%
pause
clc
%
% OPTIONAL INPUTS
%
% NESOLVE has up to four optional input arguments.  To supply any particular
% optional argument, all arguments preceding it in the input argument list
% must also be supplied, but may be empty matrices to hold the positions.
%
% DETAILS, the third input argument, is an optional 16-vector whose elements
% specify various algorithmic options and set various tolerances.  Default
% values are assumed if DETAILS is not present or if an element is set to
% zero.  If DETAILS is present, but has less than 16 elements, the remaining
% elements assume their default values.  The next screen lists the specific
% function of each element of DETAILS.
%
pause
clc
%
% DETAILS elements (defaults are in square brackets [ ]):
%
% ELEMENT   NAME        DESCRIPTION
% -------   ---------   ---------------------------------------------------
%    1      PRINTCODE   [0]=No trace;  1=Trace;  2=Trace & Statistics.
%    2      GLOBMETH    [1]=Line Search;  2=Hookstep; (3-4 reserved).
%    3      CHEAPF      [0]=No (Use secant update if function is expensive to
%                       evaluate);  1=Yes (Always use finite differences).
%    4      ANALJAC     [0]=No; 1=Yes, there is a jacobian function (.M file)
%                       (if so, supply its name in input argument JAC).
%    5      FACTSEC     [0]=No; (1 reserved for future use).
%    6      ITNLIMIT    Maximum number of iterations allowed [100].
%    7      DELTA       Initial trust radius for GLOBMETH=2,3 [cauchy step].
%    8      FVECTOL     How small the scaled norm(F) must get [eps^(1/3)].
%    9      STEPTOL     Minimum step size [eps^(2/3)].
%   10      MINTOL      For detecting non-root minima [eps^(2/3)].
%   11      MAXSTEP     Largest allowed step size [see NEINCK.M].
%   12      FDIGITS     Number of good digits returned by FVEC [-log10(eps)].
%   13      ETA         ^* Internal use only ^*
%   14      SAVEPATH    ^* Internal use only ^*
%   15      PASSPARM    ^* Internal use only ^*
%   16      SCALEFLG    [0]=No;  1=Yes, by starting point; 2=Yes, by SCALE.
%                       (use if units in different variables are mismatched)
pause
clc
%
% OPTIONAL INPUT ARGUMENTS (CONTINUED)
%
% FPARAM, the fourth input argument, is an optional input.  NESOLVE does
% not care what it contains, but if it exists and is nonempty it will be
% passed on to FVEC (and JAC) as a second argument.  This feature allows
% the function to contain parameters which are held constant during the
% solution process but may be changed for the next run.  Any number of
% parameters of any kind may be passed, but it is up to the user to pack
% them all into one matrix FPARAM and up to the user-written function FVEC
% (and JAC) to unpack them and use them appropriately.
%
% JAC, the fifth input argument, is optional.  To use a function file
% which evaluates all the first derivatives of the function, the name of
% that file must be supplied in JAC and the flag DETAILS(4) must be set
% to one.  The output of JAC should be the Jacobian matrix of the function
% specified in FVEC, evaluated at the input argument X.  The i,j element of
% the Jacobian matrix is defined as the derivative of function element i
% with respect to variable j.  JAC needs to provide a second output
% variable indicating the computational cost of one call to JAC in terms
% of the equivalent number of calls to FN (this is only used for reporting
% statistics on the number of function evaluations when PRINTCODE=2).
%
pause
clc
%
% OPTIONAL INPUT ARGUMENTS (CONTINUED)
%
% SCALE, the sixth and last input argument, is optional.  If supplied, it
% should be an (n x 2) matrix whose first column contains a 'typical' X
% vector and whose second column contains a 'typical' F vector.  SCALE must
% not contain any zeros.  It is used to improve convergence behavior in
% in situations where the units of the variables are such that some variables
% are typically several orders of magnitude larger than others.  The flag
% controling scaling is DETAILS(16) which must be set to 2 if SCALE is
% to be used.  If DETAILS(16) is 1, X0 and F(X0) are used for scaling.
%
pause
clc
%
% OUTPUT ARGUMENTS OF NESOLVE
%
% XF is the primary output of NESOLVE.  It contains the final approximation
% of the solution to the system of nonlinear equations.
%
% TERMCODE is an optional output (though highly recommended).  If a second
% output variable is supplied, It will contain the termination code giving
% the reason the iteration process was stopped.  The reasons are outlined
% below.  For additional information the user is refered to the files
% NESTOP.M, NESTOP0.M and NEINCK.M or to the book referenced on screen one
% of this help file.
%
%    TERMINATION CODES
%    -2 : Input error in Fdigits (DETAILS(12)).
%    -1 : Input error in starting point (X0).
%     1 : Normal termination, XF is probably near a root unless fvectol
%         (DETAILS(8)) is too large.
%     2 : Two steps too small ( < steptol), maybe near a root.
%     3 : Can't find a good step. Maybe near a root or jacobian inacurate.
%     4 : Iteration limit exceeded. Increase DETAILS(6) or restart from XF.
%     5 : Five steps too big ( > maxstep), looks like asymtotic behavior.
%     6 : Stuck at a minimizer which is not a root.  Restart with new X0.
pause
clc
%
% OUTPUT ARGUMENTS (CONTINUED)
%
% PATH is an optional output.  If a third output variable is supplied it will
% contain the sequence of points generated by the iterations.  The i-th row
% of PATH is the transpose of vector X at iteration (i-1), so you can use
% the command plot(path) to see how the variables converge.
%
echo off

