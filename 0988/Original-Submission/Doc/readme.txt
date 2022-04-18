AMGKQ: [RES, ERR, NSUB, FL] = amgkq(F, A, B, C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, CFLAG, VERB, P1, P2, ...)

DESCRIPTION: Adaptive, multidimensional Gauss-Kronrod quadrature for simultaneous integrands.
Computes the integral of F(X) from A to B in ND dimensions.  Improper integrals with infinite
limits A or B are supported.  Limits A and B do not have to be ordered low to high. Complex
line (contour) integrals require ND = 1 and finite limits A and B.  Multiple simultaneous
integrands are supported, e.g. F(X) = [F1(X); F2(X); ...].  The number of dimensions ND and
the number of integrands NF are not limited. 

INPUT: (F, A, B) are required; (C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, CFLAG, VERB, P1, P2, ...) are optional.
Pass an empty matrix [] to select the default value of an optional input argument.
F = integrand F(X) to be described later
[A, B] = limits of integration as column vectors with ND rows
C = matrix of size [ND, NC] giving breakpoint (boundary) locations, default to (A + B)/2
EAER = requested error with one or two elements, e.g. EAER = EA or EAER = [EA, ER]
  EA = requested absolute error >= 0, default to sqrt(eps) < 1.5e-08
  ER = requested relative error >= 0, default to 0, TOL = max(EA, ER*ABSRES0)
MAXNSUB = maximum number of subregions to evaluate, default to ND*1e2, required > NC*2^ND - NC + 1
NGK = order of Gauss-Kronrod quadrature rules, default to 7, required > 1
TTYPE = trigonometric (1) or rational (2) transformation type, e.g. TTYPE = TBTH or TTYPE = [TFIN, TINF]
  TFIN = type for finite limit transformation (both A and B), default to 1 = trigonometric
  TINF = type for infinite limit transformation (either or both A and B), default to 1 = trigonometric
  TBTH = type for both finite and infinite limit transformations
SFLAG = boolean flag for subregion culling, e.g. SFLAG = SBTH or SFLAG = [SUB1, SUB2]
  SUB1 = 0 disables culling of subregions meeting convergence criterion, default to 1
  SUB2 = 0 disables culling of subregions surpassing width criterion, default to 1
  SBTH = 0 disables all culling of subregions
CFLAG = boolean flag for reordering C, CFLAG = 0 disables reordering, default to 1
VERB = integer verbosity level from -2 to 2, default to 0
     = -2 turns off termination, culling, and tolerance warnings
     = -1 turns off termination warnings
     = 0 leaves warnings on with no iteration output
     = 1 turns on some iteration output
     = 2 turns on more iteration output
[P1, P2, ...] = additional parameters passed to F when given as a string

QUADRATURE RULES: For order NGK in {7, 10, 15, 20, 25, 30}, values tabulated to 1e-25 are loaded
which get rounded to double precision ~ 2e-16.  For other orders, a double precision routine is
called with round-off error ~ 2e-16.  The Gauss and Kronrod orders are NG = NGK and NK = 2*NGK+1.

INTEGRAND: F can be given as an inline function, a function handle, or a string for a function name.
When given as a string, F is converted to a function handle with additional parameters [P1, P2, ...].
F must accept a matrix X of size [ND, NX] as input and return a matrix Y of size [NF, NX] as output.
Expected calling format is Y = F(X), for example F = @(X) somefunction(X, P1, P2, ...).

OUTPUT: [RES, ERR, NSUB, FL] have the following descriptions.
RES = column vector of results with NF rows 
ERR = column vector of estimated errors with NF rows 
NSUB = scalar integer giving number of subregions evaluated
FL = output flag in {-2, -1, 0, 1, 2}
   = -2 means NaN encountered in integrand
   = -1 means Inf encountered in integrand
   = 0 means maximum number of subregions evaluated
   = 1 means all subregions meet convergence criterion
   = 2 means estimated errors meet convergence criterion

REFERENCES:
  P. van Dooren and L. de Ridder, Algorithm 6: 
    An adaptive algorithm for numerical integration over 
    an N-dimensional cube, J. Comput. Appl. Math., 2 (1976) 207-217.
  A. C. Genz and A. A. Malik, Algorithm 019. Remarks on algorithm 006:
    An adaptive algorithm for numerical integration over an
    N-dimensional rectangular region, J. Comput. Appl. Math., 6 (1980) 295-302.
  J. Berntsen, T. O. Espelid, and A. Genz, An adaptive Algorithm 
    for the approximate calculation of multiple integrals,
    ACM Trans. Math. Softw., 17 (1991) 437-451.
  L. F. Shampine, Vectorized adaptive quadrature in MATLAB, 
    J. Comput. Appl. Math., 211 (2008) 131-140.
  W. Gautschi, Algorithm 726: ORTHPOL -- a package of routines for generating
    orthogonal polynomials and Gauss-type quadrature rules, ACM Trans. Math.
    Softw., 20 (1994) 21-62.
  W. Gautschi, Orthogonal Polynomials: Computation and Approximation,
    Clarendon Press, Oxford, 2004.
  Dirk P. Laurie, Calculation of Gauss–Kronrod Quadrature Rules,
    Math. Comp., 66 (1997) 1133-1145.
  Pavel Holoborodko, Gauss-Kronrod Quadrature Nodes and Weights, (2011).
    http://www.advanpix.com/2011/11/07/Gauss-Kronrod-quadrature-nodes-weights/

Version 1.0.0
Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
This is free software; see GPLv3 or greater for copying conditions.
There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
Based on ADAPT.M (C) 2013 Alan Genz and QUADGK.M (C) 2013 David Bateman.
Developed under Octave 3.8.1 with no packages loaded.

INSTRUCTIONS:
Copy files from folders Driver/ and Src/ into a single folder for testing.
For either operating environment, type "test_amgkq" to verify functionality.
Some demonstrations of AMGKQ in action are given in the example*.m files.
Exact values from closed form expressions are given in the example*.html files
located in this (the Doc/) folder.  All tests and examples can be executed,
with output to a diary file, by running the command "runme".

