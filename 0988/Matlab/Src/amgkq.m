function [RES, ERR, NSUB, FL] = amgkq(F, A, B, C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, CFLAG, VERB, varargin)
% AMGKQ: [RES, ERR, NSUB, FL] = amgkq(F, A, B, C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, CFLAG, VERB, P1, P2, ...)
%
% DESCRIPTION: Adaptive, multidimensional Gauss-Kronrod quadrature for simultaneous integrands.
% Computes the integral of F(X) from A to B in ND dimensions.  Improper integrals with infinite
% limits A or B are supported.  Limits A and B do not have to be ordered low to high. Complex
% line (contour) integrals require ND = 1 and finite limits A and B.  Multiple simultaneous
% integrands are supported, e.g. F(X) = [F1(X); F2(X); ...].  The number of dimensions ND and
% the number of integrands NF are not limited. 
%
% INPUT: (F, A, B) are required; (C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, P1, P2, ...) are optional.
% Pass an empty matrix [] to select the default value of an optional input argument.
% F = integrand F(X) to be described later
% [A, B] = limits of integration as column vectors with ND rows
% C = matrix of size [ND, NC] giving breakpoint (boundary) locations, default to (A + B)/2
% EAER = requested error with one or two elements, e.g. EAER = EA or EAER = [EA, ER]
%   EA = requested absolute error >= 0, default to sqrt(eps) < 1.5e-08
%   ER = requested relative error >= 0, default to 0, TOL = max(EA, ER*ABSRES0)
% MAXNSUB = maximum number of subregions to evaluate, default to ND*1e2, required > NC*2^ND - NC + 1
% NGK = order of Gauss-Kronrod quadrature rules, default to 7, required > 1
% TTYPE = trigonometric (1) or rational (2) transformation type, e.g. TTYPE = TBTH or TTYPE = [TFIN, TINF]
%   TFIN = type for finite limit transformation (both A and B), default to 1 = trigonometric
%   TINF = type for infinite limit transformation (either or both A and B), default to 1 = trigonometric
%   TBTH = type for both finite and infinite limit transformations
% SFLAG = boolean flag for subregion culling, e.g. SFLAG = SBTH or SFLAG = [SUB1, SUB2]
%   SUB1 = 0 disables culling of subregions meeting convergence criterion, default to 1
%   SUB2 = 0 disables culling of subregions surpassing width criterion, default to 1
%   SBTH = 0 disables all culling of subregions
% CFLAG = boolean flag for reordering C, CFLAG = 0 disables reordering, default to 1
% VERB = integer verbosity level from -2 to 2, default to 0
%      = -2 turns off termination, culling, and tolerance warnings
%      = -1 turns off termination warnings
%      = 0 leaves warnings on with no iteration output
%      = 1 turns on some iteration output
%      = 2 turns on more iteration output
% [P1, P2, ...] = additional parameters passed to F when given as a string
%
% QUADRATURE RULES: For order NGK in {7, 10, 15, 20, 25, 30}, values tabulated to 1e-25 are loaded
% which get rounded to double precision ~ 2e-16.  For other orders, a double precision routine is
% called with round-off error ~ 2e-16.  The Gauss and Kronrod orders are NG = NGK and NK = 2*NGK+1.
%
% INTEGRAND: F can be given as an inline function, a function handle, or a string for a function name.
% When given as a string, F is converted to a function handle with additional parameters [P1, P2, ...].
% F must accept a matrix X of size [ND, NX] as input and return a matrix Y of size [NF, NX] as output.
% Expected calling format is Y = F(X), for example F = @(X) somefunction(X, P1, P2, ...).
%
% OUTPUT: [RES, ERR, NSUB, FL] have the following descriptions.
% RES = column vector of results with NF rows 
% ERR = column vector of estimated errors with NF rows 
% NSUB = scalar integer giving number of subregions evaluated
% FL = output flag in {-2, -1, 0, 1, 2}
%    = -2 means NaN encountered in integrand
%    = -1 means Inf encountered in integrand
%    = 0 means maximum number of subregions evaluated
%    = 1 means all subregions meet convergence criterion
%    = 2 means estimated errors meet convergence criterion
%
% REFERENCES:
%   P. van Dooren and L. de Ridder, Algorithm 6: 
%     An adaptive algorithm for numerical integration over 
%     an N-dimensional cube, J. Comput. Appl. Math., 2 (1976) 207-217.
%   A. C. Genz and A. A. Malik, Algorithm 019. Remarks on algorithm 006:
%     An adaptive algorithm for numerical integration over an
%     N-dimensional rectangular region, J. Comput. Appl. Math., 6 (1980) 295-302.
%   J. Berntsen, T. O. Espelid, and A. Genz, An adaptive Algorithm 
%     for the approximate calculation of multiple integrals,
%     ACM Trans. Math. Softw., 17 (1991) 437-451.
%   L. F. Shampine, Vectorized adaptive quadrature in MATLAB, 
%     J. Comput. Appl. Math., 211 (2008) 131-140.
%   W. Gautschi, Algorithm 726: ORTHPOL -- A package of routines for generating
%     orthogonal polynomials and Gauss-type quadrature rules, ACM Trans. Math.
%     Softw., 20 (1994) 21-62.
%   W. Gautschi, Orthogonal Polynomials: Computation and Approximation,
%     Clarendon Press, Oxford, 2004.
%   Dirk P. Laurie, Calculation of Gauss–Kronrod Quadrature Rules,
%     Math. Comp., 66 (1997) 1133-1145.
%   Pavel Holoborodko, Gauss-Kronrod Quadrature Nodes and Weights, (2011).
%     http://www.advanpix.com/2011/11/07/Gauss-Kronrod-quadrature-nodes-weights/
%
% Version 1.0.0
% Copyright (C) 2014 Robert W. Johnson, Alphawave Research.
% This is free software; see GPLv3 or greater for copying conditions.
% There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  For details, see GPLv3 or greater.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Based on ADAPT.M (C) 2013 Alan Genz and QUADGK.M (C) 2013 David Bateman.
% Developed under Octave 3.8.1 with no packages loaded.

% Sanity checks and initialization
if (nargin < 3), help amgkq, return, end
if isempty(F)
    error('amgkq: integrand F(X) must be specified'), end
if isempty(A) || isempty(B)
    error('amgkq: limits A and B must be specified'), end
[RA, CA] = size(A); [RB, CB] = size(B);
if any([CA, CB] > 1) || (RA ~= RB)
    error('amgkq: A and B must both be column vectors with ND rows')
else ND = RA; end
if (nargin < 4) || isempty(C), C = (A + B)/2; DEFC = true;
else DEFC = false; end
[RC, NC] = size(C);
if (RC ~= ND)
    error('amgkq: C must be a matrix with ND rows'), end
NCNSUB = NC*2^ND - NC + 1;
if (nargin < 5) || isempty(EAER), EAER = [sqrt(eps), 0]; end
if (numel(EAER) > 2)
    error('amgkq: EAER must be of the form EAER = EA or EAER = [EA, ER]'), end
if (numel(EAER) == 2), EA = EAER(1); ER = EAER(2);
else EA = EAER(1); ER = 0; end
if (EA < 0) || (ER < 0)
    error('amgkq: EA and ER must both be non-negative >= 0'), end
if (nargin < 6) || isempty(MAXNSUB), MAXNSUB = 2^ND*1e2; end
if ~isscalar(MAXNSUB) || (MAXNSUB ~= round(MAXNSUB)) || (MAXNSUB <= NCNSUB)
    error('amgkq: MAXNSUB must be a scalar integer > NC*2^ND - NC + 1'), end
ISCONTOUR = ~isreal([A, B, C]);
if ISCONTOUR && (ND > 1)
    error('amgkq: complex line integrals must have ND = 1'), end
if ~ISCONTOUR && any(A == B)
    error('amgkq: A and B values must not be equal for real integrals'), end
checkACB = sign(bsxfun(@minus,C,A).*bsxfun(@minus,B,C));
if ~ISCONTOUR && any(checkACB(:) < 0)
    error('amgkq: all C must be between A and B inclusive for real integrals'), end
ABINF = isinf([A, B]);
if ISCONTOUR && any(ABINF(:))
    error('amgkq: complex integrals with infinite limits are not supported'), end
if (nargin < 7) || isempty(NGK), NGK = 7; end
if ~isscalar(NGK) || (NGK ~= round(NGK)) || (NGK < 2)
    error('amgkq: NGK must be a scalar integer > 1'), end
if (nargin < 8) || isempty(TTYPE), TTYPE = [1, 1]; end
if (numel(TTYPE) > 2)
    error('amgkq: TTYPE must be of the form TTYPE = TBTH or TTYPE = [TFIN, TINF]'), end
if (numel(TTYPE) == 2), TFIN = TTYPE(1); TINF = TTYPE(2);
else TFIN = TTYPE; TINF = TTYPE; end
if ~any(TFIN == [1, 2]) || ~any(TINF == [1, 2])
    error('amgkq: TTYPE values must be 1 or 2'), end
if (nargin < 9) || isempty(SFLAG), SFLAG = [1, 1]; end
if (numel(SFLAG) > 2)
    error('amgkq: SFLAG must be of the form SFLAG = SBTH or SFLAG = [SUB1, SUB2]'), end
if (numel(SFLAG) == 1), SFLAG = [SFLAG, SFLAG]; end
SFLAG = logical(SFLAG);
if (nargin < 10) || isempty(CFLAG), CFLAG = 1; end
if (numel(CFLAG) > 1)
    error('amgkq: CFLAG must be either 0 or 1'), end
CFLAG = logical(CFLAG);
if (nargin < 11) || isempty(VERB), VERB = 0; end
VERB = round(VERB);
if (numel(VERB) > 1)
    error('amgkq: VERB must be an integer from -2 to 2'), end
if (VERB < -2) || (VERB > 2)
    error('amgkq: VERB must be an integer from -2 to 2'), end

% Set warning ids according to environment and verbosity
persistent iamoctave
if isempty(iamoctave), iamoctave = (exist('OCTAVE_VERSION','builtin') > 0); end
if iamoctave, warndiv0 = 'Octave:divide-by-zero';
else warndiv0 = 'MATLAB:divideByZero'; end
warninfnan = 'amgkq:inf-or-nan';
warnnsub = 'amgkq:insufficient-subregions';
warntol = 'amgkq:error-exceeds-tolerance';
warncull = 'amgkq:small-subregions-culled';

% Must save warning states because 'local' option is Octave specific
origwarnstate = warning;
warning('off',warndiv0);
if VERB < 0
    warning('off',warninfnan), warning('off',warnnsub), end
if VERB < -1
    warning('off',warntol), warning('off',warncull), end

% Convert function given as a string to a function handle
if ischar(F)
    F = @(X) feval(F, X, varargin{:});
elseif (nargin > 11)
    error('amgkq: additional parameters can only be specified for F given as a string')
end

% Account for limit ordering
if ISCONTOUR
    INTSGN = 1;
else
    INTSGN = prod(sign(B - A));
    if any(B < A)
        ABNEW = sort([A, B],2);
        A = ABNEW(:,1); B = ABNEW(:,2);
    end
end

% Account for presence of edge singularities
% and infinite limits (improper integrals)
% using standard variable transformations
ABXFORM = findinfdir(F, A, B, ND);
ABCHECK = ABXFORM + ABINF;
if any(ABCHECK(:) > 1)
    error('amgkq: integral is divergent'), end
ABXND = find(ABXFORM(:,1) & ABXFORM(:,2));
AXND = find(ABXFORM(:,1) & ~ABXFORM(:,2));
BXND = find(ABXFORM(:,2) & ~ABXFORM(:,1));
ABC = [A, B, C];
if ~isempty(ABXND)
    OLDF = F;
    ABCY = ABC(ABXND,:);
    BMA = B(ABXND) - A(ABXND);
    if (TFIN == 1)
        F = @(Y) bsxfun(@times,prod(bsxfun(@times,BMA,sin(Y(ABXND,:))/2),1),OLDF(xofyab1(Y, A(ABXND), BMA, ABXND)));
        ABCY = bsxfun(@rdivide,bsxfun(@minus,ABCY,A(ABXND)),BMA);
        ABC(ABXND,:) = 2*atan(sqrt(ABCY./(1 - ABCY)));
    else
        BPA = B(ABXND) + A(ABXND);
        F = @(Y) bsxfun(@times,prod(bsxfun(@times,BMA*3/4,(1 - Y(ABXND,:).^2)),1),OLDF(xofyab2(Y, BMA, BPA, ABXND)));
        ABCY = bsxfun(@rdivide,bsxfun(@minus,2*BPA,4*ABCY),BMA);
        ABCY = ((sqrt(ABCY.^2 - 4) + ABCY)/2).^(1/3);
        ABCY = (sqrt(-3)*(1 - ABCY.^2) - (1 + ABCY.^2))/2./ABCY;
        ABC(ABXND,:) = real(ABCY);
    end
end
if ~isempty(AXND)
    OLDF = F;
    F = @(Y) bsxfun(@times,prod(2*Y(AXND,:),1),OLDF(xofya(Y, ABC(AXND,1), AXND)));
    ABC(AXND,:) = sqrt(bsxfun(@minus,ABC(AXND,:),ABC(AXND,1)));
end
if ~isempty(BXND)
    OLDF = F;
    F = @(Y) bsxfun(@times,prod(-2*Y(BXND,:),1),OLDF(xofyb(Y, ABC(BXND,2), BXND)));
    ABC(BXND,:) = -sqrt(bsxfun(@minus,ABC(BXND,2),ABC(BXND,:)));
end
A = ABC(:,1); B = ABC(:,2); C = ABC(:,3:end);
if any(ABINF(:))
    ABINF = any(ABINF,2);
    OLDF = F;
    if (TINF == 1)
        F = @(Y) bsxfun(@times,prod(1./cos(Y(ABINF,:)).^2,1),OLDF(xofyinf1(Y, ABINF)));
        ABC(ABINF,:) = atan(ABC(ABINF,:));
        A = ABC(:,1); B = ABC(:,2); C = ABC(:,3:end);
    else
        F = @(Y) bsxfun(@times,prod((1 + Y(ABINF,:).^2)./(1 - Y(ABINF,:).^2).^2,1),OLDF(xofyinf2(Y, ABINF)));
        ABC(ABINF,:) = 2*ABC(ABINF,:)./(1 + sqrt(1 + 4*ABC(ABINF,:).^2));
        A(isinf(A)) = -1; B(isinf(B)) = 1; C = ABC(:,3:end);
    end
    % If we took default C then we need to fix NaN and Inf values
    if DEFC, C(ABINF) = (A(ABINF)+B(ABINF))/2; end
end

% Setup initial subregions
XAXS = amgkcoef(NGK); XRAT = (1 - XAXS(end))/2;
[CLS, HWS, NSUB] = subsetup(ND, NC, A, B, C, CFLAG, ISCONTOUR);
for NS = 1:NSUB
    [QV, E2, F4] = gkint(NGK, ND, CLS(:,NS), HWS(:,NS), F);
    RESS(:,NS) = QV; ER2S(:,NS) = E2; FD4S(:,:,NS) = F4;
    if VERB > 1, disp(' '), NS, QV, E2, F4, end
end
if SFLAG(2) && ~ISCONTOUR
    HWSTOL = eps(abs(CLS)+HWS); XRATHWS = XRAT*HWS;
elseif SFLAG(2)
    HWSTOL = eps(abs(CLS)+abs(HWS)); XRATHWS = XRAT*abs(HWS);
end
RES = sum(RESS,2); ER2 = sum(ER2S,2);
NF = size(QV,1); first = true;

% Loop while the error is too large and NSUB < MAXNSUB
FL = [];
CULLFLAG = false;
while NSUB <= MAXNSUB
    % Break if Inf encountered
    if any(isinf(RESS(:)))
        warning(warninfnan,'amgkq: Inf encountered')
        FL = -1; break
    end
    % Break if NaN encountered
    if any(isnan(RESS(:)))
        warning(warninfnan,'amgkq: NaN encountered')
        FL = -2; break
    end
    % Break if converged globally
    ABSRES = abs(RES);
    TOL = max(EA, ER*ABSRES);
    if all(ER2 < TOL.^2)
        FL = 2; break
    end
    % Accept the subintervals that meet the convergence criterion
    % Exactly what that criterion should be is not clear
    % The following seems to work well:
    if SFLAG(1), idx1 = find(all(bsxfun(@lt,ER2S,eps(ABSRES).^2),1)); else idx1 = []; end
    % We also need to cull any subintervals that are getting too small
    % Exactly what the definition of "too small" should be is not clear
    % The following seems to work well, at least for real integrals:
    if SFLAG(2), idx2 = find(any(XRATHWS < HWSTOL,1)); else idx2 = []; end
    if ~CULLFLAG && ~isempty(idx2), CULLFLAG = true; end
    % Combine indices for culling
    idx = [idx1, idx2];
    if first
        RESCULL = sum(RESS(:,idx),2);
        ER2CULL = sum(ER2S(:,idx),2);
        first = false;
    else
        RES = RESCULL + sum(RESS,2);
        ER2 = ER2CULL + sum(ER2S,2);
        RESCULL = RESCULL + sum(RESS(:,idx),2);
        ER2CULL = ER2CULL + sum(ER2S(:,idx),2);
    end
    if VERB > 1, disp(' '), NSUB, idx, RES, ER2
    elseif VERB > 0, disp(' '), NSUB, RES, ER2, end
    % Remove culled subintervals from consideration
    CLS(:,idx) = []; HWS(:,idx) = [];
    RESS(:,idx) = []; ER2S(:,idx) = []; FD4S(:,:,idx) = [];
    if SFLAG(2), HWSTOL(:,idx) = []; XRATHWS(:,idx) = []; end
    % Break if no subintervals left
    if isempty(CLS), FL = 1; break, end
    % Locate next subregion
    NSUB = NSUB + 1; SUBN = size(CLS,2) + 1;
    % Find subregion and integrand with greatest error
    [DUM, SUBNDX] = max(ER2S(:));
    [FK, SUBE] = ind2sub([NF, SUBN-1], SUBNDX);
    % Find direction for subdivision
    if (ND == 1), SDIR = 1; else
    [DUM, SDIR] = max(abs(squeeze(FD4S(FK,:,SUBE)))); end
    % Divide the subregion in two halves along SDIR
    HWS(SDIR,SUBE) = HWS(SDIR,SUBE)/2;
    CLS(:,SUBN) = CLS(:,SUBE); HWS(:,SUBN) = HWS(:,SUBE);
    CLS(SDIR,SUBE) = CLS(SDIR,SUBE) - HWS(SDIR,SUBE);
    CLS(SDIR,SUBN) = CLS(SDIR,SUBN) + HWS(SDIR,SUBN);
    % Compute integral and error over each half and store results
    [QV, E2, F4] = gkint(NGK, ND, CLS(:,SUBE), HWS(:,SUBE), F);
    RESS(:,SUBE) = QV; ER2S(:,SUBE) = E2; FD4S(:,:,SUBE) = F4;
    if VERB > 1, SUBE, QV, E2, F4, end
    [QV, E2, F4] = gkint(NGK, ND, CLS(:,SUBN), HWS(:,SUBN), F);
    RESS(:,SUBN) = QV; ER2S(:,SUBN) = E2; FD4S(:,:,SUBN) = F4;
    if VERB > 1, SUBN, QV, E2, F4, end
    % Also compute HWS tolerance
    if SFLAG(2) && ~ISCONTOUR
        HWSTOL(:,SUBE) = eps(abs(CLS(:,SUBE))+HWS(:,SUBE));
        HWSTOL(:,SUBN) = eps(abs(CLS(:,SUBN))+HWS(:,SUBN));
        XRATHWS(:,SUBE) = XRAT*HWS(:,SUBE);
        XRATHWS(:,SUBN) = XRAT*HWS(:,SUBN);
    elseif SFLAG(2)
        HWSTOL(:,SUBE) = eps(abs(CLS(:,SUBE))+abs(HWS(:,SUBE)));
        HWSTOL(:,SUBN) = eps(abs(CLS(:,SUBN))+abs(HWS(:,SUBN)));
        XRATHWS(:,SUBE) = XRAT*abs(HWS(:,SUBE));
        XRATHWS(:,SUBN) = XRAT*abs(HWS(:,SUBN));
    end
end

% Now check for warnings
if isempty(FL)
    warning(warnnsub,['amgkq: maximum number of subregions ' num2str(MAXNSUB) ' evaluated'])
    NSUB = NSUB - 1; FL = 0;
end
if CULLFLAG
    warning(warncull,'amgkq: subregions smaller than width tolerance have been culled')
end

% Account for sign and sqrt and check error tolerance
RES = INTSGN*RES;
ERR = sqrt(ER2);
TOL = max(EA, ER*abs(RES));
if any(ERR > TOL)
    warning(warntol,['amgkq: error tolerance ' num2str(max(TOL)) ' not satisfied by ' num2str(ERR')])
end

% Restore warning states because 'local' option is Octave specific
if VERB < 0
    warning('on',warninfnan), warning('on',warnnsub), end
if VERB < -1
    warning('on',warntol), warning('on',warncull), end
warning(origwarnstate)

end % function % amgkq

function ABXFORM = findinfdir(F, A, B, ND);
% Return divergent directions for F near A and B
% Checks that "near A" is also "far from B" in the infinite case
isinfA = isinf(A); isinfB = isinf(B);
AK = A; A(isinfA) = min(-1,B(isinfA))/eps;
BK = B; B(isinfB) = max(1,A(isinfB))/eps;
BMA = B - A;
NEARA = A + BMA/100;
X = NEARA(:,ones(1,ND));
for k = 1:ND
    X(k,k) = AK(k);
end
FA = F(X);
NEARB = B - BMA/100;
X = NEARB(:,ones(1,ND));
for k = 1:ND
    X(k,k) = BK(k);
end
FB = F(X);
FAINF = any(isinf(FA),1);
FBINF = any(isinf(FB),1);
ABXFORM = [FAINF; FBINF].';
end % function % findinfdir

function X = xofyab1(Y, A, BMA, ABXND)
X = Y;
X(ABXND,:) = bsxfun(@plus,A,bsxfun(@times,BMA,(1-cos(Y(ABXND,:)))/2));
end % function % xofyab1

function X = xofyab2(Y, BMA, BPA, ABXND)
X = Y;
X(ABXND,:) = bsxfun(@plus,bsxfun(@times,BMA/4,Y(ABXND,:).*(3 - Y(ABXND,:).^2)),BPA/2);
end % function % xofyab2

function X = xofya(Y, A, AXND)
X = Y;
X(AXND,:) = bsxfun(@plus,A,Y(AXND,:).^2);
end % function % xofya

function X = xofyb(Y, B, BXND)
X = Y;
X(BXND,:) = bsxfun(@minus,B,Y(BXND,:).^2);
end % function % xofyb

function X = xofyinf1(Y, ABINF)
X = Y;
X(ABINF,:) = tan(Y(ABINF,:));
end % function % xofyinf1

function X = xofyinf2(Y, ABINF)
X = Y;
X(ABINF,:) = Y(ABINF,:)./(1 - Y(ABINF,:).^2);
end % function % xofyinf2

function [CLS, HWS, NSUB] = subsetup(ND, NC, A, B, C, CFLAG, ISCONTOUR)
% Create initial subregions
if ISCONTOUR
    ACB = [A, C, B];
    CLS = (ACB(2:end)+ACB(1:end-1))/2;
    HWS = (ACB(2:end)-ACB(1:end-1))/2;
    NSUB = size(CLS,2);
    return
end
if CFLAG
    CPERMS = perms(1:NC);
    NPERMS = size(CPERMS,1);
    DACB = zeros(1,NPERMS);
    for k = 1:NPERMS
        ACB = [A, C(:,CPERMS(k,:)), B];
        DACB(k) = sum(sqrt(sum(diff(ACB,[],2).^2,1)));
    end
    [DUM, KACB] = min(DACB);
    C = C(:,CPERMS(KACB,:));
end
CLS = []; HWS = [];
for l = 1:NC
    CL = C(:,l);
    if l > 1
        absdiffCLCLS = abs(bsxfun(@minus,CLS,CL));
        FC = find(all(absdiffCLCLS <= HWS,1),1);
        if isempty(FC)  % if at first you don't succeed,
            FC = find(all((absdiffCLCLS - HWS) < eps,1),1);
        end
        if isempty(FC)  % try trying again
            FC = find(all((absdiffCLCLS - HWS) < eps*5,1),1);
        end
        A = CLS(:,FC) - HWS(:,FC); B = CLS(:,FC) + HWS(:,FC);
        CLS(:,FC) = []; HWS(:,FC) = [];
    end
    D = [A+CL, CL+B]/2; E = abs([A-CL, CL-B])/2;
    CLSL = D(:,1); HWSL = E(:,1);
    for k = 1:ND
        CLSL = [CLSL, CLSL]; HWSL = [HWSL, HWSL];
        CLSL(k,1:end/2) = D(k,1); CLSL(k,end/2+1:end) = D(k,2); 
        HWSL(k,1:end/2) = E(k,1); HWSL(k,end/2+1:end) = E(k,2); 
    end
    CLS = [CLS, CLSL]; HWS = [HWS, HWSL];
    NULLV = isinf(log(prod(HWS,1)));
    CLS(:,NULLV) = []; HWS(:,NULLV) = [];
end
NSUB = size(CLS,2);
end % function % subsetup

function [QK, E2, F4] = gkint(NGK, ND, CLS, HWS, F)
% Gauss-Kronrod quadrature of order NG = NGK and NK = 2*NG+1
persistent NGKLAST XGK WK WG XSTAR X
persistent NDLAST XIND XGCOL WKX WGX XD4ND WD4
if isempty(NGKLAST) || (NGKLAST ~= NGK)
    NGKLAST = NGK;
    [XGK, WK, WG] = amgkcoef(NGK);
    NDLAST = [];
end
if isempty(NDLAST) || (NDLAST ~= ND)
    NDLAST = ND;
    NG = NGK;
    NK = 2*NG + 1;
    NXK = NK^ND;
    NXG = NG^ND;
    X = zeros(ND,NXK); 
    XSUB = cell(ND,1);
    [XSUB{:}] = ind2sub(NK*ones(1,ND),1:NXK);
    XIND = cell2mat(XSUB);
    XGCOL = find(all(~mod(XIND,2),1));
    WKX = ones(NXK,1);
    WGX = ones(NXG,1);
    XD4ND = zeros(ND,NK);
    for KND = 1:ND
        WKX = WKX.*WK(XIND(KND,:));
        WGX = WGX.*WG(XIND(KND,XGCOL));
        % The following works even for ND = 1
        NKND = 1:ND; NKND(KND) = [];
        XD4ND(KND,:) = find(all(XIND(NKND,:) == (NG + 1),1));
    end
    WD4 = amgkfdc(4,0,XGK);
    XSTAR = zeros(ND,NK);
end
VOL = prod(HWS);
XSTAR = bsxfun(@plus, HWS * XGK, CLS);
for KND = 1:ND
    X(KND,:) = XSTAR(KND,XIND(KND,:));
end
Y = F(X);
QK = Y*WKX*VOL;
QG = Y(:,XGCOL)*WGX*VOL;
E2 = (QK - QG).^2;
F4 = 0*Y(:,1:ND);
if ND > 1
    % Should be vectorized
    for KND = 1:ND
        F4(:,KND) = Y(:,XD4ND(KND,:))*WD4;
    end
end
end % function % gkint

% AMGKQ: [RES, ERR, NSUB, FL] = amgkq(F, A, B, C, EAER, MAXNSUB, NGK, TTYPE, SFLAG, CFLAG, VERB, P1, P2, ...)

%!assert(amgkq('sin',-pi,pi,[],[],[],7),0,1e-14)
%!assert(amgkq(inline('sin'),-pi,pi,[],[],[],8),0,1e-14)
%!assert(amgkq(@sin,-pi,pi,[],[],[],9),0,1e-14)
%!assert(amgkq(@(x) sin(x),-pi,pi,[],[],[],10),0,1e-14)

%!assert(amgkq(@(x) exp(-x.^2),-inf,inf),sqrt(pi),1e-14)
%!assert(amgkq(@(x) exp(-x.^2).*cos(x),0,inf),sqrt(pi)/2/exp(1/4),1e-14)
%!assert(amgkq(@(x) exp(-x.^2)./(1+x.^2),0,1),pi/4*exp(1)*(1-erf(1)^2),1e-14)

%!assert(amgkq(@(x) exp(-x.^2),-inf,inf,[],[],[],[],2),sqrt(pi),1e-14)
%!assert(amgkq(@(x) 1./(sqrt(x.*(1-x))),0,1,[],[],[],[],2),pi,1e-13)

%!assert(amgkq(@(x) exp(-x).*x./(1-exp(-2*x)),0,inf),pi^2/8,1e-12)
%!assert(amgkq(@(x) exp(-x).*x./(1-exp(-2*x)),0,inf,[],0,1e3),pi^2/8,1e-14)
%!assert(amgkq(@(x) log(x)./(1-x.^2),0,1),-pi^2/8,1e-9)
%!assert(amgkq(@(x) log(x)./(1-x.^2),0,1,[],0,1e3),-pi^2/8,1e-14)

%!assert(amgkq(@(x) 1./sqrt(abs(x)),0,10),2*sqrt(10),1e-14)
%!assert(amgkq(@(x) 1./sqrt(abs(x)),-10,10,[],0,1e3),4*sqrt(10),1e-14)

%!assert(amgkq(@(x) 1./(sqrt (x).*(1+x)),0,inf),pi,1e-14)
%!assert(amgkq(@(x) 1./(sqrt (-x).*(1-x)),-inf,0),pi,1e-14)

%!assert(amgkq(@(x) bsxfun(@times,exp(-x),[x; x.^2; x.^3; x.^4; x.^5]),0,inf),gamma(2:6)',1e-12)

%!assert(amgkq(@(x) sin(3*x).*cosh(x).*sinh(x),10,15,[],[0 1e-14]),2.588424538641647e+10,-1e-14)
%!assert(amgkq(@(x) 1./sqrt(abs(x(1,:))).*sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:)),[0;10],[10;15],[],[0 1e-14]),2*sqrt(10)*2.588424538641647e+10,-1e-14)
%!assert(amgkq(@(x) 1./sqrt(abs(x(1,:))).*sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:)),[-10;10],[10;15],[],[0 1e-14],2e3),4*sqrt(10)*2.588424538641647e+10,-1e-14)
%!assert(amgkq(@(x) [1./sqrt(abs(x(1,:))); sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:))],[0;10],[10;15],[],[0 1e-14]),[2*sqrt(10)*5;10*2.588424538641647e+10],-1e-14)
%!assert(amgkq(@(x) [1./sqrt(abs(x(1,:))); sin(3*x(2,:)).*cosh(x(2,:)).*sinh(x(2,:))],[-10;10],[10;15],[],[],2e3),[4*sqrt(10)*5;20*2.588424538641647e+10],-5e-7)

%!assert(amgkq(@(z) [1./(1+z.^2).^2; exp(i*z)./(1+z.^2)],-1,-1,[1,2*i]),pi*[1/2; exp(-1)],1e-13)
%!error  amgkq(@(z) [1./(1+z.^2).^2; exp(i*z)./(1+z.^2)],-inf+i/2,-inf+i/2,[inf+i/2,2*i])

%!assert(amgkq(@(x) prod(exp(-x.^2),1),[-inf;inf],[inf;-inf],[],0),-pi,1e-14)
%!assert(amgkq(@(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2),[-inf;-inf],[inf;inf]),sqrt(2*pi)*pi,1e-13)

%!shared p, q
%! p = 1/2; q = 1/2;
%!assert(amgkq(@(x) x.^(p-1).*(1-x).^(q-1),0,1),beta(p,q),1e-14)
%! p = 1/3; q = 1/3;
%!assert(amgkq(@(x) x.^(p-1).*(1-x).^(q-1),0,1),beta(p,q),2e-5)
%! p = 1/4; q = 1/4;
%!assert(amgkq(@(x) x.^(p-1).*(1-x).^(q-1),0,1),beta(p,q),5e-4)

%!assert(amgkq(@(x) sin(x)./x,0,inf,[],[],1e3,[],[],0),pi/2,2e-2)
%!assert(amgkq(@(x) (sin(x)./x).^2,0,inf,[],[],1e3),pi/2,2e-6)
%!assert(amgkq(@(x) (sin(x)./x).^3,0,inf,[],[],1e3),3*pi/8,1e-7)
%!assert(amgkq(@(x) (sin(x)./x).^4,0,inf,[],[],1e3),pi/3,1e-8)

%!shared AB
%! AB = [-1 1; -1 1];
%!assert(amgkq(@(x) [(sum(x.^2,1) < 1); (sum(x.^2,1) > 1)],AB(:,1),AB(:,2),[],[],2e3),[pi; 4-pi],1e-5)
%!assert(amgkq(@(x) [(sum(x.^2,1) < 1); (sum(x.^2,1) > 1)],AB(:,1),AB(:,2),[],[],1e3,40),[pi; 4-pi],1e-5)

%!shared AB
%! AB = [-1 1; -1 1]*10;
%!assert(amgkq(@(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2),AB(:,1),AB(:,2)),sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*2*atan(AB(2,2)),1e-10)
%!assert(amgkq(@(x) [exp(-x(1,:).^2/2); 1./(1+x(2,:).^2)],AB(:,1),AB(:,2)),[sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*diff(AB(2,:));2*atan(AB(2,2))*diff(AB(1,:))],1e-13)

%!shared AB
%! AB = [[-1 1; -1 1]*10; 0 1];
%!assert(amgkq(@(x) exp(-x(1,:).^2/2)./(1+x(2,:).^2).*x(3,:).^10.*(1-x(3,:)).^10,AB(:,1),AB(:,2),[],0),sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*2*atan(AB(2,2))*beta(11,11),1e-15)
%!assert(amgkq(@(x) [exp(-x(1,:).^2/2); 1./(1+x(2,:).^2); x(3,:).^10.*(1-x(3,:)).^10],AB(:,1),AB(:,2)),[sqrt(2*pi)*erf(AB(1,2)/sqrt(2))*diff(AB(2,:));2*atan(AB(2,2))*diff(AB(1,:));beta(11,11)*prod(diff(AB(1:2,:),[],2))],1e-12)

