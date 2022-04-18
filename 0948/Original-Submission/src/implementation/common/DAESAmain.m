function sadata = DAESAmain(fcn, n, varargin)

global problemSize; problemSize = n;

t = cputime;

% Compute sigma
[sigma, exitflag, fNotEval, varNotUsed] = compSigma(fcn, n, varargin{:});

if exitflag==-3 % more equations than required, syntax error
     sadata = OutputError([], -3, fNotEval, varNotUsed, n, fcn, ...
         cputime-t, [], [], [], [], varargin{:});
     return;
end

[c, d, ccl, dcl, HVT, DOF, p, q, rC, cc, rr, tmpflag] = lapdm(sigma);

exitflag = min(exitflag, tmpflag);
if exitflag<0 
    % -1 -> solving LAP fails
    % -2 -> missing equation/variable
    sadata = OutputError(sigma, exitflag, fNotEval, varNotUsed, n, fcn, ...
        cputime-t, p, q, rr, cc, varargin{:});
    return;
end
% Debugging issues:
% exitflag = -3, more equations than needed
%          = -2, missing equation or variable
%          = -1, failure in solving LAP (under- & over-determined)
%          =  0, SWP

index = max(c) + (~isempty(find(d == 0,1)));
% Calculate the index

[cfl, dfl, sigmaPerm, p, q, rF, JNZ] = dm2blockSP(sigma, c, d, p, q);
% Returns local local offsets, permuted sigma, 
%   boundaries of blocks, permutations, essential sparsity pattern of Jacobian
[IVset, constraintSet, alpha, alphaf, lin] = ...
    QuasiLinAnal(fcn, n, sigmaPerm, c, ccl, cfl, dfl, p, q, ...
    rC, rF, JNZ, varargin{:});
             
% Do Quasi-Linear Analysis for each block
% Returns sets of IVs, sets of constraints, quasi-linearity in each block
%
% Do Linearity Analysis for each equation
% Returns linearity results for each equation, stored in lin
sadata = SAdata(n, fcn, sigma, DOF, HVT, index, c, d, ccl, dcl, cfl, dfl, ...
            IVset, constraintSet', p, q, rC, rF, alpha, alphaf, lin, ...
            JNZ, true, 0, cputime-t, [], [], ...
            varargin);
end

function sadata = OutputError(sigma, exitflag, fNotEval, varNotUsed, n, fcn, ...
    time, p, q, r, s, varargin)
    sadata = SAdata(n, fcn, sigma, NaN, NaN, NaN, [], [], [], [], [], [], ...
            [], [], p, q, r, s, [], [], [], ...
            [], false, exitflag, time, fNotEval, varNotUsed, ...
            varargin);
end